"""Task A runner: nested grouped CV for every (target, variant, split, model) task, in parallel.

Usage:  python p4_run.py [--jobs 18] [--only primary,positive,neg,variants,loo]
Outputs (results/taskA/):
  metrics.csv            one row per task (outer held-out metrics, chosen hyperparameters, timing)
  predictions/<run>.csv.gz  fold-by-fold held-out predictions for primary/positive-control runs
  coefficients.csv       linear-model coefficients per split (stability)
  importance.csv         permutation importance on held-out genes (primary runs; ridge/rf/hgb)
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from p4_data import (
    RES,
    VARIANTS,
    dl_tau2,
    feature_names,
    feature_table,
    gene_groups,
    loo_names,
    outer_splits,
    target,
)
from p4_models import BASELINES, LOO_SET, REDUCED, TUNED, make
from scipy.stats import spearmanr
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV, GroupKFold

OUT = RES / "taskA"
(OUT / "predictions").mkdir(parents=True, exist_ok=True)


def build_runs(only):
    runs = []  # (run_id, target_name, variant, models, save_pred)
    if "primary" in only:
        for t in ("PE_in_CD", "PE_in_MD"):
            runs.append((f"primary__{t}", t, "primary", BASELINES + TUNED, True))
    if "positive" in only:
        runs.append(("positive__Diet_in_Veh", "Diet_in_Veh", "primary", BASELINES + TUNED, True))
    if "neg" in only:
        for t in ("PE_in_CD", "PE_in_MD"):
            runs.append((f"neg_permuted__{t}", f"PERM:{t}", "primary", REDUCED, False))
        for d in ("CD", "MD"):
            runs.append((f"neg_depth__DEPTH_{d}", f"DEPTH_{d}", "primary", REDUCED, False))
    if "variants" in only:
        for v in VARIANTS:
            if v == "primary":
                continue
            for t in ("PE_in_CD", "PE_in_MD"):
                runs.append((f"{v}__{t}", t, v, REDUCED, False))
    if "loo" in only:
        for d in loo_names():
            for t in ("PE_in_CD", "PE_in_MD"):
                runs.append((f"loo_{d}__{t}", f"LOO:{d}:{t}", "primary", LOO_SET, False))
    return runs


def prepare(target_name, variant):
    X = feature_table(variant)
    t = target(target_name)
    j = X.join(t, how="inner").dropna(subset=["y", "s"])
    j = j[j["s"] > 0]
    return j


def run_task(run_id, target_name, variant, model, rep, fold, save_pred, data_cache={}):
    key = (target_name, variant)
    if key not in data_cache:
        data_cache[key] = prepare(target_name, variant)
    j = data_cache[key]
    feats = feature_names(variant)
    grouping = VARIANTS[variant].get("grouping", "module")
    for r, k, tr, te, grp in outer_splits(j.index, grouping):
        if r == rep and k == fold:
            break
    X, y, s = j[feats].values, j["y"].values, j["s"].values
    tau2 = dl_tau2(y[tr], s[tr])
    w = 1.0 / (s**2 + tau2)
    est, grid = make(model, feats)
    t0 = time.time()
    if grid is None:
        est.fit(X[tr], y[tr], model__sample_weight=w[tr])
        best = {}
    else:
        gs = GridSearchCV(est, grid, cv=GroupKFold(n_splits=4), scoring="neg_mean_squared_error", n_jobs=1, refit=True)
        gs.fit(X[tr], y[tr], groups=grp[tr], model__sample_weight=w[tr])
        est, best = gs.best_estimator_, gs.best_params_
    pred = est.predict(X[te])
    null = np.average(y[tr], weights=w[tr])
    yt = y[te]
    sse, sse0 = ((yt - pred) ** 2).sum(), ((yt - null) ** 2).sum()
    wt = w[te]
    wsse, wsse0 = (wt * (yt - pred) ** 2).sum(), (wt * (yt - null) ** 2).sum()
    deg = j["deg"].values[te] if "deg" in j else np.zeros(te.sum(), bool)
    row = dict(run=run_id, target=target_name, variant=variant, model=model, repeat=rep, fold=fold,
               n_train=int(tr.sum()), n_test=int(te.sum()), tau2=tau2,
               rmse=float(np.sqrt(mean_squared_error(yt, pred))), mae=float(mean_absolute_error(yt, pred)),
               r2_oos=float(1 - sse / sse0), r2_heldout_mean=float(r2_score(yt, pred)), wr2_oos=float(1 - wsse / wsse0),
               spearman=float(spearmanr(yt, pred).statistic) if np.std(pred) > 0 else np.nan,
               sign_agree=float((np.sign(yt) == np.sign(pred - 0)).mean()),
               sign_agree_deg=float((np.sign(yt[deg]) == np.sign(pred[deg])).mean()) if deg.any() else np.nan,
               n_deg_test=int(deg.sum()), best_params=json.dumps(best), seconds=round(time.time() - t0, 2))
    extra = {}
    if save_pred:
        extra["pred"] = pd.DataFrame({"gene": j.index[te], "y": yt, "pred": pred, "run": run_id, "model": model,
                                      "repeat": rep, "fold": fold})
        m = est.named_steps["model"] if hasattr(est, "named_steps") else est
        if hasattr(m, "coef_") and model in ("ols_all", "ridge", "elasticnet"):
            extra["coef"] = pd.DataFrame({"run": run_id, "model": model, "repeat": rep, "fold": fold,
                                          "feature": feats, "coef_per_SD": np.ravel(m.coef_)})
        if model in ("ols_all", "ridge", "rf", "hgb", "catboost"):
            pi = permutation_importance(est, X[te], yt, n_repeats=5, random_state=0, scoring="neg_mean_squared_error", n_jobs=1)
            extra["imp"] = pd.DataFrame({"run": run_id, "model": model, "repeat": rep, "fold": fold, "feature": feats,
                                         "perm_importance_mse": pi.importances_mean})
    return row, extra


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=18)
    ap.add_argument("--only", default="primary,positive,neg,variants,loo")
    a = ap.parse_args()
    only = a.only.split(",")
    gene_groups()  # build/cache once before forking
    tasks = []
    heavy = {"rf": 0, "catboost": 1, "hgb": 2}
    for run_id, tn, v, models, save in build_runs(only):
        # Deviation D1 (compute only; logged in the decision report): 3 repeats for primary/positive runs,
        # 1 repeat (repeat 0, 5 folds) for variants, negative controls and LOO.
        reps = [0, 1, 2] if run_id.startswith(("primary__", "positive__")) else [0]
        for m in models:
            for r in reps:
                for k in range(5):
                    tasks.append((run_id, tn, v, m, r, k, save))
    tasks.sort(key=lambda t: heavy.get(t[3], 9))
    print(f"{len(tasks)} tasks on {a.jobs} workers", flush=True)
    t0 = time.time()
    out = Parallel(n_jobs=a.jobs, backend="loky", verbose=5, batch_size=1)(delayed(run_task)(*t) for t in tasks)
    rows = pd.DataFrame([o[0] for o in out])
    tag = "_".join(only)
    rows.to_csv(OUT / f"metrics_{tag}.csv", index=False)
    for kind in ("pred", "coef", "imp"):
        fr = [o[1][kind] for o in out if kind in o[1]]
        if fr:
            df = pd.concat(fr)
            if kind == "pred":
                for rid, d in df.groupby("run"):
                    d.to_csv(OUT / "predictions" / f"{rid}.csv.gz", index=False)
            else:
                df.to_csv(OUT / f"{'coefficients' if kind == 'coef' else 'importance'}_{tag}.csv", index=False)
    print(f"done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
