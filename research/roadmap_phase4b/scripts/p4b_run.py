"""Phase 4B runner: representation x model x target x Phase-4 outer split, nested GroupKFold tuning.
Usage: python p4b_run.py --stage main|controls [--jobs 19]"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

P4B = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P4B / "scripts"))
from p4b_features import RAW_COLS, REPS, RepBuilder, raw_table  # noqa: E402

sys.path.insert(0, str(P4B.parent / "roadmap_phase4" / "scripts"))
from p4_data import dl_tau2, outer_splits, target  # noqa: E402
from p4_models import make  # noqa: E402

OUT = P4B / "results"
MODELS = {"ols": "ols_all", "ridge": "ridge", "elasticnet": "elasticnet", "hgb": "hgb", "catboost": "catboost"}
_CACHE = {}


def data(target_name, shuffle_b):
    k = (target_name, shuffle_b)
    if k not in _CACHE:
        X = raw_table(shuffle_b)
        t = target(target_name)
        j = X.join(t, how="inner").dropna(subset=["y", "s"])
        _CACHE[k] = j[j["s"] > 0]
    return _CACHE[k]


def run_task(stage, rep, model, target_name, rep_i, fold, shuffle_b=False):
    j = data(target_name, shuffle_b)
    for r, k, tr, te, grp in outer_splits(j.index, "module"):
        if r == rep_i and k == fold:
            break
    X, y, s = j[RAW_COLS].values, j["y"].values, j["s"].values
    w = 1.0 / (s**2 + dl_tau2(y[tr], s[tr]))
    base, grid = make(MODELS[model], RAW_COLS)
    est = Pipeline([("rep", RepBuilder(rep)), ("scale", StandardScaler()), ("model", base.named_steps["model"])])
    t0 = time.time()
    if grid is None:
        est.fit(X[tr], y[tr], model__sample_weight=w[tr])
        best = {}
    else:
        gs = GridSearchCV(est, grid, cv=GroupKFold(4), scoring="neg_mean_squared_error", n_jobs=1, refit=True)
        gs.fit(X[tr], y[tr], groups=grp[tr], model__sample_weight=w[tr])
        est, best = gs.best_estimator_, gs.best_params_
    pred = est.predict(X[te])
    yt, null = y[te], np.average(y[tr], weights=w[tr])
    sse, sse0 = ((yt - pred) ** 2).sum(), ((yt - null) ** 2).sum()
    wt = w[te]
    row = dict(stage=stage, rep=rep, model=model, target=target_name, shuffle_b=shuffle_b, repeat=rep_i, fold=fold,
               n_train=int(tr.sum()), n_test=int(te.sum()), n_features=len(REPS[rep]),
               rmse=float(np.sqrt(mean_squared_error(yt, pred))), mae=float(mean_absolute_error(yt, pred)),
               r2_oos=float(1 - sse / sse0), wr2_oos=float(1 - (wt * (yt - pred) ** 2).sum() / (wt * (yt - null) ** 2).sum()),
               spearman=float(spearmanr(yt, pred).statistic) if np.std(pred) > 0 else np.nan,
               best_params=json.dumps(best), seconds=round(time.time() - t0, 2))
    coef = None
    m = est.named_steps["model"]
    if hasattr(m, "coef_") and repeat_is_primary(stage):
        coef = pd.DataFrame({"rep": rep, "model": model, "target": target_name, "repeat": rep_i, "fold": fold,
                             "feature": REPS[rep], "coef_per_SD": np.ravel(m.coef_)})
    return row, coef


def repeat_is_primary(stage):
    return stage == "main"


def tasks_for(stage):
    t = []
    if stage == "main":
        for rep in REPS:
            for model in MODELS:
                for tgt in ("PE_in_CD", "PE_in_MD"):
                    for r in range(3):
                        for k in range(5):
                            t.append(("main", rep, model, tgt, r, k, False))
    elif stage == "controls":
        for model in ("ridge", "hgb"):
            for tgt in ("PE_in_CD", "PE_in_MD"):
                for k in range(5):
                    t.append(("ctrl_permuted_target", "FULL", model, f"PERM:{tgt}", 0, k, False))
                    t.append(("ctrl_shuffled_B", "F8+B", model, tgt, 0, k, True))
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="main")
    ap.add_argument("--jobs", type=int, default=19)
    a = ap.parse_args()
    tasks = tasks_for(a.stage)
    order = {"catboost": 0, "hgb": 1}
    tasks.sort(key=lambda t: order.get(t[2], 9))
    print(f"{len(tasks)} tasks", flush=True)
    t0 = time.time()
    out = Parallel(n_jobs=a.jobs, backend="loky", verbose=5, batch_size=1)(delayed(run_task)(*t) for t in tasks)
    pd.DataFrame([o[0] for o in out]).to_csv(OUT / f"metrics_{a.stage}.csv", index=False)
    co = [o[1] for o in out if o[1] is not None]
    if co:
        pd.concat(co).to_csv(OUT / f"coefficients_{a.stage}.csv", index=False)
    print(f"done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
