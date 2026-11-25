"""Task 5 (predeclared): F8-OLS, F8+C-OLS, F8+C-Ridge, F8+B+C-Ridge on the exact Phase 4 outer folds.
Ridge alpha grid = Phase 4 grid, tuned by inner GroupKFold(4) on training genes. Saves OOF predictions,
coefficients, per-fold condition numbers; paired dRMSE vs F8-OLS (descriptive; folds are reused)."""
import json

import numpy as np
import pandas as pd
from c00_common import (
    RAW_COLS,
    REPS,
    RS,
    RepBuilder,
    T,
    dl_tau2,
    nb_interval,
    outer_splits,
    raw_table,
    target,
)
from joblib import Parallel, delayed
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

# F8+B+C is identical to Phase 4B "FULL-A" (F8 + block B + block C); reported under the label F8+B+C
assert REPS["FULL-A"] == REPS["F8"] + REPS["F8+B"][8:] + REPS["F8+C"][8:]
MODELS = [("F8", "ols"), ("F8+C", "ols"), ("F8+C", "ridge"), ("FULL-A", "ridge")]
ALPHAS = [0.01, 0.1, 1, 10, 100, 1000]


def run(tgt, rep, mo, r_i, k_i):
    X0 = raw_table()
    j = X0.join(target(tgt), how="inner").dropna(subset=["y", "s"])
    j = j[j.s > 0]
    for r, k, tr, te, grp in outer_splits(j.index):
        if (r, k) == (r_i, k_i):
            break
    X, y, s = j[RAW_COLS].values, j.y.values, j.s.values
    w = 1 / (s**2 + dl_tau2(y[tr], s[tr]))
    est = Pipeline([("rep", RepBuilder(rep)), ("scale", StandardScaler()),
                    ("model", LinearRegression() if mo == "ols" else Ridge())])
    alpha = None
    if mo == "ridge":
        gs = GridSearchCV(est, {"model__alpha": ALPHAS}, cv=GroupKFold(4), scoring="neg_mean_squared_error")
        gs.fit(X[tr], y[tr], groups=grp[tr], model__sample_weight=w[tr])
        est, alpha = gs.best_estimator_, gs.best_params_["model__alpha"]
    else:
        est.fit(X[tr], y[tr], model__sample_weight=w[tr])
    Z = est[:-1].transform(X[tr])
    sv = np.linalg.svd(Z, compute_uv=False)
    pred = est.predict(X[te])
    yt, null = y[te], np.average(y[tr], weights=w[tr])
    met = dict(target=tgt, rep=rep, model=mo, repeat=r_i, fold=k_i, n_train=int(tr.sum()), n_test=int(te.sum()),
               alpha=alpha, cond_number=float(sv[0] / sv[-1]), rmse=float(np.sqrt(((yt - pred) ** 2).mean())),
               r2_oos=float(1 - ((yt - pred) ** 2).sum() / ((yt - null) ** 2).sum()),
               spearman=float(spearmanr(yt, pred).statistic), pred_sd=float(pred.std()), obs_sd=float(yt.std()))
    coef = pd.DataFrame({"target": tgt, "rep": rep, "model": mo, "repeat": r_i, "fold": k_i, "feature": REPS[rep],
                         "coef_per_SD": est.named_steps["model"].coef_})
    oof = pd.DataFrame({"gene": j.index[te], "y": yt, "pred": pred, "target": tgt, "rep": rep, "model": mo,
                        "repeat": r_i, "fold": k_i})
    return met, coef, oof


tasks = [(t, rep, mo, r, k) for t in ("PE_in_CD", "PE_in_MD") for rep, mo in MODELS for r in range(3) for k in range(5)]
out = Parallel(n_jobs=19)(delayed(run)(*x) for x in tasks)
M = pd.DataFrame([o[0] for o in out])
M["label"] = M.rep.replace({"FULL-A": "F8+B+C"})
M.to_csv(T / "fix_experiment_metrics_by_fold.csv", index=False)
pd.concat([o[1] for o in out]).to_csv(T / "fix_experiment_coefficients.csv", index=False)
pd.concat([o[2] for o in out]).to_csv(T / "fix_experiment_oof_predictions.csv.gz", index=False)

# historical comparators from Phase 4B (not rerun)
h = pd.read_csv(RS / "roadmap_phase4b" / "results" / "metrics_main.csv")
h = h[(h.rep.isin(["FULL"]) & h.model.isin(["ols", "ridge"])) | ((h.rep == "F8") & (h.model == "catboost"))]
h = h[["target", "rep", "model", "repeat", "fold", "n_train", "n_test", "rmse", "r2_oos"]].assign(historical=True)
A = pd.concat([M.assign(historical=False), h], ignore_index=True)
rows = []
base = A[(A.rep == "F8") & (A.model == "ols")]
for (tgt, rep, mo, hist), d in A.groupby(["target", "rep", "model", "historical"]):
    d = d.sort_values(["repeat", "fold"])
    b = base[base.target == tgt].sort_values(["repeat", "fold"])
    dd = d.rmse.values - b.rmse.values
    m, lo, hi = nb_interval(dd, d.n_test.values, d.n_train.values)
    rows.append(dict(target=tgt, rep=rep, model=mo, historical=hist, rmse=d.rmse.mean(), r2_oos=d.r2_oos.mean(),
                     r2_oos_sd=d.r2_oos.std(), d_rmse_vs_F8_ols=m, ci_low=lo, ci_high=hi, n_better=int((dd < 0).sum()),
                     d_r2_vs_F8_ols=float((d.r2_oos.values - b.r2_oos.values).mean()),
                     cond_number_median=float(d.cond_number.median()) if "cond_number" in d and d.cond_number.notna().any() else np.nan,
                     alpha_choices=json.dumps(d.alpha.value_counts().to_dict()) if "alpha" in d and d.alpha.notna().any() else ""))
S = pd.DataFrame(rows).sort_values(["target", "historical", "rmse"])
S.to_csv(T / "fix_experiment_summary.csv", index=False)
# reproducibility vs Phase 4B for overlapping configurations
p4b = pd.read_csv(RS / "roadmap_phase4b" / "results" / "metrics_main.csv")
rep_chk = M.merge(p4b, on=["target", "rep", "model", "repeat", "fold"], suffixes=("", "_4b"))
json.dump({"n_overlap": len(rep_chk), "max_abs_rmse_diff_vs_phase4b": float((rep_chk.rmse - rep_chk.rmse_4b).abs().max()) if len(rep_chk) else None},
          open(T / "fix_experiment_reproducibility.json", "w"), indent=2)
with pd.option_context("display.width", 220):
    print(S.round(4).to_string(index=False))
print(open(T / "fix_experiment_reproducibility.json").read())
