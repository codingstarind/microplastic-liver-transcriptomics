"""Lipid-gene prediction follow-up, step 1: assemble held-out (out-of-fold) gene-level predictions on the locked
Phase 4 co-expression-module outer splits (3 repeats x 5 folds, seeds 101-103, 20 Ward modules).

- F8-OLS, F8+C-OLS, F8+C-Ridge: read unchanged from Phase 4C `fix_experiment_oof_predictions.csv.gz`
  (Phase 4C test: reproduces Phase 4B exactly).
- FULL (28 features; F8+A+B+C) OLS and Ridge: Phase 4B saved only per-split metrics, not predictions. They are
  regenerated here with the unchanged Phase 4B code path (p4b_run.data / RepBuilder / p4_models.make / nested
  GroupKFold(4) tuning / DL weights) and checked against the stored Phase 4B per-split metrics.
- Null baseline per split: weighted training-gene mean, as in Phase 4/4B (w = 1/(s^2 + tau2_DL)).
Nothing earlier is modified; outputs go to research/lipid_prediction_followup/results/.
Run: .venv/Scripts/python.exe research/lipid_prediction_followup/scripts/lp01_oof_predictions.py
"""
import json
import os
import sys
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"
RS = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RS / "roadmap_phase4b/scripts"))
sys.path.insert(0, str(RS / "roadmap_phase4/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from sklearn.model_selection import GridSearchCV, GroupKFold  # noqa: E402
from sklearn.pipeline import Pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

import p4b_run  # noqa: E402
from p4b_features import RAW_COLS, RepBuilder  # noqa: E402
from p4_data import dl_tau2, outer_splits  # noqa: E402
from p4_models import make  # noqa: E402

OUT = RS / "lipid_prediction_followup"
RES = OUT / "results"
RES.mkdir(parents=True, exist_ok=True)
TARGETS = ("PE_in_CD", "PE_in_MD")


def split_null(tgt):
    """Per-split null predictions (weighted training mean) and fold membership, identical to p4b_run.run_task."""
    j = p4b_run.data(tgt, False)
    rows = []
    for r, k, tr, te, grp in outer_splits(j.index, "module"):
        y, s = j.y.values, j.s.values
        w = 1.0 / (s**2 + dl_tau2(y[tr], s[tr]))
        rows.append(pd.DataFrame({"gene": j.index[te], "target": tgt, "repeat": r, "fold": k, "y": y[te], "s": s[te],
                                  "null": np.average(y[tr], weights=w[tr]), "module": grp[te]}))
    return pd.concat(rows)


def full_task(tgt, model, r_i, k_i):
    """Exact copy of p4b_run.run_task fitting logic for rep FULL, returning held-out predictions."""
    j = p4b_run.data(tgt, False)
    for r, k, tr, te, grp in outer_splits(j.index, "module"):
        if r == r_i and k == k_i:
            break
    X, y, s = j[RAW_COLS].values, j["y"].values, j["s"].values
    w = 1.0 / (s**2 + dl_tau2(y[tr], s[tr]))
    base, grid = make(p4b_run.MODELS[model], RAW_COLS)
    est = Pipeline([("rep", RepBuilder("FULL")), ("scale", StandardScaler()), ("model", base.named_steps["model"])])
    if grid is None:
        est.fit(X[tr], y[tr], model__sample_weight=w[tr])
    else:
        gs = GridSearchCV(est, grid, cv=GroupKFold(4), scoring="neg_mean_squared_error", n_jobs=1, refit=True)
        gs.fit(X[tr], y[tr], groups=grp[tr], model__sample_weight=w[tr])
        est = gs.best_estimator_
    pred = est.predict(X[te])
    return pd.DataFrame({"gene": j.index[te], "y": y[te], "pred": pred, "target": tgt, "rep": "FULL", "model": model,
                         "repeat": r_i, "fold": k_i})


if __name__ == "__main__":
    N = pd.concat([split_null(t) for t in TARGETS], ignore_index=True)
    old = pd.read_csv(RS / "roadmap_phase4c_diagnostics/tables/fix_experiment_oof_predictions.csv.gz")
    old = old[old.rep.isin(["F8", "F8+C"])]
    jobs = [(t, m, r, k) for t in TARGETS for m in ("ols", "ridge") for r in range(3) for k in range(5)]
    new = pd.concat(Parallel(n_jobs=int(os.environ.get("LP_JOBS", 12)), verbose=5)(delayed(full_task)(*a) for a in jobs),
                    ignore_index=True)
    P = pd.concat([old, new], ignore_index=True)
    P["model_id"] = P.rep + "-" + P.model.map({"ols": "OLS", "ridge": "Ridge"})
    P = P.merge(N.drop(columns="y"), on=["gene", "target", "repeat", "fold"], how="left", validate="many_to_one")
    assert P.null.notna().all()
    # y in the 4C file must equal the target used here
    yy = N.set_index(["gene", "target", "repeat"]).y
    assert np.allclose(P.set_index(["gene", "target", "repeat"]).y.values, yy.loc[P.set_index(["gene", "target", "repeat"]).index].values)
    assert (P.groupby(["model_id", "target", "repeat"]).gene.nunique() == 11758).all()
    P.to_csv(RES / "oof_predictions_with_null.csv.gz", index=False)

    # reproduction check vs stored Phase 4B per-split metrics (stage main, unshuffled)
    M = pd.read_csv(RS / "roadmap_phase4b/results/metrics_main.csv")
    M = M[(M.stage == "main") & (~M.shuffle_b) & M.rep.isin(["F8", "F8+C", "FULL"]) & M.model.isin(["ols", "ridge"])]
    chk = []
    for (rep, mo, t, r, k), g in P.groupby(["rep", "model", "target", "repeat", "fold"]):
        m = M[(M.rep == rep) & (M.model == mo) & (M.target == t) & (M.repeat == r) & (M.fold == k)]
        if len(m) != 1:
            continue
        sse, sse0 = ((g.y - g.pred) ** 2).sum(), ((g.y - g.null) ** 2).sum()
        chk.append(dict(rep=rep, model=mo, target=t, repeat=r, fold=k, rmse_stored=m.rmse.iloc[0],
                        rmse_reproduced=float(np.sqrt(((g.y - g.pred) ** 2).mean())), r2_stored=m.r2_oos.iloc[0],
                        r2_reproduced=float(1 - sse / sse0)))
    C = pd.DataFrame(chk)
    C["abs_diff_rmse"] = (C.rmse_stored - C.rmse_reproduced).abs()
    C["abs_diff_r2"] = (C.r2_stored - C.r2_reproduced).abs()
    C.to_csv(RES / "reproduction_check_vs_phase4b.csv", index=False)
    s = C.groupby(["rep", "model"])[["abs_diff_rmse", "abs_diff_r2"]].max()
    print(s)
    json.dump({"n_splits_checked": len(C), "max_abs_diff_rmse": float(C.abs_diff_rmse.max()),
               "max_abs_diff_r2": float(C.abs_diff_r2.max()), "by_model": s.reset_index().to_dict("records")},
              open(RES / "run_info_lp01.json", "w"), indent=2)
