"""Shortlist robustness (comparison_protocol.md section 3): Phase 4 sensitivity variants + 19 LOO targets for the
shortlisted representation x model pairs, each compared with F8 under the same model and the same 5 splits (repeat 0).
Variant V8 (effects-only) is not applicable (the representation defines the features) and is skipped.
Robustness rule (proportions carried over from the Phase 4 protocol): mean dRMSE < 0 vs F8 in >= 6/7 variants and >= 15/19 LOO."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

P4B = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P4B / "scripts"))
from p4b_features import B_COLS, B_SETS, RAW_COLS, RepBuilder  # noqa: E402

sys.path.insert(0, str(P4B.parent / "roadmap_phase4" / "scripts"))
from p4_data import P3, RE3, VARIANTS, _csv, dl_tau2, feature_table, loo_names, outer_splits, target  # noqa: E402
from p4_models import make  # noqa: E402

MODELS = {"ols": "ols_all", "ridge": "ridge", "elasticnet": "elasticnet"}
SHORT = {"PE_in_CD": ["FULL"], "PE_in_MD": ["FULL", "FULL-A"]}
VARS = [v for v in VARIANTS if v not in ("primary", "V8_effects_only")]
_C = {}


def raw_for(variant):
    if variant in _C:
        return _C[variant]
    X = feature_table(variant)
    gs = pd.read_csv(P3 / "results" / "gene_sets.csv")
    for s, c in zip(B_SETS, B_COLS):
        X[c] = X.index.isin(gs.gene[gs.set == s]).astype(float)
    for col, f in [("FASTclean_eff", "FASTclean_GSE107787_ZT0ZT4_muscleLow"), ("PPARA5_eff", "PPARA_GSE282998_WT_5v0"),
                   ("PPARA10_eff", "PPARA_GSE282998_WT_10v0")]:
        X[col] = _csv(RE3 / f"{f}.csv.gz").set_index("gene")["log2FC"].reindex(X.index)
    X = X[RAW_COLS].dropna()  # V6: re-admitted genes lack block-C inputs and are dropped (documented)
    _C[variant] = X
    return X


def task(kind, variant, target_name, rep, model, fold):
    X0 = raw_for(variant)
    t = target(target_name)
    j = X0.join(t, how="inner").dropna(subset=["y", "s"])
    j = j[j.s > 0]
    grouping = VARIANTS[variant].get("grouping", "module")
    for r, k, tr, te, grp in outer_splits(j.index, grouping):
        if r == 0 and k == fold:
            break
    X, y, s = j[RAW_COLS].values, j.y.values, j.s.values
    w = 1 / (s**2 + dl_tau2(y[tr], s[tr]))
    base, grid = make(MODELS[model], RAW_COLS)
    est = Pipeline([("rep", RepBuilder(rep)), ("scale", StandardScaler()), ("model", base.named_steps["model"])])
    if grid is None:
        est.fit(X[tr], y[tr], model__sample_weight=w[tr])
    else:
        gs = GridSearchCV(est, grid, cv=GroupKFold(4), scoring="neg_mean_squared_error", n_jobs=1)
        gs.fit(X[tr], y[tr], groups=grp[tr], model__sample_weight=w[tr])
        est = gs.best_estimator_
    pred = est.predict(X[te])
    yt = y[te]
    null = np.average(y[tr], weights=w[tr])
    return dict(kind=kind, variant=variant, target=target_name, base_target=target_name.split(":")[-1], rep=rep,
                model=model, fold=fold, n_genes=len(j), rmse=float(np.sqrt(((yt - pred) ** 2).mean())),
                r2_oos=float(1 - ((yt - pred) ** 2).sum() / ((yt - null) ** 2).sum()))


def main():
    T = []
    for tgt, reps in SHORT.items():
        for rep in reps + ["F8"]:
            for mo in MODELS:
                for k in range(5):
                    for v in VARS:
                        T.append(("variant", v, tgt, rep, mo, k))
                    for d in loo_names():
                        T.append(("loo", "primary", f"LOO:{d}:{tgt}", rep, mo, k))
    T = list(dict.fromkeys(T))
    print(len(T), "tasks", flush=True)
    t0 = time.time()
    out = pd.DataFrame(Parallel(n_jobs=19, verbose=5, batch_size=4)(delayed(task)(*t) for t in T))
    out.to_csv(P4B / "results" / "metrics_sensitivity.csv", index=False)
    g = out.groupby(["kind", "variant", "target", "base_target", "rep", "model"]).agg(rmse=("rmse", "mean"), r2=("r2_oos", "mean")).reset_index()
    rows = []
    for (kind, bt, rep, mo), d in g[g.rep != "F8"].groupby(["kind", "base_target", "rep", "model"]):
        f8 = g[(g.kind == kind) & (g.base_target == bt) & (g.rep == "F8") & (g.model == mo)].set_index(["variant", "target"])
        dd = d.set_index(["variant", "target"])
        diff = dd.rmse - f8.rmse.reindex(dd.index)
        dr2 = dd.r2 - f8.r2.reindex(dd.index)
        rows.append(dict(kind=kind, target=bt, rep=rep, model=mo, n=len(diff), n_better=int((diff < 0).sum()),
                         mean_d_rmse=float(diff.mean()), min_d_r2=float(dr2.min()), mean_d_r2=float(dr2.mean()),
                         max_d_r2=float(dr2.max())))
    R = pd.DataFrame(rows)
    R["robust"] = np.where(R.kind == "variant", R.n_better >= 6, R.n_better >= 15)
    R.to_csv(P4B / "results" / "shortlist_robustness.csv", index=False)
    g.to_csv(P4B / "results" / "sensitivity_summary.csv", index=False)
    print(R.round(4).to_string())
    print(f"done in {time.time() - t0:.0f}s")
    json.dump({"variants": VARS, "loo": loo_names()}, open(P4B / "results" / "sensitivity_run_info.json", "w"), indent=2)


if __name__ == "__main__":
    main()
