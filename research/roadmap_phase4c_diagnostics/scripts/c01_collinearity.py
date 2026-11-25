"""Task 1: multicollinearity audit of F8, A4, B12, C4, F8+C, F8+B+C, FULL (Phase 4B representations)."""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from c00_common import MATS, RAW_COLS, RS, F, RepBuilder, T, outer_splits, raw_table  # noqa: E402

X = raw_table()
rb = RepBuilder("FULL").fit(X.values)
allc = rb._all(X.values)
D = pd.DataFrame(allc, index=X.index)

rows, vif_rows, pairs = [], [], []
for name, cols in MATS.items():
    M = D[cols]
    R = np.corrcoef(M.values, rowvar=False)
    iu = np.triu_indices(len(cols), 1)
    vif = np.diag(np.linalg.pinv(R))
    sv = np.linalg.svd((M - M.mean()) / M.std(), compute_uv=False)
    rank = int((sv > 1e-10 * sv[0]).sum())
    conds = []
    for r, k, tr, te, grp in outer_splits(X.index):  # condition number after training-fold scaling
        Mt = M.values[tr]
        Z = (Mt - Mt.mean(0)) / Mt.std(0)
        s = np.linalg.svd(Z, compute_uv=False)
        conds.append(s[0] / s[-1])
    rows.append(dict(matrix=name, n_features=len(cols), max_abs_r=float(np.abs(R[iu]).max()), max_vif=float(vif.max()),
                     n_vif_gt_10=int((vif > 10).sum()), n_vif_gt_5=int((vif > 5).sum()), rank=rank,
                     smallest_singular=float(sv[-1]), cond_number_median=float(np.median(conds)),
                     cond_number_max=float(np.max(conds))))
    for c, v in zip(cols, vif):
        vif_rows.append(dict(matrix=name, feature=c, vif=float(v)))
    for i, j in zip(*iu):
        if abs(R[i, j]) >= 0.5:
            pairs.append(dict(matrix=name, a=cols[i], b=cols[j], r=float(R[i, j])))
summ = pd.DataFrame(rows)
summ.to_csv(T / "collinearity_summary.csv", index=False)
pd.DataFrame(vif_rows).to_csv(T / "vif_by_feature.csv", index=False)
P = pd.DataFrame(pairs).sort_values("r", key=abs, ascending=False)
P.to_csv(T / "correlated_pairs_absr_ge_0.5.csv", index=False)

# biological-axis similarity (effects only) vs multicollinearity
bio = D[["MASH_eff", "FAST_eff", "PPARA_eff", "CIRC_eff"]].corr(method="spearman").round(3)
bio.to_csv(T / "biological_axis_similarity_spearman.csv")

# coefficient stability from Phase 4B (OLS + ridge, 15 folds)
co = pd.read_csv(RS / "roadmap_phase4b" / "results" / "coefficients_main.csv")
co = co[co.rep.isin(["F8", "F8+C", "FULL", "FULL-A"]) & co.model.isin(["ols", "ridge"])]
st = co.groupby(["target", "rep", "model", "feature"]).coef_per_SD.agg(["mean", "std", "min", "max"]).reset_index()
st["cv"] = st["std"] / st["mean"].abs()
st["sign_flip"] = np.sign(st["min"]) != np.sign(st["max"])
st.to_csv(T / "coefficient_stability_phase4b.csv", index=False)
agg = st.groupby(["target", "rep", "model"]).agg(median_cv=("cv", "median"), max_abs_coef=("mean", lambda s: s.abs().max()),
                                                  n_sign_flips=("sign_flip", "sum"), n_coef=("feature", "size")).reset_index()
agg.to_csv(T / "coefficient_stability_summary.csv", index=False)

fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
v = pd.DataFrame(vif_rows)
full = v[v.matrix == "FULL"].sort_values("vif", ascending=False)
ax[0].barh(full.feature, full.vif, color=["tab:red" if x > 10 else "tab:grey" for x in full.vif])
ax[0].axvline(10, ls="--", c="k", lw=.7)
ax[0].set_xscale("log")
ax[0].tick_params(axis="y", labelsize=6)
ax[0].set_title("FULL (Phase 4B): VIF per feature (red > 10)", fontsize=9)
ax[1].bar(summ.matrix, summ.cond_number_median)
ax[1].set_yscale("log")
ax[1].set_ylabel("condition number (median over 15 training folds)")
ax[1].set_title("Condition number by representation", fontsize=9)
fig.tight_layout()
fig.savefig(F / "fig1_collinearity.png", dpi=110)
print(summ.round(3).to_string())
print(P.head(12).round(3).to_string())
print(agg.round(3).to_string())
json.dump({"raw_cols": RAW_COLS}, open(T / "c01_run_info.json", "w"))
