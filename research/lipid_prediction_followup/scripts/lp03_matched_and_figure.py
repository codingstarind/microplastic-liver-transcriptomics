"""Lipid-gene prediction follow-up, step 3: (i) comparison against OTHER genes matched on null-baseline error,
(ii) overlap of the three lipid sets with the other 9 block-B sets, (iii) concise-table extension, (iv) two-panel figure.
Reads lp01/lp02 outputs only. No model is fitted, tuned or selected.

Null-error matching (descriptive stratification; uses the observed target through the null squared error):
for each gene in a group, the 10 nearest OTHER genes by repeat-averaged null squared error MSE_null (identical across
models), greedy without replacement, in descending MSE_null order. Same matched set for every model within a diet.
Uncertainty: module block bootstrap (B = 2000, 20 modules) over group + matched genes; module fixed-effects regression
of per-gene dMSE with CR1 module-cluster-robust SE (t, G-1 df).
Run: .venv/Scripts/python.exe research/lipid_prediction_followup/scripts/lp03_matched_and_figure.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
RS = Path(__file__).resolve().parents[2]
OUT = RS / "lipid_prediction_followup"
RES, FIG = OUT / "results", OUT / "figures"
FIG.mkdir(exist_ok=True)
SEED, B, K = 20261125, 2000, 10
PRIMARY = ["F8-OLS", "F8+C-Ridge", "FULL-OLS"]
ALLM = PRIMARY + ["F8+C-OLS", "FULL-Ridge"]
SETS = {"LIPO": "LIPOGENESIS_CORE", "CHOL": "CHOL_SYNTHESIS_CORE", "PPARA": "SANDERSON_PPARA_TARGETS"}
B9 = ["HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_INFLAMMATORY_RESPONSE",
      "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_XENOBIOTIC_METABOLISM", "HALLMARK_BILE_ACID_METABOLISM",
      "HALLMARK_PEROXISOME", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_MTORC1_SIGNALING"]


def cr1_fe(yv, xv, module):
    from scipy import stats
    D = pd.get_dummies(module).values.astype(float)
    X = np.column_stack([xv, D])
    b, *_ = np.linalg.lstsq(X, yv, rcond=None)
    e = yv - X @ b
    XtXi = np.linalg.pinv(X.T @ X)
    gs_, n, k = np.unique(module), len(yv), X.shape[1]
    meat = sum(np.outer(X[module == g].T @ e[module == g], X[module == g].T @ e[module == g]) for g in gs_)
    V = XtXi @ meat @ XtXi * (len(gs_) / (len(gs_) - 1)) * ((n - 1) / (n - k))
    se = np.sqrt(V[0, 0])
    qt = stats.t.ppf(.975, len(gs_) - 1)
    return b[0], b[0] - qt * se, b[0] + qt * se, float(2 * stats.t.sf(abs(b[0] / se), len(gs_) - 1))


P = pd.read_csv(RES / "oof_predictions_with_null.csv.gz")
P["se"], P["se0"] = (P.y - P.pred) ** 2, (P.y - P.null) ** 2
G = P.groupby(["model_id", "target", "gene"]).agg(module=("module", "first"), mse=("se", "mean"), mse0=("se0", "mean"),
                                                   y=("y", "first")).reset_index()
G["d"] = G.mse0 - G.mse
gs = pd.read_csv(RS / "roadmap_phase3/results/gene_sets.csv")
genes = np.array(sorted(G.gene.unique()))
mem = pd.DataFrame({k: np.isin(genes, gs.gene[gs.set == v]) for k, v in SETS.items()}, index=genes)
lipid_any = mem.any(axis=1).values
GROUPS = {"LIPO": mem.LIPO.values, "CHOL": mem.CHOL.values, "PPARA": mem.PPARA.values, "LIPID_any": lipid_any}
other = ~lipid_any

# ---- (ii) overlap of lipid sets with other block-B sets (prediction universe)
ov = pd.DataFrame({s9: [int((GROUPS[g] & np.isin(genes, gs.gene[gs.set == s9])).sum()) for g in GROUPS] for s9 in B9},
                  index=list(GROUPS))
ov.insert(0, "n_genes", [int(v.sum()) for v in GROUPS.values()])
ov["n_in_any_of_9_other_B_sets"] = [int((GROUPS[g] & np.isin(genes, gs.gene[gs.set.isin(B9)])).sum()) for g in GROUPS]
ov.to_csv(RES / "lipid_sets_overlap_with_other_blockB_sets.csv")

# ---- (i) null-error-matched OTHER genes
mods = np.array(sorted(G.module.unique()))
rng = np.random.default_rng(SEED)
BOOTW = rng.multinomial(len(mods), np.ones(len(mods)) / len(mods), size=B)
rows, match_rows = [], []
for tgt in ("PE_in_CD", "PE_in_MD"):
    base = G[(G.model_id == "F8-OLS") & (G.target == tgt)].set_index("gene").loc[genes]
    m0 = base.mse0.values
    for g, gm in GROUPS.items():
        avail = other.copy()
        ref = np.zeros(len(genes), bool)
        for i in np.argsort(-m0 * gm)[: gm.sum()]:
            cand = np.where(avail)[0]
            pick = cand[np.argsort(np.abs(m0[cand] - m0[i]))[:K]]
            ref[pick], avail[pick] = True, False
        match_rows.append(dict(diet=tgt, group=g, n_group=int(gm.sum()), n_matched=int(ref.sum()),
                               median_mse0_group=float(np.median(m0[gm])), median_mse0_matched=float(np.median(m0[ref])),
                               mean_abs_y_group=float(base.y.abs().values[gm].mean()),
                               mean_abs_y_matched=float(base.y.abs().values[ref].mean()),
                               n_modules_matched=int(pd.Series(base.module.values[ref]).nunique())))
        for mid in ALLM:
            df = G[(G.model_id == mid) & (G.target == tgt)].set_index("gene").loc[genes]

            def sums(mask):
                return df[mask].groupby("module").agg(n=("d", "size"), sse=("mse", "sum"), sse0=("mse0", "sum"),
                                                      sd=("d", "sum")).reindex(mods, fill_value=0).values
            S1, S0 = sums(gm), sums(ref)
            t1, t0 = S1.sum(0), S0.sum(0)
            b1, b0 = BOOTW @ S1, BOOTW @ S0
            with np.errstate(invalid="ignore", divide="ignore"):
                br2 = (1 - b1[:, 1] / b1[:, 2]) - (1 - b0[:, 1] / b0[:, 2])
                bdm = b1[:, 3] / b1[:, 0] - b0[:, 3] / b0[:, 0]
            keep = gm | ref
            fe = cr1_fe(df.d.values[keep], gm[keep].astype(float), df.module.values[keep])
            rows.append(dict(diet=tgt, model=mid, primary=mid in PRIMARY, group=g, reference="OTHER_nullerror_matched_k10",
                             R2_group=1 - t1[1] / t1[2], R2_matched=1 - t0[1] / t0[2], diff_R2=(1 - t1[1] / t1[2]) - (1 - t0[1] / t0[2]),
                             diff_R2_ci_low=np.nanquantile(br2, .025), diff_R2_ci_high=np.nanquantile(br2, .975),
                             dMSE_group=t1[3] / t1[0], dMSE_matched=t0[3] / t0[0], diff_dMSE=t1[3] / t1[0] - t0[3] / t0[0],
                             diff_dMSE_ci_low=np.nanquantile(bdm, .025), diff_dMSE_ci_high=np.nanquantile(bdm, .975),
                             moduleFE_diff_dMSE=fe[0], moduleFE_ci_low=fe[1], moduleFE_ci_high=fe[2], moduleFE_p_CR1=fe[3]))
MT = pd.DataFrame(rows)
MT.to_csv(RES / "group_vs_nullerror_matched_other.csv", index=False)
pd.DataFrame(match_rows).to_csv(RES / "nullerror_matching_balance.csv", index=False)

# ---- (iii) concise table: append matched-reference columns (written by lp02; extended here into a new file)
T = pd.read_csv(OUT / "TABLE_results_concise.csv")
T = T[[c for c in T.columns if not c.startswith("matched_")]]  # idempotent re-runs
mt = MT[MT.primary][["diet", "model", "group", "R2_matched", "diff_R2", "diff_R2_ci_low", "diff_R2_ci_high", "diff_dMSE",
                     "diff_dMSE_ci_low", "diff_dMSE_ci_high", "moduleFE_diff_dMSE", "moduleFE_p_CR1"]]
mt.columns = ["diet", "model", "group"] + [f"matched_{c}" for c in mt.columns[3:]]
T2 = T.merge(mt, on=["diet", "model", "group"], how="left")
T2.to_csv(OUT / "TABLE_results_concise.csv", index=False)

# ---- (iv) figure
R = pd.read_csv(RES / "group_metrics_all.csv")
main = ["LIPO", "CHOL", "PPARA", "LIPID_any", "OTHER"]
COL = {"F8-OLS": "#2a6ebb", "F8+C-Ridge": "#d0782a", "FULL-OLS": "#6b6b6b"}
LAB = {"F8-OLS": "F8 (OLS)", "F8+C-Ridge": "F8+C (ridge)",
       "FULL-OLS": "FULL, 28 features (OLS) — contains membership indicators for these three sets"}
GL = {"LIPO": "Lipogenesis\ncore", "CHOL": "Cholesterol\nsynthesis core", "PPARA": "PPARα\ntargets",
      "LIPID_any": "Any of the\nthree sets", "OTHER": "Other genes"}
YLO, YHI = -0.6, 1.05
fig, axes = plt.subplots(1, 2, figsize=(12, 5.0), sharey=True)
x = np.arange(len(main))
for ax, tgt, ttl in ((axes[0], "PE_in_CD", "PE vs Veh within CD"), (axes[1], "PE_in_MD", "PE vs Veh within MD")):
    ax.axhline(0, color="#4a4a4a", lw=.8)
    for j, mid in enumerate(PRIMARY):
        d = R[(R.diet == tgt) & (R.model == mid)].set_index("group").loc[main]
        xx = x + (j - 1) * .24
        lo = d.R2_ci_low.clip(lower=YLO)
        ax.vlines(xx, lo, d.R2_ci_high.clip(upper=YHI), color=COL[mid], lw=1.6)
        for xi, (l_, raw) in zip(xx, zip(lo, d.R2_ci_low)):
            if raw < YLO:
                ax.plot(xi, YLO, marker="v", color=COL[mid], ms=5, clip_on=False)
                ax.text(xi, YLO - .045, f"{raw:.1f}", ha="center", va="top", fontsize=6, color="#4a4a4a")
        ax.scatter(xx, d.R2_group, s=34, zorder=3, color="white" if mid == "FULL-OLS" else COL[mid], edgecolor=COL[mid],
                   lw=1.4, label=LAB[mid])
    nn = R[(R.diet == tgt) & (R.model == PRIMARY[0])].set_index("group").loc[main]
    ax.set_xticks(x)
    ax.set_xticklabels([f"{GL[g]}\n{int(nn.loc[g, 'n_genes'])} genes\n{int(nn.loc[g, 'n_modules'])} modules" for g in main], fontsize=7)
    ax.axvline(3.5, color="#bdbdbd", lw=.6, ls=(0, (2, 2)))
    ax.set_ylim(YLO - .12, YHI)
    ax.set_title(ttl, fontsize=9.5, loc="left")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.tick_params(axis="y", labelsize=7.5)
axes[0].set_ylabel("Held-out R² within gene group\n1 − SSE(model) / SSE(null: training mean)", fontsize=8)
h, l_ = axes[0].get_legend_handles_labels()
fig.legend(h, l_, loc="upper left", bbox_to_anchor=(.01, .93), ncol=3, fontsize=7.2, frameon=False)
fig.suptitle("Held-out gene-level prediction of PE log2FC, locked Phase 4 co-expression-module splits (3 repeats × 5 folds). "
             "Bars: 95% module-block bootstrap CI (20 modules, B = 2000); ▼ = interval truncated at −0.6 (value shown).",
             fontsize=7.8, x=.01, ha="left")
fig.tight_layout(rect=(0, 0, 1, .89))
for e in ("png", "pdf", "svg"):
    fig.savefig(FIG / f"fig_lipid_vs_other_prediction.{e}", dpi=300)
print(ov.to_string())
print(pd.DataFrame(match_rows).round(4).to_string(index=False))
print(MT[MT.primary].round(4).to_string(index=False))
