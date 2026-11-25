"""Step 7: figures (static PNG for the report). Reads only Phase 5 outputs."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from p5_common import FIG, P5, RES  # noqa: E402

FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#8a8984",
                     "axes.labelcolor": "#52514e", "xtick.color": "#52514e", "ytick.color": "#52514e", "figure.facecolor": "#fcfcfb",
                     "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb"})
GREY, INK = "#8a8984", "#0b0b0b"
COL = {"CON": GREY, "PE": "#2a78d6", "PS": "#eb6834", "PP": "#1baf7a", "PB": "#eda100",
       "PS_80nm": "#eb6834", "PS_0.5um": "#e87ba4", "PS_5um": "#6250d6"}
SC = pd.read_csv(P5 / "per_sample_scores.csv", index_col=0)
SL = pd.read_csv(P5 / "sample_level_external_results.csv")

# ---- fig1: per-mouse scores, one row per study, one column per frozen score ----
scores = [("FAST_singscore", "FAST singscore (4D)"), ("MASH_singscore", "MASH singscore (4D)"),
          ("R0_logistic_logodds", "R0 L2-logistic log-odds (4E)"), ("R2_logistic_logodds", "R2 L2-logistic log-odds (4E)")]
orders = {"HONG": ["CON", "PE", "PS", "PP", "PB"], "GSE245069": ["CON", "PS_80nm", "PS_0.5um", "PS_5um"], "GSE277772": ["CON", "PS"]}
fig, ax = plt.subplots(3, 4, figsize=(13, 8.5))
for i, (st, order) in enumerate(orders.items()):
    D = SC[SC.study == st]
    for j, (sc, lab) in enumerate(scores):
        a = ax[i, j]
        for k, g in enumerate(order):
            v = D[D.group == g][sc].values
            a.scatter(np.full(len(v), k) + np.linspace(-.12, .12, len(v)), v, s=36, color=COL[g], edgecolor="#fcfcfb", linewidth=1, zorder=3)
            a.hlines(v.mean(), k - .28, k + .28, color=INK, lw=1.5, zorder=4)
        con = D[D.group == "CON"][sc].mean()
        a.axhline(con, color=GREY, lw=.8, ls=":", zorder=1)
        a.set_xticks(range(len(order)))
        a.set_xticklabels([o.replace("PS_", "") for o in order], fontsize=8)
        a.grid(axis="y", color="#e6e5e0", lw=.6)
        if i == 0:
            a.set_title(lab, fontsize=9, color=INK)
        if j == 0:
            a.set_ylabel(st, fontsize=10, color=INK)
fig.suptitle("Frozen per-mouse scores in independent studies (bar = group mean; dotted = shared control mean)", color=INK)
fig.tight_layout()
fig.savefig(FIG / "fig1_per_mouse_frozen_scores.png", dpi=130)
plt.close(fig)

# ---- fig2: gene level, primary contrast ----
M = pd.read_csv(RES / "gene_level_primary_merged.csv.gz", index_col=0)
G = pd.read_csv(P5 / "gene_level_external_results.csv")
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
a = ax[0]
a.hexbin(M["F8_ols__PE_in_CD"], M["log2FC"].clip(-4, 4), gridsize=60, cmap=LinearSegmentedColormap.from_list("b", ["#eef4fc", "#2a78d6", "#0d3a73"]),
         bins="log", mincnt=1)
a.axhline(0, color=GREY, lw=.6)
a.axvline(0, color=GREY, lw=.6)
a.set_xlabel("F8-OLS predicted log2FC (frozen, trained on GSE294180 PE_in_CD)")
a.set_ylabel("Hong PE vs CON observed log2FC (DESeq2 MLE, clipped ±4)")
r = G[(G.study == "HONG") & (G.contrast == "PE_vs_CON") & (G.model == "F8_ols") & (G.training_target == "PE_in_CD")].iloc[0]
a.set_title(f"Spearman {r.spearman:.2f} [{r.spearman_ci_low:.2f}, {r.spearman_ci_high:.2f}] (gene-block CI)\n"
            f"weighted R² vs own mean {r.R2w_vs_own_mean:.3f}; mouse-label permutation p = {r.perm_p_one_sided:.2f}", fontsize=8.5)
a = ax[1]
sub = G[(G.study == "HONG") & (G.contrast.isin(["PE_vs_CON", "PS_vs_CON", "PP_vs_CON", "PB_vs_CON"])) & (G.training_target == "PE_in_CD")
        & G.model.isin(["F8_ols", "F8C_ridge", "FAST_only_ols", "MASH_only_ols"])]
mods = ["F8_ols", "F8C_ridge", "FAST_only_ols", "MASH_only_ols"]
for k, c in enumerate(["PE_vs_CON", "PS_vs_CON", "PP_vs_CON", "PB_vs_CON"]):
    for jj, m in enumerate(mods):
        rr = sub[(sub.contrast == c) & (sub.model == m)].iloc[0]
        x = k + (jj - 1.5) * .2
        a.vlines(x, rr.spearman_ci_low, rr.spearman_ci_high, color=["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"][jj], lw=2)
        a.scatter(x, rr.spearman, color=["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"][jj], s=30, zorder=3, label=m if k == 0 else None)
        a.scatter(x, rr.perm_null_q95, marker="_", color=INK, s=80, zorder=4, label="mouse-label permutation 95th pct\n(log2CPM proxy, all relabellings)" if (k == 0 and jj == 0) else None)
a.axhline(0, color=GREY, lw=.6)
a.set_xticks(range(4))
a.set_xticklabels(["PE", "PS", "PP", "PB"])
a.set_ylabel("Spearman(predicted, observed)")
a.set_title("Hong lab, each polymer vs the SAME 3 controls", fontsize=9)
a.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
fig.tight_layout()
fig.savefig(FIG / "fig2_gene_level_external.png", dpi=130)
plt.close(fig)

# ---- fig3: R2 decomposition ----
DC = pd.read_csv(RES / "r2_feature_contributions_external.csv")
DC = DC[DC.feature != "(TOTAL)"]
p = DC[(DC.study == "HONG") & (DC.contrast == "PE_vs_CON")].set_index("feature")
q = DC[(DC.study == "GSE245069") & (DC.contrast == "PSpooled_vs_CON")].set_index("feature")
feat = p.index.tolist()
order = np.argsort(-np.abs(p.gse294180_contrib_PE_in_CD.values))
feat = [feat[i] for i in order]
cols = [("GSE294180 PE_in_CD", p.gse294180_contrib_PE_in_CD, "#8a8984"), ("GSE294180 PE_in_MD", p.gse294180_contrib_PE_in_MD, "#c3c2b7"),
        ("Hong PE vs CON", p.contribution_logodds, "#2a78d6"), ("GSE245069 PS pooled vs CON", q.contribution_logodds, "#eb6834")]
fig, a = plt.subplots(figsize=(10, 5.2))
y = np.arange(len(feat))
for k, (lab, s, c) in enumerate(cols):
    a.barh(y + (k - 1.5) * .2, s.reindex(feat).values, height=.19, color=c, label=lab)
a.set_yticks(y)
a.set_yticklabels([f.replace("R1_", "").replace("HALLMARK_", "H_") for f in feat], fontsize=8)
a.invert_yaxis()
a.axvline(0, color=INK, lw=.6)
a.set_xlabel("contribution to the R2 log-odds shift (coef/SD × group-mean feature difference)")
a.legend(fontsize=8, frameon=False)
a.set_title("What moves the frozen R2 score: lipid/cholesterol/PPARα features dominate;\nthe MASH-program feature (R0_signature) is small", fontsize=9)
fig.tight_layout()
fig.savefig(FIG / "fig3_r2_decomposition.png", dpi=130)
plt.close(fig)

# ---- fig4: pathway concordance heatmap (primary panel) ----
PC = pd.read_csv(P5 / "pathway_concordance.csv")
PC = PC[PC.panel == "primary12"]
cols = ["HONG|PE_vs_CON", "HONG|PS_vs_CON", "HONG|PP_vs_CON", "HONG|PB_vs_CON", "GSE245069|PSpooled_vs_CON", "GSE245069|PS80nm_vs_CON",
        "GSE245069|PS05um_vs_CON", "GSE245069|PS5um_vs_CON", "GSE277772|PS_vs_CON"]
H = PC.assign(k=PC.study + "|" + PC.contrast).pivot(index="set", columns="k", values="median_log2FC_external")[cols]
ref = PC.drop_duplicates("set").set_index("set")[["median_log2FC_GSE294180_PE_in_CD", "median_log2FC_GSE294180_PE_in_MD", "median_log2FC_REF_FAST",
                                                   "median_log2FC_REF_MASH1"]]
H = pd.concat([ref, H], axis=1)
H.columns = ["GSE294180 PE|CD", "GSE294180 PE|MD", "REF fasting", "REF MASH1", "Hong PE", "Hong PS", "Hong PP", "Hong PB", "Li PS pooled",
             "Li 80nm", "Li 0.5um", "Li 5um", "GSE277772 PS*"]
cm = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f2f1ee", "#eb6834"])
fig, a = plt.subplots(figsize=(11, 5))
v = 1.0
im = a.imshow(H.values.clip(-v, v), cmap=cm, vmin=-v, vmax=v, aspect="auto")
a.set_xticks(range(H.shape[1]))
a.set_xticklabels(H.columns, rotation=40, ha="right", fontsize=8)
a.set_yticks(range(H.shape[0]))
a.set_yticklabels([s.replace("HALLMARK_", "H_") for s in H.index], fontsize=8)
for i in range(H.shape[0]):
    for j in range(H.shape[1]):
        a.text(j, i, f"{H.values[i, j]:.2f}", ha="center", va="center", fontsize=6.5, color=INK)
a.axvline(3.5, color=INK, lw=1)
plt.colorbar(im, ax=a, fraction=.025, label="median gene log2FC (clipped ±1)")
a.set_title("Fixed 12-set panel: median gene log2FC. Left of line = GSE294180 PE and reference axes; right = independent studies (*normalised input)", fontsize=8.5)
fig.tight_layout()
fig.savefig(FIG / "fig4_pathway_concordance.png", dpi=130)
plt.close(fig)
print("figures written")
