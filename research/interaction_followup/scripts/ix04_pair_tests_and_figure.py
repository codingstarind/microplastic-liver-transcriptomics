"""Interaction follow-up, step 4: held-out tests for the pair-only model, surface stability, concise table, figure.
Reads ix02/ix03 outputs. Run: .venv/Scripts/python.exe research/interaction_followup/scripts/ix04_pair_tests_and_figure.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

RS = Path(__file__).resolve().parents[2]
OUT = RS / "interaction_followup"
RES, FIG = OUT / "results", OUT / "figures"
FIG.mkdir(exist_ok=True)
SEED, B = 20261125, 2000

P = pd.read_csv(RES / "pair_only_heldout_predictions.csv.gz")
for c in ("full", "pair", "additive"):
    P[f"se_{c}"] = (P.y - P[f"pred_{c}"]) ** 2
P["se0"] = (P.y - P.null) ** 2
G = P.groupby(["target", "gene"]).agg(module=("module", "first"), **{f"mse_{c}": (f"se_{c}", "mean") for c in ("full", "pair", "additive")},
                                       mse0=("se0", "mean")).reset_index()
gs = pd.read_csv(RS / "roadmap_phase3/results/gene_sets.csv")
genes = np.array(sorted(G.gene.unique()))
mem = {k: np.isin(genes, gs.gene[gs.set == v]) for k, v in
       {"LIPO": "LIPOGENESIS_CORE", "PPARA": "SANDERSON_PPARA_TARGETS", "CHOL": "CHOL_SYNTHESIS_CORE"}.items()}
GROUPS = {"ALL": np.ones(len(genes), bool), "LIPO": mem["LIPO"], "PPARA": mem["PPARA"], "CHOL": mem["CHOL"],
          "OTHER": ~(mem["LIPO"] | mem["PPARA"] | mem["CHOL"])}
mods = np.array(sorted(G.module.unique()))
BOOTW = np.random.default_rng(SEED).multinomial(len(mods), np.ones(len(mods)) / len(mods), size=B)


def cluster_t(v, m):
    u = np.unique(m)
    Gc, n, mean = len(u), len(v), v.mean()
    se = np.sqrt((Gc / (Gc - 1)) * sum(((v[m == g] - mean).sum()) ** 2 for g in u)) / n
    q = stats.t.ppf(.975, Gc - 1)
    return mean, mean - q * se, mean + q * se, float(2 * stats.t.sf(abs(mean / se), Gc - 1)), Gc


def holm(p):
    p = np.asarray(p, float)
    o, adj, run = np.argsort(p), np.empty(len(p)), 0
    for r_, j in enumerate(o):
        run = max(run, (len(p) - r_) * p[j])
        adj[j] = min(1, run)
    return adj


def bh(p):
    p = np.asarray(p, float)
    n, o = len(p), np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


rows = []
for t in ("PE_in_CD", "PE_in_MD"):
    d = G[G.target == t].set_index("gene").loc[genes]
    for g, mk in GROUPS.items():
        x = d[mk]
        gain = (x.mse_additive - x.mse_pair).values
        mean, lo, hi, p, Gc = cluster_t(gain, x.module.values)
        Sm = x.groupby("module")[["mse_pair", "mse_additive", "mse0"]].sum().reindex(mods, fill_value=0).values
        bs = BOOTW @ Sm
        with np.errstate(invalid="ignore", divide="ignore"):
            bd = (bs[:, 1] - bs[:, 0]) / bs[:, 2]
        r2 = {c: 1 - x[f"mse_{c}"].sum() / x.mse0.sum() for c in ("full", "pair", "additive")}
        rows.append(dict(diet=t, group=g, n_genes=int(mk.sum()), n_modules=Gc, R2_additive=r2["additive"], R2_pair_only=r2["pair"],
                         R2_full=r2["full"], dR2_pair_minus_additive=r2["pair"] - r2["additive"], dR2_ci_low=np.nanquantile(bd, .025),
                         dR2_ci_high=np.nanquantile(bd, .975), share_of_full_gain=(r2["pair"] - r2["additive"]) / (r2["full"] - r2["additive"]),
                         mean_gain_mse=mean, gain_ci_low=lo, gain_ci_high=hi, p_cluster_t=p, n_genes_pair_better=int((gain > 0).sum())))
PT = pd.DataFrame(rows)
pm = PT.group.isin(["ALL", "LIPO", "PPARA"])
PT.loc[pm, "p_holm_6"] = holm(PT.loc[pm, "p_cluster_t"])
PT.loc[~pm, "p_BH_secondary"] = bh(PT.loc[~pm, "p_cluster_t"])
PT.to_csv(RES / "pair_only_vs_additive_heldout.csv", index=False)

# surface stability
SU = pd.read_csv(RES / "pair_interaction_surface_by_split.csv.gz")
SU["split"] = SU["repeat"] * 5 + SU["fold"]
stab = []
med = {}
for t, g in SU.groupby("target"):
    M = g.pivot_table(index=["i", "jx"], columns="split", values="pure_interaction")
    C = np.corrcoef(M.values.T)
    iu = np.triu_indices_from(C, 1)
    med[t] = M.median(1)
    stab.append(dict(target=t, median_pairwise_r_between_splits=float(np.median(C[iu])), min_r=float(C[iu].min()),
                     surface_sd_median=float(M.std(0).median()), surface_range_median_over_splits=float((M.max(0) - M.min(0)).median())))
stab.append(dict(target="CD_vs_MD", median_pairwise_r_between_splits=float(np.corrcoef(med["PE_in_CD"], med["PE_in_MD"])[0, 1])))
ST = pd.DataFrame(stab)
ST.to_csv(RES / "pair_surface_stability.csv", index=False)

# ---------------- concise table
W = pd.read_csv(RES / "pair_summary_CD_MD.csv")
HO = pd.read_csv(RES / "heldout_full_vs_additive_by_group.csv")
top = W.sort_values("mean_excess", ascending=False).groupby("model").head(3)[
    ["model", "pair", "median_inter_sd_CD", "median_null_sd_CD", "n_splits_exceed_CD", "median_inter_sd_MD", "median_null_sd_MD",
     "n_splits_exceed_MD", "reproducible"]]
top.insert(0, "section", "A_training_fold_interaction_vs_null (top 3 pairs per model)")
ho = HO[["diet", "model", "group", "n_genes", "R2_full", "R2_additive", "dR2_full_minus_additive", "dR2_ci_low", "dR2_ci_high",
         "p_cluster_t", "p_holm_primary12", "p_BH_secondary"]].copy()
ho.insert(0, "section", "B_heldout_full_vs_additive")
pt = PT[["diet", "group", "n_genes", "R2_pair_only", "R2_additive", "dR2_pair_minus_additive", "dR2_ci_low", "dR2_ci_high",
         "p_cluster_t", "p_holm_6", "p_BH_secondary"]].copy()
pt.insert(0, "section", "C_heldout_HGB_pair_only(FAST_eff x CIRC_se)_vs_additive")
pd.concat([top, ho, pt], ignore_index=True).to_csv(OUT / "TABLE_results_concise.csv", index=False)

# ---------------- figure (2 x 2)
S = pd.read_csv(RES / "pair_by_split_observed_vs_null.csv")
fig = plt.figure(figsize=(12, 9.2))
gsp = fig.add_gridspec(2, 2, hspace=.42, wspace=.28)
info = json.load(open(RES / "run_info_ix03.json"))
GA, GB = np.array(info["grid_FAST_eff"]), np.array(info["grid_CIRC_se"])
vmax = max(np.abs(med["PE_in_CD"]).max(), np.abs(med["PE_in_MD"]).max())
for c, (t, ttl) in enumerate((("PE_in_CD", "A  CD"), ("PE_in_MD", "B  MD"))):
    ax = fig.add_subplot(gsp[0, c])
    Z = med[t].values.reshape(len(GA), len(GB))
    im = ax.imshow(Z.T, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   extent=[0, len(GA), 0, len(GB)])
    ax.set_xticks(np.arange(len(GA))[::2] + .5)
    ax.set_xticklabels([f"{v:.2f}" for v in GA[::2]], fontsize=6.5, rotation=45)
    ax.set_yticks(np.arange(len(GB))[::2] + .5)
    ax.set_yticklabels([f"{v:.3f}" for v in GB[::2]], fontsize=6.5)
    ax.set_xlabel("FAST_eff (fasting reference log2FC; grid = 2.5–97.5% quantiles)", fontsize=7.2)
    ax.set_ylabel("CIRC_se (circadian reference SE)", fontsize=7.2)
    r = ST.set_index("target").loc[t]
    ax.set_title(f"{ttl}: HGB pure interaction FAST_eff × CIRC_se, median of 15 splits\n"
                 f"(between-split surface r median {r.median_pairwise_r_between_splits:.2f})", fontsize=8, loc="left")
    cb = fig.colorbar(im, ax=ax, fraction=.04)
    cb.ax.tick_params(labelsize=6.5)
    cb.set_label("PD₂ − PD(FAST_eff) − PD(CIRC_se), log2FC", fontsize=6.8)
# C: observed vs null interaction SD for top pairs
ax = fig.add_subplot(gsp[1, 0])
pairs = W[W.model == "hgb"].sort_values("mean_excess", ascending=False).pair.head(5).tolist()
for mi, (m, mk) in enumerate((("hgb", "o"), ("catboost", "s"))):
    for di, (t, col) in enumerate((("PE_in_CD", "#2a6ebb"), ("PE_in_MD", "#d0782a"))):
        for pi_, pr in enumerate(pairs):
            g = S[(S.model == m) & (S.target == t) & (S.pair == pr)]
            xx = pi_ + (mi * 2 + di - 1.5) * .17
            ax.scatter(np.full(len(g), xx), g.inter_sd, s=7, color=col, marker=mk, alpha=.55, lw=0)
            ax.plot([xx - .07, xx + .07], [g.null_median.median()] * 2, color="#222", lw=1.4)
ax.set_xticks(range(len(pairs)))
ax.set_xticklabels([p.replace(" x ", " ×\n") for p in pairs], fontsize=6.8)
ax.set_ylabel("Pure-interaction SD (log2FC), training genes", fontsize=7.5)
ax.set_title("C  Per-split interaction strength (dots) vs interaction-free null (black bar = median)\n"
             "order per pair: HGB-CD, HGB-MD, CatBoost-CD, CatBoost-MD", fontsize=8, loc="left")
for lab, col, mk in (("CD", "#2a6ebb", "o"), ("MD", "#d0782a", "o"), ("HGB", "#777", "o"), ("CatBoost", "#777", "s")):
    ax.scatter([], [], color=col, marker=mk, s=14, label=lab)
ax.legend(fontsize=6.5, frameon=False, ncol=4, loc="upper right")
ax.tick_params(axis="y", labelsize=7)
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)
# D: held-out delta R2 vs additive
ax = fig.add_subplot(gsp[1, 1])
grp = ["ALL", "LIPO", "PPARA"]
lab = {"ALL": "All genes\n(11,758)", "LIPO": "Lipogenesis\n(12)", "PPARA": "PPARα targets\n(14)"}
series = [("HGB full", HO[HO.model == "hgb"], "dR2_full_minus_additive", "p_holm_primary12", "#9c4f12", "o"),
          ("CatBoost full", HO[HO.model == "catboost"], "dR2_full_minus_additive", "p_holm_primary12", "#d0782a", "s"),
          ("HGB pair-only", PT, "dR2_pair_minus_additive", "p_holm_6", "#2a6ebb", "D")]
for di, t in enumerate(("PE_in_CD", "PE_in_MD")):
    for si, (nm, df, col_d, col_p, colr, mk) in enumerate(series):
        d = df[df.diet == t].set_index("group").loc[grp]
        xx = np.arange(len(grp)) + (di * 3 + si - 2.5) * .12
        ax.vlines(xx, d.dR2_ci_low.clip(lower=-.6), d.dR2_ci_high, color=colr, lw=1.3)
        ax.scatter(xx, d[col_d], color=colr if t == "PE_in_CD" else "white", edgecolor=colr, marker=mk, s=24, zorder=3,
                   label=f"{nm} − additive" if t == "PE_in_CD" else None)
ax.axhline(0, color="#444", lw=.8)
ax.set_xticks(range(len(grp)))
ax.set_xticklabels([lab[g] for g in grp], fontsize=7)
ax.set_ylim(-.62, .45)
ax.set_ylabel("Held-out ΔR²_group vs additive model", fontsize=7.5)
ax.set_title("D  Held-out benefit of interactions (filled = CD, open = MD; 95% module-bootstrap CI)\n"
             f"Smallest Holm-adjusted p over the 18 primary tests = {min(HO.p_holm_primary12.min(), PT.p_holm_6.min()):.3f}; CIs clipped at −0.6", fontsize=8, loc="left")
ax.legend(fontsize=6.5, frameon=False, loc="lower left")
ax.tick_params(axis="y", labelsize=7)
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)
fig.suptitle("Interactions among reference-signature features in the Phase 4 HGB / CatBoost models (20 co-expression modules, 5 folds × 3 repeats). "
             "Diagnostics on training genes; benefit on held-out genes.", fontsize=8, x=.01, ha="left")
for e in ("png", "pdf", "svg"):
    fig.savefig(FIG / f"fig_interaction_FAST_CIRCse.{e}", dpi=300, bbox_inches="tight")
pd.set_option("display.width", 250)
print(PT.round(4).to_string(index=False))
print(ST.round(3).to_string(index=False))
