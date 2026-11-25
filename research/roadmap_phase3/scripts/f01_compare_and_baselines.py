"""Stages E-I (Python part): reference robustness summaries, reference-vs-reference comparison,
collinearity, and the two baselines against GSE294180 PE effects (analysis_plan.md sections 6-10)."""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from p3lib import P3, RE, REPO, compare, read_gz_csv, target, universe  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

OUT = P3 / "results" / "comparisons"
OUT.mkdir(parents=True, exist_ok=True)
U, US = universe("main"), universe("sensitivity")


def load(label, eff="log2FC"):
    t = pd.read_csv(RE / f"{label}.csv.gz").set_index("gene")
    return t.rename(columns={eff: "log2FC"})[["log2FC", "padj"]]


REFS = {
    "MASH1": "MASH1_GSE316023_GANvsChow_TA",
    "MASH2": "MASH2_GSE246328_GANvsChow_24w",
    "FAST": "FAST_GSE107787_FASTvsFED_ZTadj",
    "FAST_clean": "FASTclean_GSE107787_ZT0ZT4_muscleLow",
    "FAST2": "FAST2_GSE225695_FastvsFed",
    "PPARA_20v0": "PPARA_GSE282998_WT_20v0",
    "CIRC_T1": "CIRC_GSE117134_WTfemale_T1_ZT6_minus_ZT2",
    "CIRC_T2": "CIRC_GSE117134_WTfemale_T2_ZT12_minus_ZT0",
}
R = {k: load(v) for k, v in REFS.items()}
cls = read_gz_csv(P3 / "results" / "ppara" / "genotype_classification.csv.gz", index_col=0)
dep = set(cls.index[cls["class"] == "PPARa_dependent"])
R["PPARA_dependent"] = R["PPARA_20v0"].loc[R["PPARA_20v0"].index.intersection(list(dep))]
res = {}

# ---------------- fasting robustness ----------------
f, fc, f2 = R["FAST"], R["FAST_clean"], R["FAST2"]
fdeg = f.index[(f.padj < .05) & (f.log2FC.abs() >= .38)]
g = f.index.intersection(fc.index).intersection(list(U))
res["fasting"] = dict(
    rho_main_vs_clean=float(spearmanr(f.loc[g, "log2FC"], fc.loc[g, "log2FC"]).statistic),
    sign_agree_clean_on_main_DEG=float((np.sign(f.loc[fdeg, "log2FC"]) == np.sign(fc.loc[fdeg, "log2FC"])).mean()),
    n_main_deg=int(len(fdeg)),
    rho_main_vs_GSE225695=float(spearmanr(f.loc[g.intersection(f2.index), "log2FC"], f2.loc[g.intersection(f2.index), "log2FC"], nan_policy="omit").statistic),
    sign_agree_GSE225695_on_main_DEG=float((np.sign(f.loc[fdeg.intersection(f2.index), "log2FC"]) == np.sign(f2.loc[fdeg.intersection(f2.index), "log2FC"])).mean()),
)
r = res["fasting"]
res["fasting"]["label"] = "CONDITIONAL" if (r["rho_main_vs_clean"] < .60 or r["sign_agree_clean_on_main_DEG"] < .80) else "STABLE"
inter = json.loads((P3 / "results" / "other_refs_run_info.json").read_text())["fasting"]
res["fasting"]["frac_deg_with_interaction"] = inter["frac_deg_with_interaction_fdr05"]
res["fasting"]["time_dependent_label"] = inter["frac_deg_with_interaction_fdr05"] > 0.10

# ---------------- PPARa robustness ----------------
p20 = R["PPARA_20v0"]
pp = {k: load(f"PPARA_GSE282998_WT_{k}") for k in ("5v0", "10v0", "20v0_muscleAdj", "trend")}
ko = load("PPARA_GSE282998_KO_20v0")
gp = p20.index.intersection(list(U))
sig20 = gp[(p20.loc[gp, "padj"] < .05)]
res["ppara"] = {f"rho_20v0_vs_{k}": float(spearmanr(p20.loc[gp, "log2FC"], v.loc[gp, "log2FC"], nan_policy="omit").statistic) for k, v in pp.items()}
res["ppara"].update({f"sign_agree_{k}_on_20v0_sig": float((np.sign(p20.loc[sig20, "log2FC"]) == np.sign(v.loc[sig20, "log2FC"])).mean()) for k, v in pp.items()})
res["ppara"]["rho_WT20_vs_KO20"] = float(spearmanr(p20.loc[gp, "log2FC"], ko.loc[gp, "log2FC"], nan_policy="omit").statistic)
res["ppara"]["n_20v0_padj05"] = int(len(sig20))
res["ppara"]["class_counts"] = cls["class"].value_counts().to_dict()
cons = (np.sign(pp["5v0"].loc[sig20, "log2FC"]) == np.sign(p20.loc[sig20, "log2FC"])) & (
    np.sign(pp["10v0"].loc[sig20, "log2FC"]) == np.sign(p20.loc[sig20, "log2FC"]))
res["ppara"]["frac_20v0_sig_dose_consistent"] = float(cons.mean())
key = ["Cyp4a14", "Cyp4a10", "Acot1", "Ehhadh", "Cpt1a", "Hmgcs2", "Pdk4", "Angptl4", "Fgf21", "Plin2", "Cd36",
       "Acox1", "Fasn", "Acly", "Hmgcr", "Sqle", "Ppara"]
kt = pd.DataFrame({k: v.reindex(key)["log2FC"] for k, v in {"5v0": pp["5v0"], "10v0": pp["10v0"], "20v0": p20, "KO_20v0": ko}.items()})
kt["class"] = cls["class"].reindex(key)
kt.to_csv(OUT / "ppara_key_genes.csv")

# ---------------- circadian ----------------
fr = read_gz_csv(P3 / "results" / "circadian" / "CIRC_GSE117134_WTfemale_rhythm.csv.gz", index_col=0)
mr = read_gz_csv(P3 / "results" / "circadian" / "CIRCsens_GSE107787_FEDmale_rhythm.csv.gz", index_col=0)
both = fr.index[(fr.fdr_harmonic < .05)].intersection(mr.index[mr.fdr_harmonic < .05])
d = np.angle(np.exp(1j * 2 * np.pi * (fr.loc[both, "peak_phase_ZT"] - mr.loc[both, "peak_phase_ZT"]) / 24)) * 24 / (2 * np.pi)
res["circadian"] = dict(n_rhythmic_female=int((fr.fdr_harmonic < .05).sum()), n_rhythmic_male_fed=int((mr.fdr_harmonic < .05).sum()),
                        n_both=int(len(both)), median_abs_phase_diff_h=float(np.median(np.abs(d))),
                        frac_phase_within_3h=float((np.abs(d) <= 3).mean()))
clock = ["Bmal1", "Clock", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2", "Nr1d1", "Nr1d2", "Dbp", "Tef", "Nfil3"]
fr.reindex(clock)[["amplitude_log2", "peak_phase_ZT", "phase_se_h", "fdr_harmonic"]].join(
    mr.reindex(clock)[["peak_phase_ZT", "fdr_harmonic"]], rsuffix="_maleFED").to_csv(OUT / "circadian_clock_genes.csv")

# ---------------- reference vs reference ----------------
names = ["MASH1", "MASH2", "FAST", "FAST2", "PPARA_20v0", "PPARA_dependent", "CIRC_T1", "CIRC_T2"]
pairs, mat = [], pd.DataFrame(index=names, columns=names, dtype=float)
for i, a in enumerate(names):
    for b in names[i:]:
        c = compare(R[a], R[b], U, a, b)
        mat.loc[a, b] = mat.loc[b, a] = c["spearman"]
        if a != b:
            pairs.append(c)
pd.DataFrame(pairs).to_csv(OUT / "reference_pairwise.csv", index=False)
mat.to_csv(OUT / "reference_spearman_matrix.csv")
# collinearity on genes present in all primary signed references
core = ["MASH1", "FAST", "PPARA_20v0", "CIRC_T1", "CIRC_T2"]
common = set(U)
for k in core:
    common &= set(R[k].dropna().index)
Z = pd.DataFrame({k: R[k].loc[sorted(common), "log2FC"].rank() for k in core})
Z = (Z - Z.mean()) / Z.std()
vif = {}
for k in core:
    others = Z.drop(columns=k).values
    beta, *_ = np.linalg.lstsq(np.c_[np.ones(len(Z)), others], Z[k].values, rcond=None)
    pred = np.c_[np.ones(len(Z)), others] @ beta
    r2 = 1 - ((Z[k] - pred) ** 2).sum() / ((Z[k] - Z[k].mean()) ** 2).sum()
    vif[k] = dict(R2_on_others=float(r2), VIF=float(1 / (1 - r2)))
res["collinearity"] = dict(n_common_genes=len(common), rank_corr=Z.corr().round(3).to_dict(), vif=vif,
                           flagged_pairs=[(a, b, round(float(Z.corr().loc[a, b]), 3)) for i, a in enumerate(core) for b in core[i + 1:]
                                          if abs(Z.corr().loc[a, b]) >= .7])
# pathway-level: mean log2FC per core set (sign only comparable) + fry where available
sets = pd.read_csv(P3 / "results" / "gene_sets.csv")
CORE = ["LIPOGENESIS_CORE", "CHOL_SYNTHESIS_CORE", "SANDERSON_PPARA_TARGETS", "HALLMARK_FATTY_ACID_METABOLISM",
        "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
        "HALLMARK_XENOBIOTIC_METABOLISM"]
tg = {"Diet_in_Veh": target("Diet_in_Veh"), "PE_in_CD": target("PE_in_CD"), "PE_in_MD": target("PE_in_MD")}
allsig = {**{k: R[k] for k in names if k != "PPARA_dependent"}, **tg}
pm = pd.DataFrame({k: [v.loc[v.index.intersection(list(set(sets.gene[sets.set == s]) & U)), "log2FC"].median() for s in CORE]
                   for k, v in allsig.items()}, index=CORE)
pm.to_csv(OUT / "core_pathway_median_log2FC.csv")

# ---------------- baselines ----------------
base, loo_rows = [], []
for tname, t in tg.items():
    for rname in ["MASH1", "MASH2", "FAST", "FAST2", "PPARA_20v0", "PPARA_dependent", "CIRC_T1", "CIRC_T2"]:
        c = compare(t, R[rname], U, tname, rname)
        cs = compare(t, R[rname], US, tname, rname)
        c["spearman_sensitivity_universe"] = cs["spearman"]
        base.append(c)
loodirs = sorted((REPO / "results" / "deseq2" / "loo").glob("provisional_deposited__without_*"))
for dd in loodirs:
    for tname in tg:
        t2 = pd.read_csv(dd / f"results_{tname}.csv.gz").set_index("gene")[["log2FC", "padj"]]
        for rname in ["MASH1", "MASH2", "FAST", "PPARA_20v0", "CIRC_T1", "CIRC_T2"]:
            gq = sorted(set(t2.index) & set(R[rname].index) & U)
            loo_rows.append(dict(target=tname, ref=rname, dropped=dd.name.split("__")[-1],
                                 spearman=spearmanr(t2.loc[gq, "log2FC"], R[rname].loc[gq, "log2FC"], nan_policy="omit").statistic))
B = pd.DataFrame(base)
L = pd.DataFrame(loo_rows).groupby(["target", "ref"]).spearman.agg(["min", "max"]).reset_index()
B = B.merge(L.rename(columns={"target": "a", "ref": "b", "min": "loo_min", "max": "loo_max"}), on=["a", "b"], how="left")
B.to_csv(OUT / "baselines.csv", index=False)
hist = json.loads((REPO / "results" / "concordance" / "provisional_deposited" / "concordance_summary.json").read_text())
res["historical_split_reference"] = hist["split_reference"]["PE_in_CD vs Diet_in_Veh"]["spearman_all_genes"]
(OUT / "comparison_summary.json").write_text(json.dumps(res, indent=2, default=float))

# ---------------- figures ----------------
fig, ax = plt.subplots(figsize=(7.5, 6.2))
m = mat.astype(float)
im = ax.imshow(m.values, cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(range(len(names)))
ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(len(names)))
ax.set_yticklabels(names, fontsize=8)
for i in range(len(names)):
    for j in range(len(names)):
        ax.text(j, i, f"{m.values[i, j]:.2f}", ha="center", va="center", fontsize=7)
fig.colorbar(im, label="Spearman (gene effects, main universe)")
ax.set_title("Reference-axis similarity (rank correlation; scales not compared)", fontsize=9)
fig.tight_layout()
fig.savefig(P3 / "figures" / "reference_similarity_heatmap.png", dpi=110)

fig, ax = plt.subplots(1, 2, figsize=(14, 4.5))
refs_plot = ["MASH1", "MASH2", "FAST", "FAST2", "PPARA_20v0", "PPARA_dependent", "CIRC_T1", "CIRC_T2"]
w = .27
for k, tname in enumerate(["Diet_in_Veh", "PE_in_CD", "PE_in_MD"]):
    sub = B[B.a == tname].set_index("b").reindex(refs_plot)
    x = np.arange(len(refs_plot)) + (k - 1) * w
    ax[0].bar(x, sub.spearman, w, label=tname)
    lo, hi = sub.loo_min, sub.loo_max
    ax[0].errorbar(x, sub.spearman, yerr=[sub.spearman - lo, hi - sub.spearman], fmt="none", c="k", lw=.8)
    ax[1].bar(x, sub.jaccard, w, label=tname)
ax[0].axhline(0, c="k", lw=.6)
ax[0].axhline(res["historical_split_reference"]["median"], ls="--", c="grey", lw=.8)
ax[0].set_xticks(range(len(refs_plot)))
ax[0].set_xticklabels(refs_plot, rotation=30, fontsize=8)
ax[0].set_ylabel("Spearman of gene effects (bars: target LOO range)")
ax[0].set_title("Baseline 2: effect correlation (dashed: historical same-data split-reference, PE_in_CD)", fontsize=9)
ax[0].legend(fontsize=7)
ax[1].set_xticks(range(len(refs_plot)))
ax[1].set_xticklabels(refs_plot, rotation=30, fontsize=8)
ax[1].set_ylabel("Jaccard of DEG lists (padj<0.05, |log2FC|>=0.38)")
ax[1].set_title("Baseline 1: DEG overlap", fontsize=9)
fig.tight_layout()
fig.savefig(P3 / "figures" / "baseline_comparisons.png", dpi=110)

fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
ax[0].scatter(f.loc[g, "log2FC"], fc.loc[g, "log2FC"], s=2, alpha=.3)
ax[0].set_xlabel("FAST vs FED, all ZT (18 v 18)")
ax[0].set_ylabel("FAST vs FED, ZT0/4 muscle-low (4 v 6)")
ax[0].set_title(f"Spearman {res['fasting']['rho_main_vs_clean']:.2f}", fontsize=9)
pz = read_gz_csv(P3 / "results" / "fasting" / "GSE107787_fasting_effect_by_ZT.csv.gz", index_col=0)
ax[1].hist(pz.interaction_F_p, bins=50, color="grey")
ax[1].set_title(f"ZT x state interaction p-values; {inter['n_interaction_fdr05']} genes FDR<0.05", fontsize=9)
gg = [x for x in ["Pdk4", "Hmgcs2", "Fasn", "Srebf1", "Hmgcr", "Cyp4a14", "Angptl4", "Fgf21"] if x in pz.index]
zts = [c for c in pz.columns if c.startswith("fast_effect_")]
for x in gg:
    ax[2].plot(range(6), pz.loc[x, zts].values, marker="o", label=x)
ax[2].set_xticks(range(6))
ax[2].set_xticklabels([z.replace("fast_effect_", "") for z in zts])
ax[2].axhline(0, c="k", lw=.5)
ax[2].set_ylabel("fasting effect (log2)")
ax[2].legend(fontsize=7)
ax[2].set_title("Fasting effect depends on ZT", fontsize=9)
fig.tight_layout()
fig.savefig(P3 / "figures" / "fasting_reference.png", dpi=110)

fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
ax[0].scatter(p20.loc[gp, "log2FC"].clip(-5, 8), ko.loc[gp, "log2FC"].clip(-5, 8), s=2, alpha=.3,
              c=np.where(cls.reindex(gp)["class"] == "PPARa_dependent", "tab:red", "tab:grey"))
ax[0].set_xlabel("WT B6129SF2/J GW7647 20 vs 0 log2FC")
ax[0].set_ylabel("Ppara-null 20 vs 0 log2FC")
ax[0].set_title("red = PPARa-dependent by formal interaction", fontsize=9)
kk = kt.dropna(subset=["20v0"]).head(12)
for x in kk.index:
    ax[1].plot([0, 5, 10, 20], [0, kk.loc[x, "5v0"], kk.loc[x, "10v0"], kk.loc[x, "20v0"]], marker="o", label=x)
ax[1].set_xlabel("GW7647 mg/kg/day")
ax[1].set_ylabel("log2FC vs vehicle (WT)")
ax[1].legend(fontsize=6, ncol=2)
fig.tight_layout()
fig.savefig(P3 / "figures" / "ppara_reference.png", dpi=110)

fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
rh = fr[fr.fdr_harmonic < .05]
ax[0].hist(rh.peak_phase_ZT, bins=24, range=(0, 24), color="tab:purple")
ax[0].set_xlabel("peak phase (ZT h), rhythmic genes FDR<0.05 (female WT)")
ax[0].set_title(f"{len(rh)} rhythmic genes", fontsize=9)
ax[1].scatter(fr.loc[both, "peak_phase_ZT"], mr.loc[both, "peak_phase_ZT"], s=3, alpha=.4)
for x in clock:
    if x in both:
        ax[1].annotate(x, (fr.loc[x, "peak_phase_ZT"], mr.loc[x, "peak_phase_ZT"]), fontsize=7)
ax[1].set_xlabel("peak ZT, GSE117134 female WT")
ax[1].set_ylabel("peak ZT, GSE107787 male FED")
ax[1].set_title(f"genes rhythmic in both: {len(both)}; |phase diff|<=3h: {res['circadian']['frac_phase_within_3h']:.2f}", fontsize=9)
fig.tight_layout()
fig.savefig(P3 / "figures" / "circadian_reference.png", dpi=110)
print(json.dumps(res, indent=1, default=float)[:5000])
print(B[["a", "b", "n_genes", "spearman", "loo_min", "loo_max", "spearman_sensitivity_universe", "sign_agree_on_b_sig",
         "n_deg_a", "n_deg_b", "n_deg_overlap", "jaccard", "hypergeom_p", "n_overlap_concordant", "n_overlap_discordant"]].round(3).to_string())
print(mat.round(2).to_string())
