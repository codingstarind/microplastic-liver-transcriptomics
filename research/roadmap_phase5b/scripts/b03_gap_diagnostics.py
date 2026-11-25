"""Stage 2: quantitative diagnosis of the internal-vs-external gap (no new model fitting except none; all frozen).
Outputs tables/gap_*.csv and figures/fig1-3. Mouse = experimental unit throughout."""
import itertools

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from b00_common import FIG, P5, PB, REPO, RS, TAB, WORK5, log2cpm, wr2  # noqa: E402
from c00_common import loo_names, target  # noqa: E402
from scipy import stats  # noqa: E402

FIG.mkdir(exist_ok=True)
P = pd.read_csv(P5 / "frozen/gene_level_frozen_predictions.csv.gz", index_col=0)
CO = pd.read_csv(P5 / "frozen/gene_level_frozen_coefficients.csv")
TM = CO[CO.feature == "(training weighted mean y)"].set_index("model").coef_per_SD
DE = P5 / "results/external_de"
man = pd.read_csv(P5 / "results/external_sample_manifest.csv")
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def ext(st, tag):
    t = pd.read_csv(DE / f"{st}__{tag}.csv.gz").set_index("gene")
    return t[np.isfinite(t.log2FC) & np.isfinite(t.lfcSE) & (t.lfcSE > 0)]


def metr(y, p, ref):
    return dict(r2_vs_train_mean=1 - np.sum((y - p) ** 2) / np.sum((y - ref) ** 2), r2_vs_own_mean=1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2),
                spearman=stats.spearmanr(y, p).statistic, rmse=np.sqrt(np.mean((y - p) ** 2)), mae=np.mean(np.abs(y - p)))


# ---------- 1. train vs outer-CV vs external, same unweighted metrics ----------
fx = pd.read_csv(RS / "roadmap_phase4c_diagnostics/tables/fix_experiment_metrics_by_fold.csv")
rows = []
for tgt in ("PE_in_CD", "PE_in_MD"):
    t = target(tgt)
    for mkey, rep, mo in (("F8_ols", "F8", "ols"), ("F8C_ridge", "F8+C", "ridge")):
        key = f"{mkey}__{tgt}"
        g = P.index.intersection(t.dropna().index)
        rows.append(dict(stage="train (all genes, in-sample)", target=tgt, model=mkey, **metr(t.loc[g, "y"].values, P.loc[g, key].values, TM[key])))
        f = fx[(fx.target == tgt) & (fx.rep == rep) & (fx.model == mo)]
        rows.append(dict(stage="outer grouped CV (15 folds, mean)", target=tgt, model=mkey, r2_vs_train_mean=f.r2_oos.mean(),
                         spearman=f.spearman.mean(), rmse=f.rmse.mean(), r2_fold_sd=f.r2_oos.std()))
        for st, tag in (("HONG", "PE_vs_CON"), ("GSE245069", "PSpooled_vs_CON")):
            e = ext(st, tag)
            g = P.index.intersection(e.index)
            rows.append(dict(stage=f"external {st} {tag}", target=tgt, model=mkey, **metr(e.loc[g, "log2FC"].values, P.loc[g, key].values, TM[key])))
GAP = pd.DataFrame(rows)
GAP.to_csv(TAB / "gap_train_cv_external.csv", index=False)

# ---------- 2. target distributions / heavy tails / SE ----------
t294 = pd.read_csv(REPO / "results/deseq2/provisional_deposited/results_PE_in_CD.csv.gz").set_index("gene")
dist = []
for lab, t in (("GSE294180 PE_in_CD", t294), ("GSE294180 PE_in_MD", pd.read_csv(REPO / "results/deseq2/provisional_deposited/results_PE_in_MD.csv.gz").set_index("gene")),
               ("HONG PE_vs_CON", ext("HONG", "PE_vs_CON")), ("HONG PS_vs_CON", ext("HONG", "PS_vs_CON")),
               ("GSE245069 PSpooled", ext("GSE245069", "PSpooled_vs_CON"))):
    g = P.index.intersection(t.dropna(subset=["log2FC", "lfcSE"]).index)
    y, s = t.loc[g, "log2FC"].values, t.loc[g, "lfcSE"].values
    dist.append(dict(target=lab, n_genes=len(g), sd=y.std(), mad_sd=1.4826 * np.median(np.abs(y - np.median(y))), excess_kurtosis=stats.kurtosis(y),
                     frac_abs_gt_2=np.mean(np.abs(y) > 2), median_se=np.median(s), signal_share=1 - np.mean(s**2) / np.var(y),
                     signal_share_robust=1 - np.median(s) ** 2 / (1.4826 * np.median(np.abs(y - np.median(y)))) ** 2))
DIST = pd.DataFrame(dist)
DIST.to_csv(TAB / "gap_target_distributions.csv", index=False)

# ---------- 3. Hong PE errors by expression / uncertainty strata ----------
e = ext("HONG", "PE_vs_CON")
g = P.index.intersection(e.index)
d = pd.DataFrame({"y": e.loc[g, "log2FC"], "s": e.loc[g, "lfcSE"], "bm": e.loc[g, "baseMean"], "p": P.loc[g, "F8_ols__PE_in_CD"],
                  "p_md": P.loc[g, "F8_ols__PE_in_MD"], "src_se": t294.reindex(g).lfcSE})
strat = []
for col, lab in (("bm", "Hong baseMean tercile"), ("s", "Hong lfcSE tercile"), ("src_se", "GSE294180 lfcSE tercile")):
    d["q"] = pd.qcut(d[col], 3, labels=["low", "mid", "high"])
    for q, dd in d.groupby("q", observed=True):
        strat.append(dict(stratum=lab, level=q, n=len(dd), spearman_F8_CD=stats.spearmanr(dd.y, dd.p).statistic,
                          spearman_F8_MD=stats.spearmanr(dd.y, dd.p_md).statistic, mean_resid=float((dd.y - dd.p).mean()),
                          rmse=float(np.sqrt(((dd.y - dd.p) ** 2).mean())), calib_slope=float(np.polyfit(dd.p, dd.y, 1)[0])))
pd.DataFrame(strat).to_csv(TAB / "gap_strata_hong_pe.csv", index=False)

# ---------- 4. external target uncertainty (mouse resampling, DESeq2 refits from b02) ----------
R = pd.read_csv(PB / "results/hong_resampling/hong_pe_target_resampling.csv.gz")
W = R.pivot(index="gene", columns="fit", values="log2FC")
full = W["full"]
gg = P.index.intersection(W.dropna().index)
res_rows = []
for c in W.columns:
    res_rows.append(dict(fit=c, spearman_with_full=stats.spearmanr(W.loc[gg, c], full[gg]).statistic,
                         transfer_spearman_F8_CD=stats.spearmanr(W.loc[gg, c], P.loc[gg, "F8_ols__PE_in_CD"]).statistic,
                         transfer_spearman_F8_MD=stats.spearmanr(W.loc[gg, c], P.loc[gg, "F8_ols__PE_in_MD"]).statistic))
RR = pd.DataFrame(res_rows)
loo = [c for c in W.columns if c.startswith("LOO_")]
n = len(loo)
jk_var = (n - 1) / n * ((W.loc[gg, loo].sub(W.loc[gg, loo].mean(1), axis=0)) ** 2).sum(1)
rel_hong_jk = 1 - jk_var.mean() / full[gg].var()
two = [c for c in W.columns if c.startswith("2v2")]
pairs = [(a, b) for a, b in itertools.combinations(two, 2) if len(set(a.split("_")[2:]) & set(b.split("_")[2:])) == 0]
rel_hong_split = np.mean([stats.spearmanr(W.loc[gg, a], W.loc[gg, b]).statistic for a, b in pairs]) if pairs else np.nan
RR.to_csv(TAB / "gap_hong_target_resampling.csv", index=False)

# GSE294180 PE_in_CD jackknife over CD mice (existing Phase 1/4 LOO refits)
cdl = [m for m in loo_names() if m.startswith("CD")]
L294 = pd.DataFrame({m: target(f"LOO:{m}:PE_in_CD").y for m in cdl})
g2 = gg.intersection(L294.dropna().index).intersection(t294.dropna(subset=["log2FC"]).index)
k = len(cdl)
jk2 = (k - 1) / k * ((L294.loc[g2].sub(L294.loc[g2].mean(1), axis=0)) ** 2).sum(1)
rel_294_jk = 1 - jk2.mean() / t294.loc[g2, "log2FC"].var()
r_obs = stats.spearmanr(t294.loc[g2, "log2FC"], full[g2]).statistic
r_obs_md = stats.spearmanr(pd.read_csv(REPO / "results/deseq2/provisional_deposited/results_PE_in_MD.csv.gz").set_index("gene").reindex(g2).log2FC, full[g2], nan_policy="omit").statistic
def robust_rel(jkv, y):
    mad2 = (1.4826 * np.median(np.abs(y - np.median(y)))) ** 2
    return 1 - np.median(jkv) / mad2


rel_hong_r, rel_294_r = robust_rel(jk_var.values, full[gg].values), robust_rel(jk2.values, t294.loc[g2, "log2FC"].values)
rel = dict(hong_jackknife_reliability_meanvar=rel_hong_jk, hong_jackknife_reliability_robust=rel_hong_r,
           gse294180_PE_in_CD_jackknife_reliability_robust=rel_294_r,
           disattenuated_robust_PE_in_CD_vs_HONG=r_obs / np.sqrt(rel_hong_r * rel_294_r) if min(rel_hong_r, rel_294_r) > 0 else np.nan,
           transfer_rho_range_over_15_mouse_resamples_F8_CD=f"{RR[RR.fit != 'full'].transfer_spearman_F8_CD.min():.3f}..{RR[RR.fit != 'full'].transfer_spearman_F8_CD.max():.3f}",
           transfer_rho_range_over_15_mouse_resamples_F8_MD=f"{RR[RR.fit != 'full'].transfer_spearman_F8_MD.min():.3f}..{RR[RR.fit != 'full'].transfer_spearman_F8_MD.max():.3f}",
           hong_jackknife_reliability=rel_hong_jk, hong_disjoint_2v2_split_correlation=rel_hong_split, n_disjoint_split_pairs=len(pairs),
           gse294180_PE_in_CD_jackknife_reliability=rel_294_jk, n_cd_loo=k,
           between_study_spearman_PE_in_CD_vs_HONG=r_obs, between_study_spearman_PE_in_MD_vs_HONG=r_obs_md,
           disattenuated_PE_in_CD_vs_HONG=r_obs / np.sqrt(max(rel_hong_jk, 1e-6) * max(rel_294_jk, 1e-6)) if rel_294_jk > 0 and rel_hong_jk > 0 else np.nan,
           ceiling_transfer_spearman_hong=np.sqrt(max(rel_hong_jk, 0)))
pd.Series(rel).to_csv(TAB / "gap_reliability_and_disattenuation.csv")

# ---------- 5. mouse-level PCA of Hong (all 15) and the permutation-null mechanism ----------
X = pd.read_csv(WORK5 / "HONG_expression.csv.gz", index_col=0)
mh = man[man.study == "HONG"].set_index("mouse_id")
X = X[mh.index]
Lc = log2cpm(X)
Lc = Lc[(X >= 10).sum(axis=1) >= 3]
top = Lc.var(axis=1).sort_values(ascending=False).index[:2000]
Z = Lc.loc[top].sub(Lc.loc[top].mean(axis=1), axis=0)
u, sv, vt = np.linalg.svd(Z.values, full_matrices=False)
ve = sv**2 / (sv**2).sum()
ref_fast = pd.read_csv(RS / "roadmap_phase3/reference_effects/FAST_GSE107787_FASTvsFED_ZTadj.csv.gz").set_index("gene").log2FC
pcs = []
for i in range(5):
    load = pd.Series(u[:, i] * sv[i], index=top)
    gi = load.index.intersection(P.index)
    pcs.append(dict(pc=i + 1, var_explained=ve[i], loading_spearman_F8_CD=stats.spearmanr(load[gi], P.loc[gi, "F8_ols__PE_in_CD"]).statistic,
                    loading_spearman_REF_FAST=stats.spearmanr(load[load.index.intersection(ref_fast.index)], ref_fast[load.index.intersection(ref_fast.index)]).statistic,
                    **{f"score_{m}": vt[i, j] for j, m in enumerate(mh.index)}))
PCS = pd.DataFrame(pcs)
PCS.to_csv(TAB / "gap_hong_mouse_pca.csv", index=False)
# 20 labellings of the 6 PE/CON mice: transfer rho vs alignment of labelling with PC scores
six = list(mh.index[mh.group.isin(["PE", "CON"])])
L6 = Lc.loc[:, six]
gp = P.index.intersection(L6.index)
pr = stats.rankdata(P.loc[gp, "F8_ols__PE_in_CD"])
lab_rows = []
for exp in itertools.combinations(range(6), 3):
    a = np.zeros(6, bool)
    a[list(exp)] = True
    dd = L6.loc[gp].values[:, a].mean(1) - L6.loc[gp].values[:, ~a].mean(1)
    sc = {f"pc{i + 1}_score_diff": vt[i, [list(mh.index).index(m) for m in np.array(six)[a]]].mean() - vt[i, [list(mh.index).index(m) for m in np.array(six)[~a]]].mean() for i in range(3)}
    lab_rows.append(dict(exposed=";".join(np.array(six)[a]), is_true_labelling=set(np.array(six)[a]) == set(mh.index[mh.group == "PE"]),
                         transfer_rho=stats.spearmanr(pr, dd).statistic, **sc))
LAB = pd.DataFrame(lab_rows)
LAB.to_csv(TAB / "gap_hong_20_labellings.csv", index=False)
_fit = np.linalg.lstsq(np.column_stack([np.ones(len(LAB)), LAB.pc1_score_diff, LAB.pc2_score_diff, LAB.pc3_score_diff]), LAB.transfer_rho, rcond=None)[0]
_pred = np.column_stack([np.ones(len(LAB)), LAB.pc1_score_diff, LAB.pc2_score_diff, LAB.pc3_score_diff]) @ _fit
LAB_R2 = 1 - np.sum((LAB.transfer_rho - _pred) ** 2) / np.sum((LAB.transfer_rho - LAB.transfer_rho.mean()) ** 2)
true_rank = int((LAB.transfer_rho >= LAB[LAB.is_true_labelling].transfer_rho.iloc[0]).sum())

# ---------- 6. study-shift table ----------
q = pd.read_csv(P5 / "results/external_sample_qc.csv")
s4d = pd.read_csv(RS / "roadmap_phase4d/results/per_mouse_scores.csv", index_col=0)
s5 = pd.read_csv(P5 / "per_sample_scores.csv", index_col=0)
shift = [
    dict(factor="polymer / particle", GSE294180="PE microspheres 10-150 um", HONG="PE (size NOT STATED)", GSE245069="PS 80 nm / 0.5 / 5 um", identifiable="no"),
    dict(factor="dose / route", GSE294180="fixed mass/day oral gavage (~70-80 mg/kg/d est.)", HONG="in food, dose NOT STATED", GSE245069="1 mg/L drinking water", identifiable="no"),
    dict(factor="duration / age", GSE294180="8 wk from 8 wk old", HONG="71 wk (age at start NOT STATED)", GSE245069="12 wk from 6 wk old", identifiable="no"),
    dict(factor="control diet", GSE294180="purified D09100304 (CD) / MASH diet (MD)", HONG="normal chow", GSE245069="standard chow (NOT STATED)", identifiable="no"),
    dict(factor="strain / sex", GSE294180="C57BL/6J male", HONG="C57BL/6 male", GSE245069="ICR outbred male", identifiable="no"),
    dict(factor="quantifier / annotation", GSE294180="Salmon/tximport 'no' counts, GENCODE M35 symbols", HONG="HISAT2/StringTie prepDE, GRCm39 RefSeq", GSE245069="deposited read_count, GRCm38 NCBI", identifiable="partly (mapping audited)"),
    dict(factor="counts per read pair (median)", GSE294180="0.775", HONG="1.641 (~2x: coverage-derived per-mate counts)", GSE245069="0.734", identifiable="yes (unit check)"),
    dict(factor="read layout (bases/spot)", GSE294180="NextSeq 2000 PE (pilot ~2x)", HONG="202 (2x101)", GSE245069="300 (2x150)", identifiable="yes"),
    dict(factor="depth balance exposed/control", GSE294180="PE libraries LOWER depth/yield (confounded)", HONG="0.98", GSE245069="0.97", identifiable="yes"),
    dict(factor="target heavy tail (top-1% share of sum sq, bm>500)", GSE294180="0.77", HONG="0.33", GSE245069="0.43", identifiable="yes"),
    dict(factor="baseline FAST singscore of controls (mean)", GSE294180=f"{s4d[s4d.group == 'CD+Veh'].FAST.mean():.3f} (CD+Veh)",
         HONG=f"{s5[(s5.study == 'HONG') & (s5.group == 'CON')].FAST_singscore.mean():.3f}", GSE245069=f"{s5[(s5.study == 'GSE245069') & (s5.group == 'CON')].FAST_singscore.mean():.3f}", identifiable="descriptive"),
    dict(factor="baseline MASH singscore of controls (mean)", GSE294180=f"{s4d[s4d.group == 'CD+Veh'].MASH.mean():.3f}",
         HONG=f"{s5[(s5.study == 'HONG') & (s5.group == 'CON')].MASH_singscore.mean():.3f}", GSE245069=f"{s5[(s5.study == 'GSE245069') & (s5.group == 'CON')].MASH_singscore.mean():.3f}", identifiable="descriptive"),
    dict(factor="muscle CPM (median)", GSE294180="1076", HONG=f"{q[q.study == 'HONG'].muscle_cpm_sum.median():.0f}", GSE245069=f"{q[q.study == 'GSE245069'].muscle_cpm_sum.median():.0f}", identifiable="yes"),
    dict(factor="mice per arm in the PE contrast", GSE294180="4 v 5 (CD); 5 v 5 (MD)", HONG="3 v 3 (shared CON)", GSE245069="12 v 4 (PS; shared CON)", identifiable="yes"),
]
pd.DataFrame(shift).to_csv(TAB / "gap_study_shift.csv", index=False)

# ---------- figures ----------
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
G = GAP[GAP.model == "F8_ols"]
for k, tgt in enumerate(("PE_in_CD", "PE_in_MD")):
    gg_ = G[G.target == tgt]
    ax[0].plot(range(len(gg_)), gg_.spearman.values, "o-", label=f"F8-OLS trained on {tgt}", color=["#2a78d6", "#eb6834"][k])
ax[0].set_xticks(range(len(gg_)))
ax[0].set_xticklabels(["train", "outer CV\n(held-out genes)", "Hong PE\n(new mice)", "GSE245069 PS\n(new mice)"], fontsize=8)
ax[0].axhline(0, color="#8a8984", lw=.6)
ax[0].set_ylabel("Spearman(pred, observed log2FC)")
ax[0].legend(fontsize=7, frameon=False)
ax[0].set_title("Train / CV / external gap (same metric)")
for lab, c in (("GSE294180 PE_in_CD", "#8a8984"), ("HONG PE_vs_CON", "#2a78d6"), ("GSE245069 PSpooled", "#eb6834")):
    t = {"GSE294180 PE_in_CD": t294, "HONG PE_vs_CON": ext("HONG", "PE_vs_CON"), "GSE245069 PSpooled": ext("GSE245069", "PSpooled_vs_CON")}[lab]
    y = t.loc[P.index.intersection(t.dropna(subset=["log2FC"]).index), "log2FC"].clip(-3, 3)
    ax[1].hist(y, bins=120, histtype="step", density=True, color=c, label=lab, lw=1.3)
ax[1].set_yscale("log")
ax[1].set_xlabel("log2FC (clipped ±3)")
ax[1].legend(fontsize=7, frameon=False)
ax[1].set_title("Target distributions (log density: tails)")
ax[2].scatter(LAB.pc1_score_diff, LAB.transfer_rho, c=np.where(LAB.is_true_labelling, "#eb6834", "#2a78d6"), s=40)
ax[2].set_xlabel("labelling's difference in Hong PC1 score")
ax[2].set_ylabel("transfer Spearman (F8-OLS, log2CPM proxy)")
ax[2].set_title("20 relabellings of 6 Hong mice\n(orange = true PE labelling)")
fig.tight_layout()
fig.savefig(FIG / "fig1_gap_targets_labellings.png", dpi=130)
plt.close(fig)
with pd.option_context("display.width", 250, "display.max_columns", 30):
    print(GAP.round(3).to_string(index=False))
    print(DIST.round(3).to_string(index=False))
    print(pd.DataFrame(strat).round(3).to_string(index=False))
    print(RR.round(3).to_string(index=False))
    print(pd.Series(rel).to_string())
    print(PCS[["pc", "var_explained", "loading_spearman_F8_CD", "loading_spearman_REF_FAST"]].round(3).to_string(index=False))
    print("R2 of transfer_rho ~ PC1..3 differences over 20 labellings:", round(LAB_R2, 3), "true labelling rank (1=highest):", true_rank)
    open(TAB / "gap_labelling_summary.txt", "w").write(f"R2(transfer_rho ~ PC1+PC2+PC3 score differences, 20 labellings) = {LAB_R2:.3f}; true labelling rank = {true_rank} of 20")
