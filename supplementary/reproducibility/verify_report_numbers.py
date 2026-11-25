"""Check every number quoted in the final report against the repository file that produced it.

Read-only: it opens saved result files and recomputes nothing. Writes
supplementary/tables/report_number_crosswalk.csv and exits with status 1 if any value disagrees.

Run from the repository root:
    python supplementary/reproducibility/verify_report_numbers.py
"""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "research"
rows = []


def check(section, item, reported, value, source, decimals=None):
    """reported: the value as printed in the report (string); value: the number in the source file."""
    if decimals is None:
        decimals = len(reported.split(".")[1]) if "." in reported else 0
    ok = round(float(value), decimals) == round(float(reported), decimals)
    rows.append(dict(report_section=section, item=item, reported=reported,
                     source_value=round(float(value), 6), source_file=source, match=ok))


# ---- Section 2.2 preprocessing
fs = json.load(open(ROOT / "results/tables/filter_summary.json"))
check("2.2", "Main PE gene entries", "55836", fs["n_genes_input"], "results/tables/filter_summary.json")
check("2.2", "Genes kept (>=10 reads in >=4 mice)", "13777", fs["n_genes_kept"], "results/tables/filter_summary.json")
check("2.2", "Percent of counts kept", "99.86", 100 * fs["fraction_of_total_counts_kept"], "results/tables/filter_summary.json")

# ---- Fig. 1 PCA after muscle-gene exclusion
ve = pd.read_csv(ROOT / "results/tables/pca_variance_explained.csv").set_index("feature_set")
src = "results/tables/pca_variance_explained.csv"
check("Fig. 1 / 4.1", "PC1 % after muscle exclusion", "20.8", ve.loc["top500_log2norm_muscle_genes_removed", "PC1_pct"], src)
check("Fig. 1 / 4.1", "PC2 % after muscle exclusion", "17.9", ve.loc["top500_log2norm_muscle_genes_removed", "PC2_pct"], src)
mg = pd.read_csv(ROOT / "results/tables/muscle_tracking_genes.csv")
check("3.1", "Muscle-associated genes excluded", "135", len(mg), "results/tables/muscle_tracking_genes.csv")
ari = pd.read_csv(ROOT / "results/tables/clustering_label_agreement.csv")
a = ari[ari.clustering.str.contains("muscle") & (ari.k == 2) & (ari.label == "diet")].adjusted_rand_index.iloc[0]
check("4.1", "ARI, k = 2 vs diet (muscle genes removed)", "1.00", a, "results/tables/clustering_label_agreement.csv")

# ---- DESeq2 gene counts
cs = pd.read_csv(ROOT / "results/deseq2/provisional_deposited/contrast_summary.csv").set_index("contrast")
src = "results/deseq2/provisional_deposited/contrast_summary.csv"
check("4.1", "Diet-associated genes (BH-FDR < 0.05)", "2262", cs.loc["Diet_in_Veh", "n_padj_lt_0.05"], src)
check("4.1", "PE-associated genes under CD", "858", cs.loc["PE_in_CD", "n_padj_lt_0.05"], src)
check("4.1", "PE-associated genes under MD", "472", cs.loc["PE_in_MD", "n_padj_lt_0.05"], src)

# ---- Fig. 2 / 4.1 fry
fr = pd.read_csv(ROOT / "results/pathways/provisional_deposited/limma_fry_camera.csv")
src = "results/pathways/provisional_deposited/limma_fry_camera.csv"
pp = fr[fr.pathway == "SANDERSON_PPARA_TARGETS"].set_index("contrast")
check("Fig. 2 / 4.1", "fry FDR, PPARa targets, CD", "0.0406", pp.loc["PE_in_CD", "fry_fdr"], src)
check("Fig. 2 / 4.1", "fry FDR, PPARa targets, MD", "0.5175", pp.loc["PE_in_MD", "fry_fdr"], src)
hm = fr[(fr.family == "Hallmark") & fr.contrast.isin(["PE_in_CD", "PE_in_MD"])]
check("4.1", "Hallmark sets with fry FDR < 0.05 (CD + MD)", "0", int((hm.fry_fdr < 0.05).sum()), src)

# ---- Table 2: eight pathway scores (BH across 8 scores within each diet)
t2 = pd.read_csv(R / "pathway_followup/TABLE_results_concise.csv").set_index(["pathway", "contrast"])
src = "research/pathway_followup/TABLE_results_concise.csv"
table2 = {"Lipogenesis": ("0.001", "0.0002"), "Chol_synthesis": ("0.012", "0.0001"), "FAST": ("0.016", "0.018"),
          "Chol_homeostasis": ("0.125", "0.016"), "DNA_repair": ("0.472", "0.016"), "PPARA_targets": ("0.064", "0.339"),
          "mTORC1": ("0.125", "0.125"), "Glycolysis": ("0.317", "0.418")}
for pw, (qcd, qmd) in table2.items():
    v = t2.loc[(pw, "PE_in_CD"), "q_BH_8scores"]
    if pw == "Lipogenesis":  # reported as "< 0.001"
        rows.append(dict(report_section="Table 2", item=f"{pw} q, CD", reported="< 0.001", source_value=round(v, 6),
                         source_file=src, match=bool(v < 0.001)))
    else:
        check("Table 2", f"{pw} q, CD", qcd, v, src)
    check("Table 2", f"{pw} q, MD", qmd, t2.loc[(pw, "PE_in_MD"), "q_BH_8scores"], src)

# ---- 4.2 classifier: cross-study validation and Table 3
em = pd.read_csv(R / "roadmap_phase4e/results/external_metrics.csv")
src = "research/roadmap_phase4e/results/external_metrics.csv"
for split, lab in (("S1_train316TA_test246", "MASH-1 -> MASH-2"), ("S2_train246_test316TA", "MASH-2 -> MASH-1")):
    r = em[(em.task == "mash") & (em.split == split) & (em.rep == "R2") & (em.model == "l2_logistic")].iloc[0]
    check("4.2", f"ROC-AUC {lab}", "1.00", r.auc, src)
    check("4.2", f"Balanced accuracy {lab}", "1.00", r.bal_acc, src)
    check("4.2", f"Test mice {lab}", "15" if split.startswith("S1") else "16", r.n, src)
    check("3.3", f"Selected C {lab}", "10", json.loads(r.best_params)["clf__C"], src)
ps = json.load(open(R / "roadmap_phase4e/results/permutation_summary.json"))
check("4.2", "Permutation p, MASH-1 -> MASH-2", "0.0995", ps["S1_train316TA_test246"]["p"], "research/roadmap_phase4e/results/permutation_summary.json")
check("4.2", "Permutation p, MASH-2 -> MASH-1", "0.0348", ps["S2_train246_test316TA"]["p"], "research/roadmap_phase4e/results/permutation_summary.json")
t3 = pd.read_csv(R / "roadmap_phase4e/results/gse294180_frozen_model_contrasts.csv")
t3 = t3[t3.iloc[:, 0] == "R2_l2_logistic"].set_index("contrast")
src = "research/roadmap_phase4e/results/gse294180_frozen_model_contrasts.csv"
for c, (shift, p, ph) in {"PE_in_CD": ("0.835", ".0362", ".0724"), "PE_in_MD": ("1.663", ".000208", ".000623"),
                          "Interaction": ("0.828", ".1175", ".1175")}.items():
    check("Table 3", f"{c} shift (log-odds)", shift, t3.loc[c, "estimate"], src)
    check("Table 3", f"{c} raw p", p, t3.loc[c, "p"], src)
    check("Table 3", f"{c} Holm p", ph, t3.loc[c, "p_holm"], src)
ab = pd.read_csv(R / "pathway_followup/results/ablation_gse294180_contrasts.csv").set_index(["ablation", "contrast"])
src = "research/pathway_followup/results/ablation_gse294180_contrasts.csv"
check("4.2", "MD shift, full model", "1.66", ab.loc[("full_R2", "PE_in_MD"), "estimate"], src)
check("4.2", "MD shift, both cholesterol inputs removed", "0.01", ab.loc[("drop_CHOL_both", "PE_in_MD"), "estimate"], src)
check("4.2", "MD shift, PPARa targets removed", "1.35", ab.loc[("drop_PPARA", "PE_in_MD"), "estimate"], src)

# ---- Table 4: held-out gene-level R2
lp = pd.read_csv(R / "lipid_prediction_followup/results/group_metrics_all.csv").set_index(["diet", "model", "group"])
src = "research/lipid_prediction_followup/results/group_metrics_all.csv"
for d, chol, oth in (("PE_in_CD", ".648", ".052"), ("PE_in_MD", ".600", ".057")):
    check("Table 4", f"Cholesterol synthesis (21 genes), F8 OLS, {d}", chol, lp.loc[(d, "F8-OLS", "CHOL"), "R2_group"], src)
    check("Table 4", f"Other genes, F8 OLS, {d}", oth, lp.loc[(d, "F8-OLS", "OTHER"), "R2_group"], src)
check("Table 4", "Cholesterol-synthesis genes", "21", lp.loc[("PE_in_CD", "F8-OLS", "CHOL"), "n_genes"], src)
check("4.2", "Gene modules spanned by cholesterol genes", "6", lp.loc[("PE_in_CD", "F8-OLS", "CHOL"), "n_modules"], src)
mc = pd.read_csv(R / "roadmap_phase4b/model_comparison_results.csv").set_index(["target", "rep", "model"])
src = "research/roadmap_phase4b/model_comparison_results.csv"
for d, f8, full in (("PE_in_CD", ".066", ".104"), ("PE_in_MD", ".070", ".104")):
    check("Table 4", f"All genes, F8 OLS, {d}", f8, mc.loc[(d, "F8", "ols"), "r2_oos_mean"], src)
    check("Table 4", f"All genes, FULL OLS (28 features), {d}", full, mc.loc[(d, "FULL", "ols"), "r2_oos_mean"], src)

# ---- Table 5: nonlinear vs linear, cholesterol synthesis
nl = pd.read_csv(R / "nonlinear_chol_followup/TABLE_results_concise.csv").set_index(["diet", "model"])
src = "research/nonlinear_chol_followup/TABLE_results_concise.csv"
t5 = {"F8-OLS": (".648", ".600"), "F8-Ridge": (".644", ".598"), "F8-CatBoost": (".659", ".424"), "F8-HGB": (".614", ".481")}
for m, (cd, md) in t5.items():
    check("Table 5", f"{m} R2, CD", cd, nl.loc[("PE_in_CD", m), "R2_group"], src)
    check("Table 5", f"{m} R2, MD", md, nl.loc[("PE_in_MD", m), "R2_group"], src)

# ---- 4.2 interaction (HGB pair-only, lipogenesis, CD)
ix = pd.read_csv(R / "interaction_followup/results/pair_only_vs_additive_heldout.csv").set_index(["diet", "group"])
src = "research/interaction_followup/results/pair_only_vs_additive_heldout.csv"
r = ix.loc[("PE_in_CD", "LIPO")]
check("4.2", "Lipogenesis genes", "12", r.n_genes, src)
check("4.2", "Lipogenesis CD R2, additive", "0.163", r.R2_additive, src)
check("4.2", "Lipogenesis CD R2, with FAST_eff x CIRC_se", "0.292", r.R2_pair_only, src)
check("4.2", "Holm-adjusted p", "0.054", r.p_holm_6, src)

# ---- 4.2 Hong-PE transfer
hg = pd.read_csv(R / "roadmap_phase5/gene_level_external_results.csv")
r = hg[(hg.study == "HONG") & (hg.model == "F8_ols") & (hg.training_target == "PE_in_CD")].iloc[0]
src = "research/roadmap_phase5/gene_level_external_results.csv"
check("4.2", "Hong-PE Spearman rho, frozen F8-OLS", "0.22", r.spearman, src)
check("4.2", "Hong-PE mouse-label permutation p", "0.25", r.perm_p_one_sided, src)

out = pd.DataFrame(rows)
dst = ROOT / "supplementary/tables/report_number_crosswalk.csv"
dst.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(dst, index=False)
print(f"{int(out.match.sum())}/{len(out)} report values match their source files -> {dst.relative_to(ROOT)}")
if not out.match.all():
    print(out[~out.match].to_string(index=False))
    sys.exit(1)
