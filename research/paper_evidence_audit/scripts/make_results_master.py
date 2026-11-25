"""results_master.csv = (a) numbers extracted programmatically from result CSVs (extract_tables.py) +
(b) hand-transcribed numbers from phase reports where no single machine-readable table exists (source_type = 'report').
Run extract_tables.py first."""
from pathlib import Path

import pandas as pd

AUD = Path(__file__).resolve().parents[1]
X = pd.read_csv(AUD / "results_master_extracted.csv")
X["source_type"] = "result table (extracted)"
M = [
    # phase, analysis, quantity, value, ci_low, ci_high, n, context, source
    ("P0", "DESeq2 factorial", "genes FDR<0.05 PE_in_CD", 858, None, None, "19 mice", "513 with |log2FC|>=0.38", "reports/final_report_draft.md §4.3"),
    ("P0", "DESeq2 factorial", "genes FDR<0.05 PE_in_MD", 472, None, None, "19 mice", "", "reports/final_report_draft.md §4.3"),
    ("P0", "DESeq2 factorial", "genes FDR<0.05 Diet_in_Veh", 2262, None, None, "19 mice", "LOO 1,699-2,537", "reports/final_report_draft.md §4.3"),
    ("P0", "DESeq2 factorial", "genes FDR<0.05 Interaction", 87, None, None, "19 mice", "23 muscle-tracking; 16 with muscle covariate; 49 without MD_1; 1 in edgeR/limma", "reports/final_report_draft.md §4.6"),
    ("P0", "Depth/yield adjustment", "PE_in_CD genes after depth / yield covariate", "116 / 66", None, None, "19 mice", "from 858", "research/roadmap_phase1/existing_work_inventory.md §3"),
    ("P0", "Library depth", "read pairs PE vs vehicle (M)", "13.5-16.5 vs 17.0-19.1", None, None, "19", "complete separation; one flowcell", "reports/final_report_draft.md §4.1"),
    ("P0", "Muscle carry-over", "muscle marker CPM MD_1 / CD_MP_1 / MD_2 / MD_3", "5190 / 539 / 316 / 161", None, None, "19", "PC1 30% = muscle", "reports/final_report_draft.md §4.1"),
    ("P0", "Pathway (fry)", "Cholesterol homeostasis PE_in_CD / PE_in_MD p", "0.004 / 0.005", None, None, "19", "FDR 0.11 / 0.09; interaction p 0.84", "reports/phase2/pathway_report.md"),
    ("P0", "Pathway (fry)", "Sanderson PPARa targets PE_in_CD FDR", 0.041, None, None, "19", "up; MD p 0.35", "reports/phase2/pathway_report.md"),
    ("P0", "Split-reference concordance", "rho PE_in_CD vs Diet_in_Veh (median)", 0.33, 0.29, 0.40, "6 splits", "naive 0.54; p 0.048-0.095 (min 0.048)", "reports/final_report_draft.md §4.5"),
    ("P0", "Split-reference concordance", "rho PE_in_MD vs Diet_in_Veh", 0.05, -0.03, 0.17, "20 splits", "", "reports/final_report_draft.md §4.5"),
    ("P0", "Clock phase score", "group mean CD+Veh vs others", "-2.38 vs +0.08 to +1.01", None, None, "19", "all CD+Veh below every other mouse", "reports/final_report_draft.md §4.7"),
    ("P0/P3", "Ppara discrepancy", "PE_in_CD Ppara log2FC ours vs Jung", "+0.16 to +0.17 vs +0.53", None, None, "19", "unexplained by DESeq2 settings", "research/roadmap_phase3/results/ppara_discrepancy/ppara_estimates_by_approach.csv"),
    ("P2", "Muscle exclusion", "genes in rule / in filtered set", "309 / 187", None, None, "", "", "research/roadmap_phase2/reports/phase2_signoff.md §4"),
    ("P2", "Gene universe", "main / strict", "11,888 / 8,360", None, None, "", "", "research/roadmap_phase2/reports/phase2_signoff.md §4"),
    ("P2", "Tests", "phase 2 data tests passed", 60, None, None, "", "", "research/roadmap_phase2/reports/phase2_signoff.md §1"),
    ("P3", "MASH positive control", "rho Diet_in_Veh vs MASH1", 0.355, None, None, "11,886 genes", "c1 FAIL: MC q99 0.382; exact p 0.036, q99 0.390", "research/roadmap_phase3/reports/mash_reference_report.md; reports/phase3_review_addendum.md"),
    ("P3", "MASH positive control", "sign agreement (both padj<0.05)", 0.787, None, None, "1,247 genes", "c2 met", "research/roadmap_phase3/reports/mash_reference_report.md"),
    ("P3", "MASH positive control", "rho strict universe", 0.252, None, None, "8,360 genes", "below 0.30 bar", "research/roadmap_phase3/reports/phase3_review_addendum.md §3"),
    ("P3", "Exact permutation audit", "balanced 4+4 maximum rho", 0.288, None, None, "4,900 allocations", "earlier sampled value 0.277 corrected", "research/roadmap_phase3/results/posthoc_audit/exact_permutation_by_overlap.csv"),
    ("P3", "MASH reproducibility", "rho MASH1 vs MASH2", 0.557, None, None, "", "90% sign agreement among 3,306 shared DEGs", "research/roadmap_phase3/reports/phase3_decision_report.md §3"),
    ("P3", "Fasting reference", "clean subset vs main rho", 0.919, None, None, "4 v 6", "45% of DEGs ZT-dependent", "research/roadmap_phase3/reports/fasting_reference_report.md"),
    ("P3", "Reference programmes", "median log2FC lipogenesis fasting / GW7647", "-2.79 / +1.15", None, None, "", "cholesterol synthesis -2.76 / +0.67; PPARa targets +0.82 / +1.50", "research/roadmap_phase3/reports/reference_comparison_report.md"),
    ("P3", "Collinearity", "VIF MASH1 / FAST / PPARA / CIRC_T1", "1.03 / 1.07 / 1.02 / 4.07", None, None, "11,749 genes", "", "research/roadmap_phase3/reports/reference_comparison_report.md"),
    ("P3", "Baseline 2", "rho PE_in_CD vs MASH1 [LOO]", 0.211, 0.03, 0.24, "11,888 genes", "range is LOO, not CI", "research/roadmap_phase3/results/comparisons/baselines.csv"),
    ("P3", "Baseline 2", "rho PE_in_MD vs FAST [LOO]", 0.254, 0.22, 0.27, "11,888 genes", "range is LOO, not CI", "research/roadmap_phase3/results/comparisons/baselines.csv"),
    ("P3", "Baseline 1", "PE_in_CD vs MASH2 concordant/discordant", "84 / 125", None, None, "209 overlap", "hypergeometric still 'significant'", "research/roadmap_phase3/reports/baseline_comparison_report.md"),
    ("P4C", "Target quality", "signal share DESeq2 CD / MD", "0.36 / 0.39", None, None, "", "jackknife 0.05 / 0.12", "research/roadmap_phase4c_diagnostics/tables/target_quality.csv"),
    ("P4C", "Residuals", "Spearman |resid| vs SE CD / MD", "0.40 / 0.44", None, None, "", "module-block bootstrap CIs [0.36,0.43]/[0.39,0.49]", "research/roadmap_phase4c_diagnostics/decision_memo.md"),
    ("P4C", "Subgroups", "mean residual LIPOGENESIS_CORE CD / MD", "-1.24 / -0.74", None, None, "", "CHOL_SYNTHESIS_CORE -0.42 / -0.53", "research/roadmap_phase4c_diagnostics/tables/subgroup_errors.csv"),
    ("P5B", "Count unit", "counts per read pair GSE294180 / GSE245069 / Hong", "0.775 / 0.734 / 1.641", None, None, "", "", "research/roadmap_phase5b/tables/count_unit_per_sample.csv"),
    ("P5B", "Target tails", "top-1% share of squared effect GSE294180 / Hong / GSE245069", "0.77 / 0.33 / 0.43", None, None, "", "r(MLE, mean-of-logs) 0.63 / 0.96 / 0.98", "research/roadmap_phase5b/01_pipeline_equivalence_audit.md"),
    ("P5B", "Labelling diagnostic", "R2 of 20-labelling rho on mouse PC1-3", 0.889, None, None, "20 labellings", "true labelling 5th of 20 (log2CPM proxy); 4th in intervention table (arith proxy)", "research/roadmap_phase5b/tables/gap_hong_20_labellings.csv; gap_labelling_summary.txt"),
    ("P5B", "Reliability", "robust jackknife reliability Hong / GSE294180", "0.66 / 0.55", None, None, "", "disattenuated between-study rho ~0.17", "research/roadmap_phase5b/tables/gap_reliability_and_disattenuation.csv"),
    ("P5B", "Baseline shift", "control FAST singscore GSE294180 / Hong / GSE245069", "-0.112 / +0.101 / +0.044", None, None, "", "", "research/roadmap_phase5b/tables/gap_study_shift.csv"),
    ("P5B", "Target instability", "Hong transfer rho across 15 mouse-resampled targets", "0.10-0.31", None, None, "15 refits", "MD-trained 0.14-0.39", "research/roadmap_phase5b/tables/gap_hong_target_resampling.csv"),
]
Mdf = pd.DataFrame(M, columns=["phase", "analysis", "quantity", "value", "ci_low", "ci_high", "n", "context", "source"])
Mdf["source_type"] = "phase report / summary table (hand-transcribed; re-check before citing)"
X["phase"] = "P" + X.phase.astype(str)
out = pd.concat([Mdf, X], ignore_index=True)
out.insert(0, "id", [f"R{i:04d}" for i in range(1, len(out) + 1)])
out.to_csv(AUD / "results_master.csv", index=False)
print(len(out))
