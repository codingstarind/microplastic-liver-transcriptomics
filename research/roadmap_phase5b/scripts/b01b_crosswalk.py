"""Writes pipeline_crosswalk.csv (row-by-row training vs external implementation). Evidence numbers come from
results/audit_evidence.json and tables/equivalence_checks.csv produced by b01_pipeline_audit.py."""
import json

import pandas as pd
from b00_common import PB, RES

ev = json.loads((RES / "audit_evidence.json").read_text())
cu = ev["count_unit"]
ht = {k: ev[k] for k in ev if k.startswith("heavy_tail")}
R = [
    ("1 Input assay & origin", "GSE294180 bulk RNA-seq, deposited gene-level matrix (data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz); 19 mice",
     "HONG: 4 GEO readcount files (GSE261893/261725/286116/262399); GSE245069 gene_expression.txt; GSE277772 RAW.tar (normalised)",
     "match (per-study provenance recorded)", "roadmap_phase5/dataset_provenance.csv; s01_prepare_external.py", "none", "n/a"),
    ("1 Biological units / shared controls", "19 distinct mice, 4 cells", "HONG 24 GSMs = 15 mice (3 CON byte-identical in 4 series; identical SRA read counts); GSE245069 16 mice, 4 shared CON",
     "match (handled as one study each)", "roadmap_phase5/results/external_sample_manifest.csv; tests test_biological_sample_uniqueness", "would inflate n if ignored", "n/a"),
    ("2 Count type", "Salmon/tximport countsFromAbundance='no', rounded (estimated fragment counts)",
     "HONG StringTie prepDE (coverage-derived) integers; GSE245069 estimated read_count (7.5% non-integer); GSE277772 normalised",
     "PARTIAL MISMATCH (count unit)", f"counts per read pair: GSE294180 {cu['GSE294180']['median']}, HONG {cu['HONG']['median']}, GSE245069 {cu['GSE245069']['median']} (tables/count_unit_per_sample.csv)",
     "HONG counts ~2x fragment-equivalent: Poisson part of DESeq2 variance understated -> lfcSE too small for low counts, inflated DEG numbers, precision weights over-trust low-count genes; point log2FC almost unaffected",
     "YES: rescale by an a-priori unit factor (0.5, both mates counted) - no outcome needed"),
    ("2 Annotation / identifiers", "GENCODE vM35 gene symbols (tximport by symbol)", "HONG RefSeq GRCm39 symbols, GSE245069 NCBI GRCm38 symbols -> Phase 2 mapper (exact, case, MGI synonym)",
     "match after mapping", f"HONG usable {ev['HONG_mapping']['usable']}/{ev['HONG_mapping']['rows']}; many-to-one rows {ev['HONG_mapping']['many_to_one_rows']} ({ev['HONG_mapping']['many_to_one_symbols_in_feature_table']} symbols in feature table); feature-table coverage 99.99% / 99.74%",
     "<=0.3% of evaluated genes summed from >1 row", "sensitivity only (drop many-to-one genes)"),
    ("2 Gene filtering", ">=10 counts in >=4 samples (smallest group), 19-mouse matrix", ">=10 counts in >= smallest group (3 HONG, 4 GSE245069)",
     "match in rule; units differ (see count unit)", "analysis/deseq2/run_deseq2_factorial.R L79-80; roadmap_phase5/scripts/s03_external_de.R",
     "HONG unit inflation lets lower-expressed genes pass", "YES with unit correction"),
    ("2 Gene universe / order", "Phase 2 main universe n reference genes = 11,758 feature rows (raw_table order)", "same frozen predictions, intersected with genes tested",
     "match", "equivalence_checks: gene order identical; re-prediction identical", "none", "n/a"),
    ("3 Normalisation for gene targets", "DESeq2 median-of-ratios size factors; Wald MLE log2FC (no length offset)", "DESeq2 median-of-ratios; Wald MLE log2FC",
     "match in method", "same DESeq2 version in mpliver-r env", "none", "n/a"),
    ("3 Target estimand", "PE_in_CD from ~diet*microplastic (dispersion from 19 mice)", "PE vs CON from ~group with 5 groups (dispersion from 15 mice)",
     "match in definition (exposed minus control, simple effect)", "equivalence_checks sign tests (Spearman 0.77-0.98)", "none", "n/a"),
    ("3 Target robustness (estimator behaviour)", "MLE on arithmetic means; single-mouse spikes dominate",
     "same estimator, but external targets are not spike-dominated",
     "SAME function, DIFFERENT behaviour", f"top-1% share of sum sq (bm>500): {ht['heavy_tail_GSE294180_PE_in_CD']['top1pct_share_of_sum_sq']:.2f} vs HONG {ht['heavy_tail_HONG_PE_vs_CON']['top1pct_share_of_sum_sq']:.2f}; r(MLE, mean-of-logs) {ht['heavy_tail_GSE294180_PE_in_CD']['r_mle_vs_mean_of_logs']:.2f} vs {ht['heavy_tail_HONG_PE_vs_CON']['r_mle_vs_mean_of_logs']:.2f}",
     "training target dominated by outlier genes (muscle/blood carry-over in single mice); model fits partly to noise that cannot transfer", "YES as a TRAINING-target change (robust estimator) - post-external exploratory"),
    ("4 F8 features", "MASH/FAST/PPARA/CIRC reference eff + se per gene (Phase 3 reference_effects)", "identical rows (features are reference-derived, study-independent)",
     "match", "equivalence_checks: feature columns & order, re-prediction identical", "none", "n/a"),
    ("4 Engineered features (block C)", "C_FASTxPPARA, C_MASHxFAST, C_FAST_morning_dev, C_PPARA_consistent; z-scaling fitted on training genes", "same fitted transformer",
     "match", "frozen joblib F8C_ridge; feature-name check", "none", "n/a"),
    ("4 Orientation / units of features and target", "log2(exposed/control)", "log2(exposed/control)", "match", "sign tests", "none", "n/a"),
    ("5 Train-dependent steps", "tau2 (DL), StandardScaler, RepBuilder, Ridge alpha: fitted on GSE294180 training genes only", "no refit; external only predicted",
     "match (no leakage)", "Phase 5 freeze reads no external data; no .fit in s04/s05; Phase 4C leakage_audit", "none", "n/a"),
    ("5 Evaluation weights", "w = 1/(s^2 + tau2) with GSE294180 SEs", "w = 1/(s_ext^2 + tau2_ext) with external SEs",
     "match in formula; inherits count-unit issue", "s04_gene_level.py", "HONG weights over-trust low-count genes", "YES with unit correction"),
    ("6 External target uncertainty", "4-5 mice per arm; SE from DESeq2", "3 v 3 (HONG): SE from DESeq2 only",
     "INADEQUATE in Phase 5 (no mouse resampling)", "Phase 5B b02 LOO / 2v2 DESeq2 refits", "external ground truth is itself noisy; transfer metrics unstable", "YES (diagnostic, no outcome-driven change)"),
    ("7 Sample scores: universe & ranks", "4D: 13,777-gene background, TPM-proxy (vM35 median length); 4E: 11,883-gene logCPM ranks", "background intersected with external genes; same functions",
     "match (coverage >= 85% 4D / 99.7% 4E)", "roadmap_phase5/results/score_gene_coverage.csv; coverage control <= 9%", "small", "n/a"),
    ("7 Sample scores: signature direction", "frozen up/down lists", "same lists, hash-verified", "match", "equivalence_checks: 4D signature hashes", "none", "n/a"),
    ("7 Sample scores: count unit", "ranks within sample", "ranks within sample (unit-invariant)", "match (ranks invariant to the 2x unit)", "singscore definition (Foroutan 2018)", "none", "n/a"),
]
C = pd.DataFrame(R, columns=["component", "training_implementation", "external_implementation", "match", "evidence", "expected_impact", "correction_permissible_without_PE_outcomes"])
C.to_csv(PB / "pipeline_crosswalk.csv", index=False)
print(C[["component", "match"]].to_string(index=False))
