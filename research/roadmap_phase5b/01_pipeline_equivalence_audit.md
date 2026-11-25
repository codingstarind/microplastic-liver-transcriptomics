# 01: Pipeline equivalence audit (Stage 1)

**Files:**
- `pipeline_crosswalk.csv`: 18 components, with training vs external implementation, match, evidence, impact and whether a correction is permissible.
- `tables/equivalence_checks.csv`: 16 automated checks.
- `tables/count_unit_per_sample.csv`.
- `results/audit_evidence.json`.

**Scripts:** `b01_pipeline_audit.py`, `b01b_crosswalk.py`.

## Where each step lives (verified on disk)

| Step | Training (GSE294180) | External |
|---|---|---|
| Counts | `data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz`: Salmon/tximport "no" counts, rounded. The `reports/quantification_validation` pilot confirmed this | `D:/mpliver_work/phase5/<study>_expression.csv.gz`, built by `roadmap_phase5/scripts/s01_prepare_external.py` from the GEO files in `D:/mpliver_work/phase2_refs/` |
| Gene target | `analysis/deseq2/run_deseq2_factorial.R` (prefilter L79–80; contrasts L95–107) → `results/deseq2/provisional_deposited/results_PE_in_CD.csv.gz` | `roadmap_phase5/scripts/s03_external_de.R` |
| Features | `roadmap_phase4/scripts/p4_data.py::feature_table` + `roadmap_phase4b/scripts/p4b_features.py::RepBuilder` | the same objects, frozen in `roadmap_phase5/frozen/*.joblib` |
| CV | `p4_data.outer_splits` (module groups; seeds 101–103) + inner GroupKFold(4) (`roadmap_phase4c_diagnostics/scripts/c03_fix_experiment.py`) | — |
| External scoring | — | `roadmap_phase5/scripts/s04_gene_level.py`, `s05_sample_level.py` |

## Result: what matches

All of the following hold:
- **Gene order and features:**
  - the frozen-prediction gene order equals the `raw_table` order;
  - serialised models re-predict identically from `RAW_COLS` in the locked column order;
  - a reversed-column guard changes predictions, which proves order is enforced;
  - the 8 F8 and 4 block-C feature names are as locked.
- **Features are identical for every study:** they are computed from independent reference studies per gene, not from GSE294180 or the external data.
- **Orientation:** every target is exposed − control. Spearman with the log ratio of group-mean CPM is 0.77 (GSE294180), 0.96 (Hong) and 0.98 (GSE245069).
- **Same estimator:** DESeq2 Wald MLE, median-of-ratios, the same prefilter rule (≥ 10 counts in ≥ the smallest group) and the same environment.
- **Gene mapping:**
  - 99.99% (Hong) and 99.74% (GSE245069) of the 11,758 feature-table genes are present;
  - 29 and 13 of those genes were summed from more than one source row (≤ 0.3%).
- **No leakage:**
  - the Phase 5 freeze step reads no external data;
  - the external scoring scripts contain no `.fit`;
  - the Phase 4E models reproduce Phase 4E scores to 4e-16;
  - the Phase 4D signature hashes are unchanged.
- **Sample scores:** within-sample ranks are invariant to count units. Background coverage is ≥ 85% (4D) and ≥ 99.7% (4E), and the coverage control showed ≤ 9% change.

## Result: verified mismatches and estimator issues

| # | Issue | Evidence | Impact (measured) | Correctable without PE outcomes? |
|---|---|---|---|---|
| M1 | **Count unit:** Hong counts are coverage-derived from both mates, about 2× fragment-equivalent | counts per read pair: GSE294180 0.775, GSE245069 0.734, **Hong 1.641** (2×101 bp) | Prefilter admits 2,180 more low-count genes; DEGs 1,133 → 951 after correction; median lfcSE ratio 1.006; log2FC r = 0.96. **Rank transfer unchanged: F8-OLS ρ 0.219 → 0.220** (`tables/i1_unit_harmonisation_effect.csv`) | **Yes:** a priori factor 0.5 (done: I1) |
| M2 | **Training target is outlier-dominated.** Same estimator, different behaviour: the DESeq2 MLE uses arithmetic group means, and single-mouse spikes (muscle Acta1/Myh1, haemoglobin, Gm28439 = 23 log2) dominate | Top-1% share of the sum of squares (baseMean > 500): **0.77** vs Hong 0.33 and GSE245069 0.43. r(MLE, mean-of-logs) 0.63 vs 0.96 / 0.98. Excess kurtosis 12.4 | Precision weights down-weight most of these (large lfcSE). The model still partly fits non-transferable noise | As a **training-target** change only (I2a). Post-external exploratory |
| M3 | **External-target uncertainty was not propagated in Phase 5** (DESeq2 SEs only, 3 v 3) | 15 mouse-resampled DESeq2 refits: target ρ with the full fit 0.87–0.96; **transfer ρ 0.10–0.31** | The metric itself moves by ±0.1 when one or two mice change | Yes (diagnostic; done in Stage 2/3) |

## What was **not** found
- No sign flip.
- No feature-column misorder.
- No identifier shift that would scramble gene identity.
- No scaler, imputer, filter or threshold fitted on external labels.
- No normalised-value-as-count misuse in the gene-level tests. GSE277772 was kept exploratory and rank-only.

**The Phase 5 gene-level result was not caused by a coding or pipeline error.** The one verified technical mismatch (M1) is real, and correcting it changes SE-dependent quantities, not rank transfer.

## Equivalence tests
`tests/test_phase5b.py` covers:
- normalisation and count-unit equivalence;
- sample-ID overlap and shared-control grouping;
- feature column order;
- train/test leakage;
- determinism.
