# Phase 5B: cross-study transfer audit, diagnosis and targeted rebuild

**Start with `05_results_and_generalization_verdict.md`, then `phase6_readiness_memo.md`.**

**Evidence status:**
- Phase 5 external outcomes were already seen, so every 5B comparison on Hong or GSE245069 is **post-external exploratory re-development**, not fresh validation.
- No untouched independent study exists.
- Phases 1–5 were read only; their lock tests pass. Nothing is committed.

## Execution order (from `scripts/`)

| # | Command | What it does | Time |
|---|---|---|---|
| 1 | `python b01_pipeline_audit.py` | Stage 1 evidence and 16 equivalence checks | < 1 min |
| 2 | `python b01b_crosswalk.py` | writes `pipeline_crosswalk.csv` | s |
| 3 | `Rscript b02_external_target_resampling.R` (WSL) | 16 DESeq2 refits of the Hong PE target (leave-one-out and 2-vs-2) | ~1 min |
| 4 | `python b03_gap_diagnostics.py` | Stage 2 tables and fig1 | ~1 min |
| 5 | `python b00_lock.py` | **hash-locks 03/04 and the intervention code before scoring** | s |
| 6 | `Rscript b04_intervention_targets.R` (WSL) | I1 unit-harmonised Hong targets; I2b diet-averaged PE contrast | ~1 min |
| 7 | `python b05_interventions.py 1000` (`P5B_JOBS=10`) | nested grouped CV (7 candidates × 15 outer folds), final refits, external evaluation, 20 labellings, paired mouse bootstrap (1,000), MASH-classifier bootstrap | ~25 min |
| 8 | `python b06_determinism.py` | reruns b05 with B = 50, compares outputs, restores the full run | ~5 min |
| 9 | `python b07_figures.py` | figs 2–4 | s |
| 10 | `python -m pytest ..\tests -q` | Phase 5B tests | s |

WSL Rscript path: `/root/mpliver-tools/mamba/envs/mpliver-r/bin/Rscript` (DESeq2, limma, BiocParallel).

## Inputs
- **GSE294180:** `data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz`, `results/deseq2/provisional_deposited/` (dds.rds, results_*.csv.gz) and the LOO refits in `results/deseq2/loo/`.
- **External:** `D:/mpliver_work/phase5/{HONG,GSE245069}_expression.csv.gz` (from Phase 5 `s01`) and `research/roadmap_phase5/results/external_de/`.
- **Features and CV:** `research/roadmap_phase3/reference_effects/`, `roadmap_phase4/scripts/p4_data.py`, `roadmap_phase4b/scripts/p4b_features.py`, `roadmap_phase4/results/gene_groups.csv`.
- **Frozen models:** `research/roadmap_phase5/frozen/`.

## Software and hardware
- **Software:**
  - Python 3.14.6, scikit-learn 1.9.1, numpy 2.5.1, pandas 3.0.5, scipy 1.18.1, joblib 1.6.0;
  - R (mpliver-r env) with DESeq2.
  - CatBoost is **not installed** in this environment, so HGB is the ensemble comparator.
- **Hardware and workers:**
  - 14 physical / 20 logical cores.
  - Windows had 4–8 GB free RAM, so b05 used 10 loky workers (about 300 MB each). The locked protocol text says 6. Worker count does not affect results: each bootstrap replicate is seeded `20261122 + i`, and the determinism test verifies this.
  - BLAS/OMP threads are pinned to 1; estimators are single-threaded (no nested parallelism).
  - WSL R used 8 (b02) and 3 (b04) BiocParallel workers, with 14 GB free.
- **Seeds:** 20261122 throughout; Phase 4 outer seeds 101–103.

## Deliverables
- **Documents:** `01_pipeline_equivalence_audit.md`, `pipeline_crosswalk.csv`, `02_domain_shift_and_error_diagnosis.md`, `03_resampling_and_nested_cv_protocol.md` (locked), `04_targeted_interventions.md` (locked), `05_results_and_generalization_verdict.md`, `phase6_readiness_memo.md`.
- **`tables/`:**
  - equivalence checks;
  - count units;
  - gap / strata / reliability / PCA / labellings / study shift;
  - CV fold metrics, hyperparameter choices, external candidates, labelling ranks;
  - bootstrap replicates and summary;
  - intervention verdicts;
  - MASH-classifier bootstrap;
  - determinism check;
  - coefficients.
- **`results/`:** audit evidence, Hong resampling, intervention targets, candidate models (joblib, hashed in `b05_run_info.json`), logs, `protocol_lock.json`.
- **`figures/`:**
  - fig1: gap, targets and labellings;
  - fig2: domain shift and uncertainty;
  - fig3: before/after;
  - fig4: coefficients.
- **`tests/test_phase5b.py`.**

## Methods literature (bibliographic details verified via Europe PMC / publisher, 2026-11-23)
- Bernau et al. 2014, *Bioinformatics* 30:i105 (cross-study validation ranks algorithms differently from CV).
- Varma & Simon 2006, *BMC Bioinformatics* 7:91 (nested CV to avoid optimistic tuned-CV error).
- Sweeney et al. 2017, *NAR* 45:e1 (multi-dataset effect-size reproducibility).
- Foroutan et al. 2018, *BMC Bioinformatics* 19:404 (singscore: single-sample, rank-based).
- Debray et al. 2015, *J Clin Epidemiol* 68:279 (relatedness of the validation sample: reproducibility vs transportability).
- Steyerberg & Harrell 2016, *J Clin Epidemiol* 69:245 (internal, internal–external and external validation).
- Riley et al. 2016, *BMJ* 353:i3140 (heterogeneity of performance across settings).

For Steyerberg and Riley, only the bibliographic record was verified; the abstracts were not available.

Not re-verified here: Huber 1964 (robust M-estimation); DerSimonian & Laird 1986 (τ² estimator).
