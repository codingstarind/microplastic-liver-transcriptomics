# Phase 4E: sample-level MASH-state models on a new representation

**Start with:** `decision_memo.md`. Phases 1–4D are read-only and unchanged (Phase 4D lock test-verified). Nothing is committed.

## Documents
- `00_task_and_estimands.md`: the three tasks (gene-level regression, sample-level state prediction, sample-level PE response), their units and capacity; includes a correction about Phase 4D permutation multiplicity.
- `01_model_feasibility_and_leakage_audit.md`: includes the GO/NO-GO matrix and the finding that the Phase 4D MASH signature used both MASH studies, so it is rebuilt inside each training study.
- `02_literature_and_representation_choices.md`.
- `03_predeclared_model_evaluation_plan.md`: locked (`results/plan_lock.json`).
- `04_model_comparison_results.md`, `05_failures_limits_and_sensitivity.md`, `06_phase5_readiness.md`.

## Run (from `scripts/`)

```powershell
$env:P4E_JOBS=12; python e01_mash_state_models.py
```
```powershell
python e01b_permutation_null.py
```
```powershell
python e01c_validity.py
```
```powershell
python e02_fasting_demo.py
```
```powershell
python e03_apply_to_gse294180.py
```
```powershell
python e04_diagnostics_figures.py
```
```powershell
python -m pytest ..\tests -q
```

What each script does:
1. `e01_mash_state_models.py`: 87 nested fits (4 representations × up to 8 models × 3 splits); about 4.3 min on 12 workers.
2. `e01b_permutation_null.py`: 400 permutation fits, 48 s.
3. `e01c_validity.py`: applies the predeclared validity rule.
4. `e02_fasting_demo.py`: fasting demonstration.
5. `e03_apply_to_gse294180.py`: applies frozen models only if valid.
6. `e04_diagnostics_figures.py`: feature contributions and figures.

## Compute
- 14 physical / 20 logical cores; 31.7 GB RAM with about 3–4 GB free. Hence 12 loky workers (about 190 MB each), with BLAS/OMP threads pinned to 1 and estimators single-threaded (no nested parallelism).
- Per-fold signature/PCA results are memoised, so the GridSearch refits reuse them.
- R1 pathway scores are precomputed once, because they are data-independent.
- Heavy models are scheduled first.

## Versions
Python 3.14, scikit-learn 1.9.1, numpy 2.5.1, pandas 3.0.5, scipy 1.18.1.

## Provenance
- **Expression:** `D:/mpliver_work/phase2_refs/<GSE>/prepared/expression_ready.csv.gz` (GSE316023, GSE246328, GSE107787, GSE225695) and the repo's `data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz`.
- **Samples:** Phase 2 `phase3_eligible_samples.csv`.
- **Gene sets:** Phase 3 `results/gene_sets.csv`.
- **Universe:** Phase 2 main universe ∩ datasets (11,883 genes for MASH).

## Outputs
- **`tables/`:** `sample_manifest.csv`, `feature_definitions.csv`, `train_test_study_splits.csv`.
- **`results/`:**
  - `external_metrics.csv`, `external_predictions.csv`;
  - `validity.csv`;
  - `permutation_null_auc_R2_logistic.csv`, `permutation_summary.json`;
  - `fasting_demo_*`;
  - `gse294180_frozen_model_scores.csv`, `_contrasts.csv`, `_loo.csv`, `_technical.csv`;
  - `r2_feature_contributions_gse294180.csv`;
  - `application_decision.json`;
  - run logs.
- **`figures/`:** `fig1_cross_study_performance`, `fig2_before_after_calibration`, `fig3_permutation_null`, `fig4_gse294180_2x2_frozen_scores`.
- **`tests/test_phase4e.py`** (8 passing): plan lock; disjoint studies and sizes; manifest unique with no target labels; feature-definition provenance; training-only signature code; metric completeness and ranges; application only for valid models with correct contrast direction; Phase 4D unchanged.
