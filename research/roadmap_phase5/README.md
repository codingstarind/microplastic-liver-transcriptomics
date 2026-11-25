# Phase 5: external validation of frozen models

**Start with `phase5_decision_report.md`.**
- The plan was locked first: `locked_plan.md`, with hashes in `results/plan_lock.json`.
- Earlier phases were not modified (the Phase 4D lock is test-verified).
- Nothing is committed.

## Run (from `scripts/`, Python 3.14 venv; R = WSL mpliver-r)

```powershell
python s01_prepare_external.py
```
```powershell
python s00_lock_plan.py
```
```powershell
python s02_freeze_models.py
```
Then run `Rscript s03_external_de.R <HONG|GSE245069|GSE277772>` (the 3 studies run in parallel in WSL, about 45 s), followed by:
```powershell
python s04_gene_level.py
```
```powershell
python s05_sample_level.py
```
```powershell
python s06_pathways_diagnostics.py
```
```powershell
python s07_figures.py
```
```powershell
python -m pytest ..\tests -q
```

- **Compute:** 20 logical cores, about 4 GB free RAM. Hence 6 loky workers (`P5_JOBS`), with BLAS/OMP pinned to 1 and no nested parallelism.
- **Run times:** gene-level permutation nulls (all relabellings, up to 1,820) take about 3 min; everything else takes under 1 min.

## Outputs
- **`dataset_provenance.csv`:** inputs, independence, shared controls, comparability.
- **`external_contrasts.csv`:** contrasts, methods, genes tested, signal share.
- **`gene_level_external_results.csv`:** Q1, every frozen model × contrast, with the locked verdict.
- **`sample_level_external_results.csv`, `per_sample_scores.csv`:** Q2.
- **`pathway_concordance.csv`:** fixed panel plus Hallmarks.
- **`validation_diagnostics.md`.**
- **`frozen/`:** serialised models and parameters with SHA-256 hashes.
- **`results/`:** external DE tables, camera results, R2 decompositions, coverage control, QC.
- **`figures/` fig1–fig4.**
- **`tests/test_phase5.py`** (10 passing):
  - plan and input lock; mouse uniqueness and shared controls; no overlap between study units;
  - outputs match metadata; contrast direction; gene-set coverage;
  - frozen-model serialisation and reproduction; no external-label leakage; reproducibility;
  - Phase 4D unchanged (skips if that folder is absent).
