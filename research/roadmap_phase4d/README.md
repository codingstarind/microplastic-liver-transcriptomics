# Phase 4D: per-mouse signature scores (different estimand from gene-level R²)

**Start with:** `decision_report.md`. Earlier phases are read-only; nothing is committed.

## Order and commands

From `research/roadmap_phase4d/scripts/`:

```powershell
python s01_build_signatures.py
```
```powershell
python s02_score_and_test.py
```
```powershell
python -m pytest ..\tests -q
```

What each step does:
1. **`s01_build_signatures.py`:** builds the reference-only signatures and writes `frozen_signature_manifest.csv` with SHA-256 hashes.
2. **`analysis_plan.md` + `results/plan_lock.json`:** locked after step 1 and before step 3.
3. **`s02_score_and_test.py`:** singscore, 2×2 contrasts, LOO, exact permutations, 1,000 matched null sets per signature, QC associations, sensitivities and figures.
   - Inputs, gene lengths and rank matrices are cached in `results/cache/`.
   - Runtime 34 s.

## Compute choices
- 14 physical / 20 logical cores, 31.7 GB RAM (3.5 GB free at run time).
- **No multiprocessing was used.** The saturated 2×2 contrasts are exact group-mean differences, so all permutations (31,752) and null sets (3,000) run as vectorised numpy on cached rank matrices in seconds. Worker processes would only have added overhead and memory pressure.
- BLAS threads were pinned to 1.

## Versions and assumptions
- Python 3.14, numpy 2.5.1, pandas 3.0.5, scipy.
- singscore implemented from Foroutan et al. 2018 (BMC Bioinformatics, Methods; verified against the full text). The R `singscore` and `GSVA` packages are not installed, so **GSVA was not run**.
- **TPM-proxy:** counts ÷ GENCODE vM35 median transcript length (no per-sample length offset).

## Files
- **Documents:** `source_inventory.md`, `representation_go_nogo.md`, `analysis_plan.md`, `frozen_signature_manifest.csv`, `decision_report.md`, `phase5_readiness.md`.
- **`signatures/`:** frozen lists and full candidate tables.
- **`results/`:**
  - `per_mouse_scores.csv`;
  - `contrast_estimates_primary.csv` and `contrast_estimates_sensitivity.csv`;
  - `loo_contrasts.csv` and `loo_summary.csv`;
  - `permutation_exact.csv`;
  - `null_gene_sets.csv`;
  - `qc_associations.csv`;
  - `run_summary.json`;
  - `plan_lock.json`.
- **`figures/`:** `fig1_scores_by_group`, `fig2_scores_vs_technical_covariates`, `fig3_contrasts_vs_null_sets`, `fig4_sensitivity`.
- **`tests/test_phase4d.py`** (7 passing):
  - plan and manifest locks;
  - signature hashes;
  - no PE data in signature construction;
  - signatures inside the fixed universe;
  - the singscore formula at its extremes (±1) and determinism;
  - contrast direction equals group-mean differences, with sample alignment;
  - outputs complete.
