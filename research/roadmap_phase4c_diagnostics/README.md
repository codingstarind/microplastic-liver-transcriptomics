# Phase 4C: diagnostic audit of Phase 4 / 4B models

**Start with:** `decision_memo.md`. The audit spec was predeclared and hashed (`predeclared_audit_spec.md`, `tables/spec_lock.json`). All inputs are read-only from `roadmap_phase2`–`roadmap_phase4b`.

## Scripts (`scripts/`, run in this order)

| Script | Task | Runtime |
|---|---|---|
| `c01_collinearity.py` | VIF, max \|r\|, condition numbers per training fold, rank, coefficient stability (F8, A4, B12, C4, F8+C, F8+B+C, FULL) | ~1 min |
| `c03_fix_experiment.py` | Predeclared 4-model comparison on the Phase 4 folds; writes out-of-fold predictions | ~1 min, 19 workers |
| `c02_errors_target_leakage.py` | Residual diagnostics, subgroup errors (module-block bootstrap), target quality and single-mouse influence, confound probes, dependence audit | ~2 min |

`c00_common.py` holds shared helpers.

```powershell
cd research\roadmap_phase4c_diagnostics\scripts; python c01_collinearity.py; python c03_fix_experiment.py; python c02_errors_target_leakage.py
```
```powershell
python -m pytest research\roadmap_phase4c_diagnostics\tests -q
```

## Outputs
- **`tables/`:**
  - collinearity: `collinearity_summary`, `vif_by_feature`, `correlated_pairs_absr_ge_0.5`, `biological_axis_similarity_spearman`, `coefficient_stability_*`;
  - errors: `error_overview`, `subgroup_errors`;
  - target quality: `target_quality`, `single_mouse_influence_on_target`, `loo_model_r2_F8_ols`, `confound_probes`;
  - leakage and dependence: `dependence_audit`, `leakage_audit`;
  - fix experiment: `fix_experiment_*` (by-fold metrics, out-of-fold predictions, coefficients, summary, reproducibility).
- **`figures/`:** `fig1_collinearity.png`, `fig2_residual_diagnostics_F8_ols.png`, `fig3_subgroups_and_influence.png`.

## Tests (`tests/test_phase4c.py`, 7 passing)
- spec lock unchanged;
- fix experiment reproduces Phase 4B exactly;
- fix-experiment scope and alpha grid as predeclared;
- every gene predicted exactly once per repeat;
- collinearity values (F8 VIF < 2, FULL > 10, top pair FAST_eff/A_FAST_shr);
- earlier-phase inputs unchanged (SHA-256);
- post-hoc |effect| bins labelled diagnostic-only.
