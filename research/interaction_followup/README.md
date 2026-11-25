# Interaction follow-up: Phase 4 CatBoost and HGB (exploratory; 2026-11-25)

This follow-up checks whether interactions among the fasting, PPARα, MASH and circadian reference-signature features add to held-out prediction in the Phase 4 F8 boosting models. Explanations are computed on training folds; benefit is tested on held-out genes. Start with `REPORT.md`.

## Files

| Path | Content |
|---|---|
| `figures/fig_interaction_FAST_CIRCse.{png,pdf,svg}` | A–B: HGB pure-interaction surface for the only reproducible pair (FAST_eff × CIRC_se), CD and MD. C: per-split interaction strength vs the interaction-free null for the top pairs. D: held-out ΔR² vs the additive model for all genes, lipogenesis and PPARα targets. |
| `TABLE_results_concise.csv` | A: top pairs vs null. B: held-out full vs additive. C: held-out pair-only vs additive. |
| `REPORT.md` | Numbers, decision rules, limits |
| `results/friedman_H_training_folds.csv.gz` | Observed and null H² and interaction SD, 28 pairs × splits |
| `results/pair_by_split_observed_vs_null.csv`, `results/pair_summary_CD_MD.csv` | Pair summaries, per split and over splits |
| `results/fold_stability_kendallW.csv`, `results/repeatability_CD_MD_and_models.csv` | Stability of pair rankings |
| `results/catboost_shap_interactions_training_folds.csv.gz`, `results/catboost_shap_interaction_summary.csv` | CatBoost SHAP interaction values |
| `results/heldout_full_vs_additive_predictions.csv.gz`, `results/heldout_full_vs_additive_by_group.csv` | Held-out full vs additive |
| `results/pair_only_heldout_predictions.csv.gz`, `results/pair_only_vs_additive_heldout.csv` | Held-out pair-only vs additive |
| `results/pair_interaction_surface_by_split.csv.gz`, `results/pair_surface_stability.csv` | Interaction surfaces and their stability |
| `results/refit_check_vs_saved_phase4.csv`, `results/run_info_*.json`, `results/fit_info.csv`, `logs/` | Reproduction and run records |

## Reproduce

Run these from the repository root. ix01 needs CatBoost, which is installed in `C:\Python314\python.exe` (same scikit-learn and numpy versions as `.venv`).

```powershell
C:\Python314\python.exe research\interaction_followup\scripts\ix01_fit_and_diagnose.py
```
```powershell
.venv\Scripts\python.exe research\interaction_followup\scripts\ix02_analyse.py
```
```powershell
.venv\Scripts\python.exe research\interaction_followup\scripts\ix03_pair_followup.py
```
```powershell
.venv\Scripts\python.exe research\interaction_followup\scripts\ix04_pair_tests_and_figure.py
```

Approximate runtimes: 11 min (14 workers), 10 s, 30 s and 10 s. All steps are deterministic (seed 20261125).
