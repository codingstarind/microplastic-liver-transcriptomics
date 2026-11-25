# Phase 4B: literature-guided feature engineering and model ablation (EXPLORATORY)

**Start with:** `decision_memo.md`. Phase 4 outputs are not modified; its splits, targets, grids and data layer are imported read-only from `../roadmap_phase4/scripts/`.

## Files

| File | Content |
|---|---|
| `literature_feature_model_matrix.csv` | Paper → representation → model → task → applicability; verification source and location for each |
| `feature_engineering_plan.md` | The 3 feature blocks (A reliability shrinkage, B pathway membership, C response shape): formulas, sources, leakage checks |
| `comparison_protocol.md` | Fold reuse, inner tuning, decision and shortlist rules, exploratory status |
| `results/protocol_lock.json` | SHA-256 of the plan, protocol, matrix and target-free diagnostics, locked before PE scoring |
| `results/input_checksums.sha256` | Checksums of every input file reused |
| `model_comparison_results.csv` | Every representation × model × target: mean/SD RMSE, MAE, R²_oos, Spearman; paired ΔRMSE vs F8 (same model) and vs OLS-F8, with NB-corrected CI and the rule flag |
| `feature_ablation_results.csv` | FULL minus each block, and F8 plus each block, per model and target |
| `results/` | `metrics_main.csv` (1,500 rows), `metrics_controls.csv`, `metrics_sensitivity.csv`, `shortlist_robustness.csv`, `sensitivity_summary.csv`, `extra_vs_tree_baselines.csv`, `control_shuffled_B.csv`, `ols_coefficient_stability.csv`, `all_paired_comparisons.csv`, `analysis_summary.json`, run logs |
| `figures/` | `fig1_representation_x_model_r2.png`, `fig2_feature_block_ablation.png`, `fig3_stability_per_split.png` |

## Reproduce

From `scripts/` (Python 3.14, scikit-learn 1.9.1, CatBoost 1.2.10; inputs from Phases 2–4 on D: and in the repo):

```powershell
python p4b_run.py --stage main --jobs 19
```
```powershell
python p4b_run.py --stage controls --jobs 19
```
```powershell
python p4b_analyze.py
```
```powershell
python p4b_sensitivity.py
```
```powershell
python -m pytest ..\tests -q
```

**Run information (2026-11-22, 20 cores, 19 workers):**
- main stage: 1,500 tasks, 54 min;
- controls: 40 tasks, 4 min;
- sensitivity: 1,950 tasks, 14 min.

## Data (unchanged references)
- **Target:** GSE294180 DESeq2 PE_in_CD / PE_in_MD (`results/deseq2/provisional_deposited`).
- **References** (Phase 3 `reference_effects/`):
  - GSE316023 MASH1 (TA wk23);
  - GSE107787 FAST (ZT-adjusted) and FASTclean (ZT0/ZT4 muscle-low);
  - GSE282998 PPARα WT 5v0/10v0/20v0;
  - GSE117134 CIRC T1.
- **Gene sets:** Phase 3 `results/gene_sets.csv` (MSigDB 2026.1.Mm Hallmark + 3 PPAR sets + 2 custom sets).
- **Modules:** Phase 4 `results/gene_groups.csv` (Ward, GSE117134).
