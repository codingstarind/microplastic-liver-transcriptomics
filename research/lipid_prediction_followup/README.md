# Lipid-gene prediction follow-up (exploratory, post hoc; 2026-11-25)

**Question (as posed):** do the existing gene-level regression models predict PE-associated changes in lipid-metabolism genes better than they predict other genes? Each diet (CD, MD) is analysed separately.

**Scope:**
- No model was refitted with new settings, tuned on test data, or selected.
- FULL-OLS and FULL-Ridge held-out predictions were regenerated with the unchanged Phase 4B code, because Phase 4B saved only metrics. The regenerated predictions match the stored Phase 4B metrics on all 15 splits per model and diet (`results/reproduction_check_vs_phase4b.csv`).
- No earlier outputs were modified.

## Inputs (read-only)

| Input | Source |
|---|---|
| Splits | Locked Phase 4 outer splits: 20 Ward co-expression modules (GSE117134), 3 repeats × 5 folds, seeds 101–103. Whole modules go to one fold (`roadmap_phase4/scripts/p4_data.py::outer_splits`, `roadmap_phase4/results/gene_groups.csv`). |
| F8-OLS, F8+C-OLS, F8+C-Ridge predictions | `roadmap_phase4c_diagnostics/tables/fix_experiment_oof_predictions.csv.gz`. The Phase 4C test confirms they reproduce Phase 4B exactly. |
| FULL (28 features = F8 + A + B + C) predictions | Regenerated through `p4b_run` / `RepBuilder("FULL")` / `p4_models.make`, with nested GroupKFold(4) tuning and DL weights, as in Phase 4B. |
| Targets | `results/deseq2/provisional_deposited/results_PE_in_{CD,MD}.csv.gz` (log2FC). Universe: 11,758 genes. |
| Gene sets | `roadmap_phase3/results/gene_sets.csv`: LIPOGENESIS_CORE, CHOL_SYNTHESIS_CORE, SANDERSON_PPARA_TARGETS. These are the same frozen sets used as block-B features. |

## Pre-specified model designations
These were fixed from earlier phase documents before any group metric was computed.

| Role | Model | Designation |
|---|---|---|
| F8 | **F8-OLS** | Phase 4 primary gene model |
| F8+C | **F8+C-Ridge** | Phase 4B decision memo recommendation; Phase 5 frozen secondary model |
| FULL | **FULL-OLS** | Phase 4B primary FULL comparison row |
| Secondary estimators (tables only) | F8+C-OLS, FULL-Ridge | – |

**Caveat on FULL:** FULL contains binary membership indicators (block B) for exactly the three lipid sets analysed here, plus 9 other sets. FULL's lipid-gene results are therefore not independent of set membership. F8 and F8+C contain no membership features.

## Definitions
- **Per-gene error:** squared error averaged over the 3 outer repeats. Each gene is held out exactly once per repeat.
- **Null baseline:** the weighted training-gene mean of each split, as in Phase 4/4B.
- **Group metrics:**
  - RMSE_model and RMSE_null;
  - dMSE = MSE_null − MSE_model (the improvement over the null);
  - **R²_group** = 1 − SSE_model / SSE_null;
  - fraction of genes improved;
  - Spearman correlation of y with prediction.
- **Groups:**
  - LIPO (12 genes), CHOL (21), PPARA (14);
  - LIPID_any (47): within the prediction universe the three sets do not overlap, so their exclusive versions are identical to them;
  - OTHER (11,711): none of the three sets;
  - OTHER_noB12: none of the 12 block-B sets.
- **Dependence and uncertainty**, using the 20 modules (the split grouping unit):
  - (a) module block bootstrap, B = 2,000;
  - (b) within-module label permutation, 2,000 permutations; the minimum attainable p is 1/2001 ≈ 0.0005;
  - (c) per-gene dMSE regressed on group indicator + module fixed effects, with CR1 module-cluster-robust SE and t on 19 df.
- **Gene-set overlap:**
  - pairwise overlap of the three sets (`results/gene_set_overlap_in_prediction_universe.csv`);
  - overlap of each lipid set with the 9 other block-B sets (`results/lipid_sets_overlap_with_other_blockB_sets.csv`);
  - OTHER_noB12 as an alternative reference;
  - joint regression with all three indicators plus module fixed effects (`results/joint_overlap_regression_moduleFE.csv`).
- **Scale-matched reference:** for each lipid gene, the 10 nearest OTHER genes on null squared error, matched without replacement (`results/group_vs_nullerror_matched_other.csv`, balance in `results/nullerror_matching_balance.csv`). This is a descriptive stratification that uses the observed target through the null error.
- **Model comparison:** gene-paired held-out gain of model A over model B, lipid group minus OTHER (`results/model_pair_gain_lipid_vs_other.csv`).

## Files

| Path | Content |
|---|---|
| `TABLE_results_concise.csv` | Primary models × diet × {LIPO, CHOL, PPARA, LIPID_any, OTHER}: group metrics with CIs; difference vs OTHER (bootstrap CI, permutation p, module-FE estimate, CI and p); difference vs null-error-matched OTHER (`matched_*`). |
| `figures/fig_lipid_vs_other_prediction.{png,pdf,svg}` | Two panels (CD, MD): R²_group with 95% module-bootstrap CI for F8-OLS, F8+C-Ridge and FULL-OLS. |
| `results/group_metrics_all.csv` | All 5 models × all groups. |
| `results/group_vs_reference_comparisons.csv` | All models × groups vs OTHER and vs OTHER_noB12. |
| `results/group_metrics_per_repeat.csv` | R²_group and dMSE for each outer repeat. |
| `results/model_pair_gain_lipid_vs_other.csv` | F8+C vs F8, FULL vs F8, FULL vs F8+C. |
| `results/oof_predictions_with_null.csv.gz` | All held-out predictions with null and module. |
| `results/reproduction_check_vs_phase4b.csv`, `results/run_info_lp0*.json`, `logs/` | Reproduction and run records. |
| `EVIDENCE_AUDIT.md` | Numbers and calculations only. |
| `CONCLUSION.md` | Short conclusion requested by the team (F8 and F8+C; FULL excluded as evidence). |

## Reproduce

Run these from the repository root.

```powershell
.venv\Scripts\python.exe research\lipid_prediction_followup\scripts\lp01_oof_predictions.py
```
```powershell
.venv\Scripts\python.exe research\lipid_prediction_followup\scripts\lp02_group_analysis.py
```
```powershell
.venv\Scripts\python.exe research\lipid_prediction_followup\scripts\lp03_matched_and_figure.py
```

Approximate runtimes are 35 s, 3 min and several minutes. The scripts are deterministic (seed 20261125).

**Notes:**
- The block bootstrap draws from only 5–13 modules for the lipid groups, so some intervals are very wide or asymmetric. In the figure, intervals are truncated at −0.6, and the truncated values are printed on the plot.
- **Order of analysis:**
  - Fixed before any group results were seen: the models, the groups, the OTHER and OTHER_noB12 references, the bootstrap, the permutation and the module-FE tests.
  - Added after the first group results were seen: the null-error-matched reference (lp03) and the figure's axis truncation.
  - No model, group or reference was dropped. All are reported in full.
