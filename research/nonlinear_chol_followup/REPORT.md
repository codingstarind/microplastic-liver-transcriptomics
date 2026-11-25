# Linear vs nonlinear held-out prediction of cholesterol-synthesis genes (exploratory; 2026-11-25)

This report contains numbers only, with no interpretation. Values are taken from `TABLE_results_concise.csv` and `results/*.csv`.

## 1. Inputs and protocol (reused; nothing refitted or tuned)

**Predictions:** the original Phase 4 task-A held-out predictions, `research/roadmap_phase4/results/taskA/predictions/primary__PE_in_{CD,MD}.csv.gz`.

| Model | Phase 4 name |
|---|---|
| F8-OLS | `ols_all` |
| F8-Ridge | `ridge` |
| F8-CatBoost | `catboost` |
| F8-HGB | `hgb` |

**Features:** F8 = MASH, FAST, PPARA and CIRC reference effects + their SEs. There are no pathway-membership features.

**Splits:**
- 20 co-expression modules; 5 outer folds × 3 repeats (seeds 101–103).
- Hyperparameters were tuned by inner GroupKFold(4) on training genes only.
- All 4 models were evaluated on identical held-out genes: 11,758 genes per repeat.

**Reproduction against the stored Phase 4 per-split metrics** (120 model × diet × split cells; `results/reproduction_check_vs_phase4.csv`):

| Metric | Max absolute difference |
|---|---|
| RMSE | 1.1e-16 |
| R²_oos | 4.2e-16 |
| Sign agreement | 0 |

F8-OLS predictions are identical to those used in the lipid-prediction follow-up (max difference 3.4e-21).

**Definitions:**
- **Per-gene error:** squared error averaged over the 3 repeats.
- **Null:** the weighted training-gene mean of each split.
- **R²_group** = 1 − SSE / SSE_null.
- **Direction agreement:** share of gene × repeat predictions with sign(pred) = sign(log2FC).
- **95% CIs:** module block bootstrap (B = 2,000).
- **Matched reference:** for each group gene, the 10 nearest other genes on null squared error.

## 2. Decision rules (fixed before group results were computed)

**Primary family (8 tests):**
- Scope: cholesterol synthesis (21 genes, 6 modules) × {CatBoost, HGB} vs {OLS, Ridge} × {CD, MD}.
- Statistic: per-gene gain = MSE_linear − MSE_nonlinear (> 0 means the nonlinear model is better).
- Test: module-cluster t with G − 1 = 5 df, Holm-corrected within the family.
- Sensitivity test: exact module sign-flip. With 6 modules, the smallest attainable p is 0.031.

**Secondary family:** BH correction across all of the following:
- lipogenesis, PPARα and other genes, nonlinear vs linear;
- group-gain minus other-gene gain, and group-gain minus matched-gene gain, each with module fixed effects.

**Feature-combination follow-up:** to be run only if a primary test had Holm p < 0.05 with gain > 0. **This condition was not met, so the follow-up was not run.**

## 3. Cholesterol synthesis: group metrics (21 genes, 6 modules)

| Diet | Model | R²_group [95% CI] | RMSE [95% CI] | Direction agreement | Spearman | Matched other genes R²_group [95% CI] |
|---|---|---|---|---|---|---|
| CD | F8-OLS | 0.648 [0.596, 0.851] | 0.510 [0.200, 0.639] | 1.000 | 0.709 | 0.073 [0.028, 0.115] |
| CD | F8-Ridge | 0.644 [0.592, 0.848] | 0.513 [0.202, 0.642] | 1.000 | 0.714 | 0.072 [0.028, 0.114] |
| CD | F8-CatBoost | 0.659 [0.624, 0.802] | 0.503 [0.266, 0.544] | 0.952 | 0.742 | 0.087 [0.040, 0.131] |
| CD | F8-HGB | 0.614 [0.562, 0.782] | 0.535 [0.267, 0.581] | 0.921 | 0.634 | 0.084 [0.034, 0.128] |
| MD | F8-OLS | 0.600 [0.449, 0.704] | 0.617 [0.297, 0.725] | 1.000 | 0.626 | 0.041 [−0.017, 0.125] |
| MD | F8-Ridge | 0.598 [0.448, 0.702] | 0.619 [0.299, 0.726] | 1.000 | 0.626 | 0.041 [−0.017, 0.124] |
| MD | F8-CatBoost | 0.424 [0.274, 0.568] | 0.741 [0.393, 0.816] | 1.000 | 0.199 | 0.039 [0.008, 0.078] |
| MD | F8-HGB | 0.481 [0.418, 0.649] | 0.703 [0.329, 0.763] | 1.000 | 0.269 | 0.034 [0.002, 0.078] |

Null RMSE is 0.860 for CD and 0.976 for MD.

## 4. Cholesterol synthesis: gene-paired nonlinear vs linear (primary family)

| Diet | Comparison | ΔR²_group [bootstrap CI] | Mean gain in MSE [cluster-t CI] | p (cluster t) | Holm p (8 tests) | Sign-flip p | Genes better with nonlinear |
|---|---|---|---|---|---|---|---|
| CD | CatBoost − OLS | +0.011 [−0.130, +0.190] | +0.008 [−0.060, +0.076] | 0.776 | 1.000 | 1.000 | 8/21 |
| CD | CatBoost − Ridge | +0.015 [−0.126, +0.193] | +0.011 [−0.058, +0.079] | 0.703 | 1.000 | 1.000 | 8/21 |
| CD | HGB − OLS | −0.034 [−0.207, +0.161] | −0.025 [−0.097, +0.047] | 0.411 | 1.000 | 0.531 | 7/21 |
| CD | HGB − Ridge | −0.030 [−0.204, +0.165] | −0.022 [−0.095, +0.050] | 0.461 | 1.000 | 0.531 | 7/21 |
| MD | CatBoost − OLS | −0.176 [−0.302, −0.039] | −0.168 [−0.243, −0.093] | 0.0022 | 0.018 | 0.094 | 2/21 |
| MD | CatBoost − Ridge | −0.174 [−0.300, −0.037] | −0.166 [−0.240, −0.092] | 0.0023 | 0.018 | 0.094 | 2/21 |
| MD | HGB − OLS | −0.119 [−0.143, +0.011] | −0.113 [−0.192, −0.035] | 0.014 | 0.084 | 0.250 | 4/21 |
| MD | HGB − Ridge | −0.117 [−0.141, +0.014] | −0.111 [−0.189, −0.033] | 0.014 | 0.084 | 0.250 | 4/21 |

**Primary tests with Holm p < 0.05 and gain > 0: 0 of 8.**

**Change in direction agreement (nonlinear − linear):**

| Diet | CatBoost | HGB |
|---|---|---|
| CD | −0.048 | −0.079 |
| MD | 0.000 | 0.000 |

**Cholesterol-synthesis gain compared with other genes** (module fixed effects; `results/paired_nonlinear_vs_linear.csv`):

| Diet | Comparison | vs all other genes | vs matched other genes |
|---|---|---|---|
| CD | CatBoost − OLS | +0.004 [−0.049, +0.058] | −0.048 [−0.187, +0.091], BH q 0.63 |
| CD | HGB − OLS | −0.028 [−0.084, +0.029] | −0.063 [−0.188, +0.062], BH q 0.47 |
| MD | CatBoost − OLS | −0.167 [−0.225, −0.109], BH q 0.0003 | −0.151 [−0.258, −0.044], BH q 0.090 |
| MD | HGB − OLS | −0.114 [−0.175, −0.053], BH q 0.019 | −0.128 [−0.226, −0.031], BH q 0.090 |

The Ridge comparisons differ from the OLS ones by ≤ 0.003.

## 5. Controls: dominance by genes and modules (`results/chol_dominance_controls.csv`)

**Module distribution:** 15 of the 21 cholesterol-synthesis genes are in one module (module 14). The remaining modules hold 2, 1, 1, 1 and 1 genes.

| Diet | Comparison | ΔR² full | Leave-one-gene-out range (same sign / 21) | Leave-one-module-out range (same sign / 6) | Share of absolute gain from top 3 genes |
|---|---|---|---|---|---|
| CD | CatBoost − OLS | +0.011 | −0.020 to +0.053 (17) | −0.020 to +0.084 (5) | 0.51 |
| CD | HGB − OLS | −0.034 | −0.067 to +0.013 (20) | −0.067 to +0.058 (5) | 0.55 |
| MD | CatBoost − OLS | −0.176 | −0.186 to −0.160 (21) | −0.184 to −0.130 (6) | 0.41 |
| MD | HGB − OLS | −0.119 | −0.135 to −0.100 (21) | −0.132 to −0.037 (6) | 0.54 |

- The most influential gene in the leave-one-gene-out check is Sqle for every comparison.
- Per-gene predictions and errors for each model are in `results/chol_per_gene_predictions.csv`.

## 6. Secondary groups, nonlinear − OLS (`results/paired_nonlinear_vs_linear.csv`, `results/group_metrics_by_model.csv`)

| Diet | Group | ΔR² CatBoost − OLS [CI] | p / BH q | ΔR² HGB − OLS [CI] | p / BH q |
|---|---|---|---|---|---|
| CD | Lipogenesis (12) | +0.222 [0.009, 0.286] | 0.067 / 0.30 | +0.098 [−0.090, 0.158] | 0.19 / 0.41 |
| CD | PPARα (14) | −0.007 [−0.353, 0.375] | 0.97 / 0.97 | +0.174 [0.008, 0.303] | 0.080 / 0.30 |
| CD | Other (11,711) | +0.018 [−0.001, 0.033] | 0.088 / 0.30 | +0.014 [−0.000, 0.026] | 0.084 / 0.30 |
| MD | Lipogenesis | +0.047 [−0.533, 0.113] | 0.58 / 0.67 | +0.042 [−1.025, 0.079] | 0.45 / 0.60 |
| MD | PPARα | −0.338 [−0.470, 0.082] | 0.24 / 0.45 | −0.303 [−0.493, −0.008] | 0.23 / 0.45 |
| MD | Other | +0.007 [−0.014, 0.034] | 0.56 / 0.67 | +0.006 [−0.015, 0.032] | 0.59 / 0.68 |

No secondary test has BH q < 0.05 except the two MD cholesterol-synthesis vs other-gene contrasts in section 4.

## 7. Not done, and why
- **Feature-combination analysis** (fasting, PPARα, MASH, circadian): its pre-set trigger was not met (section 2).
- **XGBoost:** not added, as requested.
- **FULL and other pathway-membership models:** not used, because the comparison is restricted to F8 features.

## Limits stated for the reader
- Exploratory: the Phase 4 splits have been reused across later phases.
- The cholesterol-synthesis set spans 6 modules, and 15 of its 21 genes are in one module.
- Module-cluster inference therefore has 5 df, and the sign-flip test cannot go below p = 0.031.
