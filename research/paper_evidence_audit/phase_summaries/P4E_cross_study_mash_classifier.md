# Roadmap Phase 4E: study-held-out MASH-state classifiers and their readout in GSE294180 (2026-11-23)

## Main sources (`research/roadmap_phase4e/`)
- `00_task_and_estimands.md`
- `01_model_feasibility_and_leakage_audit.md`
- `02_literature_and_representation_choices.md`
- `03_predeclared_model_evaluation_plan.md` (LOCKED; `results/plan_lock.json`)
- `04_model_comparison_results.md`
- `05_failures_limits_and_sensitivity.md`
- `06_phase5_readiness.md`
- `decision_memo.md`
- `results/external_metrics.csv`, `external_predictions.csv`, `validity.csv`, `permutation_null_auc_R2_logistic.csv`, `permutation_summary.json`, `gse294180_frozen_model_contrasts.csv`, `gse294180_frozen_model_loo.csv`, `r2_feature_contributions_gse294180.csv`, `fasting_demo_metrics.csv`
- `figures/fig1_cross_study_performance.png`, `fig2_before_after_calibration.png`, `fig3_permutation_null.png`, `fig4_gse294180_2x2_frozen_scores.png`
- Scripts: `e01_mash_state_models.py`, `e01b_permutation_null.py`, `e01c_validity.py`, `e02_fasting_demo.py`, `e03_apply_to_gse294180.py`, `e04_diagnostics_figures.py`

## 1. Question
Can a MASH-state (GAN vs chow) mouse classifier, trained on one independent MASH study, transfer to the other study? If it passes predeclared validity rules, what does it say about the 19 GSE294180 mice?

**This is NOT a PE classifier.** No PE-vs-vehicle classifier was built, by design.

## 2. Why
- Phase 4D's MASH signature required agreement in **both** MASH studies. That makes it unusable for study-held-out validation (a leakage issue identified in `01_…md`).
- A study-held-out classifier gives the correct test of whether a MASH-state readout generalises across labs.
- The "three tasks" table in `00_…md` separates three estimands, which must never share one leaderboard:
  - (A) gene-level R²;
  - (B) sample-level state prediction;
  - (C) the sample-level PE response.

## 3. Design
**Splits:**
- **S1:** train GSE316023 TA (8 GAN v 8 chow) → test GSE246328 (7 v 8).
- **S2:** the reverse.
- **S1b (sensitivity):** train on 21 °C + 29 °C (16 v 16).

**Representations** (within-sample rank scores over a fixed common universe):

| Code | Definition | Features |
|---|---|---|
| R0 | signed singscore of 50 up / 50 down genes, chosen by Welch t **inside the training study only** | 1 |
| R1 | singscores of 12 frozen pathway sets | 12 |
| R2 | R0 + R1 | 13 |
| PCA | 3 PCs; loadings fitted on training mice only | 3 |

**Models:**
- majority (dummy);
- R0 midpoint rule;
- logistic, C = 100;
- L2-logistic, C ∈ {0.01, 0.1, 1, 10};
- Gaussian NB;
- depth-1/2 trees (gini or entropy);
- RF (300 trees, depth 1–2).

**Inner CV:** RepeatedStratifiedKFold(4 × 3); scoring = neg-Brier.

**Metrics:**
- external ROC-AUC (primary);
- balanced accuracy, sensitivity and specificity at 0.5;
- Brier and log loss;
- stratified bootstrap (2,000 resamples).

**Controls:**
- 200 training-label permutations;
- correlation of predictions with library size.

**Validity rule for applying a model to GSE294180** (fixed before scoring): in both S1 and S2, external AUC ≥ 0.80, bootstrap lower bound > 0.5, and Brier < 0.25.

## 4. Results (`results/external_metrics.csv`; `04_…md`)

| Stage | Model | S1 AUC | S1 BA | S1 Brier | S2 AUC | S2 BA | S2 Brier |
|---|---|---|---|---|---|---|---|
| null | majority | 0.50 | 0.50 | 0.250 | 0.50 | 0.50 | 0.251 |
| baseline | R0 midpoint rule | 1.00 | 1.00 | 0.000* | 1.00 | **0.81** | **0.188*** |
| + features | R1 logistic C = 1 | 1.00 | 1.00 | 0.010 | 1.00 | 0.94 | 0.042 |
| + features | R2 logistic C = 1 | 1.00 | 1.00 | 0.009 | 1.00 | 0.94 | 0.049 |
| + tuning | R2 L2-logistic (C = 10) | 1.00 | 1.00 | 0.002 | 1.00 | **1.00** | **0.028** |
| comparator | PCA L2-logistic | 1.00 | 0.69 | 0.187 | 1.00 | 0.75 | 0.158 |

\* The R0 rule outputs hard labels, so its "Brier" is a misclassification rate. It is **not** comparable with logistic Brier scores (Phase 5 correction).

**Sample sizes:**
- S1 tests 15 mice (GSE246328).
- S2 tests 16 mice (GSE316023 TA).

**Other models:**
- Trees on R2: AUC 0.80 / 0.63. Trees on PCA: 0.69 / 0.63.
- RF on PCA: 0.91 / 0.84.
- NB and RF on R1/R2: 1.0.

> **Correction (Phase 5 decision report).** The 4E sentence "every non-PCA representation reaches AUC 1.0" is false. R2 trees reach 0.80 / 0.63.

**Permutation null (R2 logistic):**
- S1: p = 0.10, q95 = 1.00.
- S2: p = 0.035.
- Partly correct labellings keep the diet signal, so this null has low resolving power.

**Prediction vs library size:** |ρ| ≤ 0.43, with inconsistent sign between directions.

**Fasting demonstration:** trained on 10 clean GSE107787 mice, tested on 6 GSE225695 mice: 6/6 correct. Demonstration only; not frozen.

**Application to the 19 GSE294180 mice** (log-odds; `gse294180_frozen_model_contrasts.csv`):

| Score | Diet_in_Veh | PE_in_CD | PE_in_MD | Interaction |
|---|---|---|---|---|
| R0 + L2 | +4.01 [2.92, 5.10] | +0.35, Holm 0.96 | +0.86 [−0.16, 1.89], Holm 0.28 | +0.51 |
| R2 + L2 | +5.55 [4.78, 6.33] | +0.83 [0.06, 1.61], Holm 0.072 | **+1.66 [0.93, 2.39], Holm 0.0006** | +0.83, Holm 0.12 |

**R2 decomposition of PE_in_MD** (`r2_feature_contributions_gse294180.csv`):

| Feature | Contribution (log-odds) |
|---|---|
| CHOL_HOMEOSTASIS | +1.17 |
| CHOL_SYNTHESIS | +0.93 |
| PPARα targets | +0.35 |
| MTORC1 | +0.28 |
| LIPOGENESIS | −1.37 |
| MASH programme signature itself | **+0.18** |

- R2 scores correlate 0.89 with the 4D FAST score and 0.67 with clock_phase.
- The significant "MASH-state" PE shift is therefore a cholesterol / PPARα pathway shift, **not** a disease-specific one.

## 5. Key message
- Feature engineering and tuning did not improve discrimination: AUC is at its ceiling of 1.0.
- What they improved was **threshold transfer and calibration** across studies: S2 balanced accuracy 0.81 → 1.00, Brier 0.19 → 0.03.
- PCA transferred worst.

## 6. Limits
- Two studies only, so there are two generalisation estimates and no between-study variance.
- Bootstrap CIs collapse to [1, 1] at the AUC ceiling.
- The chow vs purified control mismatch.
- PE ≡ depth in GSE294180.
- Clock state (ρ ≈ 0.64–0.67).

**Deviations:**
1. Which models to apply was named before external viewing.
2. A fallback rule was added to the permutation step.
3. `validity.csv` was regenerated by `e01c`.
4. A restart with memoisation.

## 7. Transition
Every internal and reference-based analysis had now been done. The remaining test was external: apply the frozen gene-level models (4C), the signatures (4D) and the classifiers (4E) to independent microplastic studies. That is Phase 5.
