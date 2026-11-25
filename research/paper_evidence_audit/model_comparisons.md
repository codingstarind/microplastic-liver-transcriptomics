# Model comparisons and feature engineering: publication-style evidence tables

Every number below is copied from a saved output file. The file is named under each table, with machine-readable copies in `tables/T*.csv` (each row carries a `source` column).

**Rounding:** 3 decimals, except ΔRMSE (4 decimals).

**Nothing was re-run** to produce this document.

---

## 0. Which comparisons are statistically valid? (read first)

| Comparison type | Unit | Valid inference | What it does **not** show |
|---|---|---|---|
| Gene-level held-out CV (Phases 4, 4B, 4C) | gene (11,758 rows, 20 module groups) | Paired ΔRMSE with a Nadeau–Bengio corrected interval over the **same 15 outer splits** | Anything about new mice or new studies: every gene's target comes from the same 9–10 mice |
| Phase 4B / 4C vs Phase 4 | gene | **Exploratory.** The features were proposed after the Phase 4 results, and the same splits were reused, so the intervals are optimistic. | A confirmatory improvement |
| Study-held-out classifier (Phase 4E) | mouse (15–16 test mice, 2 directions) | Point estimates per direction; a bootstrap over test mice (Phase 5B) | A population generalisation error (only 2 studies). Bootstrap CIs collapse at an AUC ceiling. |
| External gene-level (Phase 5) | gene within 6 external mice | **Mouse-relabelling permutation (20 labellings)** and the mouse bootstrap (5B) | The gene-block bootstrap CI is anti-conservative. Do not use it as the headline interval. |
| Phase 5B interventions on Hong | 3 v 3 mice | **Post-external exploratory.** Hong outcomes were already seen. | Fresh validation |

**Different estimands must never share one leaderboard** (`roadmap_phase4e/00_task_and_estimands.md`):
- gene-level R² / Spearman;
- mouse-level classifier AUC;
- mouse-level signature contrasts.

---

## 1. Gene-level prediction: baselines and model families (Phase 4, locked protocol)

**Task.** Predict the held-out GSE294180 PE log2FC from the effects and SEs of 4 independent references.
- 11,758 genes;
- 20 co-expression module groups;
- 5 outer folds × 3 repeats;
- inner GroupKFold(4).

**Source:** `research/roadmap_phase4/results/taskA/summary_by_run_model.csv` → `tables/T01_phase4_taskA_models.csv`.

| Model | Features | R²_oos PE_in_CD | R²_oos PE_in_MD | R²_oos Diet_in_Veh (positive control) | Spearman CD | Spearman MD | RMSE CD | RMSE MD |
|---|---|---|---|---|---|---|---|---|
| null (training mean) | – | 0.000 | 0.000 | 0.000 | – | – | 0.326 | 0.309 |
| ols_MASH | MASH effect | 0.031 | −0.002 | 0.206 | 0.208 | 0.166 | 0.321 | 0.309 |
| ols_FAST | FAST effect | 0.022 | 0.045 | 0.012 | 0.135 | 0.246 | 0.323 | 0.302 |
| ols_PPARA | PPARα effect | −0.000 | −0.001 | 0.001 | 0.051 | 0.083 | 0.326 | 0.309 |
| ols_CIRC | CIRC_T1 effect | 0.022 | 0.005 | 0.040 | 0.104 | 0.080 | 0.323 | 0.308 |
| **ols_all (F8-OLS)** | 4 effects + 4 SEs | **0.066** | **0.070** | 0.233 | **0.264** | **0.309** | 0.316 | 0.298 |
| ridge | F8 | 0.066 | 0.070 | 0.233 | 0.264 | 0.309 | 0.315 | 0.298 |
| elastic net | F8 | 0.063 | 0.068 | 0.232 | 0.262 | 0.308 | 0.316 | 0.298 |
| CART tree | F8 | 0.038 | 0.053 | 0.179 | 0.221 | 0.288 | 0.320 | 0.301 |
| random forest | F8 | 0.073 | 0.081 | 0.247 | 0.281 | 0.338 | 0.314 | 0.296 |
| HistGradientBoosting | F8 | 0.078 | 0.079 | 0.285 | 0.276 | 0.333 | 0.313 | 0.297 |
| CatBoost | F8 | **0.086** | 0.079 | **0.306** | 0.277 | 0.341 | 0.312 | 0.297 |

### Paired comparisons
Sources: `paired_comparisons.csv` → `T02`; decision report §2. Values are ΔRMSE with the NB-corrected 95% CI.

| Comparison | PE_in_CD | PE_in_MD | Diet_in_Veh |
|---|---|---|---|
| ols_all vs null | −0.0107 [−0.017, −0.004] **better** | −0.0112 [−0.025, +0.003] not distinguishable | −0.061 [−0.080, −0.043] better |
| ols_all vs best single axis | −0.0058 [−0.011, −0.0003] vs MASH, **better** | −0.0040 [−0.012, +0.004] vs FAST, n.s. | −0.007 [−0.020, +0.006] vs MASH, n.s. |
| RF / HGB / CatBoost vs ols_all | −0.001 / −0.002 / −0.004; all CIs include 0 | −0.001 each; all CIs include 0 | HGB −0.015 [−0.026, −0.004]; CatBoost −0.022 [−0.037, −0.006] **better** |

### Reading
- **A small, robust signal.** About 7% of held-out PE-effect variance is predictable. The ols_all > null result holds in 8/8 sensitivity variants and 19/19 LOO targets.
- **Flexible models do not beat linear on the PE targets.** They do beat it on the **positive-control** target (the strong, structured diet effect). So the tree ensembles were capable of winning; the PE signal simply has no nonlinear structure to exploit.
- **The MASH single axis is fragile:** below null in 12/19 LOO targets for MD.
- **Negative controls:** a permuted target gives R² −0.002 to 0.000; a depth-thinned target gives −0.005 to +0.001.

---

## 2. Feature engineering and ablation (Phase 4B, EXPLORATORY)

**Feature blocks:**
- **A:** empirical-Bayes reliability-shrunken effects;
- **B:** 12 pathway-membership indicators;
- **C:** response shape (FASTxPPARA, MASHxFAST, FAST_morning_dev, PPARA_consistent).

**FULL** = F8 + A + B + C (28 features).

**Source:** `research/roadmap_phase4b/results/metrics_main.csv` → `tables/T03_phase4b_representation_x_model.csv` (means over 15 splits).

### R²_oos by representation × model

| Representation (n features) | OLS CD | Ridge CD | ENet CD | HGB CD | CatBoost CD | OLS MD | Ridge MD | ENet MD | HGB MD | CatBoost MD |
|---|---|---|---|---|---|---|---|---|---|---|
| F4 (4) | 0.064 | 0.064 | 0.064 | 0.065 | 0.070 | 0.062 | 0.062 | 0.062 | 0.068 | 0.072 |
| **F8 (8)** | **0.066** | 0.066 | 0.063 | 0.078 | 0.085 | **0.070** | 0.070 | 0.068 | 0.079 | 0.079 |
| A4 (4) | 0.067 | 0.067 | 0.067 | 0.080 | 0.081 | 0.067 | 0.067 | 0.066 | 0.075 | 0.074 |
| F8+A (12) | 0.074 | 0.071 | 0.068 | 0.089 | 0.100 | 0.074 | 0.077 | 0.075 | 0.079 | 0.080 |
| F8+B (20) | 0.085 | 0.085 | 0.082 | 0.081 | 0.100 | 0.081 | 0.081 | 0.075 | 0.080 | 0.083 |
| **F8+C (12)** | **0.080** | 0.080 | 0.076 | 0.084 | 0.098 | **0.096** | 0.096 | 0.094 | 0.091 | 0.092 |
| **FULL (28)** | **0.104** | **0.104** | 0.100 | 0.094 | **0.109** | **0.104** | **0.108** | 0.106 | 0.092 | 0.096 |
| FULL−A (24) = "F8+B+C" | 0.094 | 0.094 | 0.092 | 0.086 | 0.106 | 0.107 | 0.107 | 0.105 | 0.092 | 0.094 |
| FULL−B (16) | 0.091 | 0.091 | 0.088 | 0.097 | 0.105 | 0.096 | 0.098 | 0.096 | 0.091 | 0.095 |
| FULL−C (24) | 0.092 | 0.090 | 0.086 | 0.090 | 0.102 | 0.082 | 0.086 | 0.083 | 0.082 | 0.084 |

### Improvement over F8 within the same model
Source: `all_paired_comparisons.csv`, comparison `vs_F8_same_model` → `T04`.

**Locked rule:** NB CI < 0, ≥ 12/15 splits better, ΔRMSE ≤ −0.005 and ΔR² ≥ +0.01.

| Target | Representation × model | ΔRMSE [95% CI] | Splits better | ΔR² | Meets rule |
|---|---|---|---|---|---|
| CD | **FULL-OLS** | −0.0065 [−0.0124, −0.0006] | 15/15 | +0.038 | **yes** |
| CD | FULL-Ridge | −0.0065 [−0.0124, −0.0006] | 15/15 | +0.038 | yes |
| CD | F8+C-OLS | −0.0024 [−0.0059, +0.0010] | 12/15 | +0.014 | no |
| CD | F8+B-OLS | −0.0032 [−0.0080, +0.0015] | 13/15 | +0.019 | no |
| CD | FULL-CatBoost | −0.0041 [−0.0099, +0.0018] | 14/15 | +0.024 | no |
| MD | **FULL-OLS** | −0.0053 [−0.0100, −0.0005] | 14/15 | +0.034 | **yes** |
| MD | FULL-Ridge | −0.0059 [−0.0106, −0.0012] | 14/15 | +0.038 | yes |
| MD | F8+C-OLS | −0.0039 [−0.0092, +0.0015] | 13/15 | +0.026 | no |
| MD | FULL-CatBoost | −0.0027 [−0.0054, −0.0001] | 14/15 | +0.017 | no (practical size) |

### Strongest tree baseline: the comparison that matters
Source: `extra_vs_tree_baselines.csv` → `T05`.

| Target | Comparison | ΔRMSE [95% CI] | Splits better | ΔR² |
|---|---|---|---|---|
| CD | FULL-OLS vs **CatBoost-F8** | −0.0031 [−0.0111, +0.0049] | 9/15 | +0.018 |
| CD | FULL-OLS vs HGB-F8 | −0.0043 [−0.0127, +0.0042] | 12/15 | +0.025 |
| CD | FULL-OLS vs CatBoost-FULL | +0.0010 [−0.0034, +0.0054] | 7/15 | −0.005 |
| MD | FULL-OLS vs CatBoost-F8 | −0.0042 [−0.0092, +0.0009] | 14/15 | +0.025 |
| MD | FULL-Ridge vs HGB-F8 | −0.0047 [−0.0093, −0.0001] | 14/15 | +0.029 |
| MD | FULL-OLS vs CatBoost-FULL | −0.0014 [−0.0051, +0.0022] | 9/15 | +0.008 |

> **Key distinction for the paper.**
> - FULL engineered features improved the **linear** models clearly and reliably:
>   - CD: R² 0.066 → 0.104;
>   - MD: R² 0.070 → 0.104–0.108;
>   - the locked rule was met, with 7/7 variants and 19/19 LOO targets per the decision memo.
> - Against the **strongest tree baselines** the gain is statistically tied: the CI includes 0 in every CD comparison and in most MD comparisons.
> - The only CI fully below 0 vs a tree model is FULL-Ridge vs HGB-F8 for MD, with upper bound −0.0001. That is a marginal, fold-reusing, exploratory result.
> - Engineered features let a transparent model reach what boosting already achieved. They did not create a better predictor.

### Block ablations: removing a block from FULL
Source: `research/roadmap_phase4b/feature_ablation_results.csv` → `T05b`. Values are ΔRMSE(FULL − block) − RMSE(FULL); > 0 means the block helps.

| Block removed | OLS CD | Ridge CD | OLS MD | Ridge MD | HGB / CatBoost |
|---|---|---|---|---|---|
| A | +0.0017 [−0.0016, 0.0050] | +0.0017 | −0.0006 [−0.0017, 0.0005] | +0.0001 | n.s. |
| B | +0.0022 [−0.0021, 0.0064] | +0.0022 | +0.0016 | +0.0018 | n.s. |
| **C** | +0.0020 [−0.0004, 0.0045] | **+0.0023 [0.0001, 0.0045]** | **+0.0034 [0.0004, 0.0064]** | **+0.0033 [0.0001, 0.0065]** | n.s. (except CatBoost: adding C to F8 helps for CD) |

### Block verdicts (decision memo)
- **A is not useful.** It also creates near-duplicate columns, so FAST raw −0.32 vs shrunken +0.36 per SD inside FULL.
- **B helps linear CD**, but the gain disappears with shuffled membership (ridge CD 0.083 real vs 0.066 shuffled). Its largest coefficient is a negative LIPOGENESIS_CORE offset, i.e. it **re-learns GSE294180's own pathway shift**.
- **C is the most consistent contributor.** Removing it significantly hurts every linear model for MD.
  - Stable coefficients: FAST_morning_dev (+0.050 MD); FASTxPPARA (+0.024 CD, +0.013 MD).

---

## 3. Why FULL's coefficients are unreliable despite better prediction (Phase 4C)

Source: `research/roadmap_phase4c_diagnostics/tables/collinearity_summary.csv` → `T06`.

| Matrix | Features | max \|r\| | max VIF | VIF > 10 | Median condition number |
|---|---|---|---|---|---|
| F8 | 8 | 0.663 | 1.88 | 0 | 2.48 |
| A4 | 4 | 0.134 | 1.03 | 0 | 1.22 |
| B12 | 12 | 0.445 | 1.28 | 0 | 1.84 |
| C4 | 4 | 0.115 | 1.02 | 0 | 1.17 |
| F8+C | 12 | 0.856 | 3.78 | 0 | 3.76 |
| F8+B+C | 24 | 0.856 | 3.80 | 0 | 4.05 |
| **FULL** | 28 | **0.998** | **233.0** | **4** | **37.3** |

- The collinearity comes from raw vs shrunken effects: FAST r = 0.998; MASH 0.979; CIRC 0.926.
- FULL coefficients flip sign across folds in up to 12 of 28 features (MD). F8 coefficients never flip (`coefficient_stability_*.csv`).
- Predictions stay stable because collinearity only redistributes weight among near-duplicate columns. Prediction and identification are different properties.

### Fix experiment (`fix_experiment_summary.csv` → `T07`)
Values are descriptive: the same folds were reused.

| Model | R²_oos CD | R²_oos MD | ΔRMSE vs F8-OLS CD [CI] | ΔRMSE MD [CI] | Condition number | Ridge α chosen |
|---|---|---|---|---|---|---|
| F8-OLS | 0.066 | 0.070 | – | – | 2.48 | – |
| F8+C-OLS | 0.080 | 0.096 | −0.0024 [−0.0059, 0.0010] | −0.0039 [−0.0092, 0.0015] | 3.76 | – |
| F8+C-Ridge | 0.080 | 0.096 | −0.0024 [−0.0058, 0.0009] | −0.0039 [−0.0092, 0.0015] | 3.76 | 1000 in 15/15 |
| F8+B+C-Ridge | 0.094 | 0.107 | −0.0048 [−0.0107, 0.0011] | −0.0059 [−0.0103, −0.0014] | 4.05 | 1000 |
| FULL-Ridge (historical, 4B) | 0.104 | 0.108 | −0.0065 | −0.0059 | ≈ 37 | – |
| CatBoost-F8 (historical, 4) | 0.085 | 0.079 | −0.0034 [−0.0091, 0.0023] | −0.0011 [−0.0067, 0.0045] | – | – |

**Frozen for Phase 5:**
- **F8-OLS** (primary): collinearity-free, stable signs, no target-informed design choices.
- **F8+C-Ridge** (secondary).
- **Excluded:** F8+B+C (block B is target-informed) and FULL (collinear).

**Target ceiling** (`T08`): signal share is 0.356 (CD) and 0.391 (MD) from DESeq2 SEs, and 0.051 / 0.115 from the jackknife. The observed R² of 0.07–0.11 is low but not far from what the target allows.

---

## 4. Mouse-level MASH-state classifiers (Phase 4E; NOT PE classifiers)

**Design:**
- **S1:** train GSE316023 TA (8 GAN v 8 chow) → test GSE246328 (n = 15: 7 GAN, 8 chow).
- **S2:** train GSE246328 → test GSE316023 TA (n = 16).

**Representations:**
- **R0:** signed MASH programme singscore; genes chosen by Welch t inside the training study.
- **R1:** 12 pathway singscores.
- **R2:** R0 + R1.
- **PCA:** 3 PCs; loadings from training mice only.

**Source:** `research/roadmap_phase4e/results/external_metrics.csv` → `T11`; `04_model_comparison_results.md`.

### Baseline → feature engineering → tuning

| Stage | Representation / model | S1 AUC | S1 BA | S1 Brier | S2 AUC | S2 BA | S2 Brier |
|---|---|---|---|---|---|---|---|
| null | majority | 0.50 | 0.50 | 0.250 | 0.50 | 0.50 | 0.251 |
| baseline | R0 midpoint rule | 1.00 | 1.00 | 0.000† | 1.00 | **0.81** | **0.188†** |
| + features | R1, logistic C = 1 | 1.00 | 1.00 | 0.010 | 1.00 | 0.94 | 0.042 |
| + features | R2, logistic C = 1 | 1.00 | 1.00 | 0.009 | 1.00 | 0.94 | 0.049 |
| **+ tuning** | **R2, L2-logistic (inner CV chose C = 10)** | 1.00 | 1.00 | **0.002** | 1.00 | **1.00** | **0.028** |
| comparator | PCA, L2-logistic | 1.00 | 0.69 | 0.187 | 1.00 | 0.75 | 0.158 |
| other | trees on R2 (AUC) | 0.80 | – | – | 0.63 | – | – |
| other | RF on PCA (AUC) | 0.91 | – | – | 0.84 | – | – |

† The R0 rule outputs hard 0/1 labels, so its "Brier" equals the **misclassification rate**. It is not comparable with probabilistic Brier scores (Phase 5 correction).

### Bootstrap of saved predictions
Source: Phase 5B `tables/mash_classifier_bootstrap.csv` → `T17`; 2,000 stratified resamples of test mice.

| Split | Model | n | AUC | P(boot AUC = 1) | Balanced accuracy [95%] | Probabilistic metric [95%] |
|---|---|---|---|---|---|---|
| S1 | R0 rule | 15 | 1.00 | 100% | 1.00 [1.00, 1.00] | misclassification 0.00 |
| S1 | R2 L2-logistic | 15 | 1.00 | 100% | 1.00 [1.00, 1.00] | Brier 0.002 [0.001, 0.003] |
| S1 | PCA L2-logistic | 15 | 1.00 | 100% | 0.69 [0.56, 0.88] | Brier 0.19 [0.12, 0.27] |
| S2 | R0 rule | 16 | 1.00 | 100% | **0.81 [0.63, 0.94]** | misclassification 0.19 |
| S2 | R1 logistic | 16 | 1.00 | 100% | 0.94 [0.81, 1.00] | Brier 0.042 [0.018, 0.074] |
| S2 | R2 L2-logistic | 16 | 1.00 | 100% | **1.00 [1.00, 1.00]** | Brier 0.028 [0.008, 0.056] |

### Discrimination vs threshold transfer
- **Discrimination (AUC) is at its ceiling** for the logistic, NB and RF models on R0/R1/R2. GAN vs chow is easy, so AUC cannot rank models.
- **Feature engineering and tuning improved threshold transfer and probability quality.** In S2, the R0 rule's training-study threshold misclassified 3/16 mice (balanced accuracy 0.81). R2 L2-logistic reached 1.00 with Brier 0.028.

> **Correction.** `04_model_comparison_results.md` and `figures/fig1` say "AUC 1.0 for every non-PCA representation". The table itself shows trees on R2 at 0.80 / 0.63, and Phase 5 recorded that the R0 tree reached 0.81 in S2. Use the corrected statement.

### Permutation null
- 200 training-label permutations of R2 logistic: S1 p = 0.10 (q95 1.00); S2 p = 0.035.
- The null has low resolving power because partly correct labellings keep the diet signal.

### Application to GSE294180 mice (log-odds, 95% CI; `T12`)

| Score | Diet_in_Veh | PE_in_CD | PE_in_MD | Interaction |
|---|---|---|---|---|
| R0 + L2 | +4.01 [2.92, 5.10] | +0.35 (Holm 0.96) | +0.86 [−0.16, 1.89] (Holm 0.28) | +0.51 |
| R2 + L2 | +5.55 [4.78, 6.33] | +0.83 [0.06, 1.61] (Holm 0.072) | **+1.66 [0.93, 2.39] (Holm 0.0006)** | +0.83 |

The R2 PE_in_MD shift decomposes (`T13`) as follows:
- cholesterol homeostasis +1.17;
- cholesterol synthesis +0.93;
- PPARα targets +0.35;
- MTORC1 +0.28;
- lipogenesis −1.37;
- **MASH programme score +0.18**.

> **This is a cholesterol / PPARα pathway shift read through a MASH-vs-chow classifier.** It is not evidence that PE induces a MASH state.

---

## 5. External gene-level transfer of frozen models (Phase 5) and post-external interventions (Phase 5B)

### Phase 5 locked result
Source: `research/roadmap_phase5/gene_level_external_results.csv` → `T14`. Hong PE vs CON, 3 v 3; models trained on PE_in_CD.

| Model | Spearman [gene-block CI] | Mouse-relabelling p (20) | Weighted R² vs own mean | Verdict |
|---|---|---|---|---|
| F8-OLS | 0.219 [0.14, 0.28] | 0.25 | 0.041 | not passed |
| F8+C-Ridge | 0.142 [0.08, 0.20] | 0.30 | 0.004 | not passed |
| FAST-only OLS | 0.29 [0.23, 0.34] | 0.05 (minimum attainable) | 0.045 | PASS-direction |

**GSE245069 pooled PS:** F8-OLS ρ = −0.07, R² −0.14 → **FAIL**.

### Phase 5B: same task, before vs after
Source: `tables/intervention_verdicts.csv` → `T16`. All Hong numbers are post-external exploratory.

| Candidate | Train ρ | CV ρ (SD) | CV R² | Hong ρ | Δ vs B0 | Mouse-bootstrap 95% | P(Δ > 0) | Rank /20 | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| B0 F8-OLS | 0.274 | 0.264 (0.040) | 0.066 | 0.219 | – | −0.03 to 0.30 | – | 4 | comparator |
| B1 F8+C-Ridge | 0.281 | 0.270 (0.040) | 0.081 | 0.142 | −0.076 | −0.06 to 0.29 | 0.14 | 5 | comparator |
| B2 HGB-F8 | **0.378** | 0.276 (0.046) | 0.079 | 0.206 | −0.012 | −0.04 to 0.30 | 0.37 | 4 | variance gap |
| I2a robust estimator | 0.255 | 0.246 (0.044) | 0.063 | 0.217 | −0.002 | −0.02 to 0.29 | 0.56 | 3 | no gain |
| I2b diet-averaged | 0.242 | 0.234 (0.038) | 0.052 | **0.293** | **+0.074** | 0.01 to 0.35 | 0.87 | 4 | fails A2 and A3 |
| I2c robust + averaged | 0.224 | 0.217 (0.036) | 0.048 | 0.281 | +0.062 | 0.02 to 0.35 | 0.88 | 4 | fails A2 and A3 |
| I3 Huber on normal scores | 0.276 | 0.265 (0.037) | n/a | 0.220 | +0.001 | −0.05 to 0.30 | 0.21 | 4 | no gain |

- **Bias–variance:** for B0, train 0.274 ≈ CV 0.264, so it is high-bias and does not overfit. For B2 HGB, train 0.378 vs CV 0.276, a variance gap.
- **No candidate met the locked criteria.** The diet-averaged gain is post-hoc. Its mechanism is that the CIRC coefficient shrinks from −0.042 to about −0.004.

---

## 6. Why F8-OLS and F8+C were chosen for external validation

| Criterion | F8-OLS | F8+C-Ridge | FULL | F8+B+C | CatBoost-F8 |
|---|---|---|---|---|---|
| Held-out R² CD / MD | 0.066 / 0.070 | 0.080 / 0.096 | 0.104 / 0.108 | 0.094 / 0.107 | 0.085 / 0.079 |
| Condition number / max VIF | 2.48 / 1.88 | 3.76 / 3.78 | 37.3 / 233 | 4.05 / 3.80 | n/a |
| Coefficient sign stability | stable (0 flips) | stable | up to 12/28 flips | – | n/a (no coefficients) |
| Target-informed design | no | no (reference-derived) | yes (block B) | yes (block B) | no |
| Chosen? | **primary** | **secondary** | no (uninterpretable) | no (leaks PE pathway shifts) | no (not better than linear; opaque) |

Source: `roadmap_phase4c_diagnostics/decision_memo.md` §Recommendation.
