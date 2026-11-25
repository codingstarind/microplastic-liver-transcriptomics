# Predeclared model and evaluation plan (locked; hashes in results/plan_lock.json)

## B1: MASH-state classifier (primary)
- **Outer splits:**
  - S1: train GSE316023 TA (8 v 8) → test GSE246328 (7 v 8);
  - S2: train GSE246328 → test GSE316023 TA.
  - Sensitivity S1b: train GSE316023 TA + TN (16 v 16) → test GSE246328.
- **Inside the training study only:** R0 signature genes, PCA loadings, scaler, hyperparameters.
- **Inner CV:** RepeatedStratifiedKFold(4 splits, 3 repeats, seed 20261122); scoring = neg_brier_score.
- **Models and bounded grids:**

  | Model | Grid |
  |---|---|
  | majority (dummy) | — |
  | R0 rule (score > midpoint of training class means) | — |
  | logistic, weak penalty | C = 100 |
  | L2-logistic | C ∈ {0.01, 0.1, 1, 10} |
  | Gaussian NB | var_smoothing ∈ {1e-9, 1e-3} |
  | tree | depth ∈ {1, 2}, criterion ∈ {gini, entropy} |
  | random forest | 300 trees, max_depth ∈ {1, 2}, min_samples_leaf ∈ {1, 3} |

  These run for R0, R1, R2 and PCA. The rule runs only on R0.
- **"Before / after":**
  - baseline = R0 rule;
  - after feature engineering = R1/R2 with logistic at default C = 1;
  - after tuning = the inner-CV-selected model.
- **External metrics:**
  - ROC AUC (primary);
  - balanced accuracy, sensitivity, specificity at 0.5;
  - Brier and log loss;
  - stratified bootstrap 95% CI (2,000 resamples of test mice).
- **Controls:**
  - 200 training-label permutations (L2-logistic at C = 1 on R2), giving the external AUC null;
  - Spearman correlation of external predictions with library size.
- **Freeze rule:** the external study is scored once per split. No retuning after viewing.

## Validity criteria for applying a model to GSE294180 (fixed now)
A representation × model is "valid" only if, in **both** S1 and S2:
- external AUC ≥ 0.80 **and** bootstrap CI lower bound > 0.5;
- Brier < 0.25.

## C: application to GSE294180 (only if valid)
- **Final model:** retrained on both MASH studies' eligible mice (31). The R0 signature uses genes passing the rule in **both** studies (no GSE294180 data).
- **Score:** decision function (log-odds), not a calibrated probability.
- **Model:** score ~ diet × PE; contrasts PE_in_CD, PE_in_MD, Interaction, Diet_in_Veh (positive control: > 0 with CI excluding 0); Holm over the 3 PE contrasts; 19 leave-one-out refits.
- **Comparison:** against the Phase 4D MASH singscore on the same mice.

## B2: fasting (demonstration)
- Train on GSE107787 clean subset (4 FAST v 6 FED; R0 signature built within it; R0 rule + L2-logistic at C = 1).
- Test on GSE225695 (3 v 3).
- Report AUC and the confusion matrix only; no tuning comparisons.
