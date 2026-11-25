# Phase 4B comparison protocol (locked before any new feature is scored on PE)

**EXPLORATORY.** The feature families were chosen after the Phase 4 results were known. The same validation splits are reused, so any p-value or CI here is **not** untouched confirmation. Results are hypothesis-generating for the team.

## 1. Reused unchanged from Phase 4
- **Gene rows:** the 11,758 genes of `roadmap_phase4` `feature_table("primary")`.
- **Targets:** PE_in_CD and PE_in_MD (unchanged).
- **Training weights:** w = 1/(s² + τ²_train).
- **Outer splits:** `outer_splits()` (20 Ward co-expression modules; 5 folds × 3 repeats, seeds 101–103), imported from `roadmap_phase4/scripts/p4_data.py` without modification.
- **Inner tuning:** GridSearchCV + GroupKFold(4) on the training genes; scoring MSE.
- **Grids:** identical to `roadmap_phase4/model_search_spaces.csv` for ridge, elastic net, hgb and catboost.
- **Metrics:** RMSE (primary), MAE, R²_oos, Spearman, weighted R²_oos.
- **Nadeau–Bengio** corrected paired intervals.

## 2. New in Phase 4B
- **Representations:** the 10 listed in `feature_engineering_plan.md`. Every data-dependent step (τ² for block A, z-scaling for block C products, `StandardScaler`) is inside the fitted pipeline (`RepBuilder` transformer, fitted on training rows only).
- **Models:** `ols`, `ridge`, `elasticnet`, `hgb`, `catboost`.
  - **RF dropped:** slowest, and statistically tied with HGB/CatBoost in Phase 4.
  - **SVR not used:** no kernel-structure rationale for 4–28 heterogeneous tabular features; O(n²) cost at about 9,400 training genes.
  - **No Naive Bayes, classifiers or neural networks.**
- **Reproducibility check:** F8 is rerun through the new code. Its per-split metrics must equal the Phase 4 metrics to 1e-9 for ols/ridge/elasticnet/hgb/catboost; any mismatch is reported.

## 3. Comparisons and rules (fixed now)
- **Primary comparison:** each representation vs **F8** (and vs F4) within the same model, on the same 15 splits.
- **"Improvement over F8"** requires all of:
  - NB-corrected 95% CI of ΔRMSE < 0;
  - ≥ 12/15 splits better;
  - **practical size** ΔRMSE ≤ −0.005 log2 **and** ΔR²_oos ≥ +0.01.
- **"Improvement over the best baseline":** the improved representation × model must also beat **OLS-F8** (the Phase 4 recommended model) by the same rule.
- **Shortlist:** any representation × model meeting the F8 rule for either target. Only shortlisted ones get:
  - the 8 Phase 4 sensitivity variants (adapted: the replacement references swap the corresponding raw columns before blocks are built), with 1 repeat;
  - 19 LOO targets, 1 repeat;
  - a target-permutation control.

  **If nothing is shortlisted, sensitivity runs are not performed** (reported).
- **Feature-block ablations:** FULL vs FULL−A, FULL−B, FULL−C, compared within the same model. A block is "useful" if removing it significantly worsens RMSE (ΔRMSE > 0 with CI > 0).
- **Controls always run** (1 repeat; models ridge and hgb):
  - (i) target permuted across genes (seed 7, as Phase 4), representation FULL;
  - (ii) block-B membership rows shuffled across genes (seed 13), representation F8+B.
- **Depth caveat:** Phase 4's depth-thinning control (R²_oos ≈ 0) applies. No representation can address yield/quality confounding.

## 4. Parallelism and seeds
- joblib/loky, 19 workers, single-threaded tasks.
- Model seeds 20261122; split seeds as in Phase 4.
