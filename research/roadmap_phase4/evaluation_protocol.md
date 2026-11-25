# Phase 4 evaluation protocol (LOCKED before any fit)

Seeds, folds, features, models, grids and metrics are frozen here. Hashes are in `results/protocol_lock.json`. Changes after locking are logged in `phase4_decision_report.md` § Deviations; this file is not edited.

## 1. Data (task A)
- `scripts/p4_data.py` builds the gene table from:
  - the Phase 3 `reference_effects/` (read only);
  - GSE294180 `results/deseq2/provisional_deposited/results_{PE_in_CD,PE_in_MD,Diet_in_Veh}.csv.gz`;
  - Phase 2 universes.
- Rows with any missing primary feature or target are dropped, and the count is reported.

## 2. Grouping (no PE outcome used)
- **Primary groups: 20 co-expression modules.**
  - Agglomerative clustering (**Ward** linkage, Euclidean distance on per-gene z-scored log2(RPKM + 1)) of the 54 GSE117134 liver samples, `n_clusters = 20`.
  - Genes with zero variance in GSE117134 form their own extra group.
  - Deviation from the Phase 3 memo: Ward instead of average linkage, because average linkage gives one giant cluster. This was decided before any fit.
  - Every gene is in exactly one module, so overlapping gene-set membership cannot leak a gene between training and test.
  - Known residual dependence: co-regulation across modules, and the CIRC feature coming from the same dataset as the modules.
- **Sensitivity grouping: chromosome.** GENCODE vM35 chromosome of each gene, about 21 groups.

## 3. Nested cross-validation
- **Outer loop:**
  - groups randomly assigned to **5 folds**, with **3 repeats** (seeds 101, 102, 103), giving 15 outer splits per (task, variant);
  - the same splits are used for every model, so all models are scored on **exactly the same held-out genes**.
- **Inner loop:**
  - `GridSearchCV` with `GroupKFold(n_splits = 4)` on the outer-training genes only;
  - scoring `neg_mean_squared_error`;
  - `refit = True`.
- **Pipelines:** every scaler and model, and the τ² used for weights, are fitted on training rows only.
  - Weights: w = 1/(s² + τ²_train), passed to the estimator as `sample_weight`.
  - `StandardScaler` sits inside the pipeline.
  - No feature selection or imputation.
- **Leakage tests:** `tests/test_phase4.py` checks that no gene is in both training and test in any split; that groups never straddle training and test; that inner folds are subsets of outer training; that no forbidden feature is present; and that test rows never enter fitted scalers.

## 4. Models (task A)
- **Baselines:**
  - `null` (weighted training mean);
  - 4 single-axis OLS models (`ols_MASH`, `ols_FAST`, `ols_PPARA`, `ols_CIRC`; the effect only, no SE);
  - `ols_all` (all 8 features).
- **Tuned:**
  - `ridge`;
  - `elasticnet` (Lasso is its l1_ratio = 1 corner);
  - `tree` (shallow CART);
  - `rf` (random forest; `max_features = 1.0` is plain bagging);
  - `hgb` (histogram gradient boosting);
  - `catboost`.
- **Not run:** AdaBoost (see the justification file); a pooled CD+MD model; neural networks.
- Grids are in `model_search_spaces.csv`. At most 16 configurations per model.

## 5. Metrics (outer held-out genes)
- **Primary:** RMSE (log2FC units).
- **Also reported:**
  - MAE;
  - **R²_oos** = 1 − SSE_model / SSE_null, where null = training mean (can be negative);
  - sklearn R² relative to the held-out mean;
  - Spearman ρ;
  - sign agreement (all held-out genes, and target DEGs with padj < 0.05 and |log2FC| ≥ 0.38; evaluation only);
  - weighted R²_oos (using the same weight formula).

## 6. Comparisons and decision rules
- Paired per-split differences ΔRMSE (model − comparator) over the 15 outer splits.
- **Uncertainty:** Nadeau–Bengio corrected resampled t interval (n_test/n_train correction; 14 df).
- **"Better than comparator"** requires the corrected 95% CI of ΔRMSE < 0 **and** the improvement to hold in ≥ 12/15 splits.
- **Comparators:** `null`; the best single-axis baseline, **defined per target as the single-axis model with the lowest mean outer RMSE** (a baseline choice using the same outer scores; reported as such); `ols_all`.
- **"Robust improvement"** additionally requires the same direction (mean ΔRMSE < 0 vs `ols_all`) in at least 6 of 8 sensitivity variants **and** in ≥ 15/19 LOO target refits.
- **The overall winner is not selected by its single best score.** The report shows a rank table, the complexity-vs-gain plot and paired intervals.

## 7. Sensitivity variants
Run with the reduced model set: null, single-axis models, ols_all, ridge, rf, hgb.

| Variant | Change |
|---|---|
| V1 | MASH2 replaces MASH1 |
| V2 | FAST_clean replaces FAST |
| V3 | PPARA 10v0 (muscle-matched) replaces 20v0 |
| V4 | CIRC_T2 replaces T1 |
| V5 | strict universe (8,360) |
| V6 | alternative muscle-restricted universe (12,010; references refitted on that universe by `scripts/p4_alt_universe_refs.R`) |
| V7 | chromosome grouping |
| V8 | effects-only features (4) |

## 8. Negative and stress controls
- **N1 permuted target:** y shuffled across genes (seed 7). Expected R²_oos ≤ 0.
- **N2 depth-only target:**
  - Vehicle libraries of each diet binomially thinned to the median PE/vehicle total-count ratio of that diet (seed 11), then a paired DESeq2 fit `~ mouse + thinned`.
  - The "thinned vs original" log2FC is the target.
  - This tests whether references "predict" a pure-depth effect vector.
  - Yield/quality effects are **not** simulated by thinning; that is a stated limitation.
- **P1 positive control target:** Diet_in_Veh. The MASH-axis models are expected to beat null.
- **LOO:** each of the 19 existing leave-one-mouse-out target refits; models null, single-axis, ols_all, ridge and hgb; outer splits of repeat 1; inner tuning as above.

## 9. Task B (demonstration)
- **Features:** per-sample log2(CPM + 1) on genes common to GSE316023, GSE246328 and GSE294180, then per-sample z-scoring across genes (sample-wise, so no leakage).
- **Inside the pipeline:**
  - variance filter to the top k genes (k ∈ {200, 1000}), fitted on training samples only;
  - `StandardScaler`.
- **Training set:** GSE316023 wk23 TA + TN (16 GAN, 16 chow). **Inner CV:** stratified 4-fold grouped by mouse.
- **Evaluation:**
  - **outer (internal):** repeated stratified 4-fold, 5 repeats;
  - **external:** GSE246328 Vehicle 24 w vs Chow 24 w (7 v 8);
  - **descriptive:** GSE294180 vehicle mice (CD 4, MD 5), reporting the predicted probability of the "GAN" state only.
- **Models:** dummy (majority), L2 logistic regression, Gaussian NB, tree (gini), tree (entropy), random forest, HistGradientBoostingClassifier.
- **Metrics:** balanced accuracy, ROC-AUC, PR-AUC, sensitivity, specificity, Brier score, confusion matrices. **No claim.**

## 10. Seeds and computing
- Global seed 20261122.
- Parallelism: `joblib` with loky (≤ 18 workers); each task single-threaded (`n_jobs = 1`, CatBoost `thread_count = 1`), so results are deterministic per seed.
