# Phase 4C predeclared spec (written before any Phase 4C computation; SHA-256 in tables/spec_lock.json)

This is a **diagnostic audit**. Fold-level predictions and inputs are taken from Phases 4/4B (read only); no earlier file is modified.

## Diagnostics (descriptive; no model selection)
1. **Collinearity**, for the matrices F8, A4, B12, C4, F8+C, F8+B+C and FULL (Phase 4B `RepBuilder`):
   - max pairwise |r|;
   - VIF of each column, computed as the diagonal of the inverse correlation matrix;
   - condition number and singular values of the standardised matrix for each of the 15 outer training folds (scaler fitted on the training fold);
   - numerical rank (tolerance 1e-10 × largest singular value);
   - OLS/ridge coefficient stability across the 15 outer folds (CV of coefficient, sign flips).

   Pairs with |r| ≥ 0.9 are named as the collinearity source.
2. **Prediction error**, from Phase 4 saved out-of-fold predictions (OLS-F8 and CatBoost-F8; each gene averaged over its 3 repeats) and from the Phase 4C fix-experiment predictions. Reported:
   - predicted vs observed;
   - residual vs fitted;
   - residual vs target SE;
   - residual vs GSE294180 baseMean.

   **Subgroups:**
   - predefined: the 12 block-B pathway sets and "no set";
   - module (20);
   - baseMean quintiles;
   - target-SE terciles ("reliable" = lowest SE tercile; noisy = highest). SE terciles are not selected on PE significance;
   - |observed effect| quintiles, labelled **diagnostic-only (post hoc)**.

   **Uncertainty:** module-block bootstrap (1,000 resamples of the 20 modules, seed 20261122).

   **Heteroskedasticity and nonlinearity:**
   - Spearman ρ of |residual| with target SE and with baseMean;
   - binned mean residual by fitted-value decile, with block-bootstrap CI.
3. **Target quality:**
   - 19 existing LOO refits; per-gene jackknife SE;
   - ratio of jackknife SE to DESeq2 lfcSE;
   - reliability = 1 − mean(SE²)/var(log2FC), using both SEs. This is reported as an approximate **noise ceiling** on R², because the 19 refits overlap and share mice;
   - per-mouse influence = RMS change of the target vector when that mouse is dropped.

   **Confound probes (descriptive only):**
   - Spearman ρ of the PE effect with the Phase 4 depth-thinning effect;
   - Spearman ρ with GENCODE transcript length (length/yield-type bias), compared with the same correlation for Diet_in_Veh.
4. **Leakage/dependence:**
   - code-level checklist;
   - module intra-class correlation (ICC) of the target and of out-of-fold residuals;
   - a **diagnostic-only** random-gene 5-fold CV for OLS-F8, to quantify optimism when modules are ignored. It is not a benchmark.

## Fix experiment (only models; fixed now)
- **Models:** F8-OLS, F8+C-OLS, F8+C-Ridge, F8+B+C-Ridge.
- **Folds:** the exact Phase 4 outer folds (15).
- **Ridge alpha grid:** {0.01, 0.1, 1, 10, 100, 1000} (the Phase 4 grid), tuned by inner GroupKFold(4) on training genes only.
- **Weights:** as Phase 4.
- **Outputs:** out-of-fold predictions, paired ΔRMSE and ΔR²_oos vs F8-OLS (Nadeau–Bengio interval reported as **descriptive**, because the folds are reused), per-fold coefficients, condition numbers.
- **Historical comparators** (not rerun): Phase 4B FULL-OLS, FULL-Ridge, CatBoost-F8.
- **No other models, features or grids.**
