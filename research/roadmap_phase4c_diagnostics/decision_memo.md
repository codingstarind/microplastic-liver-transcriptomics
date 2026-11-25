# Phase 4C diagnostic audit: decision memo

**Scope:** a diagnosis of the Phase 4/4B models from saved predictions and inputs. The audit spec and the four-model fix experiment were written and hashed before any Phase 4C computation (`predeclared_audit_spec.md`, `tables/spec_lock.json`). No earlier file was changed (checksums are test-verified). Phase 4B and these fold-reusing comparisons are **exploratory**; their intervals are optimistic.

## CONFIRMED

| # | Issue | Evidence (tables/figures) | Course method | Action | Consequence |
|---|---|---|---|---|---|
| 1 | **Severe multicollinearity in Phase 4B FULL only**, caused by raw vs shrunken reference effects | FAST_eff vs A_FAST_shr r = 0.998, MASH 0.979, CIRC 0.926; FULL max VIF **233** (4 VIF > 10); condition number 37 vs **2.5** for F8; FULL coefficients flip sign across folds (up to 12 of 28 under MD; F8: 0) (`collinearity_summary.csv`, `fig1`) | VIF, condition number, coefficient stability | Block A removed (F8+B+C). Original FULL kept as the historical result. | FULL coefficients are uninterpretable. Its predictions are fine: F8+B+C ridge (condition number 4.0) reaches R²_oos 0.094 (CD) / 0.107 (MD) vs FULL 0.104 / 0.108 |
| 2 | **Errors are heteroskedastic and dominated by noisy, low-count genes** | Spearman(\|residual\|, target SE) = 0.40 [0.36, 0.43] CD, 0.44 [0.39, 0.49] MD (module-block bootstrap). R²_oos by SE tercile: reliable 0.11/0.28, noisy 0.05/0.03. Lowest-baseMean quintile R² ≈ 0 (`subgroup_errors.csv`, `fig2`) | Residual diagnostics | Already down-weighted (1/(s² + τ²)); stratified reporting added | Most unexplained variance is target noise, not model failure |
| 3 | **The target itself is noisy** (noise ceiling) | Signal share of target variance = 1 − mean(SE²)/var: **0.36 (CD), 0.39 (MD)** from DESeq2 SEs; 0.05/0.12 from the 19-refit jackknife (overlapping refits share mice, so approximate) (`target_quality.csv`) | Bias–variance, reliability | Reported | Achievable R² is bounded well below 1. Observed 0.07–0.11 is low but not far from what this target allows |
| 4 | **Systematic under-prediction of lipogenesis and cholesterol-synthesis genes** | Mean residual: LIPOGENESIS_CORE −1.24 (CD) / −0.74 (MD); CHOL_SYNTHESIS_CORE −0.42 / −0.53 log2 (predefined sets) | Subgroup residual analysis | None; reported | The PE decrease in these programmes is stronger than any independent reference axis predicts. Block B's gain is this offset, re-learnt from PE (not reference similarity) |
| 5 | **Single-mouse influence** | Dropping CD_MP_2 changes the PE_in_CD target by RMS 0.19 log2 (Spearman with full 0.85) and lowers F8-OLS grouped-CV R² from 0.066 to 0.048; MD_MP_2 likewise for MD (`single_mouse_influence_on_target.csv`, `fig3`) | Influence analysis (LOO) | Reported; LOO ranges already part of robustness | Conclusions rest on 9–10 mice; one PE mouse matters |
| 6 | **Gene errors are dependent through shared mice** | Every gene's PE effect comes from the same contrast mice, so module grouping cannot remove this. Random-gene CV R² 0.069 ≈ grouped 0.066 | CV design | Documented | Gene-level CV says nothing about new mice or studies |
| 7 | **Phase 4B intervals are optimistic** | Blocks proposed after Phase 4 results; same outer splits reused | CV / model selection | Labelled exploratory | No confirmatory claim |

## PLAUSIBLE BUT UNPROVEN
- **A length- or yield-related artefact in PE effects:** Spearman(PE_in_MD, transcript length) = 0.14 vs 0.05 for Diet_in_Veh, 0.08 for PE_in_CD. This is consistent with library-quality differences between PE and vehicle libraries, but not proven.
- **DESeq2 SEs understate between-mouse uncertainty:** the jackknife signal share is much lower than the DESeq2-based one. Overlapping refits make this approximate.

## NOT SUPPORTED
- **Multicollinearity in the original F8:** max VIF 1.88, max |r| 0.663 (FAST_SE vs CIRC_SE). This matches the independent check.
- **Underfitting / over-compression:** predicted SD is only 0.26–0.31 × observed, but the calibration slope is 0.85–1.04 (CIs include 1). The compression is what a weak but calibrated predictor should look like.
- **Nonlinear residual structure:** residual decile means are flat, and trees never beat linear models on PE.
- **Module leakage inflating CV:** ICC of residuals by module is 0.004/0.021.
- **Pure sequencing depth reproducing the PE pattern:** ρ ≈ 0.00 with the depth-thinning effect. This **does not test yield or quality**, so it does not resolve the confound.

## Fix experiment (predeclared; 15 Phase 4 folds; `fix_experiment_summary.csv`)

| Model | R²_oos CD | R²_oos MD | ΔRMSE vs F8-OLS (CD / MD), descriptive | Condition number |
|---|---|---|---|---|
| F8-OLS | 0.066 | 0.070 | — | 2.5 |
| F8+C-OLS / Ridge | 0.080 / 0.080 | 0.096 / 0.096 | −0.0024 / −0.0039 (CIs include 0) | 3.8 |
| F8+B+C-Ridge | 0.094 | 0.107 | −0.0048 / −0.0059 (MD CI < 0) | 4.0 |
| *FULL-Ridge (historical)* | 0.104 | 0.108 | −0.0065 / −0.0059 | ≈ 37 |
| *CatBoost-F8 (historical)* | 0.085 | 0.079 | −0.0034 / −0.0011 | — |

- Ridge selected α = 1000 (the grid maximum), yet its predictions equal OLS. With condition numbers ≤ 4, regularisation has nothing to fix.
- The fix experiment reproduces Phase 4B exactly (max difference 8e-17).

## Recommendation: frozen models for Phase 5
- **Primary:** **F8-OLS.** Collinearity-free, stable coefficient signs, no target-informed design choices.
- **Secondary:** **F8+C-Ridge.** It adds the 4 reference-derived response-shape features, condition number 3.8, and gives the same interpretation with slightly better fit.
- **F8+B+C is not frozen:** block B re-learns GSE294180's own pathway shifts, which would contaminate a cross-study test.
- **Two distinct claims:**
  - **(a) Held-out genes (done):** within GSE294180, about 7–10% of PE gene-effect variance is predictable.
  - **(b) New mice/studies (not yet tested):** the **one untouched validation for Phase 5**:
    1. Fit the two frozen models on *all* GSE294180 genes, with coefficients frozen.
    2. Predict the plastic-effect vector of an **independent study** whose effects we have not examined in this framework: Hong-lab PE vs control (GSE261893, shared controls handled as one study) as primary, GSE245069 PS as secondary.
    3. Score with R² (relative to that study's own mean), Spearman and sign agreement, against a null and the single-axis baselines.
  - **Do not split genes randomly to fabricate a test set.**

## Cannot be repaired computationally
- PE is perfectly aligned with lower depth and yield in GSE294180. Thinning tests depth only.
- n = 4–5 mice per group, with influential single mice.
- The MASH reference stays conditional (Phase 3 c1 FAIL).
- Only one PE liver experiment with a diet factor exists.

**Stop here:** the diagnostics do not justify further model improvement. Proceed to the Phase 5 study-level validation and reporting.
