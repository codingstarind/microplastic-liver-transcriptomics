# Roadmap Phase 4C: diagnostic audit and model freezing (2026-11-22)

## Main sources (`research/roadmap_phase4c_diagnostics/`)
- `predeclared_audit_spec.md` (hash in `tables/spec_lock.json`)
- `decision_memo.md`
- `tables/`: `collinearity_summary.csv`, `vif_by_feature.csv`, `coefficient_stability_*.csv`, `subgroup_errors.csv`, `target_quality.csv`, `single_mouse_influence_on_target.csv`, `confound_probes.csv`, `dependence_audit.csv`, `leakage_audit.csv`, `fix_experiment_summary.csv`, `loo_model_r2_F8_ols.csv`
- `figures/fig1_collinearity.png`, `fig2_residual_diagnostics_F8_ols.png`, `fig3_subgroups_and_influence.png`
- Scripts: `c00_common.py`, `c01_collinearity.py`, `c02_errors_target_leakage.py`, `c03_fix_experiment.py`

## 1. Question
- Why are the gene-level fits low?
- Are FULL's gains trustworthy?
- Which models should be frozen for the external test?

## 2. Methods
All diagnostics were predeclared.

| Diagnostic | Method |
|---|---|
| Collinearity | max \|r\|; VIF = diagonal of the inverse correlation matrix; condition number and singular values per training fold; numerical rank; coefficient sign flips across 15 folds |
| Prediction error | predicted vs observed; residual vs fitted, vs target SE, vs baseMean. Subgroups: 12 pathway sets, 20 modules, baseMean quintiles, SE terciles. Module-block bootstrap (1,000 resamples, seed 20261122). |
| Target quality | signal share = 1 − mean(SE²)/var(y), from DESeq2 SEs and from the 19-refit jackknife; per-mouse influence = RMS change of the target when that mouse is dropped |
| Confound probes | ρ with the depth-thinning effect; ρ with GENCODE transcript length |
| Dependence | module ICC of target and residuals; random-gene 5-fold CV as a diagnostic only |
| Fix experiment | F8-OLS, F8+C-OLS, F8+C-Ridge, F8+B+C-Ridge on the exact 15 Phase 4 folds; ridge α grid {0.01, …, 1000} |

## 3. Confirmed findings

**1. Severe collinearity in FULL only.**
- Raw vs shrunken effects: FAST r = 0.998; MASH 0.979; CIRC 0.926.
- FULL max VIF **233**, with 4 VIFs > 10.
- Condition number 37 for FULL vs **2.5** for F8.
- FULL coefficients flip sign across folds in up to 12 of 28 features (MD). F8 coefficients never flip.

**2. Heteroskedastic errors, dominated by noisy genes.**
- Spearman(|residual|, target SE): 0.40 [0.36, 0.43] (CD); 0.44 [0.39, 0.49] (MD).
- R²_oos in the reliable SE tercile: 0.11 / 0.28 (CD / MD).
- R²_oos in the noisy SE tercile: 0.05 / 0.03.
- The lowest-baseMean quintile has R² ≈ 0.

**3. The target itself is noisy (noise ceiling).**
- Signal share from DESeq2 SEs: 0.36 (CD), 0.39 (MD).
- From the jackknife: 0.05 / 0.12. This is approximate, because the refits overlap.

**4. Lipogenesis and cholesterol-synthesis genes are systematically under-predicted.**
- Mean residual, LIPOGENESIS_CORE: −1.24 (CD) / −0.74 (MD).
- Mean residual, CHOL_SYNTHESIS_CORE: −0.42 / −0.53 log2.

**5. Single mice are influential.**
- Dropping CD_MP_2 changes the PE_in_CD target by RMS 0.19 log2 (Spearman with the full target 0.85).
- The same drop lowers F8-OLS R² from 0.066 to 0.048.

**6. Gene errors are dependent through shared mice.** Random-gene CV R² 0.069 is about the same as grouped CV 0.066.

**7. Phase 4B intervals are optimistic** (folds reused after features were chosen).

## 4. Not supported
| Hypothesis | Evidence against it |
|---|---|
| Collinearity in F8 | max VIF 1.88 |
| Underfitting | calibration slope 0.85–1.04 |
| Nonlinear residual structure | residual decile means are flat |
| Module leakage | residual ICC 0.004 / 0.021 |
| Pure depth reproducing the PE pattern | ρ ≈ 0. This does **not** test yield or quality. |

## 5. Plausible but unproven
- **A length or yield artefact:** ρ(PE_in_MD, transcript length) = 0.14, vs 0.05 for Diet_in_Veh and 0.08 for PE_in_CD.
- **DESeq2 SEs may understate between-mouse uncertainty:** the jackknife signal share is much lower than the DESeq2-based one.

## 6. Fix experiment (`tables/fix_experiment_summary.csv`)

| Model | R²_oos CD | R²_oos MD | Condition number |
|---|---|---|---|
| F8-OLS | 0.066 | 0.070 | 2.5 |
| F8+C-OLS / Ridge | 0.080 | 0.096 | 3.8 |
| F8+B+C-Ridge | 0.094 | 0.107 | 4.0 |
| FULL-Ridge (historical) | 0.104 | 0.108 | ≈ 37 |
| CatBoost-F8 (historical) | 0.085 | 0.079 | — |

- Ridge selected α = 1000 (the grid maximum), but its predictions equal OLS: with condition numbers ≤ 4 there is nothing to regularise.
- The fix experiment reproduces Phase 4B exactly (maximum difference 8e-17).

## 7. Decision
- **Frozen for Phase 5:**
  - **F8-OLS** (primary);
  - **F8+C-Ridge** (secondary).
- **Not frozen:**
  - F8+B+C, because block B re-learns PE pathway shifts;
  - FULL, because it is collinear.

> **Key message.** FULL predicted better because its extra columns carried information, but its individual coefficients were unidentifiable (VIF 233, sign flips). Good prediction and interpretable coefficients are different properties.

## 8. Transition
Gene-level CV cannot speak to new mice. The next step (Phase 4D) was a per-mouse estimand: frozen signature scores per mouse, analysed with the 2 × 2 design.
