# Failures, limits, deviations and sensitivity

## Deviations (all logged; none chosen from external results)
1. **Which models to apply (clarification):** the plan said "if valid" but did not name a model. Fixed before any external metric was viewed: apply R0 + L2-logistic (simplest) and R2 + L2-logistic (tuned family).
2. **Permutation step crash:** with shuffled labels, the |diff| ≥ 1 rule could select no genes (NaN score). Fallback = top 50 by t when fewer than 10 genes pass. Real training sets always had 50 + 50 genes passing (verified), so real results are unaffected. The permutation null was rerun separately (`e01b_permutation_null.py`).
3. **`validity.csv` regenerated** by `e01c_validity.py` with the identical predeclared rule, because `e01` stopped before writing it.
4. **Optimisation after the user's request:** the first run (8 workers) was stopped and restarted with per-fold memoisation of signatures/PCA, precomputed R1 scores and 12 workers. No result had been produced by the stopped run.

## Failures and limits
- **Ceiling effect:** GAN vs chow is too easy. Every reasonable model reaches AUC 1.0, so model comparison rests on calibration (Brier, balanced accuracy at 0.5), not discrimination.
- **Two studies only:** two generalisation estimates, with no between-study variance estimate.
- **Bootstrap AUC CIs collapse at ceiling,** which understates uncertainty.
- **The permutation null has low resolving power** (partly correct labellings keep the diet signal; q95 = 1.0 in S1).
- **Purified vs chow control mismatch:** the classifiers learn "GAN vs grain chow". In GSE294180 they separate MD+Veh from CD+Veh strongly (positive control), but the features that move are diet-programme pathways (cholesterol, lipogenesis, PPARα) shared with fasting. They are not MASH-specific.
- **The R2 PE_in_MD result** (Holm p 0.0006) is **driven by cholesterol/PPARα pathway features.** It must not be reported as "PE induces a MASH state".
- **Fasting classifier:** external n = 6 with unknown sacrifice time. Demonstration only; **not applied** to the 19 mice.
- **PCA:** transfers poorly (study-specific axes). Kept as a failed comparator.
- **Depth-1/2 trees on R2 and PCA:** unstable across directions (AUC 0.63–0.80). Kept as a failed comparator.
- **No model can separate PE from lower depth/yield** in GSE294180. CD+Veh has an unusual clock state (scores correlate with clock_phase ρ ≈ 0.64–0.67).
- **GSVA not run** (not installed; population-dependent).

## Sensitivities run
- S1b: training on GSE316023 21 °C + 29 °C (32 mice). Same conclusions.
- Feature-contribution decomposition of the applied R2 score (see 04).
