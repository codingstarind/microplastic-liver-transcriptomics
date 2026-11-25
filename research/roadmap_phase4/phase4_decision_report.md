# Phase 4 decision report: model comparison and validation

**Date:** 2026-11-22. Protocol, grids and decision rules were locked before any fit (`results/protocol_lock.json`; test-checked). **Not started:** Phase 5.

Labels used throughout:
- **[obs]** observed number;
- **[ev]** statistical evidence under the locked rules;
- **[int]** interpretation needing team review.

## 1. What was predicted
- **Task A (primary):** rows are 11,758 genes (main universe with all features). The target is the GSE294180 PE log2FC (PE_in_CD, PE_in_MD), modelled separately.
- **Features:** the log2FC and SE of four independent references:
  - MASH (GSE316023; **conditional**: the Phase 3 positive control formally failed c1);
  - fasting (ZT-averaged);
  - PPARα agonist (WT 20v0);
  - time of day (T1).
- **Validation:** nested CV with whole co-expression modules held out (5 outer folds × 3 repeats; inner GroupKFold 4).
- **Scope of any score:** a held-out-gene score is **not** accuracy on held-out mice or laboratories.
- **Task B (demonstration only):** a GAN-vs-chow mouse classifier (train GSE316023, n = 32; external GSE246328, n = 15).

## 2. Task A results (out-of-sample, same held-out genes for every model)

R²_oos is 1 − SSE/SSE(training mean). Mean of 15 outer splits.

| Model | PE_in_CD R²_oos | PE_in_MD R²_oos | Positive control Diet_in_Veh |
|---|---|---|---|
| null | 0 | 0 | 0 |
| ols_MASH | 0.031 | −0.002 | **0.206** |
| ols_FAST | 0.022 | **0.045** | 0.012 |
| ols_PPARA | 0.000 | −0.001 | 0.001 |
| ols_CIRC | 0.022 | 0.005 | 0.040 |
| **ols_all** (linear, 8 features) | **0.066** | **0.070** | 0.233 |
| ridge | 0.066 | 0.070 | 0.233 |
| elastic net | 0.063 | 0.068 | 0.232 |
| tree (shallow) | 0.038 | 0.053 | 0.179 |
| random forest | 0.073 | 0.081 | 0.247 |
| hist. gradient boosting | 0.078 | 0.079 | 0.285 |
| CatBoost | 0.086 | 0.079 | **0.306** |

RMSE in log2FC units, null → best model:
- PE_in_CD: 0.326 → 0.312;
- PE_in_MD: 0.309 → 0.296;
- Diet_in_Veh: 0.490 → 0.408.

**Paired comparisons** (ΔRMSE, Nadeau–Bengio corrected 95% CI, 15 splits; `results/taskA/paired_comparisons.csv`):

| Comparison | PE_in_CD | PE_in_MD | Diet_in_Veh |
|---|---|---|---|
| ols_all vs null | −0.0107 [−0.017, −0.004] **better** | −0.0112 [−0.025, +0.003] **not distinguishable** | −0.061 [−0.080, −0.043] better |
| ols_all vs best single axis | vs MASH: −0.0058 [−0.011, −0.0003] **better** | vs FAST: −0.0040 [−0.012, +0.004] not distinguishable | vs MASH: −0.007 [−0.020, +0.006] not distinguishable |
| RF / HGB / CatBoost vs null | better (all three) | better (all three) | better |
| **RF / HGB / CatBoost vs ols_all** | −0.001 / −0.002 / −0.004, **all CIs include 0** | −0.001 / −0.001 / −0.001, **all CIs include 0** | HGB −0.015 [−0.026, −0.004], CatBoost −0.022 [−0.037, −0.006] **better**; RF not |
| Ridge / Elastic Net vs ols_all | ≈ 0 | ≈ 0 | ≈ 0 |

**Robustness** (`results/taskA/robustness_counts.csv`, `variants_r2_oos.csv`):
- **ols_all beats null** in **8/8 sensitivity variants and 19/19 leave-one-mouse-out targets**, for both PE targets.
  - Variants: MASH2, clean fasting subset, PPARα 10v0, CIRC T2, strict universe, alternative muscle universe, chromosome grouping, effects-only features.
  - LOO R²_oos ranges: 0.048–0.082 (CD) and 0.053–0.109 (MD).
- **HGB vs ols_all:** lower mean RMSE in 8/8 variants and 19/19 LOO (CD), and in 6/8 and 18/19 (MD). The **primary paired CI includes 0**, so under the locked rule this is **not a demonstrated improvement**. It is a consistent but tiny difference (≈ 0.001–0.002 log2 RMSE).
- **Single-axis fragility:**
  - ols_MASH on PE_in_CD ranges −0.016 to +0.051 across LOO targets; on PE_in_MD it is below null in 12/19.
  - ols_FAST is stable (CD 0.010–0.045; MD 0.034–0.056).
- **Effects-only features (V8)** keep most of the gain (CD 0.064 vs 0.066; MD 0.064 vs 0.070). Most of the gain comes from the reference effects, not from the SE features (expression-level proxies).

**Negative controls** (`negative_controls_r2_oos.csv`):
- targets permuted across genes: R²_oos ≈ 0 for every model (−0.002 to +0.000);
- depth-only thinned targets: R²_oos ≈ 0 (−0.005 to +0.001); thinning itself produced 0 DEGs.
- There is no sign of leakage or of the pipeline manufacturing signal. **Thinning does not simulate yield or library-quality effects**, so it does **not** clear the PE–depth confound.

**Linear coefficients** (per SD of feature; sign stable in all 15 splits; `coefficient_stability.csv`):

| Feature | PE_in_CD | PE_in_MD |
|---|---|---|
| MASH effect | +0.060 | +0.040 |
| FAST effect | +0.035 | +0.058 |
| CIRC T1 | −0.041 | +0.031 |
| PPARα effect | −0.012 | −0.021 |

Permutation importance agrees: MASH dominates for CD, fasting for MD. **Caveats:**
- the coefficients are attenuated by noise in the references (errors in variables);
- the reference scales come from different technologies;
- importance is diluted between correlated effect and SE features.
- **They describe the prediction, not biology.**

## 3. Task B (demonstration; no claim)
- Internal CV: balanced accuracy is 1.00 for logistic, Gaussian NB, RF and HGB, and 0.94 for both trees.
- External GSE246328: logistic 1.00, HGB 1.00, RF 0.94, NB 0.94, Gini tree = entropy tree 0.875.
- **Ceiling effect, as predicted at the gate.** GAN vs chow is near-perfectly separable, so these numbers cannot rank classifiers. Gini and entropy gave identical trees, so they are redundant.
- **Deviation D2:** the HGB classifier was degenerate with the default `min_samples_leaf = 20` at n = 24 (constant 0.5 predictions). It was fixed to 3; the original run is preserved.
- **Descriptive:** GSE294180 MD+Veh mice got a higher predicted "GAN-state" probability than CD+Veh mice from logistic (0.71 vs 0.21), RF (0.79 vs 0.40), NB and HGB. The trees predicted 0 for all.
  - n = 9, with a purified rather than chow control. **This is not evidence.**
- No PE classifier and no cross-laboratory multi-state classifier were trained.

## 4. Answers

**Which model beats which baseline, by how much?**
- **[ev]** PE_in_CD:
  - linear multi-axis (ols_all/ridge) beats null (ΔRMSE −0.011) and the best single axis, MASH (−0.006);
  - tree ensembles also beat null;
  - **no model beats the linear multi-axis model.**
- **[ev]** PE_in_MD:
  - RF, HGB and CatBoost beat null (ΔRMSE ≈ −0.012);
  - the linear model's improvement over null is not distinguishable under the corrected interval, although it is positive in 8/8 variants and 19/19 LOO;
  - nothing beats the fasting single axis under the locked rule.
- **[ev]** Positive control:
  - boosting beats linear (ΔRMSE −0.015 to −0.022);
  - this shows the flexible models *can* win when a strong structured signal exists.
  - For the PE targets they do not.

**Is any improvement robust?**
- The **presence of predictable signal** (ols_all > null) is robust across all variants and LOO refits.
- **Superiority of complex models over linear is not demonstrated.**
- **The recommended model for Phase 5 descriptions is the linear multi-axis model:** simplest, statistically tied with the best, stable coefficient signs.

**How large is it?**
- [obs] About **7% of the held-out variance** of the PE gene effects is predictable from the four independent reference effects (8% with boosting).
- That is small. It compares with about 23–31% for our own diet effect against the same references.

**What this DOES tell us**
- [int] Within this one experiment, the PE-associated gene-effect vector has a small, reproducible, partly predictable component aligned with independent biological-state references.
  - Under CD it is mostly MASH-like plus fasting-like (with an opposite time-of-day sign).
  - Under MD it is mostly fasting-like.
- The baseline ranking from Phase 3 (fasting stable, MASH fragile) is broadly preserved, now with out-of-sample metrics and fair baselines.

**What this does NOT tell us**
- That PE **causes** MASH-like, fasting-like or any other state. PE is perfectly aligned with lower depth and yield in GSE294180, and no model or reference removes that.
- Anything about other mice or laboratories. Held-out **genes**, not held-out animals: rows are correlated, and the effective sample size is far below 11,758.
- That the MASH axis is validated. It remains **conditional** (Phase 3 positive control: FAIL on c1). MASH-axis contributions inherit that status and its lipogenesis caveat.
- That coefficients or importances identify mechanism.

## 5. Deviations from the locked protocol
- **D1 (compute):** variants, negative controls and LOO used 1 repeat (5 splits) instead of 3. Random forest took about 4–11 min per tuned task. This was decided after a one-split smoke test of all models, before any full results; it applies equally to every model.
- **D2 (Task B):** HGB classifier `min_samples_leaf = 3` (the default was degenerate at n = 24). The original run is kept in `results/taskB/original_hgb_min_leaf20/`.
- **D3 (analysis bug):** pandas parsed the model name "null" as NA. Fixed in the analysis and tests; no numbers changed.
- The Ward-linkage modules (vs average linkage in the Phase 3 memo) were set in the locked protocol before any fit.

## 6. Decisions for the team
1. Accept the **linear multi-axis model** as the Phase 5 descriptive model (recommended), or keep boosting as a sensitivity comparator.
2. Decide how to present the MASH contribution, given its conditional status and its LOO fragility as a single axis.
3. Decide whether a yield/library-quality stress test (beyond depth thinning) is worth building before any write-up.
4. Confirm that Task B stays a demonstration and is not reported as a finding.
