# Phase 4D decision report: per-mouse signature scores

**Representation decision:** CONDITIONAL GO (`representation_go_nogo.md`). Signatures were built only from independent references and frozen with hashes before any mouse was scored; the plan was locked too (`results/plan_lock.json`). **This is a different estimand** (score shifts per mouse in the 2×2 design) from Phase 4 gene-level R². The two are not numerically comparable.

## Results (singscore on the TPM-proxy; n = 19 mice; 95% t-CI, df 15)

| Score | Contrast | Estimate [95% CI] | Std. effect | p | Holm p | Exact within-diet permutation p | Leave-one-out sign stable |
|---|---|---|---|---|---|---|---|
| **MASH** | **Diet_in_Veh (positive control)** | **+0.215 [0.158, 0.273]** | 5.3 | < 1e-4 | single | — | 19/19 |
| MASH | PE_in_CD | +0.025 [−0.033, 0.082] | 0.6 | 0.37 | 1.00 | 0.14 | 19/19 |
| MASH | PE_in_MD | +0.054 [−0.0004, 0.108] | 1.3 | 0.052 | 0.21 | 0.15 | 19/19 |
| MASH | Interaction | +0.029 [−0.050, 0.109] | 0.7 | 0.44 | 1.00 | 0.49 | 19/19 |
| **FAST** | **PE_in_CD** | **+0.065 [0.022, 0.108]** | 2.1 | 0.006 | **0.036** | 0.040 | 19/19 |
| **FAST** | **PE_in_MD** | **+0.056 [0.015, 0.096]** | 1.8 | 0.011 | 0.056 | 0.024 | 19/19 |
| FAST | Interaction | −0.010 [−0.069, 0.050] | −0.3 | 0.74 | 1.00 | 0.82 | 17/19 |
| FAST | Diet_in_Veh | +0.130 [0.087, 0.173] | 4.3 | < 1e-4 | descriptive | — | 19/19 |
| PPARα (secondary) | PE_in_CD / PE_in_MD | +0.031 [−0.009, 0.070] / +0.006 [−0.031, 0.043] | | 0.12 / 0.72 | 0.36 / 0.72 | 0.19 / 0.66 | 19/19 / 17/19 |

**Validation conditions:**
- **The positive control PASSES:** the independent MASH score separates our MD+Veh from CD+Veh mice. This is an *additional* concordance assay. **It does not reverse the locked Phase 3 c1 FAIL.**
- **Null gene sets:** 1,000 abundance-matched random signed sets per signature. The observed PE and diet shifts lie outside the 95% null range for both primary scores (for example FAST PE_in_CD +0.065 vs null −0.011 to +0.011). This shows signature specificity *across genes*. It is not extra replication, because the mouse is still the unit.
- **Permutation caveat:** the smallest attainable p is 2/126 = 0.016 (CD), and PE is aligned with lower depth. Relabelling assumes exchangeability that library differences may violate.

## Sensitivities (`results/contrast_estimates_sensitivity.csv`, `figures/fig4`)

**FAST PE shift: robust.**

| Variant | PE_in_CD | PE_in_MD |
|---|---|---|
| Original (n50) | +0.065 | +0.056 |
| n25 | +0.094 | +0.077 |
| Alternative muscle rule | identical | identical |
| log-CPM / VST input | +0.067 | +0.060 |
| 14 MASH-shared genes removed | +0.054 (CI > 0) | +0.031 (CI includes 0) |

**MASH PE shift: fragile.**
- n25 gives PE_in_MD +0.092 (CI > 0).
- **Removing the cholesterol-synthesis genes** (which dominate the MASH-down list in both MASH studies) gives PE_in_CD +0.003 and PE_in_MD +0.021, both CIs including 0.
- Removing the MASH–fasting shared genes gives +0.010 and +0.030.
- **The apparent MASH-like PE shift is carried by genes shared with the fasting/cholesterol direction.**

**Input:** log-CPM and VST give *identical* scores. Both preserve within-sample ranks, so singscore is invariant to them; only length normalisation changes ranks.

GSVA was not run (not installed).

## Technical associations (`results/qc_associations.csv`, `figures/fig2`)
- **Depth:** PE separates completely (log2 depth PE 23.69–23.98 vs vehicle 24.02–24.19). Yield (quality_proxy) overlaps slightly.
- **Scores vs depth / yield / muscle:** |ρ| ≤ 0.30, none significant.
- **Scores vs clock_phase:** MASH 0.64, FAST 0.55. The **CD+Veh group has the lowest fasting scores**, and it is the group known to have a distinct clock state. The CD estimate therefore partly reflects an unusual reference group. The MD estimate does not depend on CD+Veh.
- **None of this separates PE from depth/yield.** Low correlation with depth within 19 mice is not evidence of no technical effect.

## Relation to Phases 4 / 4B / 4C (qualitative consistency only)
- **Gene level** (held-out-gene R²: F8-OLS 0.066 / 0.070; F8+C ≈ 0.080 / 0.096; FULL 0.104 / 0.108 but collinear):
  - the fasting axis was the most stable single predictor (and dominated under MD);
  - MASH was fragile under leave-one-out;
  - block B re-learnt PE pathway shifts, chiefly lipogenesis and cholesterol synthesis.
- **Mouse level:** the same qualitative picture. **A fasting-like shift with PE in both diets.** The MASH-like shift is not distinguishable from the shared cholesterol/fasting component.
- The two approaches **agree qualitatively**. Neither is a better "model" than the other, and they answer different questions.

## Team-owned interpretation (not concluded here)
Whether a fasting-like hepatic expression shift reflects altered feeding, gut-mediated signalling, sampling time or library artefacts is for the team. The data cannot separate these, because PE ≡ depth, CD+Veh has an unusual clock state, and sacrifice time is unknown.

## Caveats
- n = 4–5 per cell; Holm-adjusted evidence is modest (FAST PE_in_CD 0.036; PE_in_MD 0.056).
- Gene length is a median-transcript approximation, with no per-sample length offset.
- The fasting reference is ZT-averaged, male, from a 24 h fast.
- The MASH reference remains conditional (Phase 3).
