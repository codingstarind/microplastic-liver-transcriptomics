# Phase 3 decision report

**Date:** 2026-11-22 · **Branch:** `research/phase2-datasets` (uncommitted Phase 3 files) · **Status:** for team review. Phase 4 has **not** been started.

Labels used throughout:
- **[obs]** directly observed number;
- **[ev]** statistical evidence;
- **[int]** interpretation needing team review;
- **[spec]** speculation.

## 1. Which reference signatures are reliable?

| Axis | Signature | Verdict |
|---|---|---|
| MASH | GSE316023 GAN vs chow (8 v 8), DESeq2, no offset | **Technically sound; candidate axis.** The positive control formally failed one criterion (§2). |
| MASH | GSE246328 (7 v 8, ≈ 58 wk GAN) | Technically sound; agrees with GSE316023 (ρ 0.56, 24/24 pathways). One influential chow mouse (16234). |
| Fasting | GSE107787 FAST vs FED, ZT-adjusted (18 v 18), limma-trend on log2 RPKM | **Stable** to the muscle-clean subset (ρ 0.92). **Strongly time-dependent** (45% of DEGs have a ZT interaction), so it is an *average over ZT*. RPKM only. |
| Fasting | GSE225695 (3 v 3) | Direction check only. Agrees with the primary signature (ρ 0.49, 80% sign agreement on primary DEGs). |
| Time of day | GSE117134 WT female, harmonic | 2,914 rhythmic genes; phases agree with male data (median 1.4 h). Usable as **phase** information and as the T1/T2 contrasts. Female, RPKM. |
| PPARα | GSE282998 WT GW7647 20v0 (4 v 4 within a 16-mouse fit) | **Robust** to the depth-matched, muscle-matched and muscle-adjusted analyses. 1,893 genes formally PPARα-dependent vs the knockout (different breeding stocks: provisional except for canonical targets). Targeted FFPE assay. |

## 2. Did the independent MASH positive control pass?

**Formally: FAIL** under the locked criteria.
- [obs] ρ = 0.355 meets the ρ ≥ 0.30 bar, but does not exceed the relabelling-null 99th percentile (0.382).
- All other criteria passed:
  - sign agreement 0.79;
  - cholesterol and PPARα-target pathways agree;
  - LOO 0.30–0.37;
  - offset 0.356.
- [ev] **Diagnosis:** the relabelling null contained partially correct labellings, which carry the very large diet effect. Among exactly balanced (4+4) relabellings the maximum was 0.277, below the observed value. The criterion I proposed was mis-specified for a dominant effect; the agreement itself is not in doubt at that level.
- [obs] **A genuine disagreement:** lipogenesis goes ↑ in GAN-vs-chow but ↓ in our MD-vs-purified-CD.
  - [int] This is plausibly a control-diet effect.
  - It means lipogenesis genes cannot be used to judge "MASH-likeness" against chow-controlled references.

**Team decision needed:** accept a **proposed deviation** (replace c1's null with balanced 4+4 relabellings, which would give PASS on all four criteria), or keep FAIL and treat the MASH axis as exploratory. I have **not** applied the deviation.

## 3. How well did the two MASH datasets agree?
- [obs] ρ = 0.557.
- 90% sign agreement among 3,306 genes significant in both.
- 24/24 Hallmark/custom sets significant in both agree in direction.
- Each agrees with our Diet_in_Veh at ρ ≈ 0.32–0.36.

That gap (0.56 between the references vs about 0.34 with our contrast) is the realistic ceiling for "independent MASH similarity" with these data.

## 4. Are fasting and PPARα distinguishable?
**Yes.**
- [obs] Genome-wide ρ = 0.04; VIF ≈ 1.0 for MASH, fasting and PPARα.
- Both raise PPARα-target sets.
- They move in **opposite directions** for lipogenesis (fasting −2.8 vs GW7647 +1.2 median log2) and cholesterol synthesis (−2.8 vs +0.7).
- The two circadian contrasts are collinear with each other (ρ 0.86) and must not both enter one model.

## 5. Which signatures are technically fragile?
- **Fasting (GSE107787):** RPKM-only; time-dependent; muscle-contaminated FED livers (handled; the clean subset agrees).
- **Circadian:** female, RPKM, single study.
- **PPARα:** targeted FFPE assay; dose-linked depth and muscle batch (shown not to drive the signature); genotype stocks differ.
- **MASH2:** one influential control mouse.
- **MASH1:** an undocumented length definition. 883 genes are offset-sensitive, with no sign reversals.

## 6. What do the conventional baselines show?
- [obs] **DEG overlap** is "significant" for almost every pair, **including pairs where most shared genes go in opposite directions** (PE_in_CD vs MASH2: 84 concordant / 125 discordant). It is not informative about state similarity.
- [obs] **Effect correlations with PE:**

  | Reference | PE_in_CD | PE_in_MD |
  |---|---|---|
  | MASH1 | 0.21 (LOO down to **0.03**) | 0.16 (LOO down to **0.02**) |
  | MASH2 | ≈ 0 | ≈ 0 |
  | Fasting | 0.14 (stable, 0.10–0.20) | 0.25 (stable, 0.22–0.27) |
  | PPARα | −0.05 | −0.09 |

- [int] With independent references, the PE effect's resemblance to MASH is weak and driven by single mice. Its resemblance to fasting is small but stable.
- [spec] This pattern is compatible with the fasting/feeding-state explanation raised in Phase 1. It is **equally compatible with technical depth effects**, which no reference can separate. No mechanism is implied.

## 7. What modelling is feasible for Phase 4?
- **Feasible:**
  - SE-weighted regression of PE gene effects on MASH1, FAST and PPARα 20v0 (VIF ≈ 1), plus **one** circadian term (T1, or harmonic phase-based features);
  - gene-set-blocked cross-validation;
  - comparison against the Baseline-2 numbers above, on the same universe and with the same LOO treatment.
- **Not feasible:**
  - separating PE from library depth;
  - magnitude-based comparisons across technologies;
  - using lipogenesis genes as MASH evidence with chow-controlled references.
- Expect weak absolute fits: every single-axis ρ is ≤ 0.25.

## 8. Unresolved assumptions and limitations
1. **The muscle exclusion list (fixed, Phase 2) removes liver-relevant genes:** Pdk4, Plin4, Col1a1/Col1a2, Lpl, Cpt1b, Cav1, Vldlr.

   This is a real technical problem; it was not driven by results. **Proposed alternative** (not applied): restrict exclusion to genes with ≥ 100× muscle/liver enrichment, or to sarcomere/contractile genes. Rerun as a sensitivity analysis and keep both versions. **Team decision.**
2. **Ppara discrepancy** (`results/ppara_discrepancy/ppara_estimates_by_approach.csv`):
   - [obs] Every DESeq2 configuration on the deposited counts gives PE_in_CD Ppara between **+0.16 and +0.17** (pairwise or joint, with or without prefilter, either set of size factors, no outlier replacement), and +0.05 under MD.
   - Jung reports +0.53 / +0.34.
   - **The difference is not explained by DESeq2 settings.** The similar gap in both diets (≈ 0.3–0.4 log2) suggests a different count input or normalisation in their analysis [spec].

   **No PPARα claim about PE should rely on Ppara itself.**
3. The c1 null flaw (§2).
4. **The GSE294180 limitations are unchanged:**
   - PE ≡ lower depth and yield;
   - CD+Veh clock state;
   - muscle carry-over;
   - unknown sacrifice time;
   - no length offset;
   - our target effects have no length offset, while the reference choice was also made no-offset.
5. The deadline (2026-11-25 per the draft) is unconfirmed. The AIBD GSE245069 files are still not found.

## Deviations from the plan
See `README.md` § Deviations.
