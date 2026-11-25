# Phase 3 review addendum (statistical audit)

**Date:** 2026-11-22. **Scope:** a post-hoc verification requested at review. **No original table, effect file or locked document was changed.** SHA-256 of all 60 pre-existing `reference_effects/` and `results/` files was checked before and after. New outputs are in `results/posthoc_audit/` only.

## 1. Positive-control decision: retained as FAIL on c1

The locked criterion c1 required ρ ≥ 0.30 **and** ρ above the 99th percentile of the label-permutation null. **It was not met, and it is not converted.** My earlier suggestion to swap in a "balanced 4+4 split" null is **withdrawn as a replacement**. The two procedures answer different questions:

| Procedure | Null hypothesis it tests | Status |
|---|---|---|
| **Original label permutation** (all 8 v 8 allocations equally likely) | H₀: the GAN/chow labels are exchangeable with respect to the statistic. That is, the true diet allocation of these 16 mice produces no stronger agreement with our Diet_in_Veh effect than an arbitrary allocation of the same mice. Partly correct allocations are legitimate draws under this randomisation null. | **Valid; pre-specified; primary** |
| **Restricted 4+4 splits** | Whether a contrast *orthogonal to diet* within these mice produces comparable agreement. This asks whether non-diet structure among the 16 mice (individual variation, processing) mimics our diet effect. It conditions on the labels and is not a randomisation test of the observed allocation. | **Diagnostic only; post-hoc** |

**What the results say:**
- The agreement is **supportive at the 5% level but not at the pre-specified 1% level**.
- Within these 16 mice, structure orthogonal to diet does not reproduce the agreement (balanced maximum 0.288).
- The permutation test has limited resolving power here: when the reference effect itself dominates the expression variance, allocations that keep most GAN mice together reproduce most of the agreement. That is a property of the design and the statistic. It is not an error in the null.

## 2. Exact permutation (post-hoc verification)

**Script:** `scripts/g01_exact_permutation_audit.py`. **Outputs:** `results/posthoc_audit/exact_permutation.json`, `exact_permutation_by_overlap.csv`, `exact_permutation_stats.npy`.

- All **choose(16,8) = 12,870** allocations were enumerated, using the identical estimator: Pearson correlation of average ranks of (mean log2 normalised count in 'case' − mean in 'control') and the Diet_in_Veh log2FC, on the same 11,888 genes.
- **Sample positions confirmed:** case = SD_wk23_A_R1–R8 (GAN, GSM9443427–434); control = NC_wk23_A_R1–R8 (chow, GSM9443443–450).
- The observed statistic from the enumeration equals the original (0.356339…, difference < 10⁻¹²).

| Quantity | Original Monte Carlo (1,000 draws; kept) | **Exact (12,870)** |
|---|---|---|
| Observed statistic | 0.3563 | 0.3563 |
| One-sided p | 0.033 (33/1001) | **0.036 (463/12,870)** |
| q95 | — | 0.344 |
| q99 | 0.382 | **0.390** |
| Max | 0.428 | 0.460 |
| Median | −0.011 | 0.000 |

**Null distribution by number of true GAN mice in the 'case' arm:**

| True GAN in 'case' | Allocations | Mean | Max | ≥ observed |
|---|---|---|---|---|
| 8 (observed) | 1 | 0.356 | — | 1 |
| 7 | 64 | 0.347 | 0.398 | 25 |
| 6 | 784 | 0.322 | 0.446 | 205 |
| 5 | 3,136 | 0.242 | 0.460 | 232 |
| 4 | 4,900 | 0.000 | **0.288** | 0 |
| ≤ 3 | mirror images (negative) | | | 0 |

*Correction:* the earlier report gave the balanced maximum as 0.277 (from 1,118 sampled splits); the exact value is 0.288.

**Exchangeability:**
- All 16 mice are male C57BL/6J at 21 °C (TA), wk23, with the same 5 h fast. Housing temperature is not a factor within this contrast; the 29 °C mice are a separate sensitivity contrast.
- **However, the GEO library names are consecutive by group: chow 8201–8208, GAN 8209–8216.** Library preparation order, and any processing batch that follows it, is therefore **aligned with diet**.
- Under H₀ the labels are exchangeable only if no such processing structure exists. If it does, it would inflate the agreement for allocations aligned with library order, as for the true one. This cannot be tested with the deposited data.

## 3. Positive-control caveats (genuine limitations)
- **Second reference:** MASH2 vs Diet_in_Veh ρ = 0.323 (LOO 0.24–0.36), while the two MASH references agree with each other at ρ = 0.557. Our diet effect agrees with independent MASH references clearly less well than they agree with each other.
- **Sensitivity universe (8,360 genes):** MASH1 vs Diet_in_Veh drops to **ρ = 0.252** (sign agreement 0.756), **below** the 0.30 bar. MASH2 drops to 0.211. The agreement is partly carried by the additional, mostly lower-expressed genes of the main universe, so the 11,888-gene value should not be quoted alone.
- **Lipogenesis is significantly opposite:**
  - LIPOGENESIS_CORE: fry ↑ in GAN-vs-chow (FDR 0.011) vs ↓ in our MD-vs-purified-CD (FDR 6×10⁻⁶); median log2FC +0.44 vs −1.43.
  - This is a **real limitation** of these references for our design. Lipogenesis is **not** removed from the universe.
  - The MASH signature is **not uniformly validated**: cholesterol synthesis and PPARα-target programmes agree; lipogenesis does not.

## 4. Muscle exclusion: primary unchanged; one exploratory sensitivity

- The locked 309-gene exclusion stays primary.
- A narrower **muscle-restricted** rule was written first (`posthoc_muscle_rule_prereg.md`; SHA-256 in `results/posthoc_audit/muscle_prereg_lock.json`). The rule: muscle median RPKM ≥ 10, ≥ 100× clean liver, **and** clean-liver median < 1 RPKM, plus 8 canonical markers. It uses tissue expression only.
- Run once (`g02_*`): it gives 88 genes; **122 genes are re-admitted**, including Pdk4, Plin4, Col1a1/2, Lpl, Cpt1b, Cav1 and Vldlr.

| Comparison | Locked 309 list | Exploratory 88 list |
|---|---|---|
| MASH1 vs Diet_in_Veh ρ / sign agreement | 0.355 / 0.787 | 0.358 / 0.791 |
| MASH2 vs Diet_in_Veh | 0.323 / 0.713 | 0.330 / 0.717 |
| MASH1 vs MASH2 | 0.557 / 0.899 | 0.560 / 0.898 |

- **The positive control is insensitive to the exclusion rule.**
- The re-admitted Plin4, Col1a1/2 and Lpl are up in all three diet contrasts. Diet_in_Veh Col1a1/2 can still contain muscle carry-over (GSE294180 MD_1), so their agreement is not counted as evidence.
- **Label: exploratory/post-hoc.**

## 5. PPARα and fasting: claims kept narrow
- **PPARα:**
  - The formal WT-minus-KO interaction classifies 1,893 genes as "PPARα-dependent".
  - Because WT (B6129SF2/J) and Ppara-null (B6;129S4) are different breeding stocks, the interaction estimates the **difference in agonist response between two stocks**, of which loss of PPARα is the intended but not the only difference.
  - Firm statement: GW7647 produces an **agonist-associated** response, and canonical targets (Cyp4a14, Acot1, Ehhadh, Cd36) are strongly attenuated in the knockout stock.
  - "PPARα-dependent" for the remaining genes is **provisional**.
- **Fasting:** the signed signature is the **ZT-averaged** fasted-minus-fed effect over ZT0–ZT20. 45% of fasting DEGs show ZT dependence. It is not a universal fasting state.
- **Counting:** the dose contrasts, trend, muscle-adjusted fit, clean subset and offset model are fitted on the **same mice** as their primary models. They are robustness checks, not independent validations.
  - **Independent replication** exists only for: MASH (2 studies); fasting direction (GSE225695, n = 3, weak); and time-of-day phase (female vs male, with the male data shared with the fasting study).

## 6. Verdict

**PRIMARY MASH reference (GSE316023 GAN vs chow, wk23, 21 °C): EXPLORATORY.**
- The locked positive-control decision remains **FAIL** (criterion c1: exact p = 0.036 > 0.01; ρ 0.356 < q99 0.390).
- Supportive evidence:
  - sign agreement 0.79;
  - cholesterol and PPARα-target pathway agreement;
  - LOO and offset stability;
  - a second independent MASH study agreeing at ρ 0.32;
  - insensitivity to the muscle rule.
- Against stronger labels:
  - the 1% criterion failed;
  - ρ is 0.25 in the stricter universe;
  - lipogenesis is opposite;
  - library order is aligned with diet.

The axis may enter Phase 4 **only as an exploratory axis**, with these caveats attached to every result that uses it. "Validated" is not supported. "Conditional", as defined in the locked plan, also required c1 to hold.
