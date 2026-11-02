# Direction A: does the PE response resemble the MASH-diet response?

> **Status: PROVISIONAL and EXPLORATORY.** Deposited counts with unverified quantification.
> The diet signature comes from the same experiment, so nothing here is independent
> confirmation.

**Code:** `analysis/concordance/pe_vs_diet_concordance.py`
**Outputs:** `results/concordance/<label>/`
**Figure:** `figures/phase2/provisional_deposited/fig17_pe_vs_diet_concordance.png`

## Question

Are PE-associated expression changes similar to, opposite to, or unrelated to the changes produced by the MASH-inducing diet? Is any similarity transcriptome-wide or confined to particular programs?

## The shared-reference problem and how we handled it

The PE effect under the control diet is (CD+PE − CD+Veh); the diet effect is (MD+Veh − CD+Veh). Both subtract the **same four CD+Veh mice**, so random noise in those mice moves both estimates in the same direction. Even with no real similarity this produces a positive correlation; for equal variances its expected size is about 0.56 of the noise component.

**What we did:**

1. **Split reference.**
   - The 4 CD+Veh mice are split into two pairs. One pair is the PE reference, the other the diet reference, so the two contrasts share no mice.
   - All 6 assignments are evaluated (3 partitions × 2 orders).
   - For PE under MD, the shared group is MD+Veh (n = 5), split 2/3 in all 20 ways.
2. **Reference distribution.**
   - Within each split, the 7 mice of the PE contrast (5 PE + 2 reference) are relabelled in all 21 possible ways.
   - This keeps gene–gene correlation and the independence from the diet contrast, and removes any real PE effect.
   - With 21 relabellings, the smallest attainable one-sided p is 1/21 = 0.048.
3. **Measures:**
   - Spearman correlation of gene log2 fold changes (9,941 genes, mean normalised count > 50);
   - sign agreement among the 500 genes most changed by diet, chosen using the diet half only;
   - Spearman correlation of Hallmark-set mean fold changes.

**What this does not fix:** a systematic difference shared by *all* CD+Veh mice, such as the time of sacrifice. Splitting removes shared *noise*, not a shared *bias*. See "Alternative explanations".

## Results (provisional)

| Comparison | Naive ρ (shared reference) | Split-reference ρ, median [range] | Null median | One-sided p (range over splits) | Top-500 sign agreement (null ≈ 0.50–0.54) |
|---|---|---|---|---|---|
| PE under CD vs diet | 0.54 | **0.33 [0.29–0.40]** | −0.01 | 0.048–0.095 | 0.76 [0.72–0.78] |
| PE under MD vs diet | −0.08 | 0.05 [−0.03–0.17] | −0.01 | 0.10–0.57 | 0.55 [0.44–0.69] |
| PE under CD vs diet, without MD_1 | – | 0.35 | – | 0.048–0.095 | 0.77 |

At the pathway level, the Hallmark mean-LFC correlation for PE under CD vs diet is 0.24 (median over splits; one-sided p 0.19–0.62), which is not distinguishable from the null.

**Reading:**

- About 40% of the naive PE(CD)–diet correlation was shared-reference artefact (0.54 → 0.33).
- A moderate positive gene-level resemblance remains under the control diet. It sits at or near the top of every relabelling distribution, but 21 relabellings allow no p below 0.048, so the evidence is moderate.
- Under the MASH diet, PE changes are **not** aligned with the diet axis.
- At the pathway level, the specific programs that move in the diet's direction with PE are **cholesterol homeostasis** (down) and **PPARα targets** (up); see the pathway report. They are not broad transcriptome-wide shifts.

## Alternative explanations (important)

1. **Circadian/feeding state of the CD+Veh group.**
   - A per-mouse clock-phase score was computed: the mean z of Bmal1/Npas2/Clock/Nfil3/Cry1 minus the mean z of Per1–3/Dbp/Tef/Nr1d1–2/Bhlhe40–41, gene list fixed before looking.
   - All four CD+Veh mice (−3.15 to −1.49) are **lower than every other mouse in the study**. The closest are CD_MP_5 (−1.39) and CD_MP_4 (−1.33). Group means: CD+Veh −2.38; CD+PE +0.08; MD+Veh +1.01; MD+PE +0.82.
   - The CD+Veh group therefore differs systematically in clock phase, i.e. probably in time of sacrifice or feeding state. That alone would push both PE_in_CD and Diet_in_Veh the same way for the many clock-controlled liver genes, including SREBP and PPARα targets.
   - It would also explain why the resemblance is absent under MD, where the two groups' clock phase is similar.
   - Adding the clock score as a covariate does not resolve this: it is nearly collinear with the CD+Veh indicator and removes much of the diet effect too (2,262 → 747 genes).
2. **Library order.** CD+Veh are the earliest libraries (S39–S42). A monotonic processing-order drift would also predict concordance for PE under MD (later MD+PE vs earlier MD+Veh), which is not seen. This argues against simple drift.
3. **Feeding behaviour under PE.** Reduced intake would produce a fasting-like profile (SREBP down, PPARα up) resembling some diet effects for a different reason.

## What would support a stronger inference

- An **independent MASH signature** from another study, so the disease axis is not built from these mice. MSigDB 2026.1 has no directional MASH signature, only a curated KEGG/WikiPathways NAFLD pathway.
- GEO lists recent mouse MASH-diet bulk RNA-seq series. Each candidate needs checking for:
  - a wild-type diet-vs-control arm;
  - male C57BL/6 mice;
  - an HFD/fructose/cholesterol diet;
  - deposited counts.

  Candidates found on 2026-11-02, not yet vetted: GSE329899, GSE342274, GSE328603, GSE333221, GSE242672, GSE272293, GSE317301/2.
- A human reference, e.g. bulk biopsies in GSE267031, would add cross-species questions.
- **This is a decision for the team** (see interpretation).
- Sacrifice time and food-intake records from the authors would address the circadian confound directly.

## Status of the claim

**Exploratory, not supported as a stronger inference.** A modest PE–diet resemblance under the control diet survives shared-reference correction, but it is confounded with a systematic circadian/feeding difference in the reference group. It is absent under the MASH diet. We do not claim that PE "mimics" metabolic dysfunction.

## Course connections

- **Correlated errors** violate the independence assumption behind a correlation test. Making the samples disjoint is a design-based remedy.
- A **permutation null must match the hypothesis.** Arbitrary label shuffling would have tested "no PE effect at all", not "PE and diet effects are unrelated".
- **Confounding vs noise:** resampling and splitting address noise; only design or external information can address confounding.
