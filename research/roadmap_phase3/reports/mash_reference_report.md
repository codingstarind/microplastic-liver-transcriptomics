# MASH reference signatures and positive control

**Scripts:** `scripts/b01_mash_references.R` (DESeq2), `scripts/c01_positive_control.py`.
**Effects:** `reference_effects/MASH1_*.csv.gz` and `MASH2_*.csv.gz` (every universe gene: log2FC, lfcSE, 95% CI, Wald stat, p, BH padj, baseMean, ashr log2FC, detected flag, max Cook's).
**Figures:** `figures/mash_positive_control.png`, `figures/diagnostics/MASH*_diagnostics.png` (dispersion, MA, Cook's).

Labels used below: **observed** = a computed number; **evidence** = a statistical test result; **interpretation** = needs team review.

## 1. Models (observed)

| Signature | Samples (fixed Phase 2 list) | Model | Genes | padj < 0.05 | DEG (padj < 0.05, \|log2FC\| ≥ 0.38) | Outlier-replaced genes |
|---|---|---|---|---|---|---|
| **MASH1** GSE316023 GAN vs chow, 21 °C, wk23 (primary, no offset) | 8 v 8 | DESeq2 `~ diet`, ref = chow, coef `diet_GAN_vs_chow` | 11,888 | 6,260 | 4,777 | 1 |
| MASH1 + length offset (sensitivity) | 8 v 8 | as above + `avgTxLength` | 11,888 | 5,763 | 4,486 | 2 |
| MASH1 TN (29 °C) sensitivity | 8 v 8 | as above | 11,888 | 7,155 | 6,209 | 63 |
| **MASH2** GSE246328 GAN vehicle 24 w (≈ 58 wk GAN) vs chow 24 w | 7 v 8 | DESeq2 `~ diet` | 11,883 | 5,252 | 4,297 | 139 |

- **Direction check** (tested): disease markers are up in GAN in both studies (Gpnmb +6.8/+7.7, Trem2 +3.6/+4.1, Lgals3, Cd68 and Col3a1 all up). The coefficient is GAN minus chow, not the reverse.
- **Diagnostics:**
  - Dispersion trends are well behaved.
  - GSE316023 has no influential sample (≤ 1 gene with Cook's > F₀.₉₉ per sample).
  - In GSE246328, **chow mouse 16234 is influential** (126 genes above the Cook's cutoff; it was also flagged by the Phase 2 correlation screen). It was kept, per the fixed sample list, and checked by leave-one-out.
- **Leave-one-out** (16 and 15 refits): the correlation of each signature with Diet_in_Veh stays in 0.336–0.369 (MASH1) and 0.242–0.359 (MASH2).

## 2. MASH positive control (GSE294180 Diet_in_Veh vs MASH1)

Criteria pre-specified in `analysis_plan.md` §4.

| Metric | Observed |
|---|---|
| M1 Spearman ρ, main universe (11,886 genes) | **0.355** |
| M2 sign agreement, genes padj < 0.05 in both (n = 1,247) | **0.787** |
| Sign agreement on all MASH1-significant genes (n = 6,260) | 0.666 |
| M3 DEG overlap | 747 of 1,363 target DEGs (expected 548); hypergeometric p = 4×10⁻³¹; Jaccard 0.139; 631 concordant / 116 discordant |
| M4 core pathways, fry FDR < 0.05 in both | 3 sets: CHOL_SYNTHESIS_CORE ↓↓ agree, SANDERSON_PPARA_TARGETS ↑↑ agree, **LIPOGENESIS_CORE ↑ (MASH1) vs ↓ (Diet_in_Veh) disagree** → 2/3 |
| M5 relabelling null (1,000 random 8 v 8 relabellings of the GSE316023 mice) | null median −0.01, **99th percentile 0.382**, max 0.428; observed 0.356; permutation p = 0.033 |
| M5 LOO | reference 0.336–0.369; target (19 existing refits) 0.304–0.374 |
| M5 offset model | ρ = 0.356 |
| M5 sensitivity universe (8,360) | ρ = 0.252 (sign agreement 0.756) |

**Pre-specified decision:**

| Criterion | Result |
|---|---|
| c1: ρ ≥ 0.30 **and** above the relabelling null's 99th percentile | **Not met** (0.356 < 0.382) |
| c2: sign agreement ≥ 0.70 | Met |
| c3: pathways ≥ 2/3 | Met |
| c4: LOO and offset ρ ≥ 0.20 | Met |

**Formal outcome: FAIL**, because c1 is not met.

### Diagnosis of the failure (evidence, not a redefinition)
- **The relabelling null includes partially correct labellings.** With 8 v 8 mice, a random split often puts 5–7 of the 8 true GAN mice together, and those splits still carry most of the very large diet effect (PC1 = 50% of variance, r = 0.97 with diet).
  - The null draws that exceeded the observed ρ came only from splits with ≥ 5 true GAN mice in one arm.
  - Among the **1,118 exactly balanced splits (4 GAN + 4 chow per arm)**, the maximum ρ was **0.277**, below the observed 0.356.
- **So the c1 test, as written, cannot reach p < 0.01 for any strong diet effect.** It is a flaw in the criterion I proposed, not evidence of disagreement.
- I have **not** replaced the criterion. A corrected null (balanced/orthogonal splits) is offered to the team as a **proposed deviation** (see `phase3_decision_report.md`).
- **Real disagreement exists in lipogenesis.**
  - GAN vs chow raises the lipogenesis core (fry ↑, median log2FC +0.44); our MD vs purified-CD contrast lowers it (−1.43).
  - Interpretation (requires team review): this is the expected direction if the purified control diet (sucrose-rich) itself induces lipogenesis relative to grain chow. A chow-controlled reference would then show "MASH" as lipogenesis-up while our purified-controlled contrast shows it as down.
  - This is a **control-diet** difference. It makes lipogenesis genes **unsuitable** for judging MASH similarity with these references.

## 3. Reproducibility of the MASH response across studies (observed)

| Pair | Spearman ρ | Sign agreement (both padj < 0.05) | Notes |
|---|---|---|---|
| MASH1 vs MASH2 (two independent GAN studies) | **0.557** | 0.899 (n = 3,306) | Relabelling-null q99 = 0.594, so it exceeds no balanced-null bound here either (same null flaw). Pathways: **24/24** sets significant in both agree in direction. |
| MASH1 vs Diet_in_Veh | 0.355 | 0.787 | see above |
| MASH2 vs Diet_in_Veh | 0.323 | 0.713 (n = 1,034) | LOO 0.242–0.359; null q99 0.349 |
| MASH1 TA vs MASH1 TN | 0.842 | 0.997 | housing temperature matters little for the axis |
| MASH1 no-offset vs offset | 0.981 | 1.000 | |

- **Benchmark:** two independent GAN-vs-chow studies agree at ρ ≈ 0.56. Our 8-week MD-vs-purified-CD contrast agrees with each at ρ ≈ 0.32–0.36, and the inflammatory, fibrosis, cholesterol and PPARα-target programmes all agree in direction.
- The weaker agreement is consistent with known differences: control diet, duration (8 vs 23 vs 58 wk), and our n = 4/5. These differences are **not** separately estimable.

## 4. What this means for the MASH axis
- MASH1 and MASH2 are **technically sound and mutually consistent** signatures of "GAN diet vs chow".
- As references for "MASH-like" in our purified-diet experiment, they are usable **with the lipogenesis caveat**, and only if the team accepts the diagnosis of the c1 failure.
- Until then the axis is labelled **candidate, positive control formally failed (c1 null flaw)**.
