# Fasting reference

**Script:** `scripts/d01_other_references.R` (fasting section), `scripts/f01_compare_and_baselines.py`.
**Effects:**
- `reference_effects/FAST_GSE107787_FASTvsFED_ZTadj.csv.gz` (primary);
- `FASTclean_GSE107787_ZT0ZT4_muscleLow.csv.gz` (fixed clean subset);
- `FAST2_GSE225695_FastvsFed.csv.gz` (direction check);
- per-ZT effects in `results/fasting/GSE107787_fasting_effect_by_ZT.csv.gz`.

**Figure:** `figures/fasting_reference.png`. Mean–variance diagnostic: `figures/diagnostics/FAST_GSE107787_meanvar.png`.

## Model (as pre-specified)
- GSE107787 liver, 18 fasted (24 h) vs 18 fed, 3 per state at each of ZT0, 4, 8, 12, 16 and 20 (male C57BL/6J).
- Values are deposited RPKM, so **no count model is used.** limma-trend on log2(RPKM + 1), `~ ZT (factor) + state`, reference FED, coefficient `stFAST` (fasted minus fed, log2 RPKM scale).
- Universe: 11,888 genes. The 309 muscle genes are already excluded.
- **Mean–variance:** the limma-trend fit is smooth. Low-RPKM genes have higher variance, which the trend absorbs. Genes with all-zero RPKM got limma's "zero variance offset" warning and are effectively uninformative (not removed, to keep the universe fixed).

## Results (observed)

| Quantity | Value |
|---|---|
| Genes padj < 0.05 / DEGs (\|log2FC\| ≥ 0.38) | 8,552 / 4,009 |
| `~ ZT * state` interaction: genes FDR < 0.05 | **4,296** |
| Fasting DEGs with interaction FDR < 0.05 | **45%** (pre-specified threshold 10%) |
| Clean subset (ZT0+ZT4, muscle-low, 4 v 6) vs main: Spearman ρ | **0.919** |
| Sign agreement of the clean subset on main DEGs | **0.995** |
| GSE225695 (3 v 3, sacrifice time unknown) vs main: ρ | 0.489; sign agreement on main DEGs 0.801 |

**Programmes** (fry, self-contained; camera competitive in brackets):
- lipogenesis core ↓ (FDR < 10⁻⁴; camera ↓);
- cholesterol synthesis ↓ (camera ↓);
- SANDERSON_PPARA_TARGETS ↑ (camera ↑);
- WP_PPAR_SIGNALING ↑;
- inflammatory response ↓ by fry, but **not** competitively (camera FDR 0.99).

Medians (log2): lipogenesis −2.79, cholesterol synthesis −2.76, PPARα targets +0.82.

## Pre-specified labels
- **Contamination stability: STABLE.** ρ (main vs clean) = 0.92 ≥ 0.60 and sign agreement 0.995 ≥ 0.80. The signature is not driven by the muscle carry-over in FED ZT8–20 livers once muscle genes are excluded.
- **Time dependence: YES.** 45% of fasting DEGs have a significant ZT × state interaction, far above the 10% threshold. **The single signed fasting effect is an average over ZT0–ZT20 and must be labelled that way.** For many genes the size, and sometimes the sign, of the fasting effect depends on the time of day (see the per-ZT panel in the figure).

## Limitations
- RPKM only: no sequencing-depth weighting, no NB model, and limma-trend inference is approximate for low-RPKM genes.
- A 24 h fast is extreme, so magnitudes are not comparable with PE effects. Use rank and sign only.
- Male; chow; single-housed.
- The GSE225695 check is weak: n = 3, unknown sacrifice time, and depth separates the groups.
- Pdk4, a classic fasting gene, is in the fixed muscle exclusion list and so is absent from this signature (see the decision report).
