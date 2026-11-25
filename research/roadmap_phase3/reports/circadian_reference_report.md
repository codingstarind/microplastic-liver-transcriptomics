# Circadian (time-of-day) reference

**Script:** `scripts/d01_other_references.R` (`harmonic()`).
**Outputs:**
- per-gene rhythm tables in `results/circadian/CIRC_GSE117134_WTfemale_rhythm.csv.gz` and `CIRCsens_GSE107787_FEDmale_rhythm.csv.gz`;
- signed time contrasts in `reference_effects/CIRC_*_T1_ZT6_minus_ZT2.csv.gz` and `*_T2_ZT12_minus_ZT0.csv.gz`;
- clock-gene table `results/comparisons/circadian_clock_genes.csv`.

**Figure:** `figures/circadian_reference.png`.

## Model
- GSE117134 WT **female** liver, ad libitum chow, ZT0–ZT20 every 4 h, n = 3 per ZT (18 samples). RPKM.
- limma-trend on log2(RPKM + 1) with `~ cos(2πZT/24) + sin(2πZT/24)`:
  - rhythmicity = F-test on both terms (BH);
  - amplitude = √(a² + b²) on the log2 scale (half peak-to-trough);
  - peak phase = atan2(b, a)·24/2π;
  - phase SE by the delta method.
- **Sensitivity:** a categorical `~ factor(ZT)` F-test.

## Results (observed)

| Quantity | GSE117134 female (primary) | GSE107787 FED male (descriptive) |
|---|---|---|
| Rhythmic genes, harmonic FDR < 0.05 | 2,914 | 4,389 |
| Rhythmic genes, categorical FDR < 0.05 | 2,379 | 5,518 |
| Rhythmic by both | 1,771 | 3,752 |

- Genes rhythmic in both datasets: 1,637.
- Median absolute peak-phase difference: **1.4 h**; 80% within ± 3 h. The core clock (Bmal1, Per1–3, Cry1, Nr1d1/2, Dbp) has the expected phases in both (table in `circadian_clock_genes.csv`).

## Signed time contrasts (for comparisons only)
Rhythmicity amplitude is **not** a signed effect. Two pre-specified predicted differences from the harmonic fit, with SE (limma contrasts), are the only circadian quantities used in signed comparisons:

| Contrast | Meaning | DEGs (padj < 0.05, \|log2FC\| ≥ 0.38) |
|---|---|---|
| **T1 = ZT6 − ZT2** | a 4-h shift inside a typical morning necropsy window | 205 |
| **T2 = ZT12 − ZT0** | light → dark transition | 1,369 |

- T1 and T2 are correlated (ρ = 0.86), because both are linear combinations of the same two harmonic coefficients for each gene. They are **two views of one rhythm**, not independent axes.
- **Phase is the scientifically primary output.** A PE effect that resembles T1 means "shifted as if sampled ~4 h later in the morning". It does not mean "circadian disruption".
- Because GSE294180's sacrifice time is unknown, any time contrast is a hypothesis about sampling time, not a measurement of it.

## Limitations
- Female (ours male). The male check uses FED mice from the fasting study, which are **not independent** of the fasting reference and carry muscle contamination at ZT8–20.
- RPKM only.
- Light–dark sampling under ad libitum feeding mixes clock-driven and feeding-driven rhythms.
- Six time points at 4 h resolution.
