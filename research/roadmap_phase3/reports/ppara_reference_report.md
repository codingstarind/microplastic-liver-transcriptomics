# PPARα activation reference (GW7647, GSE282998)

**Script:** `scripts/d01_other_references.R` (PPARα section).
**Effects:** `reference_effects/PPARA_GSE282998_WT_{5v0,10v0,20v0,trend,20v0_muscleAdj}.csv.gz`, `PPARA_GSE282998_KO_20v0.csv.gz`, `PPARA_GSE282998_WTminusKO_20v0_interaction.csv.gz`.
**Classification:** `results/ppara/genotype_classification.csv.gz`. **Key genes:** `results/comparisons/ppara_key_genes.csv`.
**Figure:** `figures/ppara_reference.png`.

## Models (as pre-specified)
- **WT B6129SF2/J, 4 per dose:** DESeq2 `~ dose` (factor, reference 0) on TempO-Seq probe counts (one probe per gene; Phase 2 mapping).
  - Contrasts 5v0, 10v0, 20v0 (coefficients `dose_X_vs_0`, checked).
  - Trend model `~ dose/10` (log2FC per 10 mg/kg).
  - Muscle sensitivity: `~ muscle_high + dose`.
- **Genotype:** WT + Ppara-null (B6;129S4-Ppara^tm1Gonz/J), `~ genotype + genotype:dose`. The interaction contrast is `genotypeWT.dose20 − genotypeKO.dose20`, with Wald SE.
- **Technical:** size factors follow dose (0 and 5 mg/kg ≈ 0.6–0.85; 10 and 20 mg/kg ≈ 1.2–1.8), mirroring the Phase 2 depth finding.

## Results (observed)

| Quantity | Value |
|---|---|
| WT 20v0 genes padj < 0.05 | 3,604 |
| ρ (20v0 vs 10v0) / (20v0 vs 5v0) | 0.689 / 0.480 |
| Sign agreement on 20v0-significant genes: 10v0 / 5v0 | 0.964 / 0.830 |
| ρ (20v0 vs muscle-adjusted 20v0) | 0.946; sign agreement 1.00 |
| ρ (20v0 vs dose trend) | 0.978 |
| 20v0-significant genes with the same sign at 5 and 10 mg/kg | 81% |
| ρ WT 20v0 vs KO 20v0 | 0.25 |
| Classes | **PPARa_dependent 1,893** · agonist-associated, dependence unresolved 2,613 · not responsive 7,335 |

- **Canonical targets** (WT 20v0 vs KO 20v0, log2FC): Cyp4a14 +5.7 vs +1.7; Acot1 +5.8 vs −0.3; Ehhadh +4.8 vs +0.5; Cd36 +4.3 vs +0.2. All classed PPARα-dependent by the formal interaction test, with clear dose responses.
- **Programmes:** fatty-acid metabolism ↑, PPAR target sets ↑, and also **lipogenesis core ↑ and cholesterol synthesis ↑** (fry FDR < 0.01; Fasn +1.95, classed PPARα-dependent).

## Technical sensitivity
- **5v0 (depth-matched) vs 20v0:** directions mostly agree (83%), with weaker correlation. The 5 mg/kg response is smaller.
- **10v0 (muscle-matched):** high agreement (96%).
- **Muscle-adjusted fit:** almost identical (ρ 0.95).

The dose-linked depth and muscle batch does **not** drive the main agonist signature. The strong dose response (PC1 r = 0.96 with dose in Phase 2) dominates.

## Specificity
- "PPARα-dependent" requires a significant WT-minus-KO interaction in the WT direction. It is not "significant in WT, not in KO".
- **Caveat:** WT (B6129SF2/J, F2 hybrid) and KO (B6;129S4 congenic) are different breeding stocks. An interaction therefore also absorbs background differences in drug response. Treat "dependent" as **strong evidence of PPARα involvement for the canonical targets**, and as provisional for the rest.
- Genes in the "unresolved" class should be described as **agonist-associated** only.

## Limitations
- Targeted FFPE TempO-Seq probes: magnitudes are **not** on the RNA-seq scale, so compare by rank and sign only.
- High-dose pharmacological agonist for 5 days. This is not physiological PPARα activity and not a measurement of PPARα protein activity in our mice.
- Pdk4 is absent from the probe panel (and from the muscle-filtered universe); Per1/Per2 are poorly detected.
