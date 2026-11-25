# Factorial DESeq2 model and diagnostics

> **Data status:** deposited GEO matrix used without a tximport length offset, **pilot-verified under assumptions A1–A3** (`reports/quantification_validation/quantification_validation_report.md` §5.1; team decision 2026-11-02). Mitochondrially encoded genes and processed pseudogenes are less reliable. All analyses remain **exploratory** until the team adopts hypotheses.

**Code:** `analysis/deseq2/run_deseq2_factorial.R`
**Outputs:** `results/deseq2/provisional_deposited/` (tables), `figures/deseq2/provisional_deposited/` (diagnostics)
**Software:** R 4.4.3, DESeq2 1.46.0, ashr 2.2.63 (environment: `envs/r-analysis.yml`)

## Biological question

How large is the PE-associated change in hepatic expression under each diet? How large is the diet effect? Does the PE effect differ between diets?

## Model (BT3041: regression with a 2 × 2 factorial design)

For gene g and mouse s, counts follow a negative binomial distribution with gene-specific dispersion:

```
log μ_gs = log s_s + β0 + β_diet·[MD] + β_PE·[PE] + β_int·[MD]·[PE]
```

- **Reference levels:** control diet (CD) and vehicle (Veh), set explicitly. The code stops if the coefficient names differ from `Intercept, diet_MD_vs_CD, microplastic_PE_vs_Veh, dietMD.microplasticPE`.
- **Contrast signs** were checked against raw group means. For genes with baseMean > 100, the sign agrees in 99.7–99.9% of genes for every contrast; the remainder are near-zero effects.

| Contrast | Coefficients | Meaning |
|---|---|---|
| PE_in_CD | β_PE | CD+PE vs CD+Veh |
| PE_in_MD | β_PE + β_int | MD+PE vs MD+Veh |
| Diet_in_Veh | β_diet | MD+Veh vs CD+Veh |
| Interaction | β_int | (MD+PE − MD+Veh) − (CD+PE − CD+Veh), on the log2 scale |

**Why a negative binomial:** Phase 1 (fig05) showed that count variance grows faster than the mean (overdispersion), which rules out Poisson or ordinary *t*-tests on counts. DESeq2 shrinks each gene's dispersion toward a fitted mean–dispersion trend. That stabilises estimates when only 4–5 mice per group are available.

**Pre-filter:** ≥ 10 counts in ≥ 4 samples (as in Phase 1), giving 13,777 genes.

**Testing:** each contrast is its own test family. The test is a Wald test, followed by Benjamini–Hochberg FDR after DESeq2's independent filtering. CIs are MLE ± 1.96 × SE. ashr-shrunken LFCs are reported for ranking and plotting.

## Results (provisional)

| Contrast | Genes tested | FDR < 0.05 | ... and \|log2FC\| ≥ 0.38 | π₀ estimate* |
|---|---|---|---|---|
| PE_in_CD | 13,722 | 858 | 513 | 0.78 |
| PE_in_MD | 13,722 | 472 | 295 | 0.80 |
| Diet_in_Veh | 13,722 | 2,262 | 1,582 | 0.59 |
| Interaction | 13,722 | 87 | 84 | 0.91 |

\*Storey's π₀ (λ = 0.5): a rough estimate of the share of genes with no effect.

Full per-gene tables are in `results_<contrast>.csv.gz`. Columns: baseMean, log2FC, lfcSE, ci95_low, ci95_high, Wald stat, p-value, padj, ashr log2FC and SE, and the naive group-mean difference.

**Significance vs magnitude vs relevance:**

- Genes with FDR < 0.05 have a detectable change *in these 19 mice*. They need not be large or biologically important.
- The median |log2FC| of FDR genes is about 0.5–1 (1.4–2-fold).
- A non-significant gene is **not** evidence of no effect. The median interaction CI width is 0.80 log2, so interactions smaller than about ±0.4 log2 are invisible.

## Diagnostics

| Diagnostic | Figure | Observation |
|---|---|---|
| Mean–dispersion trend | `dispersion_trend.png` | Usual decreasing trend, levelling at ~0.012 (biological CV ≈ 11%) for well-expressed genes. A tail of high-dispersion outlier genes is not shrunk; those are mostly low-count genes and contamination-affected genes. |
| P-value histograms | `pvalue_histograms.png` | All four contrasts show a flat background with a peak near 0. There is no U-shape or hump at intermediate p, so no obvious misspecification. The interaction peak is the weakest. |
| MA plots | `ma_<contrast>.png` | Effects span the expression range. PE_in_MD has a cluster of large negative effects at moderate mean; these are the muscle genes (see the interaction report). |
| Cook's distance | `cooks_distance_by_sample.png`, `sample_influence.csv` | Most influential: **MD_1 (19 genes), MD_MP_2 (18), CD_MP_2 (15)**. These are the samples Phase 1 flagged for muscle, cholangiocyte and Kupffer/blood marker excess. With 4–5 replicates per group, DESeq2 **cannot replace** outliers (that needs ≥ 7); 55 genes had their p-values set to NA instead. |
| Size factors | `sample_influence.csv` | 0.78–1.26, identical to the Phase 1 median-of-ratios values |

## Checks before interpreting large numbers of genes

These are summarised in the robustness report:

- **PE_in_MD and the interaction are partly contamination artefacts.** With a muscle-marker covariate, PE_in_MD drops from 472 to 135 genes and the interaction from 87 to 16. 23 of the 87 interaction genes are muscle-tracking, against 1% of all tested genes.
- **PE_in_CD is robust to contamination** (928 genes with the covariate) but **not separable from library depth/quality**. With a depth or quality-proxy covariate, which are collinear with PE, it falls to 116 or 66 genes.
- **Mitochondrial and top-abundant genes** are over-represented among PE_in_CD FDR genes: 11 genes, 1.3% of the FDR list, against 0.11% of tested genes for mitochondrial genes. A quality-linked artefact is possible here.
- **Circadian clock genes** shift in PE_in_CD and Diet_in_Veh in the same direction, because the CD+Veh reference mice have a distinct clock phase (see interpretation).

## Course connections

- **Factorial regression** separates main effects from the interaction instead of comparing DEG counts.
- **Reference coding** decides what each coefficient means; getting it wrong silently flips signs, as our pydeseq2 smoke test showed.
- **Multiple testing:** BH-FDR within each contrast. With ~13,700 tests, an unadjusted α = 0.05 would give about 685 false positives under the null.
- **Model diagnostics:** dispersion fit, p-value distribution and influence measures (Cook's distance) are the RNA-seq counterparts of residual and leverage plots in linear regression.
