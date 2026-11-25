# Reference-axis comparison and collinearity

**Script:** `scripts/f01_compare_and_baselines.py`. **Tables:** `results/comparisons/reference_spearman_matrix.csv`, `reference_pairwise.csv`, `core_pathway_median_log2FC.csv`, `comparison_summary.json`. **Figure:** `figures/reference_similarity_heatmap.png`.

All comparisons are rank- and sign-based on the pairwise overlap within the 11,888-gene universe. Different technologies (RNA-seq counts, RPKM, TempO-Seq) are **never compared on magnitude**.

## Pairwise Spearman ρ of gene effects (observed)

| | MASH1 | MASH2 | FAST | FAST2 | PPARA 20v0 | CIRC T1 | CIRC T2 |
|---|---|---|---|---|---|---|---|
| MASH1 | 1 | **0.56** | 0.06 | −0.17 | −0.01 | −0.12 | −0.05 |
| MASH2 | | 1 | −0.04 | −0.19 | 0.14 | −0.13 | −0.02 |
| FAST | | | 1 | **0.49** | 0.04 | −0.07 | 0.06 |
| PPARA 20v0 | | | | | 1 | −0.12 | −0.12 |
| CIRC T1 | | | | | | 1 | **0.86** |

Restricted to the 1,893 PPARα-dependent genes: FAST ρ = 0.24, MASH2 ρ = 0.33.

**Collinearity** (5 primary signed axes on 11,749 common genes, rank-standardised):

| Axis | VIF |
|---|---|
| MASH1 | 1.03 |
| FAST | 1.07 |
| PPARA 20v0 | 1.02 |
| CIRC T1 | 4.07 |
| CIRC T2 | 4.04 |

The only flagged pair (|ρ| ≥ 0.7) is CIRC T1/T2.

## Answers
- **Does MASH resemble fasting?** No, genome-wide (ρ 0.06 / −0.04). Both lower cholesterol synthesis, but they differ in lipogenesis (fasting strongly ↓, MASH1 ↑ vs chow) and inflammation/EMT (MASH ↑, fasting ↓ or flat).
- **Does fasting resemble PPARα activation?** Weakly genome-wide (ρ 0.04), moderately within PPARα-dependent genes (ρ 0.24).
  - Both raise PPARα-target sets (fasting +0.82, GW7647 +1.50 median log2).
  - **They are opposite for lipogenesis** (fasting −2.79, GW7647 +1.15) **and cholesterol synthesis** (fasting −2.76, GW7647 +0.67).
  - So "PPARα targets up with lipogenesis and cholesterol down" is a fasting-like combination, **not** what this PPARα agonist reference produces.

  This is an observation about the references. It is not a claim about PE.
- **Are lipogenesis and cholesterol effects shared across states?**
  - Cholesterol synthesis is down in fasting and in both MASH references, but up with GW7647.
  - Lipogenesis is down only in fasting; it is up in MASH-vs-chow and with GW7647.
- **Collinearity:** MASH, fasting and PPARα are **distinguishable** axes. VIFs near 1 mean separate coefficients would be estimable in Phase 4. The two circadian contrasts must not both enter one model; choose one, or use the harmonic coefficients.

## Caveats
- Genome-wide ρ near 0 does not mean the programmes are unrelated. It reflects many genes with small or noisy effects; the programme-level table is more informative for the shared metabolic genes.
- The references differ in sex, strain, assay, comparator diet and the extremity of the state. Similarity is about direction patterns only.
