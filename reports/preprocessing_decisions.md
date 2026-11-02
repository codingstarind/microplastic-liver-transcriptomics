# Preprocessing decisions

For each operation that might be applied to the count matrix, this file records:

- its purpose;
- whether it is appropriate here;
- what we actually did.

The raw file is never changed, and every derived matrix is a separate, named file (see `metadata/data_dictionary.md`).

| Operation | Purpose | Appropriate here? | Decision |
|---|---|---|---|
| **Rounding** | Count models (DESeq2, edgeR) need integers | Already integers | **Not done** (nothing to do) |
| **Filtering low-expression genes** | Remove genes with too little information. They add noise to PCA/clustering, inflate multiple testing later, and their log values are dominated by the pseudocount | Yes. 53% of genes are all-zero; most others have a handful of reads | **Done:** keep genes with ≥ 10 counts in ≥ 4 samples (4 = smallest group, so a gene expressed in only one group survives). 13,777 genes kept, 99.86% of counts retained |
| **Library-size normalisation** | Raw counts scale with sequencing depth; libraries here differ 1.5-fold (10.4–15.5 M). Without it, depth differences would look like expression differences | Essential for any between-sample comparison | **Done:** DESeq2 **median-of-ratios** size factors (0.78–1.26) on the filtered genes. Robust to a few very abundant or strongly changed genes, unlike total-count scaling. Size factors track library size closely (Spearman ρ = 0.95) |
| **CPM** | Simple depth scaling | Fine for display, less robust when a few genes (*Alb*, *Mup*s) take a large share | **Computed** for a sensitivity PCA and for marker scores only |
| **log2(x + 1)** | Count variance grows with the mean (fig05A), so untransformed PCA would be dominated by a few highly expressed genes. The log makes fold-changes additive and compresses the 0 – 900,000 range | Yes, for EDA, PCA and clustering. Not as input to count-based tests | **Done** on normalised counts. Pseudocount 1. Drawback: low-count genes get inflated variance (fig05B). The filter limits this, and a VST (below) would handle it better |
| **Variance-stabilising transformation (VST / rlog)** | Like log2 but shrinks the noisy low-count end. It is DESeq2's recommendation for PCA | Would be the better choice | **Not done yet:** DESeq2 is an R package and R is not installed on the analysis machine. Planned for Phase 2 (R, or `pydeseq2`), as a comparison against log2(x + 1) |
| **Gene-wise standardisation (z-scores)** | Gives each gene equal weight | Not by default: it would let thousands of noisy, barely-varying genes weigh as much as real signals | **Not done** for PCA (genes centred only). Correlation heatmaps are scale-free anyway |
| **Selecting the most variable genes** | Focus PCA and clustering on informative genes | Common practice, but the choice of N shapes the result | **Done:** top 500 (primary). **Sensitivity:** 1,000, 2,000, all 13,777 |
| **Removing contaminated samples or genes** | Muscle contamination dominates PC1 (see the EDA report) | Removing data is a strong decision and should be agreed by the team | **Not done** in the primary matrices. **Sensitivity analyses only:** (a) remove 135 genes whose expression tracks the muscle score (\|Spearman ρ\| ≥ 0.7), (b) drop MD_1. Note that (a) is crude: with 19 samples, some non-muscle genes pass the threshold by chance or through diet confounding |
| **Batch correction** | Remove technical batch effects | No batch variable is recorded (all "B1"). Library order is fully confounded with treatment, so "correcting" it would remove the treatment signal | **Not done.** The confounding is reported as a limitation |
| **Gene-length normalisation (TPM/FPKM)** | Compare different genes within a sample | Not needed to compare the same gene across samples | **Not done** |

## Open choices for the team (Phase 2)

1. Should MD_1 be excluded, down-weighted, or kept with a muscle covariate in formal models?
2. Should log2(x + 1) be replaced by VST, which needs R/DESeq2?
3. Should we reprocess the FASTQs (SRA) to recover transcript-level information and Ensembl-ID gene counts? This resolves audit issues 4–5, but costs time (19 × ~15 M read pairs).
