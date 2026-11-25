# GSE316023 count input and FPKM-derived length offset: validation

**Script:** `scripts/a01_offset_audit.py`. **Outputs:** `results/offset_audit/audit.json`, `offset_gene_table.csv.gz`, `figures/offset_audit.png`. **Scope:** the 16 primary-contrast samples (GAN wk23 TA ×8, chow wk23 TA ×8) × the 11,888-gene main universe (11,872 present in GSE316023).

The decision rule was fixed in `analysis_plan.md` §2 (locked, SHA-256 in `results/plan_lock.json`) before any model was fitted.

## Verdict
**The FPKM-derived offset cannot be justified as a tximport-equivalent length correction. The primary MASH model is DESeq2 on the rounded deposited estimated counts without a length offset.** The offset model is retained as an **approximate sensitivity model**. Its effect on the conclusions is small (§4).

There are two independent reasons:
1. **The provenance is undocumented.** The quantification software, the length definition (effective vs annotated, transcript vs union-gene) and the FPKM denominator are not stated by the authors (see check 2).
2. **Check (ii) failed narrowly.** Implied lengths correlate only moderately with annotated transcript lengths (ρ = 0.497).

   Caution: the threshold of 0.50 was written into the audit script before it ran, but it is **not** in the locked plan, which said only "correlates with GENCODE transcript length". Reason 1 alone already sets the decision under the plan, so this threshold does not change the outcome.

## Checks

| # | Question | Result |
|---|---|---|
| 1 | Same pipeline, samples and genes? | Identical sample columns and gene rows in identical order. The zero pattern matches except for 144 cells (all with count > 0 and FPKM = 0 because of 4-decimal rounding; none in the 16 × universe block). The column scaling is exact (next row), so the two files are consistent with one pipeline. |
| 2 | FPKM calculation documented? | **No.** The GEO `data_processing` says only "FPKM was generated" (mm10). The paper's Methods (PMC13392390) give library and sequencing details but no quantification method. |
| 3 | Lengths identifiable up to a sample constant? | Implied L = count·10⁹/(FPKM·N) with N = the count column total. The per-sample scale factors are 1.0000–1.0006, so the FPKM denominator is exactly the column total and the column constant is ~0% of the between-sample variance. What remains is gene-by-sample variation (residual SD 0.094 on the log scale). The lengths are identifiable; the sample constant would in any case be absorbed by size factors. |
| 4 | Unstable estimates from small FPKM? | 4-decimal FPKM gives a relative length error > 1% in 53 cells (25 detected genes). 195 universe genes lack a length in ≥ 1 sample (zero count); they get the gene's mean length. 16 genes get no offset at all. **This is minor.** |
| 5 | Consistent with what DESeq2/tximport assume? | tximport's `avgTxLength` is the abundance-weighted average of **transcript effective lengths** per sample (tximport vignette; `DESeqDataSetFromTximport` stores it and DESeq2's `estimateSizeFactors` turns it into normalisation factors via `estimateNormFactors`). Here the median implied length is **1.52× the GENCODE median transcript length**, with ρ = 0.497. That pattern looks more like a gene-level length (for example the union of exons, or another annotation release) than an abundance-weighted effective transcript length. **It cannot be confirmed to be the quantity tximport expects.** |
| 6 | Normalisation factors built correctly, without double normalisation? | The offset model passes the raw lengths as the `avgTxLength` assay of a `DESeqDataSetFromMatrix` object. DESeq2 itself does the row-centring and the size-factor estimation on counts/normMatrix (the tximport-recommended route). No pre-normalised counts or FPKM enter the model, so nothing is normalised twice. |
| 7 | Factors finite, positive, centred? | Yes. Row geometric mean = 1 (max deviation 3×10⁻¹⁵); range 0.096–4.28; 1st–99th percentile 0.66–1.34. |
| 8 | Does rounding matter? | 6.2% of cells are non-integer. The maximum change to any library total is **0.0012%**. The fractional parts are genuine (spread over 0–0.5), consistent with EM-allocated counts. Rounding is negligible. |
| 9 | tximport-style or approximate? | **Approximate.** The structure is right (per-sample, per-gene lengths with an exact column scale), but the length definition is unverifiable. |

**How much would length matter?** After removing the sample scale, 304 universe genes (all counts ≥ 10) differ by more than 1.2× and 1,642 by more than 7% in mean implied length between GAN and chow. Without a valid offset, these genes' diet effects could partly reflect isoform or length changes. They are flagged.

## Consequences (observed after the decision)
- No-offset vs offset log2FC: Spearman 0.981. Among genes significant in both models, sign agreement is 100% (5,571 genes).
- **883 genes are `offset_sensitive`** (DEG status or sign at padj < 0.05 differs). They are listed in `results/positive_control/offset_sensitive_genes.csv.gz`. **No gene reverses significant direction.**
- Positive control ρ with Diet_in_Veh: 0.355 without the offset vs 0.356 with it. The MASH positive-control outcome does not depend on the offset choice.
- Assumptions required for the primary model:
  - the deposited values are estimated counts;
  - rounding is harmless (verified);
  - the length bias between conditions is ignorable at the signature level (supported by the 0.98 agreement, but **not** gene by gene for the 883 flagged genes);
  - independent mice;
  - no hidden batch (unverifiable).
