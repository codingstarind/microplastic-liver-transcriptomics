# Phase 2 methodology proposal: feasibility of three research directions

**Status:** proposal for team review. **No hypothesis tests have been run and no gene-level significance results exist.** This document describes what the data can and cannot support, so the team can finalise aims and hypotheses. Phase 1 code, reports and results are unchanged; new material is limited to:

- this file;
- `scripts/08_methodology_feasibility.py`;
- `results/methodology/`;
- `metadata/sra_run_info.csv`.

---

## 1. New facts established for this proposal

All of these are descriptive (`results/methodology/`). None involves testing genes.

| Check | Result | Consequence |
|---|---|---|
| **Sequencing run** (first FASTQ header of 9/19 runs, both PE and vehicle) | All on instrument `VH00999`, run 78, flowcell `2227J7HNX` | There is no separate sequencing-run batch. The library gap is not "PE sequenced on another day". |
| **Raw read pairs** (ENA) | PE 13.5–16.5 M vs vehicle 17.0–19.1 M. **Perfect separation** (max PE < min vehicle) | The depth gap exists *before* quantification, so it arose at library quantification and pooling. |
| **Fraction of read pairs assigned to genes** | CD+Veh 0.806, CD+PE 0.765, MD+Veh 0.797, MD+PE 0.776 | PE libraries also differ in **quality or composition**, not only depth: ~3–4 points fewer reads land in annotated genes. Possible causes: RNA integrity, more intronic or rRNA reads, or a biological effect. |
| **How much the design explains each technical variable** (R² on diet + PE + diet:PE) | log2 read pairs 0.78; assigned fraction 0.80; library index 0.94; log2 muscle CPM 0.39 | Library index is effectively the design itself. It cannot be adjusted for. Depth or assignment covariates are strongly collinear with PE (VIF ≈ 4.5–5). |
| **Muscle contamination (> 100 CPM) by group** | CD+Veh 0, CD+PE 1, MD+Veh **3**, MD+PE 0 | Concentrated in one cell. MD+Veh is the reference for "PE within MD" and one of the four cells of the interaction. |
| **Biological variability** (method-of-moments NB dispersion, within groups) | Median 0.016 overall; ~0.07 for genes with 10–50 counts; ~0.007–0.01 for 100–5,000 counts (biological CV ≈ 8–10%) | Inbred, well-controlled mice. Modest effects are detectable for well-expressed genes. Dropping MD_1 hardly changes the dispersion. |
| **Design matrix** (intercept, diet, PE, diet:PE) | Full rank (4) | The factorial model is estimable. |

### Minimum detectable effects (planning approximation)

Smallest |log2 fold change| detectable with 80% power, for a gene with a typical dispersion at that expression level. Wald approximation (`results/methodology/minimum_detectable_effects.csv`).

| Comparison | n per arm | Mean 200 counts, α = 0.05 | Mean 200 counts, α = 0.001* | Mean 1,000 counts, α = 0.001* |
|---|---|---|---|---|
| PE within CD | 4 vs 5 | 0.33 | 0.48 | 0.36 |
| PE within MD | 5 vs 5 | 0.31 | 0.46 | 0.34 |
| Diet within vehicle | 4 vs 5 | 0.33 | 0.48 | 0.36 |
| PE main effect (average over diets) | 9 vs 10 | 0.23 | 0.33 | 0.25 |
| **Diet × PE interaction** | all 19 | **0.45** | **0.67** | **0.49** |

\*α = 0.001 roughly represents the per-gene threshold that a genome-wide FDR of 5% implies when a few hundred genes are truly changed. The real threshold depends on the data.

**How to read this:**

- The paper's fold-change cut-off (|log2FC| ≥ 0.38, i.e. 1.3-fold) is **at or below** the detection limit for single-diet PE contrasts at genome-wide stringency. Many published PE DEGs will be near the edge, and their lists will be unstable between analyses.
- An interaction must be a *difference of fold changes* of ~0.5–0.7 log2 (1.4–1.6-fold) to be detectable for typical genes. **Gene-level interaction testing is underpowered, except for large effects.**
- These numbers are optimistic. They ignore dispersion-estimation error with n = 4–5, outliers and contamination, and they assume the model is correct.

---

## 2. Direction A: do PE-associated changes resemble metabolic liver dysfunction (MASH-diet changes)?

### Comparisons required

- **PE signature:** CD+PE vs CD+Veh. PE in an otherwise healthy liver; the cleanest PE contrast, with no diet effect mixed in.
- **Diet (MASH) signature:** MD+Veh vs CD+Veh.
- (Optional) PE within MD: MD+PE vs MD+Veh.
- (Optional) an **external** MASH signature: a published mouse diet-induced MASH dataset, or a human NAFLD/MASH signature.

### Statistical models (options)

1. **Gene-level fold changes from one negative-binomial GLM** (DESeq2: `~ group`, 4 levels), with shrunken log2FCs (apeglm/ashr). Then compare the two fold-change vectors:
   - rank correlation (Spearman) of shrunken LFCs over expressed genes;
   - direction concordance among genes that change under the diet;
   - **rank-rank hypergeometric overlap (RRHO)** to see where along the ranked lists overlap occurs.
2. **Gene-set level.** Run GSEA (fgsea) on each contrast with Hallmark/GO/KEGG sets and compare the normalised enrichment scores. This is closer to "metabolic dysfunction" as a biological concept than single genes are.
3. **Signature score per sample.** Define a "MASH axis" (genes changed by diet) and score each mouse, then compare CD+PE with CD+Veh on that score (a *t*-test or permutation test on one number per mouse).

### Assumptions

- NB GLM: counts are negative-binomial given the mean, dispersion depends smoothly on the mean, samples are independent, and size factors correct for depth and composition.
- Comparison step: the two fold-change vectors have **independent errors**. **This is violated**, because both contrasts use the same CD+Veh mice. Noise in CD+Veh makes the two LFC vectors positively correlated even with no biological similarity.
- Signature score: the signature is defined independently of the samples being scored.

### What it could produce

- A quantitative statement of the form "PE-induced changes are (or are not) concordant with diet-induced changes, with effect size X and uncertainty Y".
- The list of pathways shared or distinct between the two signatures.

### Limitations and required safeguards

- **Shared reference bias** (the main issue). Possible remedies:
  - estimate the null concordance by permuting labels within CD+Veh/CD+PE and MD+Veh/CD+Veh;
  - split CD+Veh between the two contrasts (n = 4 is very small for this);
  - use an **external** MASH signature, which removes the bias entirely and is preferred if a suitable dataset exists.
- **Technical confound.** PE libraries are shallower and have lower assignment. If low-quality libraries shift expression in a direction that happens to resemble MASH (e.g. mitochondrial or high-abundance transcripts), concordance would be inflated. Direction C checks this.
- **Circadian and feeding confound.** Lipogenic genes (Fasn, Acly, Me1, Pnpla3) dominate both signatures and are strongly time-of-day and feeding dependent; Dbp already differs between groups. Concordance driven by these genes may reflect sacrifice timing rather than metabolic disease.
- **Muscle contamination in MD+Veh** (3/5 mice) would add muscle genes to the "diet signature". Exclude muscle-tracking genes from the signature definition, or adjust for the muscle score.
- "Resembles" is not "causes". Concordance in direction does not show a shared mechanism.

### Reproduction vs extension

- **Reproduction:** the paper reports PE DEGs per diet and their Venn overlap, GSEA of a custom lipid-metabolism/toxicity NR gene set per contrast, and a comparison with a PS-microplastic dataset (GSE245069) **for nuclear receptors only**.
- **Extension:**
  - a genome-wide, quantitative PE-vs-diet concordance analysis with a correct null for the shared reference;
  - RRHO;
  - external-signature comparison.

---

## 3. Direction B: does the PE response depend on diet (diet × PE interaction)?

### Comparisons required

All four groups:

- PE effect within CD (CD+PE − CD+Veh);
- PE effect within MD (MD+PE − MD+Veh);
- their difference, which is the interaction.

### Statistical model

**Two-way factorial negative-binomial GLM**, per gene:

```
log(μ_gs) = log(s_s) + β0 + β_diet·MD + β_PE·PE + β_int·(MD×PE)
counts ~ NegativeBinomial(μ, dispersion φ_g)
```

- β_PE is the PE effect within CD; β_PE + β_int is the PE effect within MD; **β_int is the interaction**.
- Test β_int with a Wald test, or a likelihood-ratio test of the full model against `~ diet + PE`. Benjamini–Hochberg FDR across genes.
- Equivalent alternatives: limma-voom with the same design, or edgeR quasi-likelihood. Comparing two of them is a cheap robustness check.
- **Gene-set level interaction:** competitive tests (camera/fry) or fgsea on the interaction statistic. These have more power than single genes because they pool small concordant interactions.
- **Course-friendly view:** two-way ANOVA on per-sample summary scores (pathway scores, module eigengenes or PC scores). One number per mouse lets classical two-way ANOVA assumptions (normality, equal variance) be checked directly.

### Assumptions

- NB counts with gene-wise dispersion shared across all four groups. (DESeq2 assumes equal dispersion across groups; if the MD groups are more variable this is violated.)
- Independence of mice. Group housing violates this: cage = group, so cage effects are inseparable from treatment.
- **Interaction is scale-dependent.** It is defined on the log scale, so a purely additive effect on the log scale is "no interaction" even if the raw-count changes differ.

### What it could produce

- Genes and pathways whose PE response differs between diets, with effect sizes (difference of log2FCs) and confidence intervals.
- A distinction between **additive** (PE adds the same log-change on both diets), **synergistic** and **antagonistic** patterns.

### Limitations

- **Power** (Section 1): only interactions of roughly ≥ 0.5–0.7 log2 are detectable for typical genes. **"No significant interaction" must not be read as "PE acts the same on both diets".** Equivalence testing (TOST) with a pre-chosen margin could address this if it matters to the hypotheses.
- **Contamination concentrated in MD+Veh** directly biases β_int: β_int = (MD+PE − MD+Veh) − (CD+PE − CD+Veh) uses MD+Veh with a negative sign. Muscle genes and anything correlated with them could appear as spurious interactions.
- **The depth/quality confound is balanced across diets** (PE is shallower in both). It mostly cancels in β_int but not in β_PE. Of the three directions, the interaction is the *least* affected by the library problem and the *most* affected by contamination.
- The design is unbalanced (4/5/5/5). This is fine for the GLM, but CD+Veh is the smallest cell and appears in every contrast.

### Reproduction vs extension

- **Reproduction:** the paper reports pairwise contrasts (PE within CD, PE within MD, MD+PE vs CD+Veh) and Venn overlaps of DEG lists. It describes a "compounded" or "synergistic" effect from MD+PE vs CD+Veh, which **cannot distinguish additivity from interaction**.
- **Extension:** a formal interaction model, effect sizes for the difference in PE response, and gene-set-level interaction tests are **not in the paper**. This is the clearest genuine analytical extension.

---

## 4. Direction C: are the observations robust to the technical problems?

This direction is a set of sensitivity analyses applied to whatever A and B estimate. It is not a separate biological question.

| Technical issue | Analysis | What it would show | Limitation |
|---|---|---|---|
| **Library depth / assignment (confounded with PE)** | (i) Compare PE estimates with and without `log(read_pairs)` or `assigned_fraction` as a covariate. (ii) Check whether PE-associated genes are enriched for features linked to library quality: very short or long genes, mitochondrial genes, high-abundance genes, pseudogenes (see the ribosomal pseudogenes on PC2). (iii) Down-sample vehicle libraries to PE depth (binomial thinning of counts) and check that their pattern does not drift toward PE. | Whether "PE effects" are explainable by library properties | The covariate is collinear with PE (VIF ≈ 4.5–5, R² 0.78–0.80). Adjusting roughly doubles the SEs and may remove real PE biology. **The confound cannot be fully resolved statistically; only argued about.** Down-sampling tests depth, not quality. |
| **Muscle contamination** | (i) Muscle score (log2 muscle CPM) as a covariate. (ii) Drop MD_1, or all 4 samples > 100 CPM. (iii) Exclude muscle-tracking genes. (iv) Compare effect estimates across these versions. | Whether conclusions depend on contaminated samples | Dropping 3 MD+Veh mice leaves n = 2 in that cell, which is not viable. A covariate fitted mainly from one cell is weakly identified. |
| **Gene-level Salmon/tximport counts** | See Section 6. Treat as integer counts without a length offset, and optionally re-quantify from FASTQ to compare | Whether results depend on the missing transcript-length information | Re-quantification costs time and disk (~19 runs × ~5 Gbases, roughly 50–80 GB of FASTQ) |
| **Small n / influential mice** | Leave-one-mouse-out refits. Cook's distance (DESeq2 flags outliers automatically, but needs ≥ 3 replicates per cell, which is satisfied). Stability of top-gene rankings (Jaccard overlap across refits). | Whether conclusions rest on single animals | Descriptive only; there is no formal threshold |
| **Circadian / feeding state** | Clock-gene panel (Dbp, Per1–3, Nr1d1/2, Bmal1, Cry1) compared across groups, then a "clock score" covariate | Whether lipogenic findings could be timing artefacts | If PE truly shifts the clock, adjustment removes real biology |
| **Cell composition** | Marker scores or deconvolution; compare PE effects with and without composition covariates | Whether bulk effects are composition shifts | Reference mismatch; approximate |
| **Method dependence** | DESeq2 vs edgeR-QL vs limma-voom on the same design | Whether conclusions depend on the software | Agreement does not prove correctness |

### Reproduction vs extension

None of this is in the paper. The paper reports no sample-level QC, contamination check or sensitivity analysis for the bulk data. **The whole direction is an extension.**

---

## 5. What the publication already did (to avoid uncredited repetition)

| Analysis | In paper? | Our status |
|---|---|---|
| PCA of bulk samples (DESeq2 `plotPCA`) | Yes (Fig. 2A) | Done in Phase 1, extended with sensitivity analyses and a contamination diagnosis |
| Pairwise DEGs per diet, |log2FC| ≥ 0.38, FDR < 0.05 (Fig. 2B–C, tables S5–S7) | Yes | Would be a **reproduction**. Useful as validation of our pipeline if tables S5–S7 can be obtained |
| Venn overlap of PE DEGs under CD and MD | Yes | Reproduction; superseded by the interaction model |
| GSEA of a custom NR lipid-metabolism/toxicity gene set | Yes (Fig. 2H–J) | Reproduction if repeated; the gene set is not public as far as we know |
| NR heatmap; comparison with a PS dataset (GSE245069) for NRs | Yes | Reproduction; genome-wide comparison would extend it |
| Formal diet × PE interaction | **No** | **Extension** |
| Quantitative PE-vs-diet concordance with a correct null | **No** | **Extension** |
| QC, contamination, library confound, sensitivity analyses | **No** | **Extension** |
| Co-expression modules on whole-transcriptome bulk data | **No** (WGCNA only on spatial pseudo-bulk) | **Extension** |

**Validation gap:** the paper's supplementary tables S5–S7 (DEG lists) would let us check that our DESeq2 implementation reproduces their pairwise results before we extend them. Europe PMC's supplementary-file endpoint returned a server error on 2026-11-02, and science.org blocks scripted access. **A team member should download the supplementary materials manually** from the article page.

---

## 6. Can DESeq2 be implemented correctly with the files we have?

### What we have

| Input DESeq2 needs | Status |
|---|---|
| Integer count matrix (genes × samples) | ✅ The GEO file is integer, with no missing values |
| Sample table with design factors | ✅ `metadata/sample_metadata.csv` (diet, microplastic, group; verified against GEO) |
| Full-rank design | ✅ `~ diet * microplastic` has rank 4 |
| ≥ 3 replicates per cell (needed for outlier replacement) | ✅ 4–5 per cell |
| tximport **length matrix** (`avgTxLength`) for `DESeqDataSetFromTximport` | ❌ **Not available.** Only rounded gene-level counts were deposited |
| Transcript-level `quant.sf` files | ❌ Not deposited; only via re-quantification from SRA FASTQ |
| Unique Ensembl gene IDs | ⚠️ Symbols only. 228 symbols merge more than one Ensembl gene; 1,061 symbols are absent |

### Consequence of the missing length matrix

Normally Salmon → tximport → DESeq2 passes an average-transcript-length offset, which corrects for **changes in transcript isoform usage between samples** that alter effective gene length. Without it, `DESeqDataSetFromMatrix` on these counts is the standard fallback; tximport's documentation describes it as acceptable when length offsets are unavailable. The residual bias matters only for genes whose isoform mix changes between conditions. For most genes the impact is expected to be small, but it is **unquantified** for this dataset.

### Software options

| Option | Status on our machine | Notes |
|---|---|---|
| **R + Bioconductor DESeq2** (reference implementation; with apeglm, fgsea, edgeR, limma) | R is **not installed** | The authors' own tool. Needed for exact reproduction, `DESeqDataSetFromTximport`, apeglm shrinkage, edgeR/limma comparisons and RRHO. **Recommended.** |
| **pydeseq2** (Python re-implementation) | See Section 6a | Keeps everything in our Python pipeline. Covers size factors, dispersion, Wald tests, VST and LFC shrinkage. Results are close to but not identical with R DESeq2. Fewer companion tools. |

### Validation required before trusting any DESeq2 results

1. **Size factors:** DESeq2's `estimateSizeFactors` must equal our Phase 1 median-of-ratios values (`results/tables/size_factors.csv`), up to the filtering choice.
2. **Pipeline reproduction:** run the paper's pairwise contrasts with their thresholds and compare with tables S5–S7. Overlap of DEG lists and correlation of LFCs, where the paper reports them. Expect agreement but not identity (DESeq2 versions, filtering, and possibly rounding).
3. **Model diagnostics:** dispersion-trend plot, outlier (Cook's) flags, p-value histograms, MA plots, and residual checks for the interaction model.
4. **Re-quantification cross-check (optional):** Salmon on 2–3 samples from SRA, against GENCODE M35, to confirm we reproduce the deposited counts. This would also yield the length matrix if we decided to re-quantify everything.

### 6a. pydeseq2 smoke test (`scripts/09_pydeseq2_smoke_test.py`, `results/methodology/pydeseq2_smoke_test.md`)

pydeseq2 0.5.4 installs cleanly on our Python 3.14 environment. It was installed in a throwaway environment and is not yet in `requirements.txt`.

1. **Synthetic data with our exact group sizes (4/5/5/5)**, 3,000 genes, 10% with planted effects. The design `~ diet + microplastic + diet:microplastic` fits.

   | Planted effect (log2) | Mean estimate on planted genes | Mean on null genes |
   |---|---|---|
   | Diet 1.5 | 1.51 | −0.01 |
   | PE 1.0 | 1.00 | −0.01 |
   | Interaction 1.5 | 1.48 | −0.02 |

2. **Real data, size factors only:** pydeseq2's size factors are **identical** to our Phase 1 values (maximum relative difference 0). No dispersion fit or tests were run on real genes.

**Pitfall found:** the formula engine (formulaic) orders factor levels **alphabetically**, so `PE` became the reference instead of `Veh`. That flipped the meaning of every PE and interaction coefficient. The first smoke-test run showed it; the fix is to set `Veh` and `CD` as reference levels explicitly (pandas Categoricals). R's DESeq2 has the same behaviour (it needs `relevel`). Any implementation must set reference levels explicitly and include a test for it.

The first synthetic run also showed that when 50% of genes change, median-of-ratios normalisation is biased (null genes shifted by about −0.4 log2). This is expected behaviour, and it reminds us that the method assumes most genes do not change.

**Conclusion:** DESeq2-style modelling is **technically feasible now** with the deposited counts, in R or in pydeseq2. The scientifically important gaps are the missing length offset, the symbol-level aggregation and the lack of the paper's DEG tables for validation; the software is not a gap.

---

## 7. Decisions the team needs to make

1. **Which direction(s) become formal hypotheses**, and whether each is a single gene, a gene set, or a whole-profile question. This determines the multiple-testing burden and the power.
2. **Primary contrast for "PE effect":** within CD only (cleanest), within each diet, or the PE main effect (most powerful but assumes no interaction).
3. **Handling MD_1 and the muscle contamination:** keep everyone with a muscle covariate, drop MD_1, or exclude muscle-tracking genes. **Pre-specify one primary approach** and treat the others as sensitivity analyses.
4. **Library-depth confound:** accept it as a stated limitation, adjust with a covariate (losing roughly half the precision), or use it only as a sensitivity analysis. It cannot be removed completely.
5. **Software:** install R/Bioconductor (exact reproduction, full toolset), or use pydeseq2 within the current Python repository.
6. **Re-quantify from FASTQ or not.** That would resolve the length-offset and symbol-aggregation issues at a substantial time and storage cost.
7. **Reference for "metabolic dysfunction":** the internal diet contrast (shared-reference bias) or an external MASH signature (which dataset?).
8. **Significance and effect-size thresholds:** FDR level; whether to use the paper's |log2FC| ≥ 0.38 (comparability) or a threshold-based test (DESeq2 `lfcThreshold`, more defensible); whether to add equivalence testing for "no interaction".
9. **Gene-set collections** for pathway analysis (Hallmark, GO BP, KEGG, Reactome) and minimum and maximum set sizes, chosen before seeing results.

## 8. Suggested implementation order (after decisions)

1. Set up DESeq2 (R or pydeseq2); validate size factors; reproduce the paper's pairwise contrasts (needs tables S5–S7).
2. Fit the pre-specified primary factorial model `~ diet * microplastic` (+ any agreed covariates). Run diagnostics.
3. Direction B: interaction effect sizes; gene-set interaction tests.
4. Direction A: shrunken LFCs → concordance with a permutation or external null; GSEA comparison.
5. Direction C: the sensitivity-analysis table re-run for every headline result.
6. Report effect sizes with uncertainty first and significance second, and state every limitation listed here.
