# Phase 2 methodology proposal: feasibility of three research directions

**Status:** proposal for team review (revision 2, after review of PR #2).

**No hypothesis tests have been run and no gene-level significance results exist.** This document describes what the data can and cannot support, so the team can finalise aims and hypotheses. Phase 1 code, reports and results are unchanged.

New material:

- this file;
- `src/mpliver/feasibility.py` (with tests in `tests/test_feasibility.py`);
- `scripts/08_methodology_feasibility.py`;
- `scripts/09_pydeseq2_smoke_test.py`;
- `scripts/validate_methodology_outputs.py`;
- `results/methodology/`;
- `metadata/sra_run_info.csv`;
- `requirements-methodology.txt`;
- two extra CI checks.

### What changed in revision 2

| Review point | Correction |
|---|---|
| Flowcell evidence based on 9/19 runs | All 19 runs now read. Sequencing-run, library-preparation and pooling effects are discussed separately (§1.1). |
| `assigned_fraction` presented as a real proportion | Renamed `gene_counts_per_read_pair_proxy` and qualified as a proxy. Salmon mapping statistics cannot be recovered from deposited files (§1.2). |
| Count-matrix fallback presented as acceptable | Rewritten. The tximport vignette advises against using original gene-level counts without an offset. The `countsFromAbundance` setting cannot be determined from the deposited files, and the fallback is **not validated** (§6). |
| Cook's-distance documentation wrong | Corrected: automatic outlier **replacement** needs ≥ 7 replicates per group (default `minReplicatesForReplace = 7`). With 4–5 per group, DESeq2 only **flags** outliers by setting p-values to NA (§4, §6). |
| Dispersion and power presented too confidently | Labelled as planning approximations, with the estimator's limitations, the fixed-α approximation, gene-to-gene spread and leave-one-mouse-out uncertainty (§1.3). |
| Smoke test conflated recovery with validation | Separated: it checks coefficient recovery and reference levels only, not calibration or FDR control (§6a). |
| Supplementary tables | Retrieval attempted; PMC served a bot challenge. Still unresolved (§5). |
| CI only checked Phase 1 | Added unit tests for the feasibility functions, a validator for the committed methodology outputs, and a CI job running the pydeseq2 synthetic smoke test (§8). |

---

## 1. Facts established for this proposal (descriptive only)

### 1.1 Sequencing run vs library preparation vs pooling

These are three different sources of technical variation:

| Level | What would differ | Evidence available | Conclusion |
|---|---|---|---|
| **Sequencing run / flowcell** (instrument, run chemistry, cluster density) | Systematic shifts between samples on different runs | First read header of the R1 FASTQ for **all 19 runs** (ENA, 2 Nov 2026): every one reads `VH00999:78:2227J7HNX:1` (same NextSeq 2000, run 78, flowcell 2227J7HNX, lane 1) | **No evidence of a sequencing-run batch.** Caveat: only the *first* read of each file was checked; a run whose reads came from more than one flowcell cannot be excluded without scanning whole files (unlikely for a single NextSeq run). |
| **Library preparation** (RNA extraction day, poly(A) selection, fragmentation, PCR, operator) | Library complexity, insert size, rRNA or intronic content, 3′ bias | None. GEO and the paper give no extraction or prep dates. Library indices S39–S57 follow group order (CD → MD → CD+PE → MD+PE). | **Unknown.** It cannot be ruled in or out. Processing order is confounded with treatment, so any prep-order effect would align with PE. |
| **Pooling / loading** (how much of each library went into the pool) | Number of reads per sample | Raw read pairs (ENA): PE 13.5–16.5 M vs vehicle 17.0–19.1 M; **perfect separation** | **Unequal pooling is the simplest explanation** for the depth gap. Depth alone is handled by size factors. |

### 1.2 Gene counts per read pair (a proxy, not a mapping rate)

`gene_counts_per_read_pair_proxy` = (sum of deposited gene counts) ÷ (raw read pairs in SRA).

| Group | Mean | Range |
|---|---|---|
| CD+Veh | 0.806 | 0.798–0.810 |
| CD+PE | 0.765 | 0.754–0.768 |
| MD+Veh | 0.797 | 0.775–0.811 |
| MD+PE | 0.776 | 0.768–0.784 |

**Why it is only a proxy:**

- The denominator is raw read pairs, *before* fastp filtering.
- The numerator is Salmon's fragment estimates, summed to symbols and rounded. It also depends on the tximport setting (§6), on how multi-mapping fragments were distributed, and on the 1,061 GENCODE symbols missing from the matrix.

A lower value in PE libraries is *consistent with* lower library quality, more intronic/intergenic/rRNA reads, more reads removed by fastp, or a composition change. It does not show which.

**Can real Salmon statistics be recovered?**

- **No, not from the deposited files.** GEO holds only the gene-level matrix. The Salmon outputs that contain mapping rates and fragment-length distributions (`quant.sf`, `aux_info/meta_info.json`, `lib_format_counts.json`) and the fastp reports were not deposited, and neither paper nor GEO reports mapping rates.
- **They can be regenerated** by re-running fastp v0.23.4 and Salmon v1.10.3 against the GENCODE M35 transcriptome on the SRA FASTQs (§6).

### 1.3 Confounding and collinearity

R² is how much of each technical variable the factorial design explains (`~ diet + PE + diet:PE`):

| Technical variable | R² | Implied VIF if added as covariate |
|---|---|---|
| log2 read pairs | 0.78 | 4.5 |
| gene counts per read pair (proxy) | 0.80 | 4.9 |
| library index | 0.94 | 16.3 |
| log2 muscle-marker CPM | 0.39 | 1.6 |

- Library index is effectively the design itself and cannot be adjusted for.
- Depth or proxy covariates would inflate the standard errors of the PE coefficients by roughly √VIF ≈ 2.1–2.2×.
- Muscle contamination (> 100 CPM) by group: CD+Veh 0, CD+PE 1, **MD+Veh 3**, MD+PE 0.
- The design matrix is full rank (4).

### 1.4 Dispersion and minimum detectable effects: rough planning approximations

**Estimator.** Pooled within-group method of moments on size-factor-normalised counts: φ ≈ (var − μ·mean(1/s)) / μ². It is used only to get an order of magnitude. Its limitations:

- **Noise per gene.** Each group has only 3–4 residual degrees of freedom (15 in total per gene), so single-gene estimates are very noisy. Only summaries over thousands of genes (medians, quartiles) are reported.
- **Bias.** Negative estimates are clipped to 10⁻⁴, which biases low-dispersion genes upward. Normalised counts are treated as if they were NB with a Poisson term μ·mean(1/s), which is an approximation.
- **Not DESeq2's estimator.** There is no Cox–Reid adjustment and no empirical-Bayes shrinkage toward a mean–dispersion trend. DESeq2's final dispersions will differ, typically being shrunk toward the trend.
- **Contamination and outliers inflate it.** For low-count genes, CD_MP_2 is the single mouse that moves the median most; MD_MP_4 does so for high-count genes.

**Results** (`results/methodology/dispersion_by_expression.csv`):

| Mean normalised count | Genes | Dispersion median [q25–q75 across genes] | Leave-one-mouse-out range of the median |
|---|---|---|---|
| 10–50 | 2,945 | 0.069 [0.030–0.160] | 0.058–0.071 |
| 100–500 | 5,093 | 0.0094 [0.0049–0.0195] | 0.0089–0.0095 |
| 500–1k | 1,565 | 0.0069 [0.0040–0.0144] | 0.0066–0.0072 |

- The **median** is stable to dropping any one mouse (narrow leave-one-out range).
- The **gene-to-gene spread** is wide: individual genes differ several-fold from the median.

An earlier draft of this revision used a bootstrap over mice. It was dropped because resampling 4–5 mice with replacement duplicates animals and biases within-group variance downward; the bootstrap intervals fell entirely below the point estimate.

**Minimum detectable |log2FC|** (80% power, two-sided Wald test, normal approximation). Shown for typical genes at mean 200 and 1,000 normalised counts. Columns use the median dispersion, with the [q25–q75] gene-to-gene range in brackets:

| Comparison | n per arm | Mean 200, α = 0.05 | Mean 200, α = 0.001 | Mean 1,000, α = 0.001 |
|---|---|---|---|---|
| PE within CD | 4 vs 5 | 0.32 [0.27–0.42] | 0.48 [0.40–0.63] | 0.36 [0.28–0.50] |
| PE within MD | 5 vs 5 | 0.31 [0.25–0.40] | 0.45 [0.37–0.59] | 0.34 [0.27–0.47] |
| Diet within vehicle | 4 vs 5 | 0.32 [0.27–0.42] | 0.48 [0.40–0.63] | 0.36 [0.28–0.50] |
| PE main effect (average over diets) | 9 vs 10 | 0.22 [0.19–0.29] | 0.33 [0.27–0.43] | 0.24 [0.19–0.34] |
| **Diet × PE interaction** | all 19 | **0.45 [0.37–0.58]** | **0.66 [0.55–0.86]** | **0.49 [0.39–0.68]** |

For low-count genes (mean ~50), all values roughly double: 1.1–1.6 log2 at α = 0.001.

**Why these are rough approximations:**

1. **Fixed α.** α = 0.001 is only an *illustrative* stand-in for a genome-wide FDR threshold. The per-gene p-value cut-off that BH-FDR actually implies depends on how many genes are truly changed and on the p-value distribution, both unknown until the analysis is run.
2. **Normal approximation with tiny n.** Real Wald statistics with 4–5 mice per group and estimated dispersions behave more like *t* than normal. That widens the detectable effect.
3. **Ignores dispersion-estimation error, LFC shrinkage, independent filtering, and Cook's-distance NA flags.**
4. **19 mice.** Every number here is computed from the same 19 animals whose effects we want to detect. They describe this dataset's noise level, not a population.

**How to use them:**

- The paper's cut-off (|log2FC| ≥ 0.38 ≈ 1.3-fold) is close to the detection limit for single-diet PE contrasts at genome-wide stringency, so published DEG lists near that edge are likely to be unstable.
- Gene-level interaction effects must be roughly ≥ 0.5–0.9 log2 to be detectable.

---

## 2. Direction A: do PE-associated changes resemble metabolic liver dysfunction (MASH-diet changes)?

**Comparisons.**

- PE signature: CD+PE vs CD+Veh (cleanest, no diet mixed in).
- Diet signature: MD+Veh vs CD+Veh.
- Optional: PE within MD.
- Optional: an **external** MASH signature (published mouse diet-induced MASH dataset or a human NAFLD/MASH signature).

**Models (options).**

1. A single NB GLM (`~ group`) with shrunken LFCs (apeglm/ashr). Then compare the two LFC vectors by Spearman correlation, direction concordance among diet-responsive genes, and rank-rank hypergeometric overlap (RRHO).
2. Gene-set level: run GSEA on each contrast and compare normalised enrichment scores.
3. Per-mouse signature score (one number per mouse), then a two-group comparison.

**Assumptions.**

- NB distribution; dispersion is a smooth function of the mean; mice are independent; size factors are adequate.
- For the comparison step, the two LFC vectors must have independent errors. **This is violated:** both contrasts share the CD+Veh mice, so noise in CD+Veh induces positive correlation even with no biological similarity.

**Possible outputs.** An effect-size statement of PE–diet concordance with uncertainty; shared vs distinct pathways.

**Limitations and safeguards.**

- **Shared reference bias.** Use a permutation null that respects the shared reference, or an external signature (preferred if a suitable dataset exists).
- **Technical confound.** PE libraries differ in depth and in the quality proxy (§1.1–1.2). If library quality shifts expression toward a MASH-like pattern, concordance is inflated. Direction C checks this.
- **Circadian and feeding confound.** Lipogenic genes dominate both signatures, and *Dbp* differs between groups (Phase 1).
- **Muscle contamination in MD+Veh** contaminates the diet signature. Exclude or adjust for it.
- Concordant direction does not establish a shared mechanism.

**Reproduction vs extension.**

- **Reproduction:** PE DEGs per diet, Venn overlap, custom NR-gene-set GSEA, and a polystyrene (GSE245069) comparison restricted to nuclear receptors.
- **Extension:** a genome-wide quantitative concordance analysis with a correct null, RRHO, and an external-signature comparison.

---

## 3. Direction B: does the PE response depend on diet (diet × PE interaction)?

**Comparisons.** All four groups:

- PE effect within CD;
- PE effect within MD;
- their difference (the interaction).

**Model.** Two-way factorial NB GLM, per gene:

```
log(μ_gs) = log(s_s) + β0 + β_diet·MD + β_PE·PE + β_int·(MD×PE)
counts ~ NegativeBinomial(μ, φ_g)
```

- β_PE is the PE effect within CD; β_PE + β_int is the PE effect within MD; β_int is the interaction.
- Test β_int by Wald test, or by LRT of the full model against `~ diet + PE`. Apply BH-FDR across genes.
- Alternatives for robustness: edgeR quasi-likelihood, limma-voom.
- Gene-set interaction tests (camera/fry or fgsea on the interaction statistic).
- Two-way ANOVA on per-mouse summary scores, which lets the classical assumptions be checked directly.
- **Reference levels must be set explicitly** (CD, Veh). See §6a.

**Assumptions.**

- NB counts with **one dispersion per gene shared across all four groups.** If MD groups are more variable, this is violated.
- Independent mice. Violated by group housing, since cage equals group.
- The interaction is defined on the **log scale**.

**Possible outputs.** The difference in PE response with confidence intervals, separating additive, synergistic and antagonistic patterns; gene sets with diet-dependent PE responses.

**Limitations.**

- **Power** (§1.4): only large interactions are detectable. "No significant interaction" ≠ "same response". Equivalence testing (TOST) with a pre-set margin would be needed to claim similarity.
- **Contamination is concentrated in MD+Veh**, which enters β_int directly. This could create spurious interactions in muscle-correlated genes.
- **The depth and quality confound is balanced across diets**, so it largely cancels in β_int (not in β_PE).
- **Cook's distance.** With 4–5 mice per group, DESeq2 will **not replace** outlier counts (replacement needs ≥ 7 replicates). Genes where one mouse has a large Cook's distance get **NA p-values**. MD_1 (heavy contamination) could therefore remove genes from the interaction test rather than distort them. Count and report those NA genes.

**Reproduction vs extension.**

- **Reproduction:** pairwise contrasts and Venn overlaps. The paper's "compounded" MD+PE effect is based on MD+PE vs CD+Veh, which cannot distinguish additivity from interaction.
- **Extension:** a formal interaction model, interaction effect sizes and gene-set interaction tests. This is **the clearest genuine extension**.

---

## 4. Direction C: robustness to technical problems

These are sensitivity analyses applied to whatever A and B estimate:

| Issue | Analysis | What it shows | Limitation |
|---|---|---|---|
| **Depth / pooling** | Compare PE estimates with and without `log(read_pairs)`; thin vehicle libraries to PE depth (binomial thinning) and check that their profile does not drift toward PE | Whether depth alone can mimic PE effects | The covariate is collinear with PE (VIF ≈ 4.5). Thinning tests depth, not quality. |
| **Library quality / composition** | Re-quantify to get real mapping statistics (§6). Check whether PE-associated genes are enriched for quality-sensitive features (gene length, GC, mitochondrial, very abundant, pseudogenes). Add a quality covariate. | Whether "PE effects" track library quality | **Library-prep order is confounded with PE.** This can be argued about but not statistically resolved. |
| **Muscle contamination** | Muscle score as covariate; drop MD_1; exclude muscle-tracking genes; compare estimates | Dependence on contaminated samples | Dropping all 3 contaminated MD+Veh mice leaves n = 2, which is not viable. A covariate fitted mainly from one cell is weakly identified. |
| **Count format (tximport)** | Re-quantify a subset and compare against the deposited counts (§6) | Whether the missing length offset matters | Requires FASTQ processing |
| **Small n / influential mice** | Leave-one-mouse-out refits; stability of top-gene rankings; report Cook's-NA genes | Whether conclusions rest on single animals | Descriptive only |
| **Circadian / feeding state** | Clock-gene panel, then a clock-score covariate | Whether lipogenic findings are timing artefacts | Adjustment may remove real biology if PE shifts the clock |
| **Cell composition** | Marker scores or deconvolution as covariates | Whether bulk effects are composition shifts | Reference mismatch |
| **Method dependence** | DESeq2 vs edgeR-QL vs limma-voom | Software dependence | Agreement ≠ correctness |

**Reproduction vs extension.** None of this is in the paper; the whole direction is an extension.

---

## 5. What the publication already did

| Analysis | In paper? | Our status |
|---|---|---|
| PCA (DESeq2 `plotPCA`) | Yes, Fig. 2A | Done in Phase 1, with sensitivity analyses and a contamination diagnosis |
| Pairwise DEGs, \|log2FC\| ≥ 0.38, FDR < 0.05 (tables S5–S7) | Yes | Reproduction; useful for validating our pipeline |
| Venn overlap of PE DEGs under CD and MD | Yes | Reproduction; superseded by the interaction model |
| GSEA on a custom NR lipid-metabolism/toxicity gene set | Yes | Reproduction; the gene set may be in the supplement |
| NR heatmap; PS-dataset comparison for NRs | Yes | Reproduction; a genome-wide comparison would extend it |
| Formal diet × PE interaction | **No** | **Extension** |
| Quantitative PE-vs-diet concordance with a correct null | **No** | **Extension** |
| QC, contamination, library confound, sensitivity analyses | **No** | **Extension** |
| Whole-transcriptome bulk co-expression | **No** (WGCNA on spatial pseudo-bulk only) | **Extension** |

**Supplementary tables: unresolved.**

- `https://pmc.ncbi.nlm.nih.gov/articles/instance/13274612/bin/sciadv.aec8681_tables_s1_to_s18.zip` returned a "Preparing to download" proof-of-work bot challenge (HTTP 200, 1.8 KB HTML), not the zip. Europe PMC's supplementary endpoint returned HTTP 500.
- No table contents are described or used anywhere in this proposal.
- A team member needs to download the zip in a browser. It is needed to validate any DESeq2 reproduction against tables S5–S7, and possibly to obtain the custom gene set.

---

## 6. Can DESeq2 be implemented correctly with the files we have?

### Inputs

| DESeq2 needs | Status |
|---|---|
| Integer count matrix | ✅ The GEO matrix is integer, with no missing values |
| Sample table with factors | ✅ `metadata/sample_metadata.csv` |
| Full-rank design | ✅ Rank 4 for `~ diet * microplastic` |
| Cook's-distance outlier **replacement** | ❌ Not available: needs ≥ 7 replicates per group (`minReplicatesForReplace = 7`); we have 4–5. DESeq2 still computes Cook's distances and **flags** outliers by setting p-values to NA (flagging applies to samples in groups with ≥ 3 replicates). pydeseq2 behaves the same (`min_replicates = 7`). |
| tximport length matrix (`avgTxLength`) for `DESeqDataSetFromTximport` | ❌ Not deposited |
| Transcript-level `quant.sf`, Salmon/fastp logs | ❌ Not deposited; only via re-quantification |
| Unique Ensembl gene IDs | ⚠️ Symbols only (228 symbols merge more than one Ensembl gene; 1,061 symbols absent) |

### The tximport length-offset problem

**What the tximport vignette says** (Bioconductor, read 2 Nov 2026). There are two supported routes for count-based DGE tools:

- *original estimated counts plus an offset* calculated from the transcript-level abundances (`DESeqDataSetFromTximport` does this automatically); or
- *bias-corrected counts without an offset*, by running tximport with `countsFromAbundance = "scaledTPM"` or `"lengthScaledTPM"` and using that matrix as an ordinary count matrix.

It explicitly advises **against** passing the original gene-level counts to downstream methods without an offset. The only exceptions it gives are data without length bias, such as 3′-tagged RNA-seq or protocols with little or no fragmentation. **Our data are full-length, fragmented poly(A) RNA-seq, so that exception does not apply.**

**What the deposited matrix could be:**

| tximport setting | What the counts are | Correct downstream use |
|---|---|---|
| `"no"` (the tximport default) | Σ Salmon NumReads per gene; depends on sample-specific effective lengths | Needs the length offset, **which we do not have** |
| `"scaledTPM"` | TPM scaled to library size; no length dependence | May be used directly as counts |
| `"lengthScaledTPM"` | TPM × average transcript length (across samples) scaled to library size | May be used directly as counts |
| (any of the above, then rounded) | Integers, as deposited | Same as above |

**What can and cannot be established from the deposited files:**

- **Established:**
  - the values are integers, so they were rounded after tximport;
  - GEO calls them "Raw counts";
  - GEO says tximport was used, without naming the setting;
  - the paper says DESeq2 size factors were estimated with `estimateSizeFactors`.
- **Not established:** which `countsFromAbundance` setting was used, and whether the authors' own DESeq2 run used the length offset.
  - Column totals cannot discriminate between the settings, because all three preserve the library total.
  - Gene-level ratios cannot either: the settings differ only through per-sample transcript lengths, which were not deposited.
  - The "Raw counts" label and the tximport default (`"no"`) make the *offset-requiring* case **plausible, but this is not confirmed**.
- **Consequence:** feeding this matrix to `DESeqDataSetFromMatrix` without an offset is **an unvalidated fallback, not an established equivalent.**
  - If the matrix is `"no"`-type counts, it is the configuration the vignette advises against.
  - The resulting bias affects genes whose effective length changes between conditions (isoform switching, changes in fragment-length or 3′ bias). PE libraries differ in the quality proxy, so a library-quality-linked length bias that correlates with PE is a real possibility, not a negligible one.
  - The size of this bias in our data is **unknown**.

**How to resolve it (empirically).**

1. Re-run fastp v0.23.4 and Salmon v1.10.3 (selective alignment, `-l A --seqBias --useVBOpt`, GENCODE M35 transcriptome) on a **subset of 3–4 samples**.
2. Run tximport three times (`"no"`, `"scaledTPM"`, `"lengthScaledTPM"`), summarising to GENCODE symbols, and round.
3. Compare each version with the deposited columns for those samples. An exact or near-exact match identifies the setting. This also yields real mapping rates (§1.2).
4. If the deposited matrix matches `"no"`, either re-quantify all 19 samples and use `DESeqDataSetFromTximport` (the proper offset), or run our own tximport with `lengthScaledTPM`.

**Cost:**

- 38 FASTQ files, 44.6 GB compressed in total (~2.3 GB per sample).
- Salmon has no native Windows build. Docker is installed on the analysis machine, but its engine was not running at the time of checking; alternatives are WSL or a Linux machine.
- A 3–4 sample pilot is ~10 GB of downloads.

### 6a. pydeseq2 smoke test: what it does and does not show

`scripts/09_pydeseq2_smoke_test.py`; output in `results/methodology/pydeseq2_smoke_test.md`. pydeseq2 0.5.4 is pinned in `requirements-methodology.txt` (not yet adopted).

**What it checks:**

- **Coefficient recovery on synthetic data.** Real group sizes (4/5/5/5), 3,000 genes, 10% with planted effects, constant dispersion 0.015. The design `~ diet + microplastic + diet:microplastic` fits, and the *average* estimates over planted genes are 1.51, 1.00 and 1.48 log2 against planted 1.5, 1.0 and 1.5. Null genes average about 0.
- **Reference levels.** The design columns must be `diet[T.MD]`, `microplastic[T.PE]` and their interaction. Without explicit Categoricals, the formula engine orders levels **alphabetically**, making PE the reference; the first run hit exactly this.
- **Size factors on real data:** identical to Phase 1 (maximum relative difference 0). Normalisation only; no dispersion fit or tests on real genes.

**What it does not check:**

- **Statistical calibration**: whether null p-values are uniform.
- **FDR control** at the nominal level.
- **Power**: how many planted genes reach significance.
- Behaviour with the real dispersion structure (gene-specific, mean-dependent), with outliers or contamination, with Cook's-NA flagging, or with the missing length offset.

Average coefficient recovery on one simulated dataset says the parameterisation is right, not that the inference is right. Validating calibration would need a simulation study:

- many replicate datasets under a null and under realistic alternatives;
- dispersions drawn from our estimated distribution;
- an empirical FDR and power computed at each nominal level.

That is listed as open work. It is not needed to *choose* hypotheses, but it is needed before trusting gene-level calls at n = 4–5.

The first smoke-test draft used a synthetic set where 50% of genes changed; that biased median-of-ratios normalisation by about −0.4 log2 on null genes. This is expected behaviour, and a reminder of the "most genes unchanged" assumption.

---

## 7. Decisions the team needs to make

1. **Which direction(s) become formal hypotheses**, and at what level (gene, gene set, whole profile). This sets the multiple-testing burden and the power.
2. **Primary PE contrast:** within CD only, within each diet, or the main effect (most powerful, but assumes no interaction).
3. **Count-matrix route:** accept the unvalidated no-offset fallback with its limitation stated, or run the re-quantification pilot first (§6). Should only the pilot be run, or all 19 samples?
4. **MD_1 / muscle contamination:** pre-specify one primary approach (keep with covariate, drop MD_1, exclude muscle-tracking genes); the others become sensitivity analyses.
5. **Library depth and quality:** state as a limitation, adjust with a covariate (losing ~half the precision), or use only as a sensitivity analysis.
6. **Software:** R/Bioconductor DESeq2 (exact reproduction, tximport offset, apeglm, edgeR/limma, RRHO) or pydeseq2 (stays in the current Python repo).
7. **"Metabolic dysfunction" reference:** internal diet contrast (shared-reference bias) or an external MASH signature (which dataset?).
8. **Thresholds:** FDR level; the paper's |log2FC| ≥ 0.38 (comparability) or `lfcThreshold`-based testing (more defensible); equivalence margins if "no interaction" matters.
9. **Calibration study:** run a simulation of p-value calibration and FDR before trusting gene-level calls, or rely on published DESeq2 benchmarks.
10. **Gene-set collections and size limits**, chosen before seeing results.

## 8. Validation now in CI

- `tests/test_feasibility.py`:
  - design rank;
  - collinearity detection;
  - dispersion estimator recovers a known φ on simulated NB data;
  - minimum-detectable-effect formula against a hand calculation;
  - MDE ordering (main effect < simple effect < interaction);
  - leave-one-out splitting;
  - jackknife SE against the textbook value.
- `scripts/validate_methodology_outputs.py`:
  - committed tables exist with the qualified column names;
  - 19 unique runs matching GEO;
  - summed counts never exceed read pairs;
  - interval bounds ordered;
  - the expected MDE ordering holds.
- CI job `pydeseq2-smoke-test`: installs pydeseq2 0.5.4 and runs the synthetic smoke test with `--strict`. The job fails if the reference levels are wrong or the planted coefficients are not recovered within 0.15 log2.

## 9. Suggested implementation order (after decisions)

1. Decide the count-matrix route; run the re-quantification pilot if chosen.
2. Set up DESeq2 (R or pydeseq2) with explicit reference levels. Validate size factors and reproduce the paper's pairwise contrasts against tables S5–S7 once obtained.
3. Fit the pre-specified primary model. Run diagnostics: dispersion trend, Cook's-NA counts, p-value histograms, MA plots.
4. Direction B, then A, then the Direction C sensitivity table for every headline result.
5. Report effect sizes with uncertainty first and significance second, with every limitation listed here.
