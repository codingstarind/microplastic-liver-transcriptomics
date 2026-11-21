# Methods literature for modelling a small, confounded liver RNA-seq factorial and a few external microplastic studies

Compiled 2026-11-21. How each paper was checked:
- **FT**: full text read (Europe PMC / NCBI E-utilities XML, or the author-hosted PDF). The sections and figures cited were read directly.
- **ABS**: only the abstract and bibliographic record were read (the full text was paywalled, or the fetch was blocked). Any detail beyond the abstract is marked **UNVERIFIED**.

All DOIs, journals, volumes and pages were checked against Europe PMC records. Citation counts are Europe PMC counts on the date above and are only a rough guide.

**Our design, for reference.** 2x2 diet x PE, n = 4/5/5/5 (N = 19). Four cell means give 15 residual df. PE is perfectly confounded with lower sequencing depth, and cage = group. External studies are small (for example 4x4 = 16 samples) and differ in polymer, particle size, strain and design.

---

## Direction A: factorial and gene-set inference within a small RNA-seq design

### A1. CAMERA: competitive gene-set test that adjusts for inter-gene correlation (FT)
- **Citation:** Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Research* 2012; 40(17): e133. PMID 22638577, PMCID PMC3458527.
- **DOI:** https://doi.org/10.1093/nar/gks461
- **Problem:** Competitive gene-set tests ask whether a set is more DE than other genes. Tests that permute genes (geneSetTest, PAGE, sigPathway, and GSEA in its gene-permutation mode for small n) assume genes are independent. Positive correlation within a set inflates their type I error.
- **Method:** Fit the limma linear model and compute genewise moderated t, converted to z. From the d = n - p residual df, extract independent residuals for the genes in the set, standardise them, and estimate the variance inflation factor VIF = 1 + (m - 1) * rho_bar, where rho_bar is the mean inter-gene correlation. The set's mean z is then compared with the background genes by a two-sample t-test (parametric version) or a correlation-adjusted Wilcoxon test (rank version), with the standard error inflated by the VIF and referred to a t distribution on d df.
- **Assumptions:** A linear model on log-expression (for RNA-seq, voom or logCPM). The estimated VIF is treated as chi-square with d df. A single average correlation is assumed per set. The null is competitive (the set is no more DE than the rest).
- **Evaluation:**
  - Simulations with 10,000 genes, sets of m = 100, and two groups of four (d = 6) or a third group added (d = 27). Results are in "CAMERA controls type I error correctly" and Table 3 ("CAMERA retains good power").
  - With rho_bar = 0.05, CAMERA keeps nominal size where gene-permutation tests do not. Power at d = 6 is about 0.45 to 0.65 for small consistent log2FC (0.05 to 0.25).
  - Correlation estimates are more precise than the nominal chi-square approximation implies (Table 2: m = 40, d = 27), so CAMERA is slightly conservative at small n.
  - Real data: breast cancer subtypes, mammary epithelial cells, and mouse HSCs from genetically identical mice.
- **Relevant limitation:**
  - "Even well-controlled data sets show inter-gene correlations": correlation in curated sets was clearly positive even for sorted cells from inbred mice. That is the reason to use CAMERA rather than an ORA/hypergeometric or gene-permutation test in our study.
  - Precision rests on residual df. We have d = 15, between the paper's 6 and 27.
  - CAMERA cannot separate a PE effect from a depth effect. Both are the same column of the design. It only tests the confounded contrast more honestly.
  - A competitive test is somewhat insensitive to global shifts shared by all genes. Depth artefacts that are gene-specific, such as low-count or short or GC-biased genes, can still enrich particular sets. (This last point is our inference, not a claim in the paper.)
- **Sections:** Results ("Variance inflation factors", "CAMERA", "CAMERA controls type I error correctly", "CAMERA retains good power", "Even well-controlled data sets show inter-gene correlations"); Table 2; Table 3; Discussion.

### A2. ROAST: self-contained rotation gene-set test for any linear-model design (FT)
- **Citation:** Wu D, Lim E, Vaillant F, Asselin-Labat ML, Visvader JE, Smyth GK. ROAST: rotation gene set tests for complex microarray experiments. *Bioinformatics* 2010; 26(17): 2176-2182. PMID 20610611, PMCID PMC2922896. (`fry` is a later fast closed-form approximation in limma. It is not covered by this paper.)
- **DOI:** https://doi.org/10.1093/bioinformatics/btq401
- **Problem:** Self-contained gene-set tests (H0: no gene in the set is DE for a given contrast). Gene permutation destroys correlation between genes. Sample permutation needs many replicates, works only for simple two-group comparisons, and cannot test general linear-model contrasts such as interactions.
- **Method:**
  - For each gene, project the data onto the (d + 1)-dimensional residual space orthogonal to the nuisance coefficients. This removes the other factorial effects.
  - Randomly rotate the residual vectors, using the same rotation for every gene so that gene-gene correlation is preserved. Recompute the moderated t and z for each gene, then the set statistic (mean, msq, mean50 or floor-mean), with up, down or mixed alternatives.
  - P is the proportion of rotated statistics at least as extreme as the observed one, usually from B = 10,000 rotations. Because rotation is exact under normality, the number of rotations does not depend on sample size.
- **Assumptions:** Multivariate normal residuals (shown to be robust to strong non-normality for the mean-type statistics, Section 5.2). The linear model must be correctly specified. Sample weights and correlations can be included.
- **Evaluation:**
  - Section 5.1 simulations included **two-way 2x2 factorial designs with 3 or 5 samples per cell (D2-3, D2-5)** and three-group designs with 3 or 5 replicates.
  - Size is held exactly (Section 5.2, Supplementary Tables S1, S2, S4).
  - Power was higher than GSEAlm (Section 5.4). GSEAlm failed to hold size with 3 per group because of granular permutation p-values.
  - Application: conservation of a mouse-defined mammary stem-cell set in human data (Section 6).
- **Relevant limitation:**
  - It is the closest match in the literature to our design: a 2x2 with 3 to 5 per cell, testing any contrast, including the interaction.
  - A self-contained test answers "is anything in this set DE for this contrast". With a contrast perfectly aliased with depth, any depth-driven change in the set will be declared significant. Self-contained tests are therefore more exposed to the PE/depth confound than competitive tests.
  - P-values are per set. Testing many sets requires FDR adjustment, and the Monte Carlo floor is 1/B.
- **Sections:** Sections 3.1, 3.3, 3.4, 4.2, 5.1-5.4, 6, 7.

### A3. Schurch et al. 2016: how many biological replicates for RNA-seq DE (FT)
- **Citation:** Schurch NJ, Schofield P, Gierlinski M, Cole C, Sherstnev A, Singh V, Wrobel N, Gharbi K, Simpson GG, Owen-Hughes T, Blaxter M, Barton GJ. How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? *RNA* 2016; 22(6): 839-851. PMID 27022035, PMCID PMC4878611.
- **DOI:** https://doi.org/10.1261/rna.053959.115
- **Problem:** An empirical gold standard for power, FDR control and tool choice as a function of the number of replicates.
- **Method:** 48 vs 48 replicates of *S. cerevisiae* (WT vs Delta-snf2), of which 42 and 44 were "clean". Bootstrap subsamples of n_r replicates are compared, for 11 DE tools, against each tool's own full-data calls. FPR is checked with WT-vs-WT null comparisons.
- **Assumptions:** The full-replicate call set is treated as truth. The material is a clonal yeast strain under lab conditions, which has low biological variability.
- **Evaluation:**
  - With 3 replicates, 9 of 11 tools find 20-40% of the significant DE genes, rising to more than 85% for genes with more than fourfold change.
  - Finding more than 85% of all DE genes needs more than 20 replicates.
  - edgeR and DESeq2 control FDR at all n_r (Fig. 1 for edgeR exact; per-tool Supplementary Figs S2-S12).
  - DEGSeq, NOISeq and SAMSeq fail null FPR tests (Fig. 4).
  - Recommendation: at least 6 replicates, and at least 12 when most DE genes matter. Use edgeR or DESeq2 below 12 ("Recommendations for RNA-seq experiment design").
- **Relevant limitation:**
  - With 4 to 5 per cell, we are near the n = 3 to 6 regime, where most moderate-fold-change DE genes are missed while FDR stays controlled.
  - A null result at the gene level is therefore largely uninformative. This argues for gene-set or effect-size-level inference rather than counting DEGs.
  - Mouse liver from different cages has more variability than clonal yeast, so the required n is probably higher for us. The yeast design also has no confounding or interaction terms.
- **Sections:** Abstract; Results ("Tool performance", Fig. 1; "Testing tool false positive rates", Fig. 4); Discussion ("Recommendations").

---

## Direction B: cross-study synthesis, batch confounding, and mouse-to-human benchmarking

### B1. Ramasamy et al. 2008: key issues in meta-analysis of expression datasets (FT; a methods guide, not a new method)
- **Citation:** Ramasamy A, Mondry A, Holmes CC, Altman DG. Key issues in conducting a meta-analysis of gene expression microarray datasets. *PLoS Medicine* 2008; 5(9): e184. PMID 18767902, PMCID PMC2528050.
- **DOI:** https://doi.org/10.1371/journal.pmed.0050184
- **Problem:** A practical checklist for combining heterogeneous expression studies, from study identification and data extraction through cross-platform probe-to-gene mapping, annotation and choice of combination method.
- **Method:** Reviews four combination families: vote counting; combining ranks (RankProd, METRADISC); combining p-values (Fisher's sum of logs, chi-square on 2k df); and combining effect sizes (Hedges' g with small-sample correction, inverse-variance weighting, fixed or random effects). Illustrated with an inverse-variance random-effects meta-analysis of 21 cancer studies (Tables 3-5, forest plot Fig. 4).
- **Assumptions:** Studies measure the same contrast. Per-study effect sizes and variances can be estimated. Gene identifiers can be harmonised.
- **Evaluation:** A worked example only. There is no benchmark. Table 6 lists advantages and disadvantages of each choice.
- **Relevant limitation:**
  - Vote counting is the weakest option and ignores effect size and direction strength.
  - Fisher's method can be driven by one study and does not give a pooled effect.
  - Hedges' g is biased upward at small n, hence the correction.
  - Most relevant to us: the guide presumes each study estimates the *same* biological contrast. Different polymers, sizes, doses, strains and diets break this. A random-effects tau^2 is the only formal accommodation, and it is poorly estimated with k of about 3 to 5 studies.
- **Sections:** Issue 6, "Choosing a Meta-Analysis Technique" (vote counting, combining ranks, p-values, effect sizes); Table 1 checklist; Table 6; Discussion.

### B2. Sweeney et al. 2017: empirical benchmark of random-effects meta-analysis for reproducibility (FT)
- **Citation:** Sweeney TE, Haynes WA, Vallania F, Ioannidis JP, Khatri P. Methods to increase reproducibility in differential gene expression via meta-analysis. *Nucleic Acids Research* 2017; 45(1): e1. PMID 27634930, PMCID PMC5224496.
- **DOI:** https://doi.org/10.1093/nar/gkw797
- **Problem:** Which meta-analysis choices (random-effects estimator, significance or effect-size thresholds, heterogeneity filters, number of datasets K vs total sample size N) give reproducible DE genes?
- **Method:**
  - Three very large clinical meta-analyses: cardiomyopathy, kidney transplant rejection and lung adenocarcinoma, each with more than 1,000 samples and at least 13 datasets.
  - The "silver standard" is the set of genes significant at q < 0.01 under the Sidik-Jonkman (SJ), DerSimonian-Laird (DL) and Hunter-Schmidt (HS) estimators on all data.
  - Every subset of 5 datasets, and every subset with at most N/2 samples and at most K/2 datasets, is meta-analysed with Hedges' g and random effects, and its calls are scored against the silver standard.
- **Assumptions:** The full-cohort result is the best available truth (relative accuracy). Datasets with batch-confounded designs were **excluded** at curation ("excluding studies where, for instance, healthy controls came from one batch and diseased samples were processed separately").
- **Evaluation:**
  - Effect-size thresholds (about 1.2 to 1.3-fold) cut false positives more than true positives. Filtering on heterogeneity removed more true positives than false positives (Fig. 3, "Significance thresholds", "Effect sizes and significance").
  - BH q-values substantially underestimated the actual FPR. For example, cardiomyopathy SJ had an FPR of 33% at q < 0.01.
  - For fixed total N, accuracy rose monotonically with K (Fig. 4; Table 3 regressions).
- **Relevant limitation:**
  - Every K here comes from clinical cohorts with about 1,500 samples in total. Nothing was tested at k = 3 to 5 with about 16 to 20 samples each.
  - The authors removed confounded studies, which is exactly our PE/depth situation, so their framework gives no guidance for a confounded study except to exclude it or down-weight it.
  - The "more datasets beats more samples" result does support pooling several small external studies over relying on any single one.
- **Sections:** Methods ("Construction of meta-analysis cohorts", "Silver standard true positives"); Results ("Significance thresholds", "Effect sizes and significance", "Number of datasets and size of datasets", "Accuracy as a function of sample size and K"); Figs 1, 3, 4; Tables 1, 3; Discussion.

### B2b. Haynes et al. 2017: MetaIntegrator multi-cohort framework (FT; software and framework paper, supporting)
- **Citation:** Haynes WA, Vallania F, Liu C, Bongen E, Tomczak A, Andres-Terre M, Lofgren S, Tam A, Deisseroth CA, Li MD, Sweeney TE, Khatri P. Empowering multi-cohort gene expression analysis to increase reproducibility. *Pacific Symposium on Biocomputing* 2017; 22: 144-153. PMCID PMC5167529.
- **DOI:** https://doi.org/10.1142/9789813207813_0015
- **Method as written (Section 2.2):**
  - Per-dataset Hedges' g with correction J.
  - Pooled g_s by DerSimonian-Laird random effects, with weights W_i = 1/(V_i + tau^2); BH FDR across genes.
  - Cochran's Q for heterogeneity.
  - Fisher's combined p separately for up and down.
  - Signature selection by thresholds on effect size, FDR, heterogeneity and the number of studies in which a gene is present, plus leave-one-dataset-out robustness (Section 2.3).
  - Sample score = geometric mean of up genes minus geometric mean of down genes, z-scored within each study (Section 2.4).
- **Relevant limitation:**
  - There is no new benchmark here; that is in B2. The Discussion notes that rank-based methods such as RankProd "fail to produce a summary effect size".
  - At the time it was microarray-only (RNA-seq counts were planned).
  - Leave-one-dataset-out is a cheap robustness check that fits our setting with k of about 4 to 5.

### B2c. RankProd (ABS; supporting)
- **Citation:** Hong F, Breitling R, McEntee CW, Wittner BS, Nemhauser JL, Chory J. RankProd: a Bioconductor package for detecting differentially expressed genes in meta-analysis. *Bioinformatics* 2006; 22(22): 2825-2827. PMID 16982708.
- **DOI:** https://doi.org/10.1093/bioinformatics/btl476
- **Verified from the abstract:** Extends the Breitling 2004 rank product to integrate studies from different labs and platforms for two-condition DE. Significance comes from a non-parametric permutation test with p and FDR. It accepts any pre-processed data and needs no cross-platform normalisation.
- **Relevance:** Rank-based pooling is robust to scale differences between labs and to outliers. It gives no pooled effect size or heterogeneity estimate (see Haynes 2017 Discussion), and with 4 per group the within-study rank product is coarse. Performance details beyond the abstract are **UNVERIFIED**.

### B3. Leek et al. 2010: batch effects and confounding (FT; a Nature Reviews Genetics *Perspective*, not a primary methods paper)
- **Citation:** Leek JT, Scharpf RB, Bravo HC, Simcha D, Langmead B, Johnson WE, Geman D, Baggerly K, Irizarry RA. Tackling the widespread and critical impact of batch effects in high-throughput data. *Nature Reviews Genetics* 2010; 11(10): 733-739. PMID 20838408, PMCID PMC3880143. Read from the author-hosted PDF.
- **DOI:** https://doi.org/10.1038/nrg2825
- **Problem:** Technical variation (processing date, lab, reagent lot, personnel) is widespread and survives normalisation. It is most damaging when correlated with the outcome.
- **Content and analyses:**
  - Fig. 1: a bladder cancer study in which CIS status was confounded with processing date. Normalisation removed global differences (Fig. 1a-b), but gene-subset batch effects remained (Fig. 1c), and controls clustered perfectly by date (Fig. 1d).
  - Table 1: in 8 more public datasets across technologies, batch surrogates (date or processing group) correlated with a top PC more strongly than the biological outcome did.
  - Fig. 2: in 1000 Genomes sequencing data, 32% of features were associated with date, up to 73% with PC2, and only 17% with biology.
  - Fig. 3: batch changes gene-gene correlation structure, which affects rank-based and network methods.
  - Fig. 4: a workflow of exploratory analysis, then adjustment (known covariates, ComBat or SVA), then diagnostics.
- **Key statements:**
  - Adjustment approaches "are most effective when batch effects are not highly confounded, or correlated, with the biological variables of interest."
  - In a confounded experiment, "results due to batch will be impossible to distinguish from real biological effects."
  - Recommendations: distribute biological groups equally across processing groups and times, and report processing date and group.
- **Relevant limitation:** Directly applicable. PE is perfectly confounded with depth, and cage is aliased with group, so no statistical method can partition the PE main effect into biology vs depth or cage. Only design (balance) or external replication can.
- **Sections:** "Group and date are only surrogates"; "Examination of public data"; "Downstream consequences"; the paragraph on SVA and ComBat ("most effective when … not highly confounded"); Conclusion; Figs 1-4; Table 1.

### B4. ComBat-seq: negative-binomial batch adjustment for counts (FT)
- **Citation:** Zhang Y, Parmigiani G, Johnson WE. ComBat-seq: batch effect adjustment for RNA-seq count data. *NAR Genomics and Bioinformatics* 2020; 2(3): lqaa078. PMID 33015620, PMCID PMC7518324.
- **DOI:** https://doi.org/10.1093/nargab/lqaa078
- **Problem:** Gaussian ComBat applied to logCPM is a poor fit for over-dispersed counts. Downstream edgeR or DESeq2 need integer counts.
- **Method:**
  - Genewise NB GLM: log mu_gij = alpha_g + gamma_gi + X_j beta_g + log N_j, with a batch-specific dispersion phi_gi.
  - Batch mean and dispersion parameters are estimated edgeR-style (Fisher scoring; Cox-Reid adjusted profile likelihood).
  - Each count is mapped by quantile matching to a "batch-free" NB that keeps the biological covariate term X_j beta_g. The output is integers, and zeros stay zero.
- **Assumptions:** Batch is known. The biological covariates are supplied and are **estimable separately from batch**. NB distribution.
- **Evaluation:**
  - polyester simulations with two conditions x two batches in a **balanced** case-control design, varying mean and dispersion batch effects. Without batch effects, FPR 0.048 and TPR 0.96.
  - ComBat-seq gave the best power and controlled FPR when dispersion differed by batch. FPR was 0.059-0.067 when only the mean differed, where adding batch as a covariate suffices (Results, "Simulations", Fig. 3).
  - Real data: the GFRN signature dataset (Fig. 4, PCA and variance explained). The authors describe its design as "challenging": controls balanced across batches, each treatment nested in one batch.
- **Relevant limitation:** The Discussion states that the model "may not work well on data with severely or even completely confounded study designs" and that "batch correction in confounded designs is challenging for most if not all the state-of-the-art batch adjustment methods". Consequences for us:
  - (i) Within our study, ComBat-seq cannot remove a depth effect aliased with PE. Specifying batch = depth/PE and keeping PE as a covariate is non-estimable, and dropping PE would remove the PE signal.
  - (ii) Across studies, lab = study is nested with study-specific design. Correction works only if each study contains common reference groups (for example, each has its own control). Even then, the treatments differ, so any "batch-corrected pooled matrix" must keep within-study contrasts.
  - Vacca et al. (B6) used ComBat on quantile-normalised counts across facilities, relying on each model having its own matched controls.
- **Sections:** Methods ("ComBat-seq model", "ComBat-seq adjustment", "Simulations"); Results ("Simulations", Fig. 3; "Application to the GFRN signature dataset", Fig. 4); Discussion (the confounding caveat).

### B5. Teufel et al. 2016: mouse NAFLD models vs human liver transcriptomes (ABS)
- **Citation:** Teufel A, Itzel T, Erhart W, Brosch M, Wang XY, Kim YO, von Schonfels W, Herrmann A, Bruckner S, Stickel F, Dufour JF, Chavakis T, Hellerbrand C, Spang R, Maass T, Becker T, Schreiber S, Schafmayer C, Schuppan D, Hampe J. Comparison of gene expression patterns between mouse models of nonalcoholic fatty liver disease and liver tissues from patients. *Gastroenterology* 2016; 151(3): 513-525.e0. PMID 27318147. Not in PMC.
- **DOI:** https://doi.org/10.1053/j.gastro.2016.05.051
- **Verified from the abstract:**
  - Chow controls plus 9 mouse NAFLD models (HFD ± fructose, Western-type diet, MCD, HFD + streptozotocin, hepatocyte Pten knockout), compared with human liver from 25 NASH, 27 NAFLD, 15 healthy obese and 39 healthy non-obese people.
  - Within-species limma DE; functional overlap by GSEA.
  - Human-mouse differences were larger than differences between stages or models. Of 65 NASH and 177 NAFLD human DEGs, only 1-18 were also DE in any mouse model.
  - At the pathway level (lipid metabolism), overlap was better, and HFD was the closest model. Data are at nash-profiler.com.
- **Relevance:** Gene-level overlap of DEG lists between species and labs is very low, while pathway-level concordance is better. That supports comparing our diet signature to human MASLD at the gene-set or NES level rather than by Venn overlap. Section and figure numbers, ortholog mapping and exact statistics are **UNVERIFIED** (full text not accessed).

### B6. Vacca et al. 2024: unbiased ranking of murine MASLD models against human MASLD (FT)
- **Citation:** Vacca M, Kamzolas I, Harder LM, Oakley F, Trautwein C, Hatting M, Ross T, Bernardo B, Oldenburger A, Hjuler ST, Ksiazek I, Linden D, Schuppan D, Rodriguez-Cuenca S, et al. An unbiased ranking of murine dietary models based on their proximity to human metabolic dysfunction-associated steatotic liver disease (MASLD). *Nature Metabolism* 2024; 6(6): 1178-1196. PMID 38867022, PMCID PMC11199145.
- **DOI:** https://doi.org/10.1038/s42255-024-01043-6
- **Problem:** Which rodent models best resemble human MASLD and MASH, metabolically, histologically and transcriptomically?
- **Method:**
  - Data: 598 animals (509 mice, 89 rats) in 39 model designs, each with matched controls (Fig. 1).
  - Human reference: UCAM-VCU (n = 136; E-MTAB-9815, GSE130970) and EPoS (n = 168; GSE135251). DESeq2 with design ~batch + condition, contrasts mild vs control and moderate-severe vs mild fibrosis.
  - Human genes and pathways were grouped as "early disease development", "all disease stages" and "disease progression". KEGG pathways were scored by fgsea NES.
  - Rodent data came from different facilities. They were quantile-normalised and then ComBat-corrected (sva) to merge. DE was always model vs its own control.
  - Agreement was defined as the % of human DEGs or DRPs significant and in the same direction in each model, tested with a one-sided hypergeometric test (Fig. 4). Pathways were compared as an NES heatmap (Fig. 3), and key genes as a log2FC heatmap (Fig. 5).
  - A composite "MHPS" score was built from phenotype, histology and transcriptome layers (Fig. 6), with PLSR on design factors (Fig. 7).
- **Assumptions:** Models share comparable matched controls. Batch correction across facilities works because every batch contains both controls and treated animals. Human-mouse gene mapping by HGNC symbols via biomaRt.
- **Evaluation and results:**
  - Whole-transcriptome direction agreement was only "approximately 50%" (chance level).
  - Restricting to human-significant DEGs and DRPs gave improved, significant agreement for most models.
  - Pathways were more concordant than single transcripts.
  - Western diets ranked closest. Choline-deficient models induced fibrosis fast but translated poorly.
- **Relevant limitation:**
  - The template is directly usable: compare contrast against matched controls, not pooled raw expression. Score agreement only on genes or pathways with significant human signal, use direction-aware metrics, and test overlap against a defined background.
  - Each model typically had about 5 to 15 animals, and the paper notes histologic fibrosis was not predicted well by fibro-inflammatory transcription.
  - The reference is human MASLD, not microplastics. No microplastic reference signature of comparable quality exists.
- **Sections:** Results ("Study design", "Murine MASLD transcriptomes are close to humans but do not predict fibrosis efficiently"); Figs 1, 3-6; Methods ("Human next-generation sequencing datasets", "Processing RNA-seq and microarray data", "Statistical analysis and reproducibility", "Enrichment analysis", "MHPS").

---

## Direction C: predictive ML on small cohorts and cross-study validation

### C1. Bernau et al. 2014: cross-study validation (CSV) (FT)
- **Citation:** Bernau C, Riester M, Boulesteix AL, Parmigiani G, Huttenhower C, Waldron L, Trippa L. Cross-study validation for the assessment of prediction algorithms. *Bioinformatics* 2014; 30(12): i105-i112. PMID 24931973, PMCID PMC4058929.
- **DOI:** https://doi.org/10.1093/bioinformatics/btu279
- **Problem:** Within-study CV does not measure how a learner performs on independent studies from other settings.
- **Method:**
  - For each algorithm k, build an I x I "CSV matrix" Z^k whose entry (i, j) is the performance (here the C-index) of a model trained on study i and validated on study j. The diagonal holds 4-fold CV estimates (Section 2.3).
  - Summarise the off-diagonal entries by mean or median (Section 2.4).
  - Simulation by two-level bootstrap: resample studies, then patients (Section 2.7).
- **Assumptions:** Several independent studies measure the same outcome with the same features. The super-population view is that studies are draws from a population of studies.
- **Evaluation:**
  - 8 ER+ breast cancer microarray cohorts (Table 1: n = 99 to 512 per study) predicting distant metastasis-free survival with 6 algorithms (Ridge, Lasso, CoxBoost, SuperPC, UniCox, plusminus).
  - Simulations: mean CV C-index about 0.65 vs CSV about 0.55 (Fig. 2; Section 3.1).
  - Real data: CV is about 0.06 higher than CSV (Section 3.2; Fig. 4).
  - Algorithm rankings differ between CV and CSV. Ridge, Lasso and CoxBoost act as "specialists" that adapt to the training study (Section 3.4; Table 2; Fig. 3).
- **Relevant limitation:**
  - Even with 100 to 500 samples per study and harmonised platforms, within-study CV was optimistic relative to transfer.
  - With our studies (about 16 to 19 samples, different polymer, size, strain and diet), a CSV matrix would be 4x4 or 5x5, with each cell estimated from 16 to 20 test samples. The labels are not even the same task (PE vs PS exposure; with or without diet). CSV is the right evaluation in principle, but its estimates would be extremely noisy here.
- **Sections:** Sections 2.3, 2.4, 2.7, 3.1, 3.2, 3.4, 4; Table 1; Figs 1-4.

### C2. Patil & Parmigiani 2018: training replicable predictors across multiple studies (ABS)
- **Citation:** Patil P, Parmigiani G. Training replicable predictors in multiple studies. *PNAS* 2018; 115(11): 2578-2583. PMID 29531060, PMCID PMC5856504.
- **DOI:** https://doi.org/10.1073/pnas.1708283115
- **Verified from the abstract:** Proposes ensembles of prediction models, each trained on a different study, and combines them with weighting criteria designed for replicability across studies. Quantifies how training on a single study contributes to failures of replication.
- **UNVERIFIED (full text blocked, so recalled and not checked):**
  - The data are curatedOvarianData / breast cancer compendia.
  - Weights come from stacking (regressing outcomes on study-specific predictions) and CSV-based schemes, compared with merging all studies.
  - A later paper by Zhang, Johnson and Parmigiani is described in search snippets as concluding that merging beats ensembling when studies are homogeneous and ensembling wins under heterogeneity. That is a secondary description, also unverified.
- **Relevance:** Conceptually the right way to build a "microplastic exposure" predictor across labs (study-specific learners plus a replicability-weighted ensemble, evaluated by leave-one-study-out). It needs several studies, each large enough to train a model. 16 to 19 samples per study with different exposures is far below that.

### C3. Ambroise & McLachlan 2002: selection bias in gene extraction (ABS)
- **Citation:** Ambroise C, McLachlan GJ. Selection bias in gene extraction on the basis of microarray gene-expression data. *PNAS* 2002; 99(10): 6562-6566. PMID 11983868, PMCID PMC124442.
- **DOI:** https://doi.org/10.1073/pnas.102102699
- **Verified from the abstract:**
  - Published near-zero error rates from a few genes came from CV or test error computed **without** repeating gene selection inside each training fold.
  - Correction: CV or bootstrap *external* to selection. Prefer 10-fold over LOOCV, and use the .632+ bootstrap for overfitted rules.
  - On two published datasets, corrected CV error is no longer zero for small gene subsets.
- **UNVERIFIED (full text not accessed):** the datasets are colon cancer (Alon, 62 samples) and leukaemia (Golub, 72 samples), with SVM-RFE and forward selection, and figures show internal-CV error near 0 vs external-CV error well above 0.
- **Relevance:** With 19 samples and about 15,000 genes, any feature filtering (DE-based, variance or "top genes") done before CV will make apparent accuracy near-perfect even with no signal. Every selection step must sit inside the resampling loop.

### C4. Varma & Simon 2006: bias of CV when tuning by CV; nested CV (FT)
- **Citation:** Varma S, Simon R. Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics* 2006; 7: 91. PMID 16504092, PMCID PMC1397873.
- **DOI:** https://doi.org/10.1186/1471-2105-7-91
- **Problem:** Is the minimum CV error over tuning parameters an unbiased estimate of true error?
- **Method:**
  - Simulated datasets of **40 samples (20 vs 20) x 6,000 features**, either null or with 10 DE genes at mean shift 1.
  - Tune shrunken centroids (Delta, 10-fold CV) or SVM (LOOCV). Compare the minimum CV error with true error on 20,000 independent samples.
  - Nested CV: an inner loop tunes, an outer loop estimates error.
- **Evaluation (Table 1; Figs 1-4):**
  - On null data, true error is 50%, but the mean CV error of the optimised classifier was 37.8% (shrunken centroids) and 41.7% (SVM).
  - CV error was below 30% in 18.5% of shrunken-centroid runs. The abstract reports 38% of SVM runs below 30%, while Table 1 shows 22.2%; the discrepancy is in the source.
  - Bias exceeded 20 percentage points in about 22 to 25% of runs.
  - Nested CV was nearly unbiased in both the null and non-null settings (Figs 3, 4).
- **Relevant limitation:** At n = 40, optimistic bias of 8 to 12 points on average, with wide variance, arose even with no signal. At n = 19 split across four cells, CV variance would be larger still. Nested CV removes the bias, not the variance, and with only about 5 samples per class per fold the error estimate is too coarse to be useful.
- **Sections:** Background; Results; Implementation; Conclusion; Table 1; Figs 1-4.

### C5 (optional). Toxicogenomics classifiers: the sample sizes they needed (ABS)
- **Fielden MR, Brennan R, Gollub J.** A gene expression biomarker provides early prediction and mechanistic assessment of hepatic tumor induction by nongenotoxic chemicals. *Toxicological Sciences* 2007; 99(1): 90-100. PMID 17557906. DOI https://doi.org/10.1093/toxsci/kfm156.
  - Verified from the abstract: rat liver expression after 5 days of dosing with **100** structurally and mechanistically diverse nongenotoxic hepatocarcinogens and non-hepatocarcinogens (DrugMatrix-type data) trained a multigene signature.
  - Validation on **47 independent test chemicals** gave sensitivity 86% and specificity 81%, beating liver weight, histopathology, ALT and Cyp induction.
  - Number of animals per chemical, classifier type and signature length: **UNVERIFIED**.
- **Natsoulis G, et al.** Classification of a large microarray data set: algorithm comparison and analysis of drug signatures. *Genome Research* 2005; 15(5): 724-736. PMCID PMC1088301. DOI https://doi.org/10.1101/gr.2807605.
  - Verified from the abstract: **597 microarrays** from a rat in vivo database of hundreds of drugs and toxicants. Linear SVMs and logistic regression with feature selection gave interpretable short "reward and penalty" gene signatures, and multiple disjoint signatures could reach the same endpoint.
- **Relevance:** Exposure-class prediction from liver transcriptomes that held up in independent validation needed tens to hundreds of *compounds* (hundreds of arrays) in a single-platform, single-lab compendium. That is 1 to 2 orders of magnitude beyond what pooled microplastic studies offer. The non-uniqueness of signatures (Natsoulis) warns against interpreting the "top genes" of a small-n classifier.

---

## Synthesis (feasibility with n of about 19 plus a few small external studies)

**A (within-study factorial and gene-set inference): feasible, and the most defensible primary analysis.**
- ROAST was explicitly validated on 2x2 factorials with 3 to 5 per cell and tests any contrast, including the interaction.
- CAMERA holds size with d = 6 to 27 and is needed because inter-gene correlation persists even in inbred mice. We have d = 15.
- Schurch shows gene-level DE at 4 to 5 per group misses most moderate changes while FDR stays valid. Report gene-set results and effect sizes, and treat gene-level non-significance as uninformative.
- Neither test can resolve PE vs depth (Leek; ComBat-seq Discussion). The PE main effect must be labelled "PE-or-depth". The diet effect, and the diet x PE interaction, are less exposed only if the depth artefact is additive.

**B (cross-study synthesis): feasible only as a descriptive, contrast-level concordance analysis, not as pooled inference.**
- Ramasamy and Sweeney/Haynes support Hedges' g random-effects synthesis, with leave-one-study-out checks and effect-size thresholds. Sweeney shows more studies helps, but its evidence base is more than 1,000 samples and it *excluded* confounded studies. With k of about 4 small, heterogeneous studies, tau^2 is not estimable with any precision.
- Vacca and Teufel show the workable template: contrast vs matched control within each study. Then compare direction agreement on reference-significant genes, plus pathway NES correlation, with hypergeometric or rotation tests against a defined background. Pathway-level concordance is far more stable than DEG overlap.
- Raw-count merging with ComBat-seq is invalid where batch is nested with treatment.

**C (predictive ML): not feasible as an inferential claim.**
- Ambroise and Varma show optimistic CV error even at n = 40 with no signal.
- Bernau shows within-study CV overstates transfer even with 100 to 500 samples per study.
- Toxicogenomic classifiers needed about 100 compounds and hundreds of arrays.
- At most, run a nested-CV / leave-one-study-out demonstration, framed as illustrating these biases (permutation-null comparison), not as a predictor.

---

## Verification status

| Paper | Status |
|---|---|
| Wu & Smyth 2012 (camera) | FT |
| Wu et al. 2010 (ROAST) | FT |
| Schurch et al. 2016 | FT |
| Ramasamy et al. 2008 | FT (guide, not original method) |
| Sweeney et al. 2017 | FT |
| Haynes et al. 2017 | FT (software and framework) |
| Hong et al. 2006 (RankProd) | ABS; performance details UNVERIFIED |
| Leek et al. 2010 | FT via author PDF (Perspective, not original research) |
| Zhang et al. 2020 (ComBat-seq) | FT |
| Teufel et al. 2016 | ABS; figure and section refs UNVERIFIED |
| Vacca et al. 2024 | FT |
| Bernau et al. 2014 | FT |
| Patil & Parmigiani 2018 | ABS; datasets, weighting schemes, figures UNVERIFIED |
| Ambroise & McLachlan 2002 | ABS; datasets and figures UNVERIFIED |
| Varma & Simon 2006 | FT |
| Fielden et al. 2007 | ABS; per-group n and classifier details UNVERIFIED |
| Natsoulis et al. 2005 | ABS |

Not located or not used: Hansen et al. 2020 *BMC Genomics* (a Europe PMC search did not find it, so it was not included), Sweeney et al. 2015/2016 (Sweeney 2017 NAR was used instead), and `fry` (no dedicated paper; it is a limma function documented in the limma user guide).
