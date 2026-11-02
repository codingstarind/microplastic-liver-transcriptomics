# Preliminary observations and candidate questions (Phase 1 → Phase 2)

**Status:** these are observations from exploratory analysis of 19 mice, plus candidate questions. They are **not** results, hypotheses or claims. No significance tests have been run. The team should decide which questions to pursue and fix hypotheses **before** running tests.

## What the exploration showed (facts about these data)

1. **Diet is the dominant designed source of variation.**
   - After contamination is handled, the first PC separates CD from MD (η² diet ≈ 0.80) and an unsupervised tree splits exactly by diet.
   - About 980 well-expressed genes have more than half of their variance tied to diet; with shuffled labels the figure is 3.
2. **PE has a smaller but visible footprint.**
   - 252 genes have more than half of their variance tied to PE, against 11 with shuffled labels.
   - CD+PE and CD+Veh separate cleanly on the design axis; MD+PE and MD+Veh separate less cleanly.
3. **On the main design axis, PE moves samples in the same direction as the MASH diet**, on both diets. The order is CD+Veh → CD+PE → MD+Veh → MD+PE.
4. **Gene-level patterns match several published PE effects.**
   - Lower lipogenic genes (*Fasn*, *Acly*, *Me1*, *Pnpla3*); higher *Plin4* and *Cyp4a14*.
   - In most of these, PE shifts expression the same way as the diet, and MD+PE is the most extreme group.
   - *Ppara* itself barely differs in bulk.
5. **Technical and design issues that could mimic or hide effects.**
   - PE status is confounded with library order and with smaller libraries (every PE library < every vehicle library).
   - Five samples carry skeletal-muscle tissue: MD_1 heavily; CD_MP_1, MD_2 and MD_3 moderately; CD_2 lightly.
   - Some samples have higher macrophage, cholangiocyte or stellate marker reads.
   - The clock gene *Dbp* differs between groups, which suggests different circadian or feeding states at sacrifice.

## Candidate questions

### Q1. Do PE-associated expression changes resemble the changes caused by the MASH diet?

- **Why it is interesting:** observation 3 suggests PE pushes the liver along the same transcriptional axis as a steatohepatitis-inducing diet. If true, PE could act like a mild "metabolic stressor". If false, the apparent similarity could come from shared technical structure or from a few genes.
- **Data:**
  - PE effect within CD: CD+PE vs CD+Veh.
  - Diet effect within vehicle: MD+Veh vs CD+Veh.
  - These two contrasts share the CD+Veh reference. That creates correlated noise, so a naive correlation between them is biased upwards.
- **Possible approach:**
  - Estimate per-gene log fold changes with a negative-binomial GLM (DESeq2/edgeR), with shrinkage (apeglm/ashr).
  - Compare the two fold-change vectors: rank correlation, concordance of direction among genes affected by both, and gene-set-level comparison (GSEA on each contrast, then compare normalised enrichment scores).
  - To avoid the shared-reference bias: split the CD+Veh mice between the two contrasts, or use a permutation null that re-labels samples.
  - Projection idea: define the "diet axis" from MD vs CD within vehicle, then score PE samples on it.
- **Limitations:**
  - n = 4 for CD+Veh.
  - Shared reference group.
  - PE confounded with library order.
  - Circadian differences may produce diet-like lipogenic changes.
  - Similarity of direction does not show a shared mechanism.
- **Already done in the paper:** pairwise DEG lists, a Venn diagram of PE DEGs under CD and MD, and GSEA of a lipid-toxicity NR gene set. **Not done:** a direct quantitative comparison of PE-induced vs diet-induced profiles.

### Q2. Does an already metabolically stressed liver respond differently to PE? (diet × PE interaction)

- **Why it is interesting:** the paper suggests MD+PE is "compounded", but it used pairwise contrasts, which cannot distinguish additive from synergistic effects. Our descriptive look shows mostly additive-looking shifts (*Fasn*, *Plin4*). A few genes look diet-dependent: *Cd36* rises with PE under CD (+0.7) but not under MD (+0.1); *Cyp4a14* rises +2.6 under CD vs +1.7 under MD.
- **Data:** all four groups.
- **Possible approach:**
  - **Two-way factorial model:** `~ diet + PE + diet:PE` with a negative-binomial GLM per gene. The interaction term is the formal test of "responds differently". FDR control across genes.
  - **Gene-set level:** test the interaction for pathways (camera/roast or fgsea on interaction statistics).
  - Counts of DEGs per condition are a weak substitute (power-dependent) and should not be used to claim an interaction.
- **Limitations:**
  - Interaction tests need much more power than main effects; with n = 4–5 per cell, only large interactions are detectable. A null result is not evidence of no interaction.
  - Scale matters: an effect additive on the log scale can look interactive on the raw scale.
  - Contamination (MD_1, MD_2, MD_3 are all MD+Veh) could create spurious interactions unless handled.
- **Already done in the paper:** pairwise comparisons only. A formal interaction model would be new.

### Q3. Are there coordinated gene-expression programs (co-expression modules)?

- **Why it is interesting:** individual-gene tests ignore that genes move together. Modules could summarise the lipogenic, cholesterol-synthesis, inflammation, PPARα-target, mitochondrial and contamination programs, and show which programs respond to diet, PE or both.
- **Data:** the log2 normalised (or VST) matrix of the filtered, well-expressed genes.
- **Possible approaches:**
  - Hierarchical clustering of **genes** on correlation (the course method);
  - WGCNA-style modules, or k-means on gene profiles;
  - PCA loadings / sparse PCA;
  - non-negative matrix factorisation.

  Then relate module eigengenes (first PC of each module) to diet, PE, muscle score and library size using a two-way ANOVA on the eigengene.
- **Limitations:**
  - 19 samples is small for co-expression (WGCNA recommends ≥ 15–20; correlations are noisy).
  - Modules will partly reflect contamination and cell composition, not regulation.
  - The modules must be defined *before* relating them to treatment, or the analysis becomes circular.
- **Already done in the paper:** WGCNA on **spatial pseudo-bulk central hepatocytes** (Xenium panel genes only), focused on *Ppara*. **Not done:** whole-transcriptome bulk co-expression.

### Q4. Additional questions the data raise

| Question | Data / approach | Limitations |
|---|---|---|
| **a. Are PE effects robust to the technical confounds?** | Refit models with library size, muscle score and chrM fraction as covariates; leave-one-out (drop each mouse) stability of top genes | Covariates correlated with PE can absorb real effects. Over-adjustment trades bias for variance |
| **b. Do circadian or feeding genes differ between groups?** | Clock-gene panel (*Dbp*, *Per1–3*, *Nr1d1/2*, *Bmal1*, *Cry1*); compare groups descriptively, then as a covariate | If clock genes differ, lipogenic "effects" may be timing artefacts. If the clock is genuinely disrupted by PE, that is itself interesting |
| **c. Does cell-type composition differ (inflammation, Kupffer cells, ductular reaction)?** | Marker scores (started in fig09) or deconvolution with a public mouse-liver single-cell reference (e.g. MuSiC, CIBERSORTx-like NNLS) | Reference mismatch (diet, strain); bulk deconvolution is approximate |
| **d. Are PE-responsive genes enriched for xenobiotic, oxidative-stress or mitochondrial pathways?** | GSEA / over-representation on PE contrasts (Hallmark, GO, KEGG) | Results are only as good as the DE input; mitochondrial signals are partly technical here |
| **e. Is the PE signal consistent with polystyrene MPs?** | Compare with GSE245069 (PS liver), as the paper did for NRs, but genome-wide | Different polymer, size, dose, lab, platform: cross-study batch effects |
| **f. How much of the bulk "PE effect" is carried by a few samples?** | Influence and outlier diagnostics (Cook's distance), robust regression | Small n makes every sample influential |

## Recommended order for Phase 2 (for team discussion)

1. **Decide** how to handle contamination (MD_1) and whether to switch to VST in R/DESeq2.
2. **Q2 (factorial model)** answers both "PE effect within each diet" and "interaction" in one coherent model. It is the natural core analysis and builds on the course's hypothesis-testing material.
3. **Q1** follows from the Q2 model's fold changes.
4. **Q3** for the "coordinated patterns" aim, using the clustering material from the course.
5. **Q4a–b** as robustness checks that make any conclusion credible.
