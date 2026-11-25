# Polyethylene microplastic exposure and hepatic transcriptional programs under control and MASH-inducing diets: a re-analysis of GEO GSE294180

**BT3041: Analysis and Interpretation of Biological Data, IIT Madras**
**Authors:** [TEAM: names and roll numbers]
**Draft:** 2026-11-02. Final report due 2026-11-25.

> **Draft status.** Computational results were generated with Claude Code and have **not yet been independently checked by the team**. Items marked **[TEAM]** need a team decision or text. All analyses are **exploratory**: the team has not yet formally adopted hypotheses (Section 1.3). The count data are the deposited GEO matrix, used under assumptions A1–A3 (Section 3.4).

---

## Abstract

**Background.** Polyethylene (PE) is among the most abundant microplastics found in human tissues. A recent study (Jung *et al.*, *Sci. Adv.* 2026) reported that oral PE causes hepatotoxicity in mice and worsens diet-induced steatohepatitis, based partly on bulk liver RNA-seq (GEO GSE294180).

**Aim.** To re-analyse the bulk data with an explicit 2 × 2 factorial framework and ask three things:

- (A) whether PE-associated expression changes resemble those of a MASH-inducing diet;
- (B) whether the PE response depends on diet;
- (C) how robust these observations are to the dataset's technical limitations.

**Data and methods.**

- 19 male C57BL/6 mice: control diet (CD) or MASH diet (MD), × daily gavage with vehicle or 2 mg PE, for 8 weeks; n = 4/5/5/5.
- We audited the deposited count matrix and verified its quantification with a 5-sample Salmon re-quantification pilot.
- We fitted a negative-binomial factorial model (DESeq2), ran mouse-level gene-set tests (limma fry/camera; MSigDB 2026.1 mouse Hallmark/PPAR sets), and used a split-reference concordance analysis that removes shared-control bias.
- Robustness: leave-one-mouse-out refits, contamination and technical covariates, and alternative methods.

**Results.**

- **Diet** had a strong, robust effect: 2,262 genes at FDR < 0.05.
- **PE under CD** was associated with lower lipogenic and cholesterol-synthesis gene expression and higher PPARα-target expression. This reproduces the published gene-level directions.
- **Cholesterol homeostasis** was lower with PE under **both** diets, with no detectable interaction.
- **Gene-level interactions** were largely artefacts: 23 of 87 were skeletal-muscle carry-over. Only 1 interaction gene survived in edgeR or limma-voom, and no pathway-level interaction was detected.
- **PE vs diet:** PE changes under CD correlated modestly with diet changes (split-reference Spearman ρ ≈ 0.33), but not under MD (ρ ≈ 0.05).

**Two confounds limit causal interpretation:**

1. Every PE library was sequenced more shallowly and had a lower gene-count yield than every vehicle library. Adjusting for this removed most of the PE signal.
2. The CD+Veh reference mice had a distinct circadian clock-gene state.

**Conclusions.** A reproducible PE-associated hepatic signature exists in these data, centred on cholesterol and lipid synthesis programs and similar under both diets. Its attribution to PE is limited by library and circadian confounding. "PE mimics metabolic dysfunction" and "PE synergises with MASH" are not supported by the bulk transcriptome.

---

## 1. Introduction

### 1.1 Background

- **Microplastics** (particles 1 µm – 5 mm) are found in human blood, placenta and liver. PE is reported to be the dominant polymer accumulating in human liver.
- **The liver** integrates lipid metabolism and xenobiotic detoxification. **Metabolic dysfunction-associated steatohepatitis (MASH)** is a common diet-driven liver disease.
- Jung *et al.* (2026) exposed mice to PE under a control or a MASH-inducing diet and reported:
  - higher ALT and steatosis with PE;
  - PE-responsive genes in bulk RNA-seq;
  - from spatial transcriptomics (10x Xenium) plus ChIP and knockdown experiments, a PPARα–Anxa2 axis in pericentral hepatocytes.

### 1.2 Why re-analyse

The original bulk analysis used pairwise comparisons and lists of differentially expressed genes (DEGs). Such lists cannot:

- separate a PE effect that **depends on diet** (an interaction) from one that simply adds to the diet effect;
- quantify whether PE "resembles" the disease diet without bias from shared controls.

The original paper also reports no sample-level quality control for the bulk data.

### 1.3 Research directions

The team's working directions are listed below. **[TEAM: replace this section with the formally adopted hypotheses.]**

- **A. Similarity.** Do PE-associated expression changes resemble, oppose or differ from MASH-diet-associated changes, across the transcriptome or only in specific programs?
- **B. Diet dependence.** Is the PE response under the MASH diet different from that under the control diet (diet × PE interaction)?
- **C. Robustness.** Do observations from A and B survive the dataset's technical limitations?

Candidate testable hypotheses suggested by the exploratory work (not yet adopted):

- *H1:* PE lowers expression of the hepatic cholesterol-homeostasis program under both diets.
- *H2:* the PE effect on PPARα target genes differs between diets.

---

## 2. Data

| Item | Detail |
|---|---|
| Source | GEO GSE294180 (BioProject PRJNA1248512); SRA SRR33050945–SRR33050963 |
| Paper | Jung W. *et al.*, *Sci. Adv.* 12(25): eaec8681 (2026), doi:10.1126/sciadv.aec8681 |
| Animals | Male C57BL/6 mice, 8 weeks old at start; group-housed by treatment |
| Diets (8 weeks) | CD: Research Diets D09100304 (10 kcal% fat). MD: D09100310 (40 kcal% fat, 20 kcal% fructose, 2% cholesterol). |
| Exposure | Daily oral gavage: 2 mg PE microspheres (10–150 µm) in 0.1% Tween-20, or vehicle |
| Groups | CD+Veh n = 4; CD+PE n = 5; MD+Veh n = 5; MD+PE n = 5 |
| Assay | Poly(A) stranded bulk RNA-seq, NextSeq 2000, 2 × 150 bp; fastp 0.23.4; Salmon 1.10.3 (GENCODE M35); tximport |
| Deposited data | Integer gene-level count matrix: 55,836 gene symbols × 19 samples |

**Single-cell vs spatial vs bulk.** Only the bulk data are analysed here. The paper's Xenium spatial data, which resolve cell types and liver zones, are on Zenodo and were not used. Bulk RNA-seq averages over all cells in a tissue piece, so a change confined to one cell type or zone is diluted (Section 5.3).

---

## 3. Methods

All code is in the project repository. Phase 1 (audit and exploratory analysis) is in `scripts/01–07`; Phase 2 is in `analysis/` and `pipelines/`. Reproduction instructions are in Section 7. Analytic choices were recorded **before** results were seen (`reports/phase2/analysis_plan.md`), and every deviation is logged there.

### 3.1 Data audit (Phase 1)

- **Structure:** integer values, no missing or negative values, sample IDs matched to GEO records.
- **Sex:** all male (*Xist* ≈ 0, Y-chromosome genes expressed).
- **Gene identifiers:** GENCODE M35 gene symbols.
- **Low-expression filter:** ≥ 10 counts in ≥ 4 samples, keeping 13,777 genes.

### 3.2 Exploratory analysis (Phase 1)

- Library sizes and median-of-ratios normalisation.
- log2(normalised + 1) transformation, chosen after inspecting the mean–variance relationship.
- PCA on the 500 most variable genes.
- Hierarchical clustering (1 − Pearson distance, average linkage).
- Marker-gene screens for tissue contamination.

### 3.3 Sequencing and technical metadata (Phase 2 feasibility)

- Run metadata from ENA; the first FASTQ header of all 19 runs identifies the instrument, run and flowcell.
- A gene-counts-per-read-pair ratio, used as a **proxy** for library yield (not a mapping rate).

### 3.4 Quantification validation

**The problem.** The tximport vignette advises against using original gene-level Salmon counts without an average-transcript-length offset for full-length RNA-seq. The deposited matrix came without that offset, and the tximport setting was not reported.

**The pilot.**

- 5 samples (CD_1, CD_MP_3, CD_MP_5, MD_4, MD_MP_3), MD5-verified ENA FASTQs; CD_1 quantified twice.
- Re-quantified with the authors' reported settings: fastp → Salmon selective alignment, `--seqBias --useVBOpt`, GENCODE M35.
- Imported with all three tximport `countsFromAbundance` settings.
- Compared with the deposited matrix using **pre-registered metrics and decision rules** (`reports/quantification_validation/pilot_protocol.md`).

**Team decision (2026-11-02).** The analysis proceeds on the deposited matrix under three stated assumptions:

- **A1:** the values are estimated counts of the original (`"no"`) type, not TPM-scaled.
- **A2:** the missing length offset adds gene-level noise but no systematic PE or group bias.
- **A3:** differences between our Salmon run and the authors' do not affect protein-coding conclusions. Mitochondrially encoded genes and processed pseudogenes are flagged as less reliable.

### 3.5 Factorial differential-expression model

For each gene, a negative-binomial GLM (DESeq2 1.46.0, R 4.4.3):

```
log μ = log s + β0 + β_diet·[MD] + β_PE·[PE] + β_int·[MD]·[PE]
```

- Reference levels: CD and vehicle. These are set explicitly and checked in code, because a wrong reference silently flips every sign.
- Contrasts:

| Contrast | Coefficients | Meaning |
|---|---|---|
| PE_in_CD | β_PE | PE effect under control diet |
| PE_in_MD | β_PE + β_int | PE effect under MASH diet |
| Diet_in_Veh | β_diet | Diet effect without PE |
| Interaction | β_int | Difference of the two PE effects, on the log2 scale |

- Wald tests; Benjamini–Hochberg FDR within each contrast (one test family per contrast). 95% CIs = MLE ± 1.96 × SE. ashr-shrunken fold changes are used for ranking.
- Contrast signs were checked against raw group-mean differences (≥ 99.7% sign agreement for well-expressed genes).

### 3.6 Pathway analysis

- **Gene sets:** MSigDB 2026.1.Mm (mouse-native): 50 Hallmark sets (primary family) and 3 PPAR sets (separate family). Reactome is exploratory only.
- **Primary test:** limma **fry**, a self-contained rotation test on the VST matrix with the same factorial design. The unit of analysis is the mouse, and correlation between genes is respected.
- **Supporting tests:**
  - **camera** (competitive, correlation-adjusted);
  - **fgsea** (preranked on the Wald statistic; gene-level, anti-conservative);
  - per-mouse mean-z scores fitted with a two-way linear model.
- Each pathway's PE response was then classified across diets as similar, amplified, attenuated, opposing or one-diet-only, using pre-specified rules.

### 3.7 PE vs diet concordance (Direction A)

- **The bias.** PE_in_CD and Diet_in_Veh share the CD+Veh mice, which induces a spurious positive correlation (expected ≈ 0.56 of the noise component).
- **The fix: split reference.** The CD+Veh mice were split into disjoint pairs, one pair as the PE reference and the other as the diet reference, in all 6 assignments. The same was done for MD+Veh in the PE_in_MD comparison (20 splits).
- **Statistic:** Spearman correlation of gene log2 fold changes (9,941 genes).
- **Reference distribution:** exhaustive relabelling of the PE contrast's mice (21 relabellings per split).

### 3.8 Robustness (Direction C)

Every check was compared with the primary model:

- 19 leave-one-mouse-out refits;
- dropping MD_1 (heaviest contamination);
- covariates: skeletal-muscle marker score, sequencing depth, yield proxy, clock-phase score;
- edgeR quasi-likelihood and limma-voom with the same design;
- gene-class checks: muscle-tracking, mitochondrial, clock and highly abundant genes.

The primary model was pre-specified **without** covariates, and the team kept it as primary (2026-11-02).

---

## 4. Results

### 4.1 Data quality and exploratory structure

- **Muscle carry-over** (diaphragm/skeletal muscle) in some samples:

  | Sample | Muscle-marker CPM |
  |---|---|
  | MD_1 | 5,190 |
  | CD_MP_1 | 539 |
  | MD_2 | 316 |
  | MD_3 | 161 |
  | Others | ≤ 58 |

  Muscle genes formed PC1 (30% of variance) of the primary PCA (`figures/pca_clustering/fig10–fig11`).
- **After removing muscle-tracking genes,** PC1 (21%) separated diet, and hierarchical clustering at k = 2 recovered diet exactly (adjusted Rand index 1.0). PE was a sub-structure within diet (`fig13`, `fig15`).
- **Library depth and yield:** all 19 runs came from one flowcell. PE libraries had fewer read pairs (13.5–16.5 M vs 17.0–19.1 M; complete separation) and a lower gene-counts-per-read-pair yield (≈ 0.77 vs ≈ 0.80).

### 4.2 Quantification validation

| Check | Result |
|---|---|
| Gene set | Our pipeline drops exactly the same 1,061 GENCODE symbols as the deposited matrix |
| Run-to-run noise | Zero (CD_1 quantified twice: 99.5% of genes give identical integers) |
| `scaledTPM` | **Excluded.** Deposited counts scale with transcript length (slope ≈ 0.96) |
| `"no"` vs `lengthScaledTPM` | `"no"` fits best (median \|log2 ratio\| 0.029–0.037 vs 0.066–0.072) |
| Formal verdict | **Not identifiable**: our totals are 4.5–7.2% higher, concentrated in mitochondrial genes and pseudogenes, so the authors' Salmon run differed slightly from ours |
| Real mapping rates | 83.0–86.1%. PE samples are 0.6–3.0 points below the vehicle sample of the same diet. |
| Length offset | Effective-length differences between mice look like random per-mouse noise, not a PE effect. A group-level offset built from 1–2 pilot mice was **invalid** (it injected that noise as group effects). |

### 4.3 Differential expression

Primary model; FDR < 0.05.

| Contrast | Genes significant | … and \|log2FC\| ≥ 0.38 |
|---|---|---|
| PE under CD | 858 | 513 |
| PE under MD | 472 | 295 |
| Diet (vehicle) | 2,262 | 1,582 |
| Interaction | 87 | 84 |

**Diagnostics** (`figures/deseq2/provisional_deposited/`):

- well-behaved p-value histograms and dispersion trend;
- the most influential samples by Cook's distance were MD_1, MD_MP_2 and CD_MP_2, the samples flagged for contamination or composition differences;
- with 4–5 mice per group DESeq2 cannot replace outliers, so 55 genes received NA p-values.

**Reproduction of published genes** (directions; within-diet PE contrasts):

| Pattern | Genes |
|---|---|
| Lower with PE in **both** diets | Fasn, Acly, Me1, Pnpla3 (−0.8 to −3.1 log2) |
| Higher with PE in both diets | Plin4, mt-Nd6 |
| Lower with PE under CD | Angptl8, Acp3, Hmgcr |
| Higher with PE under CD | Cyp4a12a, Cyp26b1 |
| Significant only under MD | Tm7sf2, Ldlr, Mylip, Oit3 |
| **Not reproduced** as within-diet PE effects | Ppara (+0.16 under CD, FDR 0.62); Nr1i3 |

### 4.4 Pathways

| Program | PE under CD | PE under MD | Diet | Interaction |
|---|---|---|---|---|
| Cholesterol homeostasis | Down (fry p = 0.004, FDR 0.11) | Down (fry p = 0.005, FDR 0.09) | Down (FDR < 0.001) | None (fry p = 0.84; mean-z estimate −0.003, 95% CI −0.41 to 0.40) |
| PPARα targets (mouse liver) | **Up (fry FDR 0.041)** | Weaker (fry p = 0.35) | Up (FDR < 0.001) | camera FDR 0.003; fry p = 0.21 |
| Inflammation, xenobiotic, fibrosis/EMT | No effect | No effect | Mixed | None by fry |

- Across 53 Hallmark and PPAR sets, fry found **no pathway-level interaction** at FDR < 0.05.
- Competitive tests (camera, fgsea) suggested PPAR, peroxisome and bile-acid programs might be less PE-induced under MD. This is exploratory.
- Figure: `figures/phase2/provisional_deposited/fig16_pathway_pe_response_forest.png`.

### 4.5 Direction A: PE vs diet resemblance

| Comparison | Naive ρ | Split-reference ρ, median [range] | One-sided p |
|---|---|---|---|
| PE under CD vs diet | 0.54 | **0.33 [0.29–0.40]** | 0.048–0.095 (minimum attainable 0.048) |
| PE under MD vs diet | −0.08 | 0.05 [−0.03–0.17] | 0.10–0.57 |

- About 40% of the naive PE(CD)–diet correlation was shared-reference artefact.
- The remaining modest resemblance is stable when any mouse is left out (median ρ 0.25–0.36).
- At the pathway level it is concentrated in cholesterol (down) and PPARα (up) programs.
- Figure: `fig17_pe_vs_diet_concordance.png`.

### 4.6 Direction B: diet × PE interaction

| Fit | Interaction genes (FDR < 0.05) |
|---|---|
| Primary | 87, of which 23 (26%) muscle-tracking, against 1% of all genes |
| + muscle covariate | 16 (none muscle-tracking) |
| Without MD_1 | 49 |
| edgeR-QL | **1** |
| limma-voom | **1** |

- Of 640 genes with clear, same-direction PE effects in both diets, only 3 had interaction FDR < 0.05.
- The median interaction CI width was 0.80 log2, so differences up to about 1.3-fold cannot be excluded.
- Figure: `fig18_interaction_vs_muscle.png`.

### 4.7 Robustness

Full matrix: `reports/phase2/robustness_report.md`.

| Finding | Leave-one-out | Muscle | Clock | Depth / yield covariates | edgeR / limma |
|---|---|---|---|---|---|
| Diet effect | ✔ | ~ (1,483 genes) | ~ | ✔ | ✔ |
| PE under CD (genes) | ✔ | ✔ | ~ | ✖ (116 / 66 genes) | ~ (544 / 517) |
| PPARα up, PE under CD | ✔ 19/19 | ✔ | ✔ | ✖ | – |
| Cholesterol down, PE in both diets | ✔ 19/19 | ✔ | ✔ | ~ / ✖ | – |
| Gene-level interaction | ~ | ✖ | ~ | ✔ / ~ | ✖ |
| No pathway interaction | ✔ 18/19 | ✔ | ✔ | ✔ | – |
| PE(CD)–diet resemblance | ✔ | – | not interpretable | – | – |

✔ survives; ~ weakened; ✖ lost; – not applicable.

**Confounders that cannot be separated from treatment:**

- **Library depth and yield** (aligned with PE in both diets).
- **Clock state** (aligned with the CD+Veh group). All four CD+Veh mice scored below every other mouse on a pre-defined clock-gene phase score: group mean −2.38 vs +0.08 to +1.01 in the other groups.
- **Library-preparation order** and **cage** (both aligned with group).

---

## 5. Discussion

### 5.1 What the bulk data support

- A **reproducible PE-associated hepatic signature**: lower de novo lipogenesis and cholesterol-synthesis genes, and higher PPARα targets under control diet. It reproduces the paper's gene-level observations.
- The most coherent pathway-level change, **lower cholesterol homeostasis, is similar under both diets**. The bulk data give no evidence that metabolic stress amplifies the PE response, and no robust evidence of opposing or diet-specific responses.
- The modest PE–diet resemblance appears only under the control diet. It cannot be distinguished from the CD+Veh group's distinct clock state, so the hypothesis that **"PE mimics metabolic dysfunction" is not supported**.

### 5.2 Alternative explanations

1. **Library properties.** Adjusting for depth or yield removes most PE signal, but those variables are almost perfectly aligned with PE, so this neither confirms nor refutes a biological effect. PE-associated genes are over-represented among mitochondrial and very abundant transcripts, the classes most sensitive to library quality.
2. **Feeding and circadian state.** SREBP-driven cholesterol and lipogenic genes and PPARα targets are among the most fasting- and time-of-day-responsive programs in the liver. A fasting-like state, whether from PE affecting food intake or from different sacrifice timing, would produce the observed pattern.
3. **Tissue composition.** Skeletal-muscle carry-over created most gene-level "interactions". Kupffer, cholangiocyte and stellate marker differences between samples mean bulk effects may partly reflect cell composition.

### 5.3 Comparison with the original study

| | Status in our analysis |
|---|---|
| Lipid-gene directions, number of PE-responsive genes under CD, diet separation | **Reproduced**, which is reproducibility evidence (same mice and reads), not independent validation |
| "Compounded" MD+PE effect | **Not supported as an interaction**. The paper's comparisons against CD+Veh mix diet and PE effects; our factorial model separates them. |
| Ppara induction | **Not seen** as a within-diet PE effect in bulk |
| PPARα–Anxa2 axis, inflammatory hotspots | Rest on spatial and experimental data; dilution in bulk (Anxa2 shows no PE effect) is expected, not contradictory |

**Extensions beyond the paper:**

- the formal interaction model;
- the bias-corrected concordance analysis;
- systematic QC, confounding and robustness analysis;
- quantification verification.

**[TEAM: a literature search for other microplastic–liver transcriptomic studies is still needed before calling these novel.]**

### 5.4 Limitations

- n = 4–5 mice per group. Interaction power is low: about 0.5–0.7 log2 is needed for typical genes.
- PE is confounded with library properties; CD+Veh with clock state; group with cage and library-preparation order.
- The quantification was verified on 5 samples only (assumptions A1–A3). Full re-quantification was not done.
- Bulk RNA cannot resolve cell types or zones, and mRNA changes do not measure protein, flux or pathway activity.
- The phenotype data (ALT, triglycerides, histology) were not deposited per mouse, so expression cannot be linked to phenotype.
- Supplementary DEG tables S1–S18 have not yet been obtained for a direct list comparison.

### 5.5 Future work

- Full re-quantification (the pipeline and download manifest are ready).
- Sacrifice times and food-intake data from the authors.
- An independent MASH signature from an external dataset (candidate GEO series listed in the similarity report).
- Integration with the authors' Xenium data to test cell-type attribution.

---

## 6. Conclusions

1. Oral PE exposure is associated with a reproducible bulk hepatic signature: lower lipogenic and cholesterol-synthesis gene expression and higher PPARα-target expression. The cholesterol-homeostasis component is similar under control and MASH diets.
2. There is no robust evidence that the MASH diet amplifies or reverses the PE response. Apparent gene-level interactions were mostly skeletal-muscle contamination.
3. PE-associated changes resemble diet changes only modestly, and only under the control diet, where a circadian/feeding confound in the reference group offers an alternative explanation.
4. Because PE is perfectly aligned with library depth and yield, the data alone cannot establish that these changes are caused by PE.

**[TEAM: restate conclusions in terms of the adopted hypotheses.]**

---

## 7. Reproducibility

- **Repository:** https://github.com/codingstarind/microplastic-liver-transcriptomics. Phase 1 is on `main`; methodology is in PR #2; Phase 2 is in PR #3.
- **Environments:** Python 3.14 (`requirements.txt`); R 4.4.3 / Bioconductor and Salmon 1.10.3 / fastp 0.23.4 via micromamba (`envs/`).
- **Order of execution:** `scripts/01–07` (Phase 1); `scripts/08–09` (feasibility); `analysis/deseq2`, `analysis/pathways`, `analysis/concordance`, `analysis/interaction`, `analysis/robustness`, `analysis/figures` (Phase 2). Commands are in `reports/phase2/README.md`.
- **Provenance:** checksums and download dates are in `metadata/data_provenance.json` and `metadata/sra_run_info.csv`. Pilot FASTQs were MD5-verified against ENA.
- **Checks:** CI runs ruff, 33 unit tests, metadata and methodology validators, and a pydeseq2 synthetic smoke test.

## 8. How the course methods were used

| Course topic | Where it appears and why |
|---|---|
| Descriptive statistics and distributions | Library sizes, count distributions, zero inflation: deciding filters and transforms |
| Normalisation and transformation | Median-of-ratios size factors; log2/VST, chosen from the mean–variance relationship |
| PCA | Revealed the contamination axis before the biological axis |
| Hierarchical clustering | Recovered diet without supervision; showed PE as a within-diet structure |
| Hypothesis formulation and inference | Pre-registered rules for the quantification pilot; pre-specified contrasts |
| Multiple testing | BH-FDR per contrast and per gene-set family; why ~685 false positives would be expected at unadjusted α = 0.05 |
| Regression and factorial models | NB GLM with interaction; reference coding; contrasts as linear combinations of coefficients |
| Model assumptions and diagnostics | Dispersion fit, p-value histograms, Cook's distance, leave-one-out influence |
| Correlated errors and permutation | Shared-reference bias; split-reference design; relabelling null matched to the hypothesis |
| High-dimensional interpretation | Gene-set tests at the mouse level vs the gene level; enrichment ≠ pathway activity |

## AI usage declaration

We used AI tools extensively in this project, mainly Claude Code (Anthropic, Claude Opus 5.5). Claude Code wrote and ran nearly all of the project's Python and R code. That covered data acquisition and auditing, preprocessing, DESeq2 factorial modelling, pathway, concordance and robustness analyses, automated tests, and all figures. It also drafted the project reports and this report's first draft. [TEAM TO CONFIRM: other tools, e.g. ChatGPT for reviewing the methodology proposal and AI assistance for slides.] AI-generated results, interpretations and candidate hypotheses were treated as drafts, not accepted findings. The research problem, project scope and key methodological decisions were set by the team, and the team takes responsibility for the hypotheses, interpretations and conclusions submitted. Full details and items still under review: `reports/ai_usage_declaration.md` in the project repository.

## References

1. Jung W., Abushukair H., … Oh T.G. Spatial transcriptome mapping identifies Ppara-Anxa2 cross-talk in microplastic-induced hepatotoxicity. *Sci. Adv.* 12(25): eaec8681 (2026). doi:10.1126/sciadv.aec8681
2. Love M.I., Huber W., Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol.* 15:550 (2014).
3. Soneson C., Love M.I., Robinson M.D. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research* 4:1521 (2015).
4. Patro R. *et al.* Salmon provides fast and bias-aware quantification of transcript expression. *Nat. Methods* 14:417–419 (2017).
5. Anders S., Huber W. Differential expression analysis for sequence count data. *Genome Biol.* 11:R106 (2010).

**[TEAM: add references for the limma fry/camera, fgsea, edgeR, ashr and MSigDB methods papers, and the microplastic and MASH background literature, after checking each citation.]**

## Appendix: key files

| Topic | File |
|---|---|
| Data audit | `reports/dataset_audit_report.md` |
| Exploratory analysis | `reports/eda_report.md` |
| Quantification validation | `reports/quantification_validation/` |
| Analysis plan, deviations, team decisions | `reports/phase2/analysis_plan.md` |
| Model and diagnostics | `reports/phase2/model_and_diagnostics.md` |
| Pathways | `reports/phase2/pathway_report.md` |
| Similarity (A) | `reports/phase2/pe_mash_similarity_report.md` |
| Interaction (B) | `reports/phase2/interaction_report.md` |
| Robustness (C) | `reports/phase2/robustness_report.md` |
| Literature comparison | `reports/phase2/literature_comparison.md` |
| Interpretation and decisions | `reports/phase2/interpretation.md` |
