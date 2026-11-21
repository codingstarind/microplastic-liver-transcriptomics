# Existing work inventory (Phase 1 roadmap, Task 1)

**Prepared:** 2026-11-21, for team review. **Nothing was re-run to make this inventory.** It was built by reading the reports, result files and scripts listed below. Claude Code produced every number quoted here, and the team has not yet independently checked them (see `reports/final_report_draft.md`, status box).

**Where things live, and what is actually on GitHub**

| Location | On GitHub? | Content |
|---|---|---|
| `main` (93 files) | Yes | Phase 1: data audit, preprocessing, EDA, PCA/clustering (`scripts/01–07`, `reports/*.md`, `figures/eda`, `figures/pca_clustering`) |
| `phase2/factorial-analysis` (current branch, in sync with `origin`) | Yes (PR #3) | Methodology proposal, re-quantification pipeline and pilot, DESeq2 factorial model, pathways, concordance, interaction, robustness, final report draft, AI-usage declaration |
| `research/literature/` | **No, local and untracked** | Earlier literature review (19 papers), `paper_matrix.csv`, `dataset_audit.md`, `research_gaps.md`, per-paper notes in `papers/` |
| `research/cross_study_replication/` | **No, local and untracked** | 8-gene cross-paper replication (Djouina, Li, Chiu) with additive/factorial models, LOO, bootstrap, clustering, PCA, report and figures. **The files exist and are complete**, with scripts, intermediate outputs and a run record (`intermediate/run_info_01.json`: R 4.4.3, DESeq2 1.46.0, 2000 bootstrap draws, seed 20261121). |
| `presentation/` | **No, untracked** | Final slide deck (.pptx/.pdf) |

> **Action for the team:** `research/` and `presentation/` would be lost if this laptop failed. Decide whether to commit them, for example on a new branch.

---

## 1. What has been completed computationally

| Block | What was done | Main files |
|---|---|---|
| Data audit | Deposited GSE294180 matrix checked: integers, no missing values, IDs matched to GSM/SRA, all male, GENCODE M35 symbols, 55,836 genes × 19 samples | `reports/dataset_audit_report.md`, `reports/dataset_understanding.md` |
| Preprocessing | Filter of ≥ 10 counts in ≥ 4 samples (13,777 genes), median-of-ratios normalisation, log2/VST chosen from the mean–variance plot | `reports/preprocessing_decisions.md` |
| EDA, PCA, clustering | Library sizes, distributions, tissue-marker screen, PCA (top 500, all genes, muscle genes removed), hierarchical clustering under several distances | `reports/eda_report.md`, `figures/eda`, `figures/pca_clustering` |
| Technical metadata | ENA run metadata, FASTQ headers (one flowcell), read-pair depth, gene-counts-per-read-pair yield proxy | `metadata/sra_run_info.csv`, methodology proposal |
| Feasibility and power | Dispersion-based minimum detectable effects, pydeseq2 smoke test | `results/methodology/` |
| Quantification pilot | 5 samples re-quantified (fastp → Salmon, GENCODE M35) under pre-registered decision rules; comparison of all three tximport modes | `reports/quantification_validation/` |
| Factorial DE | DESeq2 `~ diet * PE`, 4 contrasts, ashr shrinkage, diagnostics, sign check against raw means | `analysis/deseq2/`, `results/deseq2/provisional_deposited/` |
| Alternative DE methods | edgeR-QL and limma-voom with the same design | `analysis/deseq2/run_alternative_methods.R` |
| Pathways | MSigDB 2026.1.Mm Hallmark + 3 PPAR sets. fry is primary; camera, fgsea and mean-z are supporting. Each pathway's PE response classified across diets. | `analysis/pathways/`, `reports/phase2/pathway_report.md` |
| PE vs diet concordance | Split-reference Spearman correlation that removes shared-control bias, with an exhaustive relabelling null | `analysis/concordance/`, `reports/phase2/pe_mash_similarity_report.md` |
| Interaction | Gene-level interaction list, classification against muscle-tracking genes | `analysis/interaction/`, `reports/phase2/interaction_report.md` |
| Robustness | 19 LOO refits; drop MD_1; muscle, depth, yield and clock covariates; gene-class checks | `analysis/robustness/`, `reports/phase2/robustness_report.md` |
| Cross-study gene panel (local only) | 8 published genes tested against our contrasts, with a 5-category agreement rule, genome-wide BH vs panel Holm, additive vs factorial (LRT, AIC, adjusted R²), stratified bootstrap, gene clustering and panel PCA | `research/cross_study_replication/` |
| Literature (local only) | 19 papers (16 full text), a dataset-availability audit and a research-gaps list | `research/literature/` |
| Reporting | Final report draft (2026-11-02), slides, AI-usage declaration, CI (ruff, 33 tests, validators) | `reports/final_report_draft.md`, `presentation/` |

## 2. Results that are reasonably stable

"Stable" means the result survives LOO, the muscle covariate, dropping MD_1, the clock covariate and the software comparison. **None of these checks can separate PE from library depth or yield** (Section 6).

| Result | Evidence of stability | Source |
|---|---|---|
| Strong diet effect (2,262 genes FDR < 0.05) and unsupervised diet separation once muscle genes are removed (ARI = 1.0) | LOO 1,699–2,537 genes; edgeR/limma 2,080/2,101; robust to depth | robustness matrix G3 |
| Fasn, Acly and Acaca lower with PE under **both** diets (−0.7 to −1.8 log2) | 19/19 LOO keep direction and FDR; 2000/2000 bootstrap draws keep sign; muscle/clock/−MD_1 fits keep it. Point estimates barely move under depth/yield adjustment, but FDR is lost for Fasn. | `replication_report.md` §4, §6 |
| Cholesterol-homeostasis program lower with PE under both diets, interaction ≈ 0 | fry p < 0.05 in 19/19 LOO under both diets; lost under the yield covariate | robustness P2, P3 |
| PPARα targets higher with PE under CD | 19/19 LOO same direction (11/19 at FDR < 0.05); lost under depth/yield | robustness P1 |
| No detectable **pathway-level** interaction (fry, 53 sets) | 18/19 LOO; all covariate fits | robustness P4 |
| No within-diet Ppara induction ≥ ~1.4-fold | CI −0.14 to +0.45 log2 under CD; stable across LOO | replication §4.4 |
| Library depth and yield completely separate PE from vehicle libraries | Read-pair counts 13.5–16.5 M vs 17.0–19.1 M; one flowcell | final report §4.1 |
| CD+Veh mice have a distinct clock-gene phase | All 4 CD+Veh mice are below every other mouse on a pre-defined score | similarity report |

## 3. Results that are fragile or exploratory

| Result | Why it is fragile |
|---|---|
| 87 gene-level interaction genes | 23 are muscle-tracking. 16 remain with a muscle covariate, 49 without MD_1, and **1** in edgeR or limma-voom. |
| Acaca interaction (+0.91 log2, Wald FDR 0.038) | Significant in only 6/19 LOO refits. FDR 0.48 with the muscle covariate. May be a floor effect, since diet already lowers Acaca 2.9-fold. |
| PE(CD)–diet resemblance (split ρ ≈ 0.33) | Minimum attainable p = 0.048 (p range 0.048–0.095). Shares the CD+Veh reference, which has a distinct clock state. Absent under MD (ρ ≈ 0.05). |
| PE_in_MD gene count (472) | Falls to 135 with the muscle covariate; 6.6% of hits are muscle-tracking |
| Every PE gene-level and pathway result | Depth/yield covariates have R² ≈ 0.72–0.76 on PE alone. Adjusting cuts PE_in_CD from 858 to 116 or 66 genes and removes the P1/P2 pathway signals. **Unresolvable with this design.** |
| Abcc3 and Nr1i3 up with PE (Li 2025 direction) | Within-diet FDR 0.05–0.30; additive model FDR 0.005/0.042; lost under depth/yield |
| Competitive-test hint that PPAR/peroxisome induction is weaker under MD | camera/fgsea only (gene-level, anti-conservative); fry p = 0.21 |
| Panel clustering and PCA | 8 hand-picked genes; descriptive only; same mice used for testing |

## 4. Published findings already reproduced using the **same** data

This is **reproducibility** (same mice, same reads), not independent validation.

| Jung 2026 claim | Our status |
|---|---|
| PE-responsive lipid genes (Fasn, Acly, Me1, Pnpla3 down; Plin4 up) | Reproduced as within-diet PE effects |
| Many PE-responsive genes under CD, fewer under MD | Reproduced in direction (858 vs 472; their thresholds give 503 CD-only, 182 MD-only, 97 shared) |
| Diet separates samples | Reproduced |
| Ppara induced by PE | **Not reproduced.** Our PE_in_CD estimate is +0.16 log2 (FDR 0.62). **Correction (2026-11-21, from Jung's supplement):** Jung's own tables give +0.53 for PE within CD (padj 0.011) and +0.34 within MD (padj 0.23, not significant). The figure legend compares only against CD+Veh. **The gap between their +0.53 and our +0.16 on the same counts is unexplained.** Candidate causes are pairwise subsetting versus a joint 4-group dispersion fit, filtering, and shrinkage. *This needs checking before any Ppara claim.* |
| "Compounded" MD+PE effect | **Not supported as an interaction.** **Correction:** Jung's Table S7 *does* contain a DESeq2 interaction analysis (76 genes, padj < 0.05). It is never mentioned in their Methods or Results, and it is dominated by skeletal-muscle genes (Acta1, Myh1/4, Ttn, Ckm), the same artefact we found. Our earlier claim that "no study fitted an interaction" must be reworded: **an interaction was computed but neither reported nor interpreted.** |
| Anxa2 / PPARα–Anxa2 axis, inflammatory hotspots | Not testable in bulk (Anxa2 shows no PE effect). These rest on Xenium and on in-vitro and ChIP data. |

Cross-study gene comparisons (`research/cross_study_replication/`) are **replication attempts with different data**, and none counted as success. Djouina's Fasn/Acly/Acaca increases came out in the opposite direction, and Ppara showed no change. Li's Abcc3/Nr1i3 matched in direction but were inconclusive. Acot3 and Nr4a1 were uninformative.

## 5. What we have **not** independently validated

- **The counts.** Only 5 of 19 runs were re-quantified, and the verdict was "not identifiable", with our totals 4.5–7.2% higher. Assumptions A1–A3 are a team decision, not a verification. The full re-quantification pipeline and manifest exist but were never run.
- **Any PE effect against an independent dataset.** No external PE liver transcriptome exists in our analysis. GSE245069 (polystyrene) has never been re-analysed.
- **The MASH axis.** Direction A uses a diet signature built from the same 19 mice. No external MASH reference was used, and the candidates in the similarity report (GSE329899, GSE342274 and others) were **never vetted**.
- **Biological validation.** There is no qPCR, protein or flux data for these mice, and no per-mouse phenotypes (ALT, TG, histology) were deposited.
- **Sacrifice time, food intake, cage and library-prep batch.** These are unknown or aligned with group.
- **Literature fields.** Chiu 2025 was read as abstract only. Djouina fold changes were read by eye from bar charts. Li 2025 Fig. 8 asterisks and supplement were not read. Supplementary DEG tables S1–S18 of Jung 2026 were not obtained.
- **Xenium questions.** Is the slide confounded with PE (barcode prefixes `G12_`/`G34_` split exactly by treatment)? Are the Xenium samples the same mice as the bulk libraries? Both are open.
- **The team has not checked any of the computational output.**

## 6. Limitations that must stay visible in all future work

1. **PE ≡ low depth and yield.** Complete separation, so no model can isolate a PE effect from a library effect in GSE294180.
2. **CD+Veh ≡ distinct clock or feeding state.** This affects every PE_in_CD and Diet_in_Veh comparison for SREBP and PPARα programs.
3. **Skeletal-muscle carry-over** (MD_1 5,190 CPM; CD_MP_1 539; MD_2 316; MD_3 161). It drives most gene-level interactions.
4. **Transcript-length offset missing** from the deposited matrix (assumptions A1–A3). Mitochondrial genes and pseudogenes are less reliable.
5. **Group sizes of 4/5/5/5, with cage ≡ group and library-prep order ≡ group.** Interaction power is low: median CI width 0.80 log2, and about 0.5–0.7 log2 is needed for typical genes.
6. **Bulk tissue.** Cell-composition shifts (Kupffer, stellate, cholangiocyte markers) cannot be separated from per-cell expression.
7. **Selection effects.** The cross-study gene panel was examined after preliminary inspection of these data.

## 7. Analyses missing or unfinished

| Item | Status |
|---|---|
| Formal adoption of hypotheses (the report still says "[TEAM]") | Not done. Every analysis is labelled exploratory. |
| Full re-quantification (14 runs) | Pipeline ready, not run |
| External MASH reference signature | Candidates listed, not vetted |
| Re-analysis of any independent microplastic dataset with our model | Not done |
| Pathway-level cross-study comparison (as opposed to 8 genes) | Not done (deliberately not run in the replication work) |
| Xenium `.h5` inspection (slide confounding, mouse IDs) | Not done (download incomplete) |
| Contact authors (sacrifice times, food intake, slide map) | Not done |
| Committing `research/` and `presentation/` | Not done |
| Reference list for methods (limma, fgsea, edgeR, ashr, MSigDB) in the final report | Marked [TEAM] |
| Final report due 2026-11-25 (per the draft header) | **Check whether this restart changes the submission plan** |
