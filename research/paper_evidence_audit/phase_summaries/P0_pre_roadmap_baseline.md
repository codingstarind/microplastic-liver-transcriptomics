# P0: Pre-roadmap baseline analysis (course Phases 1–2, before `roadmap_*`)

> **Naming warning.** The original course project used "Phase 1" (audit/EDA) and "Phase 2" (factorial DE) on branches `main` and `phase2/factorial-analysis`. These are **not** `roadmap_phase1` and `roadmap_phase2`. The paper must keep the two numbering schemes apart. This document calls the earlier work **P0**.

## Main sources
- `reports/final_report_draft.md` (2026-11-02), which states that the team had not checked it.
- `reports/phase2/*.md`
- `research/roadmap_phase1/existing_work_inventory.md`
- `research/cross_study_replication/replication_report.md`

## 1. Question
- Re-analyse GEO GSE294180 (Jung et al. 2026, *Sci Adv*) with an explicit 2 × 2 factorial design: diet (CD vs MD) × gavage (vehicle vs PE).
- **Direction A:** does the PE response resemble the MASH diet?
- **Direction B:** does the PE response depend on diet (interaction)?
- **Direction C:** is it robust to technical limitations?

## 2. Why
The original paper used pairwise DEG lists. These cannot:
- separate interaction from additivity;
- measure PE–diet similarity without bias from shared controls.

The paper also reported no bulk QC.

## 3. Data
- GSE294180: 19 male C57BL/6 mice.
  - CD = D09100304 purified diet; MD = D09100310 (GAN-type MASH diet).
  - PE 10–150 µm, 2 mg/day by gavage, 8 weeks.
  - Group sizes: CD+Veh 4, CD+PE 5, MD+Veh 5, MD+PE 5.
- Deposited integer Salmon/tximport counts: 55,836 genes. After filtering at ≥ 10 counts in ≥ 4 samples: 13,777 genes.

## 4. Methods implemented
| Step | Method | Script / location |
|---|---|---|
| Data audit | integer, sex and ID checks | `scripts/01–03`, `reports/dataset_audit_report.md` |
| Preprocessing | filter; median-of-ratios normalisation; log2/VST chosen from the mean–variance plot | `reports/preprocessing_decisions.md` |
| Exploratory analysis | library sizes; PCA (top 500 genes; all genes; muscle genes removed); hierarchical clustering (1 − Pearson, average linkage); ARI | `scripts/04–07`, `figures/eda`, `figures/pca_clustering` |
| Quantification pilot | 5 samples re-quantified (fastp → Salmon 1.10.3, GENCODE M35); 3 tximport modes; pre-registered metrics M1–M6 | `pipelines/quantification/`, `reports/quantification_validation/` |
| Factorial DE | DESeq2 1.46 `~ diet * PE`; Wald test; BH per contrast; ashr for ranking; sign check against raw means | `analysis/deseq2/` |
| Alternative DE | edgeR-QL; limma-voom | `analysis/deseq2/run_alternative_methods.R` |
| Pathways | MSigDB 2026.1.Mm Hallmark (50) + 3 PPAR sets. fry primary; camera, fgsea and mean-z supporting | `analysis/pathways/` |
| PE vs diet concordance | split-reference Spearman (6 CD+Veh pair assignments) + exhaustive relabelling null | `analysis/concordance/` |
| Interaction audit | interaction genes vs muscle-tracking genes | `analysis/interaction/` |
| Robustness | 19 leave-one-out refits; drop MD_1; muscle, depth, yield and clock covariates | `analysis/robustness/` |
| Cross-paper gene panel | 8 published genes (Djouina, Li, Chiu); additive vs factorial (LRT/AIC); stratified bootstrap (2000 draws) | `research/cross_study_replication/` |

## 5. Key numbers
Source: `final_report_draft.md` §4. Re-verify against the result CSVs before citing.

**Differential expression (FDR < 0.05):**
| Contrast | Genes |
|---|---|
| PE under CD | 858 |
| PE under MD | 472 |
| Diet (vehicle) | 2,262 |
| Interaction | 87, of which 23 muscle-tracking. 16 remain with the muscle covariate, 49 without MD_1, and 1 in edgeR or limma |

**Muscle carry-over:** CPM is MD_1 5,190; CD_MP_1 539; MD_2 316; MD_3 161. Muscle genes form PC1 (30%). With muscle genes removed, PC1 (21%) separates diet (ARI = 1.0).

**Depth confound:** PE libraries 13.5–16.5 M read pairs vs vehicle 17.0–19.1 M (complete separation). Gene-count yield ≈ 0.77 vs 0.80.

**Lipogenic genes:** Fasn, Acly and Acaca are lower with PE in both diets (−0.7 to −1.8 log2; LOO 19/19).

**Pathways (fry):**
| Pathway | Result |
|---|---|
| Cholesterol homeostasis | down with PE: CD p = 0.004, MD p = 0.005; interaction p = 0.84 |
| PPARα targets | up with PE in CD (fry FDR 0.041) |
| Pathway interactions | none at FDR < 0.05 (53 sets) |

**Split-reference ρ (PE vs diet):**
| Comparison | ρ | p |
|---|---|---|
| PE_in_CD vs diet | 0.33 [0.29–0.40] (naive 0.54) | 0.048–0.095 |
| PE_in_MD vs diet | 0.05 | — |

**Depth/yield covariates** cut PE_in_CD from 858 genes to 116 or 66 and remove the PPARα and cholesterol pathway signals.

**Clock score:** CD+Veh −2.38 vs +0.08 to +1.01 in the other groups.

**Ppara:** our PE_in_CD estimate is +0.16 log2 (FDR 0.62). Jung's Table S5 gives +0.53 (padj 0.011). The discrepancy is unexplained (Phase 1 inventory §4).

## 6. Failures and limitations found
- Gene-level interactions were mostly muscle artefact.
- PE is completely aligned with depth/yield.
- CD+Veh mice are in a distinct clock state.
- Cage and library-prep order are aligned with group.
- Deposited counts lack length offsets (assumptions A1–A3).
- Mouse-level phenotypes were not deposited.

## 7. Why the project moved on
- "PE resembles MASH" could only be tested against a **same-data** diet axis (circular, sharing controls).
- The PE pattern (lipogenesis/cholesterol down, PPARα up) is also what **fasting or time of day** produces.
- That called for **independent reference axes** (roadmap Phases 1–3).
