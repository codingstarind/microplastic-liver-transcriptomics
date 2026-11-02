# Milestone progress (BT3041 mid-project presentation)

**Presentation constraints:** at most 10 slides, at most 8 minutes, then 2 minutes of Q&A.
**Required content:** introduction and research questions; data acquisition and preprocessing; at least 40% completion of methods, code, results and interpretation.

Last updated: 2026-11-02 (end of Phase 1).

## Requirement-by-requirement status

| Requirement | Status | Evidence in repo | What is missing |
|---|---|---|---|
| **Introduction / motivation** | Done (draft) | `README.md`, `reports/dataset_understanding.md` §1 | Condense to 1 slide |
| **Research questions** | Candidates drafted, **not yet chosen** | `reports/preliminary_observations.md` (Q1–Q4) | Team to pick 2–3 and write them as testable hypotheses |
| **Data acquisition** | Done | `scripts/01–02`, `metadata/data_provenance.json`, `metadata/sample_metadata.csv`, `metadata/study_design.csv` | – |
| **Experimental design understood** | Done | `reports/dataset_understanding.md` | – |
| **Data audit** | Done | `reports/dataset_audit_report.md`, `results/tables/audit_summary.json` | Unresolved items listed in §9 of the audit (author queries optional) |
| **Preprocessing** | Done for exploration | `scripts/04`, `reports/preprocessing_decisions.md` | Decide on VST and on MD_1 handling |
| **Descriptive statistics and visualisation** | Done | Figs 1–9, `reports/eda_report.md` | – |
| **Dimensionality reduction (PCA)** | Done, with sensitivity analyses | Figs 10–13 | – |
| **Clustering** | Done (sample-level) | Figs 14–15, `clustering_label_agreement.csv` | Gene-level clustering / modules (Q3) |
| **Hypothesis testing** | Not started (on purpose) | – | Factorial negative-binomial model (Q2), multiple-testing control |
| **Modelling / advanced analysis** | Not started | – | Fold-change comparison (Q1), co-expression modules (Q3), robustness (Q4a–b) |
| **Interpretation** | Preliminary | `reports/preliminary_observations.md` | Final interpretation after tests |
| **Reproducibility** | Done | numbered scripts, pinned `requirements.txt`, tests, CI, provenance and checksums | – |

## Rough completion estimate

Weights are our own split of the whole project, so the percentages are approximate:

| Component | Weight | Done | Contribution |
|---|---|---|---|
| Data acquisition and understanding | 15% | 100% | 15% |
| Preprocessing and QC | 15% | 90% | 13.5% |
| Exploratory analysis (descriptive, PCA, clustering) | 20% | 90% | 18% |
| Hypothesis testing / modelling | 30% | 0% | 0% |
| Interpretation and reporting | 20% | 25% | 5% |
| **Total** | | | **≈ 50%** |

The ≥ 40% requirement is **met by the work done**, provided the team signs off on the research questions. To show the "methods" side convincingly, it would help to also have **one first formal analysis** (e.g. the Q2 factorial model on the main genes) ready before the presentation.

## Suggested 10-slide outline (8 minutes)

| # | Slide | Source | ~Time |
|---|---|---|---|
| 1 | Title, team, one-line question | – | 0:20 |
| 2 | Motivation: microplastics, PE, liver, MASH | `dataset_understanding.md` §1 | 0:50 |
| 3 | Study design: 2 × 2 table, dose, diet, n; mouse → matrix pipeline | §2–4 diagram | 1:00 |
| 4 | Data acquisition and audit: integer Salmon/tximport counts, symbols, 55,836 → 13,777 genes, provenance | `dataset_audit_report.md`, fig02 | 0:50 |
| 5 | QC: library size confounded with PE (fig01) + muscle contamination (fig09) | fig01, fig09 | 1:00 |
| 6 | Normalisation and transformation choices | fig03, fig05 | 0:40 |
| 7 | PCA: primary (PC1 = muscle) vs contamination-removed (PC1 = diet, PE within diet) | fig10, fig13 | 1:00 |
| 8 | Clustering: diet recovered exactly; PE nested within diet | fig15 | 0:40 |
| 9 | Descriptive gene patterns + η² summary | fig07, fig08 | 0:50 |
| 10 | Research questions and next steps (Q1–Q3, robustness, limitations) | `preliminary_observations.md` | 0:50 |

## To do before the milestone

- [ ] Team meeting: choose research questions and write hypotheses (Q1–Q4)
- [ ] Decide on contamination handling (exclude MD_1? muscle covariate?) and on VST
- [ ] Install R + DESeq2 (or `pydeseq2`) for VST and count models
- [ ] (Optional) First factorial model `~ diet * PE` on selected genes, for the methods slide
- [ ] Check clock genes across groups (circadian confound)
- [ ] Look at SRA run metadata (flowcell/run dates) for the library-order confound
- [ ] Build the slides and rehearse to 8 minutes
- [ ] Each team member reviews and merges the Phase 1 pull request

---

# Phase 2 progress (updated 2026-11-02)

Final report due **2026-11-25**. Details: `reports/phase2/README.md`.

| Item | Status | Where |
|---|---|---|
| R/Bioconductor + Salmon environment (WSL, micromamba) | Done | `envs/` |
| Quantification-source investigation (GEO, paper, Zenodo) | Done: setting not recoverable from sources | `reports/quantification_validation/` |
| Re-quantification pilot (4 samples + replicate) | **Pending manual downloads**; pipeline ready | `pipelines/quantification/` |
| DESeq2 factorial model (CD/Veh references, 4 contrasts, CIs, FDR) + diagnostics | Done, **provisional** | `reports/phase2/model_and_diagnostics.md` |
| Pathway analysis (Hallmark/PPAR; fry, camera, fgsea, mouse-level scores) | Done, provisional | `reports/phase2/pathway_report.md` |
| Direction A: PE vs diet concordance (split reference, relabelling null) | Done, exploratory | `reports/phase2/pe_mash_similarity_report.md` |
| Direction B: interaction (gene + pathway) | Done, provisional | `reports/phase2/interaction_report.md` |
| Direction C: robustness matrix (LOO, covariates, edgeR/limma) | Done, provisional | `reports/phase2/robustness_report.md` |
| Literature comparison | Done from the main text; **supplement tables pending** | `reports/phase2/literature_comparison.md` |
| Interpretation + team decisions | Draft | `reports/phase2/interpretation.md` |
| Independent MASH signature (external dataset) | Candidates listed; not implemented | similarity report |
| Confirmatory re-run on validated counts with team-adopted hypotheses | Not started (depends on the pilot and team decisions) | – |

**Milestone presentation:** methods, code, results and interpretation now well exceed 40%. Suggested additions to the 10-slide plan:

- slide 7: replace "PCA" with the factorial model and contrasts;
- slide 8: the robustness matrix;
- slide 9: interaction vs muscle carry-over (fig18);
- slide 10: limitations (library and circadian confounds) and decisions.
