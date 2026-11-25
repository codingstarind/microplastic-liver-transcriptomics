# An Analysis of Microplastic Effects on the Liver (GSE294180)

BT3041 – Analysis and Interpretation of Biological Data, IIT Madras, November 2026.
Team: Soham, Aryan, Arjun, Mia, Liya.

**Status: complete.** The final report is [`reports/final_report.pdf`](reports/final_report.pdf). This repository holds the code, results and supplementary material behind it. The supplement starts at [`supplementary/README.md`](supplementary/README.md).

## Project overview

We re-analysed public mouse-liver bulk RNA-seq data to ask how oral polyethylene (PE) microplastic exposure relates to hepatic gene expression, under a control diet and under a metabolic-dysfunction-inducing diet. We used exploratory analysis, differential expression, pathway enrichment and machine learning with independent reference studies.

The four research questions, as in the report:

1. Which genes and biological pathways show changes associated with PE exposure, and how do these changes compare between diets?
2. Do PE-exposed mice show gene-expression patterns resembling those associated with a MASH-inducing diet?
3. Can information from other biological experiments help predict PE-associated changes in cholesterol-related genes?
4. Do more complex models capture relationships that improve these predictions?

## Dataset

| | |
|---|---|
| Main PE study | [GSE294180](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294180) (BioProject PRJNA1248512); Jung W. *et al.* (2026) *Sci. Adv.* 12:eaec8681, [doi:10.1126/sciadv.aec8681](https://doi.org/10.1126/sciadv.aec8681) |
| Animals | 19 male C57BL/6 mice, liver; bulk poly(A) RNA-seq (NextSeq 2000, 2 × 150 bp) |
| Exposure | PE microspheres, 2 mg/day, 10–150 µm, oral gavage, or vehicle, for 8 weeks |
| Input used | `GSE294180_MP_raw_counts_gene_level.csv.gz` (deposited Salmon → tximport gene counts) |

| Diet | Vehicle | PE |
|---|---|---|
| Control diet (CD) | 4 | 5 |
| Metabolic-dysfunction diet (MD) | 5 | 5 |

**Reference studies:**
- MASH-1 [GSE316023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316023) and MASH-2 [GSE246328](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246328): GAN diet vs chow.
- FAST [GSE107787](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107787): fasted vs fed.
- CLOCK [GSE117134](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117134): six sampling times.
- PPARα [GSE282998](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282998): GW7647 vs control.
- Hong-PE [GSE261893](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261893): independent PE-response check.

Full citations are in [`supplementary/references/`](supplementary/references/).

**Raw data are not committed.** The scripts download them and verify them against recorded checksums ([`metadata/data_provenance.json`](metadata/data_provenance.json)).

## Completed analysis phases

| Phase | What was done | Where | Report section |
|---|---|---|---|
| 1. Acquisition, audit, preprocessing | GEO download with checksums; count audit; filter ≥ 10 reads in ≥ 4 mice (13,777 genes); median-of-ratios normalisation | `scripts/01–04`, [`reports/dataset_audit_report.md`](reports/dataset_audit_report.md), [`reports/preprocessing_decisions.md`](reports/preprocessing_decisions.md) | §2.2 |
| 2. EDA, PCA, clustering | PCA on the 500 most variable genes before and after excluding 135 muscle-associated genes; hierarchical clustering with ARI | `scripts/05–06`, [`reports/eda_report.md`](reports/eda_report.md), [`research/section3_1_figures/`](research/section3_1_figures/) | §3.1, Fig. 1 |
| 3. Differential expression | DESeq2 `~ diet + PE + diet:PE`; edgeR and limma-voom comparison; re-quantification pilot | [`analysis/deseq2/`](analysis/deseq2/), [`pipelines/quantification/`](pipelines/quantification/), [`reports/phase2/`](reports/phase2/) | §3.1, §4.1 |
| 4. Pathway enrichment | Directional fry (primary) and camera on 50 Hallmark + 3 PPAR sets; BH within family; fgsea supporting; enrichment audit | [`analysis/pathways/`](analysis/pathways/), [`research/enrichment_audit/`](research/enrichment_audit/) | §3.1, Fig. 2 |
| 5. Robustness | Leave-one-mouse-out; muscle, depth, clock and quality covariates; drop MD_1 | [`analysis/robustness/`](analysis/robustness/), [`reports/phase2/robustness_report.md`](reports/phase2/robustness_report.md) | §2.2, §5.3 |
| 6. Reference datasets and signatures | Download, QC and gene universe for the reference studies; MASH, FAST, CLOCK and PPARα signatures | [`research/roadmap_phase2/`](research/roadmap_phase2/), [`research/roadmap_phase3/`](research/roadmap_phase3/) | §2.1 |
| 7. Eight-pathway scoring | Per-mouse singscore for 8 programmes; diet × PE linear model; BH over 8 scores; exact permutations; leave-one-out; covariates | [`research/roadmap_phase4d/`](research/roadmap_phase4d/), [`research/pathway_followup/`](research/pathway_followup/) | §3.2, Table 2 |
| 8. Cross-study MASH-state classification | L2 logistic regression and other models trained on MASH-1/MASH-2, frozen and applied to Main PE; Holm over 3 contrasts; pathway ablation | [`research/roadmap_phase4e/`](research/roadmap_phase4e/), [`research/pathway_followup/`](research/pathway_followup/) | §3.3, Table 3 |
| 9. Gene-level regression | F4, F8 and Blocks A/B/C (FULL = 28 features); OLS, Ridge, Elastic Net; 20 co-expression modules with nested grouped CV | [`research/roadmap_phase4/`](research/roadmap_phase4/), [`research/roadmap_phase4b/`](research/roadmap_phase4b/), [`research/roadmap_phase4c_diagnostics/`](research/roadmap_phase4c_diagnostics/) | §3.4, Table 4 |
| 10. CatBoost and HGB | Same splits and F8 inputs as the linear models | [`research/roadmap_phase4/`](research/roadmap_phase4/) | §3.5 |
| 11. Cholesterol-gene prediction and interactions | Lipid and cholesterol gene groups vs other genes; nonlinear vs linear for 21 cholesterol-synthesis genes; interaction diagnostics with permutation controls | [`research/lipid_prediction_followup/`](research/lipid_prediction_followup/), [`research/nonlinear_chol_followup/`](research/nonlinear_chol_followup/), [`research/interaction_followup/`](research/interaction_followup/) | §4.2, Tables 4–5 |
| 12. Independent PE validation | Frozen F8-OLS applied to Hong-PE; Spearman correlation, mouse-label permutation, bootstrap; transfer audit | [`research/roadmap_phase5/`](research/roadmap_phase5/), [`research/roadmap_phase5b/`](research/roadmap_phase5b/) | §3.5, §4.2 |

**Supporting work not reported as results:**
- [`research/literature/`](research/literature/): literature review.

Each folder has its own README with commands, inputs and outputs. [`research/README.md`](research/README.md) lists them in order.

## Main findings (as reported)

- **Exploration:**
  - After muscle-gene exclusion, PC1 explained 20.8% and PC2 17.9% of variation.
  - Two-group hierarchical clustering matched diet exactly (ARI = 1.00).
  - DESeq2 found 2,262 diet-associated genes, 858 PE-associated genes under CD and 472 under MD (BH-FDR < 0.05).
- **Pathways:**
  - In the primary directional fry screen, PPARα targets were enriched upwards after PE under CD (FDR = 0.0406), not under MD (FDR = 0.5175).
  - No Hallmark pathway met FDR < 0.05.
  - In the mouse-level scores, lipogenesis and cholesterol-synthesis scores decreased and a fasting-associated score increased under PE in both diets. Under MD, cholesterol-homeostasis and DNA-repair scores also decreased (Table 2).
- **Classification:**
  - The cross-study L2 logistic classifier reached ROC-AUC = 1.00 in both directions, with permutation p = 0.0995 and 0.0348.
  - On Main PE, the MASH-diet-like score rose under MD (+1.663 log-odds, Holm p = 0.000623). The CD shift and the diet × PE difference were not significant.
  - Removing both cholesterol inputs reduced the MD shift to about +0.01.
- **Regression:**
  - F8-OLS predicted the 21 cholesterol-synthesis genes with held-out R² = 0.648 (CD) and 0.600 (MD), against 0.052 and 0.057 for other genes.
  - Across all genes, R² rose from about 0.07 to 0.10 with feature engineering.
- **Nonlinear models:**
  - CatBoost and HGB gave no consistent gain over OLS (CatBoost 0.659 CD, 0.424 MD).
  - One lipogenesis interaction (R² 0.163 → 0.292, CD) had Holm p = 0.054, which did not meet the 5% threshold.
- **Independent PE data:** the frozen F8-OLS gave Spearman ρ ≈ 0.22 on Hong-PE, but the permutation test (p = 0.25) did not establish PE-specific generalisation.

All 75 numbers in the report were checked against the saved outputs ([`supplementary/tables/report_number_crosswalk.csv`](supplementary/tables/report_number_crosswalk.csv)).

## Limitations (as reported)

- **Sample size:** only 19 mice, 4–5 per group.
- **Technical differences:** sequencing depth differed with PE, some samples showed muscle-associated expression, and sampling conditions varied. Any of these may influence the results.
- **Pathway measures:** enrichment and mouse-level scores are not equivalent measurements, and some scores were sensitive to technical adjustment.
- **Exploratory follow-ups:** the cholesterol prediction, classifier ablation and interaction search were done after the earlier findings were seen. Their uncertainty estimates may be optimistic.
- **External check:** the Hong-PE dataset (3 vs 3 mice) was too small to verify PE-specific prediction.
- **No mechanism claimed:** the findings are related forms of evidence from the same mice, not independent replications or proven mechanisms.

## Reproducing the analysis

The full ordered run sheet, with environments and exact package versions, is in [`supplementary/reproducibility/README.md`](supplementary/reproducibility/README.md). In short:

**1. Python environment.** Python 3.14.6.

```powershell
python -m venv .venv
```

```powershell
.venv\Scripts\pip install -r requirements.txt
```

```powershell
.venv\Scripts\pip install -e .
```

**2. CatBoost.** Install `catboost==1.2.10` into the interpreter used for the CatBoost steps (`C:\Python314\python.exe` in our runs).

**3. R and quantification tools.** In WSL Ubuntu 24.04, create the R/Bioconductor environment from `envs/r-analysis.yml` and the quantification tools from `envs/quantification.yml`, using micromamba. Recorded versions: R 4.4.3, DESeq2 1.46.0, limma 3.62.1.

**4. Phase 1** (data, audit, preprocessing, EDA, PCA). Run each step in order:

```powershell
.venv\Scripts\python.exe scripts\01_download_data.py
```

```powershell
.venv\Scripts\python.exe scripts\02_build_sample_metadata.py
```

```powershell
.venv\Scripts\python.exe scripts\03_audit_counts.py
```

```powershell
.venv\Scripts\python.exe scripts\04_preprocess_counts.py
```

```powershell
.venv\Scripts\python.exe scripts\05_descriptive_eda.py
```

```powershell
.venv\Scripts\python.exe scripts\06_pca_clustering.py
```

**5. DESeq2, pathways and robustness.** Follow [`reports/phase2/README.md`](reports/phase2/README.md#reproducing-phase-2).

**6. Research phases and follow-ups.** Run them in the order listed in [`research/README.md`](research/README.md).

**7. Check the report's numbers.**

```powershell
.venv\Scripts\python.exe supplementary\reproducibility\verify_report_numbers.py
```

**Tests and lint.** `pytest` runs the 33 unit tests and `ruff check .` checks the code. GitHub Actions runs both.

## Repository layout

```
reports/          final report (PDF), Phase 1–2 reports, earlier report draft, AI declaration
supplementary/    supplement: AI-use record, prompts, chat histories, figures, tables, references, reproducibility
research/         research phases 1–5C and follow-up analyses (one folder each, with README)
scripts/          Phase 1 pipeline steps 01–09 (run in order)
analysis/         DESeq2, pathway, concordance, interaction, robustness, Phase 2 figures
pipelines/        re-quantification (fastp → Salmon → tximport)
src/mpliver/      shared Python functions
results/          Phase 1–2 result tables (DESeq2, pathways, robustness)
figures/          Phase 1–2 figures (EDA, PCA, DESeq2, Phase 2)
metadata/         sample table, study design, data dictionary, provenance and checksums
envs/             conda environment files (R/Bioconductor; quantification)
presentation/     milestone slides
notebooks/        Phase 1 walkthrough notebook
tests/            unit tests (pytest)
data/             raw and processed data (not committed; regenerated by scripts)
```

## Reports and supplementary material

- **Final report:** [`reports/final_report.pdf`](reports/final_report.pdf).
- **Supplement index:** [`supplementary/README.md`](supplementary/README.md).
  - [Figures S01–S16](supplementary/figures/README.md)
  - [Tables T00–T17](supplementary/tables/README.md)
  - [References](supplementary/references/README.md)
  - [Reproducibility](supplementary/reproducibility/README.md)
- **Report Figures 1–2:** [`research/section3_1_figures/`](research/section3_1_figures/).
- **Earlier written reports:**
  - [dataset understanding](reports/dataset_understanding.md)
  - [audit](reports/dataset_audit_report.md)
  - [EDA](reports/eda_report.md)
  - [methodology proposal](reports/methodology_proposal.md)
  - [Phase 2 reports](reports/phase2/README.md)
  - [quantification validation](reports/quantification_validation/quantification_validation_report.md)
- **Earlier draft:** [`reports/final_report_draft.md`](reports/final_report_draft.md) covers Phases 1–2 only and is superseded by the final report. It is kept for the record.
- **Slides:** [`presentation/`](presentation/).

## AI use

AI tools were used extensively in this project, mainly **Claude Code** (Anthropic). The session records show **Claude Opus 5.5** (`claude-opus-5-5`) as the main model; Claude Opus 5 and Claude Haiku 5.5 were also used, for the literature review and for sub-agent searches.

Claude Code was used for:
- code and debugging;
- running the analyses;
- figures;
- literature searches and summaries;
- drafting reports and documentation.

As stated in the final report, the students chose the research questions, reviewed the analyses and made the final decisions about interpretation, and the team is responsible for the report's contents.

Where to find the records:
- **Task-by-task record:** [`supplementary/ai_usage.md`](supplementary/ai_usage.md), covering tools, model versions, libraries and their versions.
- **Prompts:** [`supplementary/prompts/`](supplementary/prompts/).
- **Transcripts and conversation link:** [`supplementary/chat_histories/`](supplementary/chat_histories/).
- **Earlier declaration:** [`reports/ai_usage_declaration.md`](reports/ai_usage_declaration.md).

## Software

| | Versions used |
|---|---|
| Python | 3.14.6; numpy 2.5.1, pandas 3.0.5, scipy 1.18.1, scikit-learn 1.9.1, catboost 1.2.10, matplotlib 3.11.1, seaborn 0.13.2, joblib 1.6.0 |
| R / Bioconductor | R 4.4.3; DESeq2 1.46.0, limma 3.62.1 (fry, camera, voom), edgeR 4.4.0, fgsea 1.32.2, msigdbr 26.1.1 (MSigDB 2026.1.Mm), ashr 2.2.63, apeglm 1.28.0, tximport 1.34.0 |
| Quantification | salmon 1.10.3, fastp 0.23.4, GENCODE M35 |
| Platform | Windows 11 (Python), WSL Ubuntu 24.04 (R and quantification) |

Full installed-package records: [`supplementary/reproducibility/environment_records/`](supplementary/reproducibility/environment_records/).
