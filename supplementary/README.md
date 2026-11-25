# Supplementary material: master index

This folder is the supplement to the final report *An Analysis of Microplastic Effects on the Liver* (BT3041, IIT Madras, November 2026). It collects the files the report refers to as the "separate supplement":
- full gene-set tests;
- pathway and classifier-ablation plots;
- model grids;
- validation details;
- the question/prompt/response history.

Everything here was copied from, or generated from, the analysis outputs in this repository. No result was recomputed or changed.

| Folder / file | Contents |
|---|---|
| [`ai_usage.md`](ai_usage.md) | AI tools, model versions shown by the records, what each was used for, and the libraries and versions used in each task |
| [`chat_histories/`](chat_histories/) | Conversation link and readable transcripts (prompts and replies) of the seven Claude Code project sessions |
| [`prompts/`](prompts/) | Every prompt the team typed to Claude Code, by session |
| [`figures/`](figures/) | 16 supplementary figures (S01–S16), each mapped to its report section and original file |
| [`tables/`](tables/) | 17 complete result tables (T01–T17) and the report-number crosswalk (T00) |
| [`references/`](references/) | The report's 18 references, dataset accessions and their roles, gene-set and annotation versions |
| [`reproducibility/`](reproducibility/) | Environment records (installed package versions), run order for the whole pipeline, and a check of the report's numbers against the saved outputs |

## Where each part of the report comes from

| Report item | Supplement | Code and full outputs |
|---|---|---|
| §2.2 preprocessing (55,836 → 13,777 genes) | T00 | `scripts/03_audit_counts.py`, `scripts/04_preprocess_counts.py`, `results/tables/` |
| Fig. 1 PCA; ARI = 1.00 | S01 | `scripts/06_pca_clustering.py`, `research/section3_1_figures/` |
| DESeq2 gene counts (2,262 / 858 / 472) | T01 | `analysis/deseq2/`, `results/deseq2/provisional_deposited/` |
| Fig. 2 fry enrichment; PPARα targets CD FDR 0.0406 | S02, S03, S04, T02 | `analysis/pathways/pathway_analysis.R`, `results/pathways/provisional_deposited/`, `research/enrichment_audit/` |
| Table 2 eight-pathway scores | S05, T03, T04 | `research/pathway_followup/`, `research/roadmap_phase4d/` |
| Table 3 classifier; AUC 1.00; permutation p | S06–S08, T05, T06 | `research/roadmap_phase4e/` |
| Classifier ablation (+1.66 → +0.01) | S09, T07 | `research/pathway_followup/` |
| Table 4 regression (F4, F8, Blocks A/B/C, FULL) | S10–S12, T08–T11 | `research/roadmap_phase4/`, `research/roadmap_phase4b/`, `research/roadmap_phase4c_diagnostics/`, `research/lipid_prediction_followup/` |
| Table 5 CatBoost and HGB | S13, T12 | `research/roadmap_phase4/`, `research/nonlinear_chol_followup/` |
| Lipogenesis interaction (0.163 → 0.292; Holm p = 0.054) | S14, T13 | `research/interaction_followup/` |
| Hong-PE transfer (ρ ≈ 0.22; p = 0.25) | S15, S16, T14, T15 | `research/roadmap_phase5/`, `research/roadmap_phase5b/` |
| Robustness and sensitivity (§5.3) | T17 | `analysis/robustness/`, `results/robustness/`, `research/pathway_followup/` |

**Number check.** All 75 numbers printed in the report were checked against these source files by [`reproducibility/verify_report_numbers.py`](reproducibility/verify_report_numbers.py), and all match ([`tables/report_number_crosswalk.csv`](tables/report_number_crosswalk.csv)).
