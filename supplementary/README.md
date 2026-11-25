# Supplementary material: index

This folder supports the final report [`reports/final_report.pdf`](../reports/final_report.pdf). Plots, full tables and the prompt/response history are in the team's separate PDF supplement. This folder keeps the records that belong with the code.

| Folder / file | Contents |
|---|---|
| [`ai_usage.md`](ai_usage.md) | AI tools, model versions shown by the session records, what each was used for, and the libraries and versions used in each task |
| [`chat_histories/`](chat_histories/) | Claude conversation link and the list of Claude Code sessions with dates and models |
| [`prompts/`](prompts/) | Script that exports the typed prompts from the local session records |
| [`references/`](references/) | The report's 18 references, dataset accessions and roles, gene-set and annotation versions |
| [`reproducibility/`](reproducibility/) | Installed package versions, run order for the whole pipeline, and a read-only check of the report's numbers |
| [`tables/report_number_crosswalk.csv`](tables/report_number_crosswalk.csv) | Each number printed in the report, its source file, and whether it matches (75/75) |

## Where each part of the report comes from

| Report item | Code and full outputs |
|---|---|
| §2.2 preprocessing (55,836 → 13,777 genes) | `scripts/03_audit_counts.py`, `scripts/04_preprocess_counts.py`, `results/tables/` |
| Fig. 1 PCA; ARI = 1.00 | `scripts/06_pca_clustering.py`, `research/section3_1_figures/` |
| DESeq2 gene counts (2,262 / 858 / 472) | `analysis/deseq2/`, `results/deseq2/provisional_deposited/` |
| Fig. 2 fry enrichment; PPARα targets CD FDR 0.0406 | `analysis/pathways/pathway_analysis.R`, `results/pathways/provisional_deposited/`, `research/enrichment_audit/` |
| Table 2 eight-pathway scores | `research/pathway_followup/`, `research/roadmap_phase4d/` |
| Table 3 classifier; AUC 1.00; permutation p | `research/roadmap_phase4e/` |
| Classifier ablation (+1.66 → +0.01) | `research/pathway_followup/` |
| Table 4 regression (F4, F8, Blocks A/B/C, FULL) | `research/roadmap_phase4/`, `research/roadmap_phase4b/`, `research/roadmap_phase4c_diagnostics/`, `research/lipid_prediction_followup/` |
| Table 5 CatBoost and HGB | `research/roadmap_phase4/`, `research/nonlinear_chol_followup/` |
| Lipogenesis interaction (0.163 → 0.292; Holm p = 0.054) | `research/interaction_followup/` |
| Hong-PE transfer (ρ ≈ 0.22; p = 0.25) | `research/roadmap_phase5/`, `research/roadmap_phase5b/` |
| Robustness and sensitivity (§5.3) | `analysis/robustness/`, `results/robustness/`, `research/pathway_followup/` |
