# Supplementary tables

These are copies of complete result tables produced by the analysis code; the originals are unchanged at the paths shown. Column definitions are in the README of each original folder.

| # | File | Content | Report section | Original |
|---|---|---|---|---|
| T00 | [report_number_crosswalk.csv](report_number_crosswalk.csv) | Every number in the final report, its source file and whether it matches (75/75 match) | all | made by `supplementary/reproducibility/verify_report_numbers.py` |
| T01 | [T01_deseq2_contrast_summary.csv](T01_deseq2_contrast_summary.csv) | DESeq2 contrast summary: genes tested and BH-FDR < 0.05 counts per contrast | §4.1 | `results/deseq2/provisional_deposited/contrast_summary.csv` |
| T02 | [T02_fry_camera_all_sets_all_contrasts.csv](T02_fry_camera_all_sets_all_contrasts.csv) | fry (directional and mixed) and camera results for 53 sets × 4 contrasts, with BH-FDR within family and contrast. fry and camera are separate columns. | §3.1, Fig. 2 | `results/pathways/provisional_deposited/limma_fry_camera.csv` |
| T03 | [T03_eight_pathway_scores_tests.csv](T03_eight_pathway_scores_tests.csv) | Eight-pathway score tests: estimate, CI, p, exact permutation p, BH q over 8 scores, sensitivity | Table 2 | `research/pathway_followup/TABLE_results_concise.csv` |
| T04 | [T04_eight_pathway_scores_per_mouse.csv](T04_eight_pathway_scores_per_mouse.csv) | Per-mouse eight-pathway scores | §3.2 | `research/pathway_followup/results/per_mouse_pathway_scores.csv` |
| T05 | [T05_classifier_cross_study_metrics.csv](T05_classifier_cross_study_metrics.csv) | Classifier cross-study metrics for every model and representation (AUC, balanced accuracy, Brier, chosen hyperparameters) | §3.3, §4.2 | `research/roadmap_phase4e/results/external_metrics.csv` |
| T06 | [T06_classifier_frozen_contrasts_main_PE.csv](T06_classifier_frozen_contrasts_main_PE.csv) | Frozen classifier contrasts on Main PE (shift, p, Holm p) | Table 3 | `research/roadmap_phase4e/results/gse294180_frozen_model_contrasts.csv` |
| T07 | [T07_classifier_pathway_ablation_contrasts.csv](T07_classifier_pathway_ablation_contrasts.csv) | Classifier pathway-ablation contrasts | §4.2 | `research/pathway_followup/results/ablation_gse294180_contrasts.csv` |
| T08 | [T08_regression_model_grids.csv](T08_regression_model_grids.csv) | Hyperparameter grids for all regression models | §3.4, §3.5 | `research/roadmap_phase4/model_search_spaces.csv` |
| T09 | [T09_regression_model_comparison_F4_F8_blocks.csv](T09_regression_model_comparison_F4_F8_blocks.csv) | Held-out regression results for every representation (F4, F8, A4, blocks, FULL) and model | §3.4, Table 4 | `research/roadmap_phase4b/model_comparison_results.csv` |
| T10 | [T10_regression_feature_block_ablation.csv](T10_regression_feature_block_ablation.csv) | Feature-block ablation results | §3.4 | `research/roadmap_phase4b/feature_ablation_results.csv` |
| T11 | [T11_lipid_gene_groups_prediction.csv](T11_lipid_gene_groups_prediction.csv) | Held-out R² by gene group (lipogenesis, cholesterol synthesis, PPARα, other) and model, with bootstrap CIs and tests | Table 4 | `research/lipid_prediction_followup/results/group_metrics_all.csv` |
| T12 | [T12_cholesterol_nonlinear_vs_linear.csv](T12_cholesterol_nonlinear_vs_linear.csv) | Cholesterol-synthesis genes: nonlinear vs linear R², RMSE and paired tests | Table 5 | `research/nonlinear_chol_followup/TABLE_results_concise.csv` |
| T13 | [T13_interaction_followup.csv](T13_interaction_followup.csv) | Interaction follow-up: top pairs vs null, full vs additive, pair-only vs additive | §4.2 | `research/interaction_followup/TABLE_results_concise.csv` |
| T14 | [T14_hong_PE_gene_level_transfer.csv](T14_hong_PE_gene_level_transfer.csv) | Gene-level transfer of frozen models to Hong-PE (Spearman, permutation p, verdict) | §4.2 | `research/roadmap_phase5/gene_level_external_results.csv` |
| T15 | [T15_hong_PE_sample_level.csv](T15_hong_PE_sample_level.csv) | Sample-level results in the external studies | §3.5 | `research/roadmap_phase5/sample_level_external_results.csv` |
| T17 | [T17_pathway_robustness_checks.csv](T17_pathway_robustness_checks.csv) | Pathway fry results under each robustness variant | §5.3 | `results/robustness/provisional_deposited/pathway_robustness.csv` |
