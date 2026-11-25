# Supplementary figures

These are copies of figures produced by the analysis code; the original files are unchanged at the paths shown. Report Figures 1 and 2 are drawn by `research/section3_1_figures/scripts/make_figures.py`. All PNGs are 300 dpi; PDFs are vector where available.

| # | Figure | Report section | Files | Original |
|---|---|---|---|---|
| S01 | PCA of all 19 mice before and after excluding 135 muscle-associated genes (report Fig. 1 is panel b) | §3.1, Fig. 1 | [PNG](S01_pca_before_after_muscle_exclusion.png) · [PDF](S01_pca_before_after_muscle_exclusion.pdf) | `research/section3_1_figures/figures/fig_pca` |
| S02 | Directional fry enrichment, PE vs vehicle in CD and MD: PPAR sets plus Hallmark sets with fry FDR < 0.25 (report Fig. 2) | §3.1, Fig. 2 | [PNG](S02_fry_enrichment_main.png) · [PDF](S02_fry_enrichment_main.pdf) | `research/section3_1_figures/figures/fig_enrichment` |
| S03 | Directional fry enrichment for all 50 Hallmark and 3 PPAR sets | §3.1 | [PNG](S03_fry_enrichment_all_53_sets.png) · [PDF](S03_fry_enrichment_all_53_sets.pdf) | `research/section3_1_figures/figures/fig_enrichment_full_supplementary` |
| S04 | fry and camera results, PE vs vehicle and diet × PE interaction (enrichment audit) | §3.1 | [PNG](S04_fry_camera_PE_and_interaction_audit.png) · [PDF](S04_fry_camera_PE_and_interaction_audit.pdf) | `research/enrichment_audit/figures/pathway_PE_vs_Veh_CD_MD_interaction` |
| S05 | Eight mouse-level pathway scores for the 19 mice | §3.2, Table 2 | [PNG](S05_eight_pathway_scores_19_mice.png) · [PDF](S05_eight_pathway_scores_19_mice.pdf) | `research/pathway_followup/figures/fig1_pathway_scores_19_mice` |
| S06 | Cross-study MASH-state classifier performance (MASH-1 ↔ MASH-2) | §3.3, §4.2 | [PNG](S06_classifier_cross_study_performance.png) | `research/roadmap_phase4e/figures/fig1_cross_study_performance` |
| S07 | Training-label permutation null for the classifier AUC | §4.2 | [PNG](S07_classifier_permutation_null.png) | `research/roadmap_phase4e/figures/fig3_permutation_null` |
| S08 | Frozen classifier log-odds scores on the 19 Main PE mice | Table 3 | [PNG](S08_classifier_frozen_scores_main_PE.png) | `research/roadmap_phase4e/figures/fig4_gse294180_2x2_frozen_scores` |
| S09 | Classifier pathway-ablation: MD PE shift after removing selected features | §4.2 | [PNG](S09_classifier_pathway_ablation.png) · [PDF](S09_classifier_pathway_ablation.pdf) | `research/pathway_followup/figures/fig2_classifier_ablation` |
| S10 | Held-out R² by representation (F4, F8, blocks, FULL) and model | §3.4, Table 4 | [PNG](S10_regression_representation_by_model.png) | `research/roadmap_phase4b/figures/fig1_representation_x_model_r2` |
| S11 | Feature-block (A/B/C) ablation | §3.4 | [PNG](S11_regression_feature_block_ablation.png) | `research/roadmap_phase4b/figures/fig2_feature_block_ablation` |
| S12 | Lipid and cholesterol gene groups vs other genes, held-out R² | Table 4 | [PNG](S12_lipid_vs_other_gene_prediction.png) · [PDF](S12_lipid_vs_other_gene_prediction.pdf) | `research/lipid_prediction_followup/figures/fig_lipid_vs_other_prediction` |
| S13 | Cholesterol-synthesis genes: OLS/Ridge vs CatBoost/HGB | Table 5 | [PNG](S13_cholesterol_linear_vs_nonlinear.png) · [PDF](S13_cholesterol_linear_vs_nonlinear.pdf) | `research/nonlinear_chol_followup/figures/fig_chol_linear_vs_nonlinear` |
| S14 | FAST_eff × CIRC_se interaction surface and held-out tests | §4.2 | [PNG](S14_interaction_FAST_x_CIRCse.png) · [PDF](S14_interaction_FAST_x_CIRCse.pdf) | `research/interaction_followup/figures/fig_interaction_FAST_CIRCse` |
| S15 | Frozen F8-OLS gene-level transfer to Hong-PE | §4.2 | [PNG](S15_hong_PE_gene_level_transfer.png) | `research/roadmap_phase5/figures/fig2_gene_level_external` |
| S16 | Per-mouse frozen scores in the Hong-PE study | §3.5 | [PNG](S16_hong_PE_per_mouse_frozen_scores.png) | `research/roadmap_phase5/figures/fig1_per_mouse_frozen_scores` |
