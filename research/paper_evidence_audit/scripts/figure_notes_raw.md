Raw per-figure inspection notes (visual check of the PNGs via contact sheets). Indices refer to `scratchpad/sheets/index.txt`.

## [0–6] figures/deseq2/provisional_deposited/ (P0 factorial DESeq2 diagnostics, 19 mice)
- **cooks_distance_by_sample:** box per sample (log10 Cook's). Uniform; MD_MP_4 and CD_MP_2 slightly higher. Supplementary.
- **dispersion_trend:** gene-wise / fitted / final dispersion vs mean. Well behaved. Title still says [PROVISIONAL]. Supplementary.
- **ma_Diet_in_Veh / ma_Interaction / ma_PE_in_CD / ma_PE_in_MD:** log2FC vs mean; red = FDR < 0.05.
  - Many interaction "hits" sit at the ±4 cap at low counts (muscle genes).
  - All titles carry "[PROVISIONAL: deposited counts]".
  - Supplementary (one 4-panel figure).
- **pvalue_histograms:** 4 contrasts. Anti-conservative peaks; the interaction peak is the smallest. Supplementary.
- **The 26 other DESeq2 diagnostic folders** (19 LOO, cov_clock, cov_depth, cov_muscle, cov_quality_proxy, drop_MD_1, pilot_length_offset) repeat these same 7 plots for sensitivity fits. Spot-check only. EXCLUDE from the paper; cite the folder.

## [7] figures/eda/fig01_library_sizes
- A: total estimated counts per sample (millions), coloured by group. PE bars are visibly lower.
- B: size factor vs relative library size.
- **MAIN candidate:** it shows the PE ≡ depth confound.
- Caveat: these are counted reads (~10–15.5 M), not sequenced read pairs (13.5–19.1 M). The axis label is correct.

## [8] eda fig02_count_distributions
Gene mean-count histogram (55,836 genes) and per-sample density. Supplementary (preprocessing).

## [9] eda fig03_sample_boxplots_before_after_normalisation
13,777 genes. PE groups have lower medians before normalisation; aligned after. Supplementary.

## [10] eda fig04_expression_dominance
- Cumulative read fraction by gene rank: top 10 genes = 19%, top 100 = 47%, top 1000 = 74%.
- Top 15 transcripts (Alb, mt-Co1, Apoe…).
- Supplementary. Relevant to yield and length concerns.

## [11] eda fig05_mean_variance
Count-scale variance vs mean, and log2-scale SD vs mean. Justifies the transform. Supplementary.

## [12] eda fig06_sample_qc_by_group
- Panels: genes ≥ 10 counts; chrM fraction; top-gene fraction.
- PE groups have **higher chrM fraction** (CD+PE ≈ 0.115 vs CD+Veh ≈ 0.088) and fewer detected genes.
- **MAIN/SUPP:** a strong visual of library-quality differences aligned with PE.

## [13] eda fig07_selected_genes_by_group
- 12 genes (Ppara, Nr1i3, Anxa2, Plin4, Fasn, Pnpla3, Cd36, Cidea, Lcn2, Col1a1, Cyp4a14, Cyp7a1). Descriptive; no tests.
- Annotations quote the paper's claims ("paper: induced by PE").
- Supplementary. Caveat: Col1a1 is in the muscle-exclusion list.

## [14] eda fig08_variance_share_by_factor
- eta² per gene for Diet, PE and Group vs one label shuffle (13,164 genes).
- Shows diet > PE.
- Caveat: the null is a **single** shuffle (noisy). Supplementary.

## [15] eda fig09_tissue_marker_screen
- Marker CPM per sample: muscle, adipocyte, erythrocyte, pancreas, cholangiocyte, Kupffer, stellate.
- MD_1 muscle ≈ 5,000; CD_MP_1 ≈ 500.
- **MAIN/SUPP (contamination).** Log y-axes per panel.

## [16] pca fig10_pca_primary
- Top 500 genes. Scree; PC1 30.3% / PC2 17.0%; PC3/PC4.
- MD_1 is the extreme on PC1 (muscle). CD groups cluster on PC2.
- Supplementary (paired with fig13).

## [17] pca fig11_pca_top_loadings
- PC1 loadings = Myh1, Mb, Acta1, Ttn, Atp2a1, Myh2: **muscle**.
- PC2 = Chrna4, Sqle, Fabp5 vs Cfd, Cidea, Cidec.
- **MAIN/SUPP:** proves PC1 is the contamination axis.

## [18] pca fig12_pca_all_filtered_genes
All 13,777 genes; PC1 20.7%. Supplementary.

## [19] pca fig13_pca_muscle_genes_removed
- 135 muscle-tracking genes removed; PC1 20.8% now separates diet (CD left, MD right).
- **MAIN candidate.** Pair with fig10 in one figure: "contamination axis vs diet axis".

## [20] pca fig14_sample_correlation_heatmap
1 − Pearson, average linkage, top 500 genes. CD_MP_1 and MD_1–3 form an outlying muscle cluster. Supplementary.

## [21] pca fig15_dendrograms_sensitivity
- Four distance/linkage choices.
- The muscle-removed version separates diet perfectly (ARI 1.0 stated in the report).
- Supplementary.

## [22] figures/phase2/.../fig16_pathway_pe_response_forest
- 21 focus pathways. Mean-z mouse-level effect with 95% CI for PE_in_CD, PE_in_MD and interaction.
- Only SANDERSON_PPARA marked * (fry FDR < 0.05).
- Effect scale = mean z (not fry). Caveat (stated in its footnote): fry significance and the plotted mean-z estimate are different statistics.
- **MAIN candidate** for P0 pathway results.

## [23] fig17_pe_vs_diet_concordance
- Null histograms of relabelled-PE Spearman; blue lines = split-reference estimates; orange dashed = naive.
- **Title risk:** "Does the PE response resemble the MASH-diet response?" The "MASH-diet response" here is our **same-data** Diet_in_Veh, not an independent MASH reference. Retitle or caption accordingly.
- Supplementary.

## [24] fig18_interaction_vs_muscle
- Shrunken PE_in_MD vs PE_in_CD scatter: 87 interaction genes (23 muscle) in the primary model; 16 (0 muscle) with the muscle covariate.
- **MAIN/SUPP:** shows the interaction artefact.

## [25] cross_study_replication/gene_expression_heatmap
8 published genes, VST z-scores; mice in group order. Descriptive. Supplementary or exclude.

## [26] cross_study_replication/gene_panel_pca
8-gene PCA. The footer states "not a validation". Exclude or Supplementary.

## [27] cross_study_replication/gene_replication_forest
- 8 genes × (PE_in_CD, PE_in_MD, interaction); 95% Wald CI; filled = genome-wide FDR < 0.05.
- Djouina "published ↑" for Fasn/Acly/Acaca vs ours ↓: a direct visual of the direction conflict.
- Supplementary (Discussion).

## [28–33] roadmap_phase2/qc/<GSE>_qc.png (6 reference datasets)
- Each: PCA of contrast samples (top 500 genes) + technical depth by group.
- GSE107787: PC1 77% = FAST vs FED.
- GSE246328: PC1 85% = diet.
- GSE225695: n = 6.
- Supplementary (one 6-panel reference-QC figure).
- [32] GSE282998: PC1 50% = dose. Depth is clustered by dose: GW0 ~12M, GW5 ~11M, GW10/20 ~21-24M.
- [33] GSE316023: PC1 68% = diet; depth similar.

## [34] roadmap_phase3/figures/baseline_comparisons.png
- Left: Baseline 2 (Spearman of gene effects vs 8 references for Diet_in_Veh, PE_in_CD and PE_in_MD), with LOO-range bars and a dashed line at the historical split-reference value.
- Right: Baseline 1 (Jaccard of DEG lists).
- **MAIN candidate** (Phase 3 baselines). Caveat: Jaccard does not show direction, so add concordant/discordant counts in the caption.
- Rotated tick labels are small.

## [35] circadian_reference.png
- Peak-phase histogram of 2,914 rhythmic genes (female WT).
- Phase in female vs male: 1,637 genes, 80% within 3 h.
- Supplementary.

## [36–38] meanvar diagnostics
- limma-trend sqrt(sigma) vs mean for CIRC female, CIRC male and FAST; RPKM distribution.
- Supplementary only.

## [39–43] DESeq2 diagnostics of MASH1 (TA), MASH1 + length offset, MASH1 TN, MASH2 and PPARA WT
- Dispersion, MA and Cook's per sample.
- MASH2 has one sample with high Cook's (mouse 16234).
- Supplementary only. Tick labels are unreadable at print size.

## [44] fasting_reference.png
1. Main vs clean subset, ρ 0.92.
2. ZT × state interaction p-value histogram (4,296 FDR < 0.05).
3. Per-ZT fasting effect for 7 marker genes (Hmgcs2, Fasn, Srebf1, Hmgcr, Cyp4a14, Angptl4, Fgf21).
- **MAIN/SUPP.** Panel 3 is didactic: it shows the fasting effect depends on time of day.

## [45] mash_positive_control.png
1. Diet_in_Veh vs MASH1, ρ 0.36.
2. vs MASH2, 0.32.
3. MASH1 vs MASH2, 0.56.
4. Monte Carlo relabelling null (1,000 draws) with the observed value in red.
- **MAIN candidate** (positive control).
- Caveats:
  - The exact enumeration (12,870 allocations; p = 0.036, q99 = 0.390) is in `results/posthoc_audit/`, not in this figure. The caption must say "formal FAIL on c1".
  - The axis titles are tiny.

## [46] offset_audit.png
Implied vs GENCODE length (ρ 0.50) and the GAN/chow implied-length difference. Supplementary.

## [47] ppara_reference.png
- WT 20v0 vs KO 20v0 scatter (red = PPARα-dependent by the formal interaction).
- Dose–response curves for 12 genes.
- Supplementary (or a main-panel inset).

## [48] reference_similarity_heatmap.png
- 8 × 8 Spearman matrix of reference gene effects.
- MASH1–MASH2 0.56, FAST–FAST2 0.49, CIRC_T1–T2 0.86; everything else |ρ| ≤ 0.19.
- **MAIN candidate.** It shows the reference axes are nearly orthogonal.
- Caveat: the PPARA_20v0 vs PPARA_dependent cell shows 1.00 because the second is a gene subset, i.e. a restricted overlap. Say this in the caption.

## [49] taskA_complexity_vs_gain.png
- Mean R²_oos gain over the best single axis vs model complexity (0 = null … 5 = boosting), for CD and MD.
- Point labels overlap and are hard to read.
- Supplementary, or rebuild as a cleaner "model comparison matrix".

## [50] taskA_loo_stability.png
- R²_oos boxplots across 19 LOO targets for null, single-axis, ols_all, ridge and hgb.
- **MAIN/SUPP.** It shows ols_MASH is fragile in MD and ols_FAST is stable.

## [51] taskA_performance_r2.png
- Bars of mean R²_oos with NB-corrected 95% CI for 12 models × 3 targets (CD, MD, Diet_in_Veh).
- **MAIN candidate** (Phase 4 model comparison). Wide CIs are honest.
- Tick labels are tiny; rebuild as vector.

## [52] taskA_variant_stability.png
- Heatmap of mean R²_oos (5 grouped splits) for 8 variants × 9 models.
- Notable: V5 strict universe gives PE_in_MD ols_all **0.180**, rf 0.168, hgb 0.170, i.e. higher R² on the 8,360 better-detected genes.
- Supplementary.

## [53] taskB_classification_demo.png
- Balanced accuracy, internal (n = 32) vs external (n = 15). Also the predicted P(GAN) for GSE294180 vehicle mice.
- Its own title says "not a model ranking".
- Supplementary.

## [54] phase4b fig1_representation_x_model_r2.png
- Heatmaps of R²_oos for 10 representations × 5 models, CD and MD.
- \* = meets the locked rule vs F8 (same model); + = statistical only.
- **MAIN candidate** (feature-engineering matrix).
- Caveat: the cells give no CI. Pair with the paired-comparison table.

## [55] phase4b fig2_feature_block_ablation.png
- RMSE(FULL − block) − RMSE(FULL) with NB 95% CI, per model.
- Block C is the largest for MD linear models.
- **MAIN/SUPP** (ablation).

## [56] phase4b fig3_stability_per_split.png
- Per-split ΔRMSE vs F8 (box plots) for selected representations; red line = practical threshold −0.005.
- Supplementary.

## [57] phase4c fig1_collinearity.png
- VIF per FULL feature (red = VIF > 10: FAST_eff, A_FAST_shr, MASH_eff, A_MASH_shr).
- Median condition number per representation (log scale): F8 2.5 … FULL 37.
- **MAIN candidate** (why FULL coefficients are unreliable).

## [58] phase4c fig2_residual_diagnostics_F8_ols.png
- 2 × 4 panels: observed vs OOF-predicted; residual vs fitted (red decile means); residual vs target lfcSE; residual vs baseMean (CD and MD).
- Shows compression and heteroskedasticity.
- **SUPP** (or MAIN panel B for "noise dominates").

## [59] phase4c fig3_subgroups_and_influence.png
1. R²_oos by SE tercile and baseMean quintile.
2. Single-mouse RMS influence on the target (CD_MP_2 ≈ 0.19).
3. Grouped-CV R² per LOO target.
- **MAIN/SUPP.**

## [60] phase4d fig1_scores_by_group.png
- MASH, FAST and PPARA singscores per mouse by group (dots + mean bar).
- FAST: CD+Veh lowest (≈ −0.11), MD+PE highest.
- **MAIN candidate** (mouse-level result).
- Caveat: no CI on the panels.

## [61] phase4d fig2_scores_vs_technical_covariates.png
- 3 scores × 4 covariates (depth, quality_proxy, muscle, clock_phase), with ρ in the titles.
- Clock_phase ρ = 0.64 (MASH), 0.55 (FAST), 0.45 (PPARA).
- **The group colour scheme differs from fig1** (MD+Veh red here, green in fig1). Unify before publication.
- SUPP.

## [62] phase4d fig3_contrasts_vs_null_sets.png
- Contrast estimates (95% CI) vs the 95% band of random matched gene sets, for MASH, FAST and PPARA × 4 contrasts.
- **MAIN candidate.**

## [63] phase4d fig4_sensitivity.png
- Standardised effect for 28 variant × contrast labels.
- No CI, very long rotated labels, truncated y-label.
- Rebuild as a forest plot, or Supplementary only.

## [64] phase4e fig1_cross_study_performance.png
- External AUC (left) and Brier (right) per representation × model; S1 = dots, S2 = crosses.
- **The title "ceiling for every non-PCA representation" is FALSE in its own plot.** R2|tree S1 = 0.80, S2 = 0.63; R0|tree S2 = 0.81. The Phase 5 decision report corrected this. Retitle before reuse.
- The R0 "Brier" for the rule is a misclassification rate.
- **MAIN candidate after retitling.** Tiny y labels.

## [65] phase4e fig2_before_after_calibration.png
- Predicted P(GAN) of test mice by true class, for the R0 rule, R1 L2-logistic and PCA L2-logistic (S1 dots, S2 crosses).
- The R0 rule gives hard 0/1 outputs: in S2, 3 GAN mice predicted 0.
- **MAIN candidate** (the threshold-transfer story).
- Caveat: panel 2 is **R1** (not the tuned R2 model in the table). Say so in the caption.

## [66] phase4e fig3_permutation_null.png
- External AUC histogram under 200 training-label permutations for S1 and S2; observed AUC = 1.0.
- The null mass piles up near 1.0, so the null has low resolving power.
- SUPP.

## [67] phase4e fig4_gse294180_2x2_frozen_scores.png
- R0 log-odds, R2 log-odds and 4D MASH singscore per GSE294180 mouse by group.
- **MAIN/SUPP.** Must be paired with the R2 decomposition: the shift is cholesterol/PPARα-driven, not MASH-specific.

## [68] phase5 fig1_per_mouse_frozen_scores.png
- Grid of rows (HONG, GSE245069, GSE277772) × columns (FAST, MASH, R0, R2): per-mouse dots, group mean bars, dotted line = shared-control mean.
- Hong FAST: CON ≈ 0.10, PE ≈ 0.16, PS ≈ 0.125, PP ≈ 0.19, PB ≈ 0.14.
- **MAIN candidate** (external mouse-level).
- Caveats: the GSE245069 y-label is truncated, and GSE277772 is normalised input (exploratory). Say both in the caption.

## [69] phase5 fig2_gene_level_external.png
- Left: hexbin of Hong PE vs CON observed log2FC vs F8-OLS prediction. Title: ρ 0.22 [0.14, 0.28] (gene-block CI), weighted R² 0.041, permutation p = 0.25.
- Right: Spearman per polymer for 4 frozen models, with black ticks = 95th percentile of the mouse-relabelling null.
- **MAIN candidate.**
- Caveat: the error bars are gene-block CIs, which are anti-conservative. Phase 5B's mouse bootstrap gives −0.03 to 0.30. The caption must use the 5B interval or flag this.

## [70] phase5 fig3_r2_decomposition.png
- Contribution of each R2 feature to the log-odds shift for GSE294180 PE_in_CD / PE_in_MD, Hong PE and GSE245069 PS.
- **MAIN/SUPP.**

## [71] phase5 fig4_pathway_concordance.png
- Heatmap of median log2FC for 12 sets × (GSE294180 PE_in_CD/MD, REF fasting, REF MASH1 | Hong PE/PS/PP/PB, Li PS pooled/80 nm/0.5 µm/5 µm, GSE277772*).
- **MAIN candidate** (cross-study biology).
- Caveat: these are medians, not tests. Camera p-values are in `pathway_concordance.csv`.

## [72] phase5b fig1_gap_targets_labellings.png
1. Train / CV / Hong / GSE245069 Spearman for F8-OLS trained on CD vs MD.
2. Target log-density histograms (heavy tails).
3. Transfer ρ of the 20 relabellings vs each labelling's PC1 difference (orange = the true labelling).
- **MAIN candidate** (the decisive specificity diagnostic).

## [73] phase5b fig2_domain_shift_uncertainty.png
1. Frozen FAST singscore of controls and exposed mice across studies (baseline shift).
2. Transfer ρ across 15 mouse-resampled Hong targets (diamond = all mice).
3. Counts per read pair (Hong ≈ 2×).
- **MAIN/SUPP.**

## [74] phase5b fig3_before_after_interventions.png
- Left: train / CV / Hong (mouse-bootstrap bars) for B0, B1, B2, I2a, I2b, I2c, I3.
- Right: control-free ρ with bootstrap 95% CI and the true-labelling rank.
- **MAIN candidate.** Must be labelled post-external EXPLORATORY.
- Rank discrepancy: B0 is "rank 4/20" here (arith proxy), but `02_*.md` says the true labelling is "5th of 20" (log2CPM proxy, Phase 5 style). State both.

## [75] phase5b fig4_coefficients.png
- Coefficient per SD for B0, I2a, I2b and I2c.
- The CIRC_eff coefficient shrinks in I2b/I2c.
- SUPP.
