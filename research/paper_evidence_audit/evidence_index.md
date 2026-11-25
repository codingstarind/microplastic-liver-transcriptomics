# Evidence index: claim → number → file → script → figure → citation

Paths are relative to `research/` unless they start with `reports/`, `figures/`, `analysis/` or `results/` (repository root). `R####` ids refer to rows of `results_master.csv`.

| # | Claim (evidence wording, not final prose) | Number(s) | Result file | Script | Figure | Memo | Cite |
|---|---|---|---|---|---|---|---|
| 1 | PE libraries are shallower than vehicle libraries (complete separation) | 13.5–16.5 vs 17.0–19.1 M read pairs | `reports/final_report_draft.md` §4.1; `results/deseq2/covariates.csv` | `scripts/08–09` | F001, F006 | P1 inventory §6 | leek2010 |
| 2 | Muscle carry-over defines PC1; diet separates once muscle genes are removed | MD_1 5,190 CPM; PC1 30% → diet PC1 20.8%; ARI 1.0 | `reports/eda_report.md` | `scripts/04–07` | F009–F011, F013 | – | – |
| 3 | Factorial DE gene counts | 858 / 472 / 2,262 / 87 | `results/deseq2/provisional_deposited/` | `analysis/deseq2/run_deseq2_factorial.R` | F022 | final report §4.3 | love2014 |
| 4 | Interaction genes are mostly a muscle artefact | 87 → 16 (muscle covariate); 1 in edgeR/limma | `reports/phase2/interaction_report.md` | `analysis/interaction` | F018 | – | – |
| 5 | Cholesterol programme lower with PE in both diets; no pathway interaction | fry p 0.004 / 0.005; interaction p 0.84 | `reports/phase2/pathway_report.md` | `analysis/pathways` | F016 | – | wu2010roast |
| 6 | Same-data PE–diet resemblance is modest and CD-only | ρ 0.33 [0.29–0.40] vs 0.05 | `reports/phase2/pe_mash_similarity_report.md` | `analysis/concordance` | F017 | – | – |
| 7 | Ppara induction is not reproduced on the same counts | +0.16–0.17 vs Jung +0.53 | `roadmap_phase3/results/ppara_discrepancy/ppara_estimates_by_approach.csv` | `roadmap_phase3/scripts/e01_ppara_discrepancy.R` | – | P3 decision §8 | jung2026 |
| 8 | Independent MASH references agree with each other | ρ 0.557 | `roadmap_phase3/results/comparisons/reference_spearman_matrix.csv` | `f01_compare_and_baselines.py` | F033 | P3 decision §3 | gse316023paper, gse246328paper |
| 9 | MASH positive control fails the locked criterion | ρ 0.355; exact p 0.036; q99 0.390 | `roadmap_phase3/results/positive_control/`, `results/posthoc_audit/exact_permutation.json` | `c01_positive_control.py`, `g01_exact_permutation_audit.py` | F030 | P3 addendum §1–2 | – |
| 10 | Fasting and PPARα agonism are opposite for lipogenesis | −2.79 vs +1.15 | `roadmap_phase3/results/comparisons/core_pathway_median_log2FC.csv` | `f01` | F032 | P3 reference comparison | kinouchi2018, heintz2025 |
| 11 | PE–MASH correlation is fragile; PE–FAST is stable | 0.21 [LOO 0.03–0.24]; 0.25 [0.22–0.27] | `roadmap_phase3/results/comparisons/baselines.csv` | `f01` | F025 | P3 baseline report | – |
| 12 | DEG overlap is "significant" despite discordant directions | 84 / 125 concordant / discordant | same | `f01` | F025 | – | – |
| 13 | About 7% of held-out PE gene-effect variance is predictable | R²_oos 0.066 / 0.070; Spearman 0.264 / 0.309 | `roadmap_phase4/results/taskA/summary_by_run_model.csv` (`tables/T01`) | `p4_run.py`, `p4_analyze.py` | F034 | P4 decision §2 | varma2006, roberts2017 |
| 14 | Trees are not better than linear on PE targets (better on diet) | CatBoost vs ols_all CD −0.0035 [−0.0093, 0.0023]; diet −0.0215 [−0.0369, −0.0061] | `paired_comparisons.csv` (`T02`) | `p4_analyze.py` | F034 | – | nadeau2003 |
| 15 | Robust to variants and LOO | 8/8 variants; 19/19 LOO | `robustness_counts.csv` | `p4_analyze.py` | F035, F036 | – | – |
| 16 | FULL features raise linear R² but tie CatBoost | 0.066 → 0.104 (CD); vs CatBoost-F8 −0.0031 [−0.0111, 0.0049] | `roadmap_phase4b/results/metrics_main.csv`, `extra_vs_tree_baselines.csv` (`T03`, `T05`) | `p4b_run.py`, `p4b_analyze.py` | F039, F041 | P4B memo | – |
| 17 | Block C is useful; A is not; B is target-informed | ridge MD ablation +0.0033 [0.0001, 0.0065] | `roadmap_phase4b/feature_ablation_results.csv` (`T05b`), `control_shuffled_B.csv` | `p4b_analyze.py` | F040 | P4B memo | sweeney2017, schubert2018 |
| 18 | FULL coefficients are unidentifiable | VIF 233; condition number 37 vs 2.5; up to 12/28 sign flips | `roadmap_phase4c_diagnostics/tables/collinearity_summary.csv` (`T06`) | `c01_collinearity.py` | F042 | P4C memo | – |
| 19 | Target noise limits achievable R² | signal share 0.36 / 0.39 | `target_quality.csv` (`T08`) | `c02_errors_target_leakage.py` | F043, F044 | P4C memo | – |
| 20 | PE shifts a fasting-like signature per mouse in both diets | +0.065 [0.022, 0.108]; +0.056 [0.015, 0.096] | `roadmap_phase4d/results/contrast_estimates_primary.csv` (`T10`) | `roadmap_phase4d/scripts` | F045, F047 | P4D decision | foroutan2018 |
| 21 | The MASH-like PE shift depends on cholesterol genes | +0.003 / +0.021 without CHOL genes | `contrast_estimates_sensitivity.csv` | same | F048 | P4D decision | – |
| 22 | MASH-state classifiers: AUC at ceiling; threshold transfer improved | S2 balanced accuracy 0.81 → 1.00; Brier 0.19 → 0.028 | `roadmap_phase4e/results/external_metrics.csv` (`T11`); `roadmap_phase5b/tables/mash_classifier_bootstrap.csv` (`T17`) | `e01_mash_state_models.py` | F049 (retitle), F050 | P4E memo | bernau2014 |
| 23 | The classifier's PE shift is cholesterol / PPARα-driven | PE_in_MD +1.66 (Holm 0.0006); MASH feature +0.18 | `gse294180_frozen_model_contrasts.csv`, `r2_feature_contributions_gse294180.csv` (`T12`, `T13`) | `e03_apply_to_gse294180.py` | F052, F055 | P4E 04 | – |
| 24 | External gene-level transfer not passed | Hong ρ 0.22; perm p 0.25; GSE245069 −0.07 | `roadmap_phase5/gene_level_external_results.csv` (`T14`) | `s04_gene_level.py` | F054 | P5 decision | – |
| 25 | Hong FAST shift: same direction and size as GSE294180 | +0.060 vs +0.065; exact p 0.10; LOO 6/6 | `sample_level_external_results.csv` (`T15`) | `s05_sample_level.py` | F053 | P5 decision | – |
| 26 | …but not PE-specific | PP +0.088, PB +0.038, PS +0.024; non-PE pooled exact p 0.04 | same | same | F053 | P5 decision | – |
| 27 | Pathway directions across studies | Hong PE lipogenesis −0.82 (camera p 0.006); PPARα targets inconsistent | `roadmap_phase5/pathway_concordance.csv` | `s06_pathways_diagnostics.py` | F056 | – | wu2012camera |
| 28 | No pipeline error; count-unit mismatch with negligible effect | 1.641 vs 0.775; ρ 0.219 → 0.220 | `roadmap_phase5b/tables/equivalence_checks.csv`, `i1_unit_harmonisation_effect.csv` | `b01_pipeline_audit.py`, `b04_intervention_targets.R` | F058 | P5B 01 | – |
| 29 | External agreement tracks the mouse-state axis | R² 0.889; true labelling 5th (4th) of 20 | `gap_hong_20_labellings.csv` | `b03_gap_diagnostics.py` | F057 | P5B 02 | – |
| 30 | Mouse-level uncertainty includes 0 | −0.03 to 0.30 | `bootstrap_summary.csv`, `intervention_verdicts.csv` (`T16`) | `b05_interventions.py` | F059 | P5B 05 | – |
| 31 | No intervention met the locked criteria (post-external) | I2b +0.074, P 0.87; fails A2/A3 | `intervention_verdicts.csv` | `b05_interventions.py` | F059, F060 | P5B 05 | huber1964, dersimonian1986 |
| 32 | No untouched validation exists; a design is proposed | ≥ 8 v 8; own controls; fixed ZT | – | – | – | `roadmap_phase5b/phase6_readiness_memo.md` | steyerberg2016, riley2016, debray2015 |
