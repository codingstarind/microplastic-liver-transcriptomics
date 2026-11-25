# Evidence audit: pathway follow-up and R2 classifier ablation (exploratory, post hoc)

**Status:** an exploratory follow-up motivated by results already seen. It is not an independent confirmation.

**Scope:**
- Numbers only. No biological interpretation and no hypotheses (per instruction).
- No existing result file was modified.
- The Phase 4E classifier code was imported unchanged.

## 0. Verification of the quoted starting numbers
| Quoted | Stored value | Source |
|---|---|---|
| CD PPARα targets up, q 0.0406 | Up, 0.0406 (fry, PPAR family of 3 sets) | `research/enrichment_audit/tables/T2_*` |
| CD mTORC1 down, q 0.0796 | Down, 0.0796 (fry, Hallmark family of 50) | same |
| MD cholesterol homeostasis, DNA repair, glycolysis and mTORC1 down, q 0.0916 | all Down, 0.0916 | same |
| R2 +0.835 (CD), +1.663 (MD); Holm 0.0724 / 0.000623 | 0.8347 / 1.6630; 0.0724 / 0.000623 | `roadmap_phase4e/results/gse294180_frozen_model_contrasts.csv` |

**Reproduction of Phase 4E R2 L2-logistic** (`results/ablation_reproduction_check.csv`):
- S1/S2 AUC, balanced accuracy and Brier, and the three GSE294180 contrasts with their Holm p-values, all match the stored values.
- Maximum absolute difference: 1.3 × 10⁻¹⁵.

## 1. Pathway scores (19 mice): primary singscore, OLS score ~ diet × PE
- q8 = BH across the 8 scores within a contrast. This is a new family.
- The original fry families are unchanged and shown separately.
- Exact p = within-diet relabelling: 126 labellings for CD and 252 for MD (smallest attainable p is 1/126 and 1/252).

| Pathway | PE_in_CD Δ [95% CI] | q8 | exact p | PE_in_MD Δ [95% CI] | q8 | exact p | Interaction p (q8) | Original fry q, CD / MD |
|---|---|---|---|---|---|---|---|---|
| PPARα targets | +0.0229 [0.0022, 0.0436] | 0.064 | 0.087 | +0.0099 [−0.0096, 0.0294] | 0.339 | 0.206 | 0.345 (0.634) | 0.0406 / 0.518 |
| mTORC1 | −0.0098 [−0.0215, 0.0019] | 0.125 | 0.167 | −0.0093 [−0.0203, 0.0018] | 0.125 | 0.056 | 0.941 (0.941) | 0.0796 / 0.0916 |
| Cholesterol homeostasis | −0.0203 [−0.0444, 0.0037] | 0.125 | 0.095 | −0.0339 [−0.0566, −0.0112] | 0.016 | 0.024 | 0.396 (0.634) | 0.107 / 0.0916 |
| Cholesterol synthesis (core) | −0.0688 [−0.1105, −0.0272] | 0.012 | 0.024 | −0.1207 [−0.1600, −0.0814] | 0.0001 | 0.008 | 0.073 (0.551) | n/a |
| Glycolysis | −0.0043 [−0.0123, 0.0038] | 0.317 | 0.294 | −0.0030 [−0.0106, 0.0046] | 0.418 | 0.452 | 0.806 (0.922) | 0.623 / 0.0916 |
| DNA repair | −0.0020 [−0.0078, 0.0038] | 0.472 | 0.405 | −0.0078 [−0.0133, −0.0024] | 0.016 | 0.040 | 0.138 (0.551) | 0.724 / 0.0916 |
| Lipogenesis (core) | −0.0494 [−0.0645, −0.0342] | < 0.001 | 0.008 | −0.0375 [−0.0518, −0.0232] | 0.0002 | 0.008 | 0.244 (0.634) | n/a |
| FAST (Phase 4D, signed) | +0.0650 [0.0216, 0.1083] | 0.016 | 0.040 | +0.0555 [0.0146, 0.0963] | 0.018 | 0.024 | 0.738 (0.922) | n/a |

**Threshold counts on the q8 family:**
- **CD:** 3 of 8 scores have q8 < 0.05 (cholesterol synthesis, lipogenesis, FAST). PPARα targets is added at q8 < 0.10, making 4 of 8.
- **MD:** 5 of 8 have q8 < 0.05 (cholesterol homeostasis, cholesterol synthesis, DNA repair, lipogenesis, FAST). No further score is added at q8 < 0.10.
- **Interaction:** 0 of 8 at either threshold.

**Agreement with the secondary mean-z (VST) score:** same sign in all 24 rows. Values are in `TABLE_results_concise.csv`.

**Scores that met q8 < 0.05 with singscore but have mean-z p > 0.05:**
- cholesterol homeostasis in MD (p 0.104);
- FAST in MD (p 0.094).

**Glycolysis in MD:**
- fry q = 0.0916;
- singscore p = 0.418;
- mean-z p = 0.031.

## 2. Gene-set overlap and coordination
**Overlap within the 11,883-gene background** (`results/gene_set_overlap.csv`):

| Pair | Shared genes | Jaccard |
|---|---|---|
| Cholesterol homeostasis ∩ cholesterol synthesis | 17 (of 21 synthesis genes) | 0.233 |
| Cholesterol synthesis ∩ FAST-down | 11 | 0.183 |
| Cholesterol homeostasis ∩ FAST-down | 12 | 0.112 |
| Cholesterol homeostasis ∩ mTORC1 | 18 | 0.076 |
| Glycolysis ∩ mTORC1 | 23 | 0.072 |
| FAST-down ∩ lipogenesis | 4 | 0.069 |
| FAST-up ∩ PPARα targets | 3 | 0.049 |

**Unique-gene scores** (each set minus genes in any of the other sets, including FAST up/down; `results/pathway_contrasts_unique_genes.csv`):

| Score | Genes left | Contrast | Full set | Unique genes |
|---|---|---|---|---|
| Cholesterol synthesis | 1 | – | – | – |
| Cholesterol homeostasis | 38 | PE_in_MD | −0.0339, p 0.006 | +0.0021, p 0.83 |
| Lipogenesis | 5 | PE_in_CD | – | −0.0823, p < 0.0001 |
| Lipogenesis | 5 | PE_in_MD | – | −0.0321, p 0.009 |
| PPARα targets | 11 | PE_in_CD | – | +0.0272, p 0.037 |
| mTORC1 | 133 | PE_in_CD | – | −0.0059, p 0.23 |
| mTORC1 | 133 | PE_in_MD | – | −0.0021, p 0.65 |
| DNA repair | 137 | PE_in_MD | – | −0.0083, p 0.006 |
| Glycolysis | 132 | – | p ≥ 0.71 (both diets) | – |

**Score correlations** (Spearman; BH within basis; `results/pathway_score_correlations.csv`):

| Pair | All 19 mice | Within-group residuals |
|---|---|---|
| FAST vs PPARα | +0.911 | +0.723 (q 0.003) |
| FAST vs cholesterol homeostasis | −0.884 | −0.156 (q 0.65) |
| FAST vs cholesterol synthesis | −0.925 | −0.361 (q 0.24) |
| FAST vs lipogenesis | −0.877 | −0.381 (q 0.22) |
| FAST vs mTORC1 | −0.744 | −0.412 (q 0.22) |
| FAST vs glycolysis | +0.196 | +0.037 |
| FAST vs DNA repair | −0.551 | +0.025 |
| mTORC1 vs cholesterol homeostasis | +0.837 | +0.798 (q 0.001) |
| mTORC1 vs glycolysis | +0.339 | +0.632 (q 0.015) |
| Cholesterol homeostasis vs glycolysis | +0.032 | +0.728 (q 0.003) |

**Per-PE-mouse consistency:** number of PE mice on the same side as the group effect, of 5 per diet (`results/per_PE_mouse_deviation_summary.csv`).
- **All 5 of 5 in both diets:** lipogenesis (all 10 below the vehicle minimum), cholesterol synthesis, FAST.
- **PPARα targets:** CD 5/5, MD 4/5.
- **Cholesterol homeostasis and mTORC1:** 4/5 in each diet.
- **Glycolysis:** 3/5 in each diet.
- **DNA repair:** CD 3/5, MD 5/5.

**Leave-one-mouse-out** (19 refits): sign kept in 19/19 for every PE_in_CD and PE_in_MD contrast. The maximum LOO p-values are:

| Score | PE_in_CD | PE_in_MD |
|---|---|---|
| Cholesterol synthesis | 0.013 | < 0.001 |
| Lipogenesis | < 0.001 | < 0.001 |
| FAST | 0.020 | 0.035 |
| Cholesterol homeostasis | – | 0.020 |
| DNA repair | – | 0.025 |
| PPARα targets | 0.092 | – |

**Cook's D > 4/n** (`results/cooks_distance.csv`):

| Mouse | Scores flagged |
|---|---|
| CD_MP_2_S49 | mTORC1, cholesterol homeostasis, lipogenesis, FAST |
| CD_MP_3_S50 | PPARα, FAST |
| MD_MP_4_S56 | PPARα, cholesterol homeostasis |
| CD_3_S41 | lipogenesis |
| CD_MP_1_S48 | mTORC1 |
| MD_MP_2_S54 | glycolysis |

## 3. R2 classifier ablation
**Setup:**
- The Phase 4E pipeline is unchanged: same tuning, the S1/S2 study-held-out splits, and a final model trained on both MASH studies.
- GSE294180 is never used for fitting.
- Files: `results/ablation_*.csv` and `figures/fig2_classifier_ablation.*`.
- **Scale:** each PE contrast is divided by that model's own training MASH-minus-control log-odds separation ("scaled effect"; mouse bootstrap with B = 2,000).
- Glycolysis and DNA repair are **not** features of the R2 classifier, so they cannot be ablated.

**Validity:**
- Every ablated model kept AUC = 1.00 in S1 and S2.
- Balanced accuracy in S2 fell below 1.00 for:
  - dropping cholesterol homeostasis or cholesterol synthesis alone (0.938);
  - dropping both cholesterol features (0.812);
  - dropping cholesterol + lipogenesis (0.750);
  - dropping the 5-feature metabolic group (0.750);
  - dropping all lipid features (0.750);
  - keeping only R0 + inflammatory + EMT (0.812).
- S2 Brier is 0.028 for the full model and up to 0.159 for ablations.

| Model | CD log-odds (Holm p) | CD scaled [95% CI] | MD log-odds (Holm p) | MD scaled [95% CI] | MD scaled minus full [95% CI] |
|---|---|---|---|---|---|
| Full R2 | +0.83 (0.072) | 0.070 [0.017, 0.117] | +1.66 (0.0006) | 0.140 [0.094, 0.192] | – |
| − cholesterol synthesis | +0.60 (0.38) | 0.052 | +1.04 (0.074) | 0.090 [0.031, 0.158] | −0.050 [−0.071, −0.028] |
| − cholesterol homeostasis | +0.44 (0.99) | 0.038 | +0.95 (0.38) | 0.082 [−0.002, 0.190] | −0.058 [−0.105, 0.004] |
| − both cholesterol features | +0.13 (1.0) | 0.012 | +0.01 (1.0) | 0.001 [−0.120, 0.153] | −0.140 [−0.229, −0.031] |
| − mTORC1 | +0.59 (0.26) | 0.049 | +1.50 (0.0045) | 0.125 | −0.016 [−0.032, −0.001] |
| − cholesterol (both) + mTORC1 | −0.46 (1.0) | −0.040 | −0.59 (1.0) | −0.051 | −0.191 [−0.326, −0.029] |
| − PPARα targets | +0.21 (0.56) | 0.018 | +1.35 (0.003) | 0.115 | −0.025 [−0.067, 0.006] |
| − PPARα group (PPARα + FA metabolism + peroxisome) | +0.03 (0.95) | 0.003 | +1.65 (0.004) | 0.141 | 0.000 [−0.057, 0.045] |
| − lipogenesis | +2.88 (< 0.001) | 0.240 | +3.21 (< 0.001) | 0.268 | +0.127 [0.092, 0.161] |
| − cholesterol (both) + lipogenesis | +2.53 (0.037) | 0.224 | +1.68 (0.13) | 0.149 | +0.009 [−0.087, 0.140] |
| − 5 metabolic (cholesterol both + lipogenesis + PPARα + mTORC1) | +0.39 (1.0) | 0.033 | −0.56 (1.0) | −0.048 | −0.188 [−0.337, −0.014] |
| − all lipid features | +1.19 (0.22) | 0.113 | +1.45 (0.14) | 0.138 | −0.003 [−0.091, 0.091] |
| − all lipid features + mTORC1 | +0.05 (1.0) | 0.004 | +0.44 (1.0) | 0.040 | −0.100 [−0.232, 0.048] |
| − R0 (the MASH programme score) | +0.97 (0.035) | 0.084 | +1.73 (0.0004) | 0.151 | +0.010 [−0.003, 0.022] |
| Keep only R0 + inflammatory + EMT | +0.53 (1.0) | 0.046 | +0.62 (1.0) | 0.054 | −0.086 [−0.240, 0.070] |

Single removals of fatty-acid metabolism, peroxisome, bile acid, inflammatory, EMT, xenobiotic or OXPHOS features changed MD scaled effects by −0.002 to +0.018. All of these models kept MD Holm p < 0.001.

The full model's interaction contrast is +0.83, Holm p 0.118. No ablated model has an interaction with Holm p < 0.05 except "− PPARα group", at 0.039.

## 4. Sensitivity checks
**Covariate adjustment** (depth, quality, muscle, clock phase, FAST; one at a time) and **muscle-gene exclusion**: `results/covariate_and_muscle_sensitivity.csv` and the `sensitivity_*` columns of `TABLE_results_concise.csv`.
- The sign is kept in all 6 variants for every PE_in_CD and PE_in_MD contrast, except:
  - cholesterol homeostasis in CD (5/6);
  - glycolysis in both diets (5/6);
  - DNA repair in CD (5/6);
  - PPARα in MD (5/6).
- FAST has 4 applicable variants (no FAST covariate and no muscle-gene exclusion): 4/4.
- Maximum p across the variants:
  - lipogenesis: < 0.003 (CD and MD);
  - cholesterol synthesis in MD: < 0.001;
  - cholesterol synthesis in CD: 0.107;
  - FAST: 0.138 (CD) and 0.117 (MD);
  - cholesterol homeostasis in MD: 0.154;
  - PPARα in CD: 0.775.
- Classifier ablations adjusted for the same covariates: `results/ablation_covariate_adjusted.csv`.

**Independent plastic studies** (existing Phase 5 data, same scoring; `results/external_studies_pathway_contrasts.csv`):
- **Hong PE vs control (3 vs 3; 20 labellings, smallest attainable p = 0.10):**
  - the direction matches GSE294180 PE_in_CD for 8 of 8 scores;
  - exact p = 0.10 for mTORC1, cholesterol homeostasis, glycolysis, DNA repair and FAST;
  - PP, PB and PS show the same signs for most scores.
- **GSE245069 pooled PS vs control (12 vs 4; 1,820 labellings):**

| Score | Exact p |
|---|---|
| Glycolysis | 0.0011 |
| mTORC1 | 0.020 |
| PPARα | 0.039 |
| Cholesterol homeostasis | 0.27 |
| FAST | 0.24 |

  Cholesterol synthesis and lipogenesis have the opposite sign to GSE294180.

## 5. Statistical status summary
**Criteria met** (q8 < 0.05, exact p < 0.05, LOO sign 19/19):
- **Both diets:** cholesterol synthesis, lipogenesis, FAST.
- **MD only:** cholesterol homeostasis, DNA repair.

**Met only at q8 < 0.10:**
- PPARα targets in CD (q8 0.064; exact p 0.087).

**Did not meet q8 < 0.10:**
- mTORC1 in both diets;
- glycolysis in both diets;
- PPARα targets in MD;
- every interaction contrast.

**Classifier (MD):**
- The R2 PE_in_MD contrast fell to Holm p 1.0 (scaled 0.001) when both cholesterol features were removed.
- It also fell to Holm p 1.0 when the 5-feature metabolic group was removed.
- It stayed at Holm p < 0.01 when mTORC1, PPARα targets, the PPARα group or R0 was removed singly.
- The cholesterol-removed models had lower S2 balanced accuracy (0.81–0.94).

**Hypotheses and biological interpretation:** not provided, per instruction.

## 6. Negative and null results (reported as found)
- Interaction q8 ≥ 0.55 for all 8 scores. Classifier interaction Holm p 0.118.
- Glycolysis: no PE contrast with p < 0.27 on singscore, despite fry q 0.0916 in MD.
- mTORC1: exact p 0.17 (CD) and 0.056 (MD); unique-gene score p ≥ 0.23.
- The cholesterol homeostasis effect in MD is not present in its 38 unique genes (p 0.83).
- Cholesterol synthesis has 1 gene not shared with other sets.
- Hong (3 vs 3) cannot reach p < 0.10 by exact permutation.

## Reproduce
```powershell
.venv\Scripts\python.exe research\pathway_followup\scripts\pf01_pathways.py
```
```powershell
.venv\Scripts\python.exe research\pathway_followup\scripts\pf02_classifier_ablation.py
```
```powershell
.venv\Scripts\python.exe research\pathway_followup\scripts\pf03_report.py
```
