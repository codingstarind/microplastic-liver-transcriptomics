# Phase 5 locked evaluation plan (written before any external effect or score was computed)

**Locked:** 2026-11-23. SHA-256 of this file and every frozen input are in `results/plan_lock.json`.

**Already done before locking (no effect estimates or scores computed):**
- read GEO SOFT design metadata;
- mapped identifiers to vM35;
- computed per-sample library size, detected genes and muscle/sex marker CPM (`s01`).

## 0. Questions
- **Q1 (gene level):** do frozen gene-level models trained on GSE294180 predict an independent study's plastic-vs-control log2FC vector?
- **Q2 (mouse level):** do the frozen Phase 4D signatures and Phase 4E classifiers move in the same direction in independent plastic-exposed mice as in GSE294180? Is any shift cholesterol/fasting-like or broadly MASH/disease-like?

Neither question is "can we get high CV accuracy".

## 1. Study units (biological independence)

| Unit | Mice | Role | Contrast(s) |
|---|---|---|---|
| **HONG**: GSE261893 + GSE261725 + GSE286116 + GSE262399 | 15 (3 CON shared by all 4 series + 3 PE + 3 PS + 3 PP + 3 PB) | **PRIMARY PE test** | **PE vs CON (primary, the only PE contrast).** PS/PP/PB vs CON are a within-lab specificity check (exploratory) |
| **GSE245069** | 16 (4 water + 4 × 3 PS sizes) | **SECONDARY, cross-polymer (PS)** | **pooled PS (12) vs water (4)** = the one study-level contrast. The 3 sizes vs water are exploratory and correlated through the shared control |
| **GSE277772** | 6 (3 v 3) | **EXPLORATORY only** (no counts on GEO; normalised values) | PS vs CON: rank-based scores only. A gene-level limma on log values is descriptive. **Not a validation** |

- **Not independent:** the 4 Hong series (shared controls; one lab), and the 3 GSE245069 sizes (shared controls).
- Under the null with equal group sizes, two contrasts sharing a control have an expected correlation of **0.5** for gene-level effects. Cross-polymer agreement within one lab is judged against that, not against 0.

## 2. Frozen artefacts (no retuning, reselection or reorientation on external data)

| Artefact | Source | Frozen detail |
|---|---|---|
| F8-OLS (primary gene model) | Phase 4 / 4C | Features MASH/FAST/PPARA/CIRC eff + se. Target-weighted OLS, w = 1/(s² + τ²_DL), fitted on **all** GSE294180 universe genes. Two fits: target PE_in_CD (**primary**, control-diet match) and PE_in_MD (secondary) |
| F8+C-Ridge (secondary gene model) | Phase 4C | Same, plus block C. α = **1000** (the Phase 4C inner-CV choice in 15/15 CD folds and 14/15 MD folds; fixed, not re-tuned) |
| Phase 4D signatures | `roadmap_phase4d/signatures/*_main_n50.csv` (hash-checked) | **FAST_main_n50 (primary mouse-level score)**, MASH_main_n50, PPARA_main_n50 (secondary). singscore on TPM-proxy (GENCODE vM35 median transcript length) |
| Phase 4E R0 / R2 L2-logistic | Phase 4E `e03` (both MASH studies, 31 mice, C = 10, seed 20261122) | Refitted deterministically with the frozen C = 10 and serialised. Must reproduce the saved GSE294180 log-odds (max abs diff < 1e-6) before external use |
| Gene universe | Phase 2 main universe (11,888) | Gene-level: universe ∩ genes testable in the external study. 4E models: their 11,883-gene rank background ∩ external genes. 4D: the 13,777-gene GSE294180 background ∩ external genes |
| Pathway panel | Phase 4E R1 12 sets (Phase 3 `gene_sets.csv`) | Primary pathway panel; the 50 Hallmarks are exploratory |

## 3. External effect estimation (targets for Q1)
- **Count data (HONG, GSE245069):** DESeq2, Wald test, **MLE log2FC + lfcSE** (the same quantity as the GSE294180 training target), independent filtering off.
  - Prefilter: ≥ 10 counts in ≥ n_min samples, where n_min is the smallest group size (3 for HONG, 4 for GSE245069).
  - HONG: one model `~ group` with 5 groups (shared controls modelled once).
  - GSE245069: `~ group` with 4 groups; the pooled contrast weights the 3 PS groups 1/3 each.
  - GSE245069 estimated counts are rounded **only** inside DESeq2 (standard for estimated counts; these are not normalised values).
- **Sensitivities:**
  - HONG PE vs CON on its 6 mice only;
  - GSE245069 limma-voom on unrounded counts.
- **GSE277772:** limma-trend on log2(value + 1). Descriptive only.

## 4. Q1 metrics (gene level), per study contrast
- **Genes:** frozen feature table ∩ external tested genes (with finite lfcSE).
- **Primary metric:** precision-weighted R² against **(a) the zero-effect baseline** and **(b) the training-mean baseline** (a constant equal to the weighted mean GSE294180 target).
  - Weights: w = 1/(s²_ext + τ²_ext), with τ² estimated by DerSimonian–Laird from the external target alone.
- **Also reported:**
  - R² vs the external study's own mean (an oracle constant; the hardest baseline; the Phase 4C/4D readiness rule);
  - unweighted R²;
  - MAE, RMSE;
  - Spearman ρ with a module-block bootstrap 95% CI (20 Phase 4 co-expression modules, 2,000 resamples);
  - sign agreement (all genes; genes with |pred| in the top decile);
  - calibration slope (descriptive only; **never** used to rescale).
- **Comparators:**
  - the zero model;
  - single-axis frozen OLS (FAST_eff only; MASH_eff only), fitted on GSE294180 like F8;
  - the observed GSE294180 PE log2FC (ashr) used directly as the predictor (does the original PE signature replicate?).
- **Label-permutation null (accounts for gene–gene correlation):**
  - recompute a fast effect proxy (difference of mean log2CPM over prefiltered genes) for **every** relabelling of the external contrast's mice;
  - compare the observed Spearman(pred, proxy) with that distribution;
  - HONG PE vs CON: C(6,3) = 20 labellings, so the smallest attainable one-sided p is 0.05;
  - GSE245069 pooled: C(16,4) = 1,820 labellings.
- **Noise ceiling:** signal share = 1 − mean(SE²)/var(y) of the external target. If it is ≤ 0.2, R² is declared uninterpretable as a model failure (the target is mostly noise) and reported as such.
- **Pass / fail (fixed in the Phase 4C/4D readiness docs, made precise here):**
  - **FAIL** if weighted R² vs the study's own mean is ≤ 0, **or** the Spearman bootstrap CI includes 0.
  - **PASS-direction** if Spearman > 0 with CI excluding 0 and permutation p ≤ 0.10.
  - R² vs zero/training mean > 0 is reported, but cannot rescue a FAIL.

## 5. Q2 scoring and tests (mouse level)
- **4D scores:**
  - TPM-proxy (counts ÷ vM35 median transcript length → TPM);
  - within-sample ranks over the frozen background ∩ external genes;
  - signature genes absent from the background are dropped and counted;
  - GSE245069 sensitivity: deposited TPM;
  - GSE277772: normalised values ranked as deposited.
- **4E scores:**
  - log-CPM ranks over the 4E background ∩ external genes;
  - R0 = frozen both-study signature; R1 = 12 fixed pathways; frozen scaler + logistic; output = log-odds (decision function), not a calibrated probability.
- **Per contrast:**
  - difference in group means with a Welch 95% CI;
  - Hedges g;
  - exact two-sided permutation p over all relabellings (3 v 3 → 20 labellings, min p = 0.10; 4 v 4 → 70; 12 v 4 → 1,820);
  - leave-one-mouse-out signs;
  - Spearman of score vs counted reads / detected genes / muscle CPM within study.
- **Primary Q2 test:** HONG PE − CON on **FAST_main_n50**. The GSE294180 reference is PE_in_CD = +0.065.
  - **Directionally consistent:** same sign **and** all leave-one-out signs agree **and** the exact permutation is at its minimum attainable p (0.10).
  - **FAIL:** opposite sign, or |t_Welch| < 1 (an "estimate centred on 0").
  - Anything else is inconclusive.
  - p < 0.05 is **not attainable** with 3 v 3 and is not promised.
- **Secondary:** MASH_main_n50, R0-logistic, R2-logistic (Holm across these 3 for the Welch p, reported but not decisive).
- **R2:** every external R2 shift is decomposed into per-feature contributions (coef/SD × group-mean feature difference). This classifies the shift as "cholesterol/fasting-like" (> 50% of the |sum| from the CHOL_SYNTHESIS_CORE, HALLMARK_CHOLESTEROL_HOMEOSTASIS, LIPOGENESIS_CORE, SANDERSON_PPARA_TARGETS and MTORC1 features) or "broad MASH/disease-like" (the R0 + INFLAMMATORY + EMT share > 50%).
- **Coverage control:** re-score the 19 GSE294180 mice with each external study's reduced gene background and report the change in the GSE294180 PE_in_CD estimate. A denominator effect must be negligible (< 10% of the estimate).

## 6. Pathways (Q2 context)
- **Method:** limma camera (accounts for inter-gene correlation) and the mean gene log2FC on external log-CPM, same designs as §3.
- **Concordance:** direction and magnitude vs the existing GSE294180 PE_in_CD / PE_in_MD fry/camera results and the Phase 3 reference axes (FAST, MASH), on the 12-set primary panel.
- **Coverage:** genes per set in each dataset.
- **No set is chosen after seeing concordance.**

## 7. What is never done
- No refit, reselection or re-orientation on external plastic labels.
- No PE/PS classifier.
- No pooling of polymers into "replicates".
- No single leaderboard mixing gene R², AUC and score p-values.
- No causal or "PE induces MASH" language.

## 8. Phase 6 go / no-go rule (fixed now)
- **GO (narrow):** only if the primary Q2 test is directionally consistent **and** the R2 decomposition is reported. The go is for a write-up of a reproducible, fasting/cholesterol-type PE-associated shift with external directional support.
- **GO for gene-level claims:** only if Q1 passes in HONG.
- **Otherwise NO-GO** for any external-validation claim: Phase 6 becomes a reporting phase built on GSE294180-internal findings, with the external failure stated.
