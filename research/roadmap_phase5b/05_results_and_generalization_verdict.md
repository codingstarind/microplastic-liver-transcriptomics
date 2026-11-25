# 05: Results and generalisation verdict (Phase 5B)

**Labels used throughout:**
- **INTERNAL** = GSE294180 held-out genes (nested grouped CV).
- **EXPLORATORY-EXTERNAL** = Hong / GSE245069. These outcomes were already seen in Phase 5, so they are post-external re-development, **not** validation.
- No untouched study exists.

## 1. Before/after pipeline transformations

| Transformation | Before | After | Test |
|---|---|---|---|
| I1 Hong count unit (× 0.5, a priori) | 1.64 counts per read pair; 19,491 genes pass the prefilter; 1,133 DEGs (padj < 0.05) | 0.82 per read pair; 17,311 genes; 951 DEGs; log2FC r = 0.96; matched-gene median SE ratio 1.006 | `test_count_unit_mismatch_detected_and_harmonised` |
| Effect on frozen F8-OLS transfer (Hong PE) | ρ 0.219; weighted R² vs own mean 0.041 | ρ 0.220 / 0.221 (factor 0.472); R² 0.042 | `tables/external_candidates.csv` |

**Verdict:** technical mismatch verified and fixed, with **negligible effect on transfer**. The Phase 5 result was not a pipeline artefact.

## 2. Before/after models on the same gene-effect task (`tables/intervention_verdicts.csv`, `cv_fold_metrics.csv`, fig3)

**Common outcome:** Spearman vs PE log2FC. INTERNAL = outer grouped CV, 15 folds on the Phase 4 module splits, with the PE_in_CD MLE target for every candidate. Hong = exploratory, with a paired mouse bootstrap (B = 1,000) that resamples both GSE294180 and Hong mice within arm.

| Candidate | Train ρ | CV ρ (SD) | CV R² | Hong ρ | Δ vs B0 | Mouse-bootstrap Hong ρ 95% | P(Δ > 0) | True labelling rank /20 | Control-free ρ (P > 0) | Locked verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| **B0 F8-OLS (frozen)** | 0.274 | 0.264 (0.040) | 0.066 | 0.219 | — | **−0.03 to 0.30** | — | 4 | 0.002 (0.67) | comparator |
| B1 F8+C-Ridge (frozen) | 0.281 | 0.270 (0.040) | 0.081 | 0.142 | −0.076 | −0.06 to 0.29 | 0.14 | 5 | −0.057 (0.53) | comparator |
| B2 HGB-F8 | **0.378** | 0.276 (0.046) | 0.079 | 0.206 | −0.012 | −0.04 to 0.30 | 0.37 | 4 | 0.016 (0.74) | comparator: **variance gap** (train ≫ CV) |
| I2a robust estimator (CD) | 0.255 | 0.246 (0.044) | 0.063 | 0.217 | −0.002 | −0.02 to 0.29 | 0.56 | 3 | 0.048 (0.72) | no external improvement |
| I2b diet-averaged PE (MLE, 19 mice) | 0.242 | 0.234 (0.038) | 0.052 | **0.293** | **+0.074** | 0.01 to 0.35 | 0.87 | 4 | 0.078 (0.86) | fails A2 (0.87 < 0.90) and A3 (internal −0.030) |
| I2c robust + diet-averaged | 0.224 | 0.217 (0.036) | 0.048 | 0.281 | +0.062 | 0.02 to 0.35 | 0.88 | 4 | 0.125 (0.90) | fails A2 (0.88) and A3 (−0.047) |
| I3 Huber on normal scores | 0.276 | 0.265 (0.037) | n/a (different scale) | 0.220 | +0.001 | −0.05 to 0.30 | 0.21 | 4 | −0.010 (0.56) | no external improvement |

**Hyperparameters** (`tables/hyperparameter_choices.csv`; inner GroupKFold(4) on training genes):
- HGB: learning rate 0.03, min_samples_leaf 50 (14/15 folds).
- Ridge α = 1000 (I2a 12/15; I2b and I2c 15/15). Ridge again sits at the grid maximum, equivalent to near-OLS predictions with condition number ≤ 4. There is nothing for regularisation to fix.
- Huber ε = 1.35, α = 1e-4 (6/15; the choice is unstable, which means flat performance).

**Secondary results (exploratory):**
- GSE245069 pooled PS ρ: B0 −0.07; I2b +0.07; I2c +0.09.
- Weighted R² vs own mean on Hong PE: B0 0.041; I2b 0.068; I2c 0.053.

**Own-target CV:** I2b reaches CV ρ 0.31 on its own diet-averaged target, against 0.26 for B0 on PE_in_CD. The 19-mouse average is a more predictable (more reliable) target.

**Interpretation (fig4):** I2b/I2c's gain comes from **dropping the circadian feature weight** (CIRC_eff −0.042 → about −0.004) while keeping FAST and MASH. This fits the known **CD+Veh clock-state anomaly** contaminating the CD-only target: averaging over diets dilutes it. It is a plausible, mechanistically coherent improvement. But it was suggested by a Phase 5 external result (the MD model transferred better), so it is **post-hoc**.

### Was any improvement stable, reproducible and large?
- **No candidate met all locked criteria** (A1–A3; A4 for specificity).
- The diet-averaged models give the largest and most consistent exploratory gain (+0.06 to +0.07; the mouse-bootstrap interval excludes 0 only for them), but:
  - the paired bootstrap gives P(Δ > 0) of 0.87–0.88;
  - internal performance on the original target falls;
  - the true PE labelling still ranks 4th of 20.
- **Verdict:** "inconclusive external power", plus "internal trade-off". **Not** a true transfer improvement.

## 3. Resampling: what it showed and what it cannot show
- **The mouse bootstrap is the key BT3041 correction.** The frozen F8-OLS Hong ρ has a 95% interval of **−0.03 to 0.30** when mice are resampled in both studies. The module-block gene bootstrap from Phase 5 gave 0.14–0.28, a falsely precise interval, because 11,000 genes share 6 mice.
- **Leave-one-mouse-out and 2-vs-2 DESeq2 refits:** transfer ρ 0.10–0.31.
- **Exhaustive relabelling:** 20 labellings, minimum one-sided p = 0.05. No candidate reaches rank 1. With 3 v 3, permutation **cannot** remove the shared-control or husbandry confounding.
- **Bootstrap caveat:** with 3 mice per arm, within-arm resampling duplicates mice and shifts the medians (for example, the B0 control-free ρ is 0.002 observed vs a bootstrap median of 0.049). Intervals are descriptive.

## 4. Mouse-level MASH classification (Phase 4E saved predictions; `tables/mash_classifier_bootstrap.csv`)

| Split | Model | AUC | Bootstrap AUC = 1 | Balanced accuracy [95%] | Probabilistic metric [95%] |
|---|---|---|---|---|---|
| S1 (train 316 → test 246, n = 15) | R0 rule | 1.00 | 100% | 1.00 [1.00, 1.00] | misclassification rate 0.00 (hard rule, **not Brier**) |
| S1 | R2 L2-logistic (C = 10) | 1.00 | 100% | 1.00 [1.00, 1.00] | Brier 0.002 [0.001, 0.003] |
| S1 | PCA L2-logistic | 1.00 | 100% | 0.69 [0.56, 0.88] | Brier 0.19 [0.12, 0.27] |
| S2 (train 246 → test 316, n = 16) | R0 rule | 1.00 | 100% | **0.81 [0.63, 0.94]** | misclassification rate 0.19 |
| S2 | R1 logistic (C = 1) | 1.00 | 100% | 0.94 [0.81, 1.00] | Brier 0.042 [0.018, 0.074] |
| S2 | R2 L2-logistic | 1.00 | 100% | 1.00 [1.00, 1.00] | Brier 0.028 [0.008, 0.056] |

- AUC = 1 in **every** bootstrap replicate for all shortlisted models, so AUC cannot rank them. The resampling interval collapses at the ceiling and understates uncertainty with only 15–16 mice and 2 studies.
- The real difference is **threshold transfer** (balanced accuracy), shown in the reverse direction (S2).
- **Correction to earlier wording:** not every non-PCA representation reached AUC 1.0 in Phase 4E. Trees, forests and Naive Bayes ranged from 0.63 to 0.94.

## 5. Multi-study redevelopment
**Not performed.** There are 2 PE labs (3 with any polymer), and the Hong polymers share controls. Leave-one-lab-out would give 2–3 heterogeneous folds, and holding out a Hong polymer while keeping its control mice in training is invalid (`03_…md` §C).

## 6. Honest verdicts

| Question | Verdict |
|---|---|
| Pipeline mismatch? | **Technical mismatch fixed** (count unit). No sign, ID, order or leakage error. Negligible effect on transfer |
| Did the fixes help? | **Internal improvement only / inconclusive external power.** The diet-averaged target gives a consistent but non-significant (P 0.87–0.88), non-specific, post-hoc gain. The robust estimator, Huber/rank and HGB give no gain |
| Generalisation of the gene-level model to new animals? | **Inconclusive external power.** Rank transfer about equals internal CV, but it cannot be distinguished from non-specific mouse-state variation, and its mouse-level interval includes 0 |
| Biology | **No PE-specific biological transfer detected.** This is not proof of absence: the shared fasting/feeding-state axis dominates and the domain shift is unidentifiable |

## 7. What resampling and preprocessing could not fix
- **GSE294180:** PE ≡ lower depth/yield; the CD+Veh clock state; 4–5 mice per arm.
- **Hong:**
  - 3 v 3 with shared controls;
  - 71 wk in food;
  - unknown dose and size;
  - grain chow vs purified diet;
  - controls in a different baseline metabolic state (FAST −0.11 vs +0.10).
- **No resampling scheme creates the missing independent experiments.**
