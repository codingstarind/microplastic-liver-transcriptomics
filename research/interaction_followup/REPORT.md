# Interactions among reference-signature features in the Phase 4 CatBoost and HGB models (exploratory; 2026-11-25)

This report gives numbers only. All values come from `results/*.csv` and `TABLE_results_concise.csv`.

## 1. Models, data and reproduction
- **Models:** the Phase 4 F8-CatBoost and F8-HGB, refitted on each locked outer split with that split's stored training-only hyperparameters and DL weights. Phase 4 saved predictions but not fitted models.
- **Splits:** 20 co-expression modules, 5 folds × 3 repeats.
- **Diets:** PE_in_CD and PE_in_MD.
- **Reproduction:** refitted held-out predictions match the saved Phase 4 predictions (max |difference| 2.2e-16, 4 × 35,274 predictions; `results/refit_check_vs_saved_phase4.csv`).
- **Features (F8):**
  - Reference effects: MASH_eff, FAST_eff, PPARA_eff, CIRC_eff.
  - Their standard errors: MASH_se, FAST_se, PPARA_se, CIRC_se.
  - This gives 28 feature pairs: 6 effect × effect, 16 effect × SE, 6 SE × SE.
- **Interpreter:** `C:\Python314\python.exe`, which has CatBoost 1.2.10 and the same scikit-learn 1.9.1 / numpy 2.5.1 as `.venv`. `.venv` has no CatBoost, and no SHAP package is installed in either.

## 2. Methods

**Training-fold diagnostics** (training genes only):
- Friedman pure-interaction SD and H² for all 28 pairs, from 2-D minus 1-D partial dependence on 200 random training genes per split.
- CatBoost native SHAP interaction values on 1,000 training genes per split, and separately on the lipogenesis and PPARα-target genes present in each training fold.

**Permutation control** (10 replicates per split):
- An interaction-free target is built as y* = additive-HGB fit + permuted training residuals.
- The same model class with the same parameters is refitted on y*, and the same diagnostics are recomputed.

**Held-out comparison:** each full model against an additive version with the same parameters (HGB with `interaction_cst='no_interactions'`; CatBoost with depth = 1).

**Criteria (fixed before results were inspected):**
- A pair is **reproducible** if its observed interaction SD exceeds all 10 null replicates in ≥ 12/15 splits in CD and in ≥ 12/15 splits in MD.
- **Held-out primary family:** {ALL, LIPO, PPARA} × {CatBoost, HGB} × {CD, MD} = 12 tests, module-cluster t, Holm.
- **Pair-only follow-up** (only for a reproducible pair): HGB allowed only that pair's interaction vs additive, {ALL, LIPO, PPARA} × {CD, MD} = 6 tests, Holm.
- Cholesterol synthesis and other genes are reported with BH as secondary.

## 3. Training-fold interaction diagnostics vs null (`results/pair_summary_CD_MD.csv`)

Values are the median over 15 splits.

| Model | Pair | Interaction SD, CD (null) | Splits exceeding null, CD | Interaction SD, MD (null) | Splits exceeding null, MD | Reproducible |
|---|---|---|---|---|---|---|
| HGB | FAST_eff × CIRC_se | 0.0259 (0.0044) | 15/15 | 0.0151 (0.0046) | 13/15 | **yes** |
| HGB | FAST_eff × PPARA_eff | 0.0148 (0.0038) | 13/15 | 0.0118 (0.0050) | 11/15 | no |
| HGB | FAST_eff × FAST_se | 0.0184 (0.0047) | 14/15 | 0.0028 (0.0042) | 0/15 | no |
| HGB | MASH_se × FAST_eff | 0.0036 (0.0027) | 2/15 | 0.0127 (0.0053) | 14/15 | no |
| CatBoost | FAST_eff × CIRC_se | 0.0291 (0.0071) | 14/15 | 0.0114 (0.0051) | 9/15 | no |
| CatBoost | FAST_eff × FAST_se | 0.0234 (0.0075) | 14/15 | 0.0050 (0.0059) | 0/15 | no |
| CatBoost | FAST_eff × PPARA_eff | 0.0138 (0.0050) | 12/15 | 0.0056 (0.0034) | 6/15 | no |
| CatBoost | MASH_se × FAST_eff | 0.0073 (0.0051) | 4/15 | 0.0098 (0.0040) | 12/15 | no |

- **Reproducible pairs:** HGB 1 (FAST_eff × CIRC_se, an effect × SE pair); CatBoost 0.
- **Effect × effect pairs:** none meets the criterion for either model.
- **Top-ranked pair per split** (FAST_eff × CIRC_se): HGB 11/15 splits in CD and 9/15 in MD; CatBoost 10/15 and 9/15.
- **Fold stability** of pair rankings, Kendall's W over 15 splits: CatBoost 0.587 (CD), 0.644 (MD); HGB 0.607 (CD), 0.709 (MD).
- **CD–MD repeatability**, Spearman over 28 pairs:
  - median interaction SD: CatBoost 0.495, HGB 0.297;
  - excess over null: CatBoost 0.253, HGB −0.016.
- **CatBoost–HGB agreement**, Spearman of median SD over 28 pairs: 0.709 (CD), 0.673 (MD).
- **Pure-interaction surface of FAST_eff × CIRC_se**, HGB (`results/pair_surface_stability.csv`):

| Measure | CD | MD |
|---|---|---|
| Median between-split correlation | 0.918 (min 0.810) | 0.790 (min 0.551) |
| Median surface range | 0.121 log2FC | 0.062 log2FC |

  The CD vs MD median-surface correlation is 0.872.

**CatBoost SHAP interaction values** (training genes; median over splits of mean |φ_jk + φ_kj|; share = relative to the sum of main-effect |φ|; `results/catboost_shap_interaction_summary.csv`). There is no permutation null for the SHAP values; the H-statistic null above is the control.

| Genes | Diet | Pair | Mean \|SHAP\| | Share of main effects |
|---|---|---|---|---|
| All (1,000 sampled) | CD | FAST_eff × FAST_se | 0.0095 | 0.069 |
| All (1,000 sampled) | CD | FAST_eff × CIRC_se | 0.0087 | 0.064 |
| All (1,000 sampled) | MD | MASH_se × FAST_eff | 0.0082 | 0.057 |
| All (1,000 sampled) | MD | FAST_eff × CIRC_se | 0.0057 | 0.039 |
| Lipogenesis (median 11 in training fold) | CD | FAST_eff × CIRC_se | 0.353 | 0.629 |
| Lipogenesis | MD | FAST_eff × CIRC_se | 0.076 | 0.186 |
| PPARα targets (median 12 in training fold) | CD | FAST_eff × PPARA_eff | 0.048 | 0.213 |
| PPARα targets | MD | FAST_eff × CIRC_se | 0.017 | 0.084 |

## 4. Held-out: full model vs additive model (`results/heldout_full_vs_additive_by_group.csv`)

| Diet | Model | Group | R² full / additive | ΔR² [95% module-bootstrap CI] | p (cluster t) | Holm p (12 tests) |
|---|---|---|---|---|---|---|
| CD | CatBoost | ALL (11,758) | 0.085 / 0.064 | +0.021 [0.007, 0.038] | 0.019 | 0.20 |
| CD | CatBoost | LIPO (12) | 0.417 / 0.236 | +0.181 [0.135, 0.244] | 0.020 | 0.20 |
| CD | CatBoost | PPARA (14) | 0.164 / 0.088 | +0.075 [−0.110, 0.296] | 0.54 | 0.70 |
| CD | HGB | ALL | 0.079 / 0.064 | +0.014 [0.003, 0.027] | 0.031 | 0.21 |
| CD | HGB | LIPO | 0.294 / 0.163 | +0.131 [0.089, 0.157] | 0.021 | 0.20 |
| CD | HGB | PPARA | 0.345 / 0.213 | +0.132 [0.023, 0.214] | 0.078 | 0.39 |
| MD | CatBoost | ALL | 0.075 / 0.054 | +0.021 [0.008, 0.039] | 0.012 | 0.14 |
| MD | CatBoost | LIPO | 0.368 / 0.200 | +0.168 [−0.519, 0.285] | 0.031 | 0.21 |
| MD | CatBoost | PPARA | 0.020 / 0.039 | −0.019 [−0.250, 0.139] | 0.23 | 0.70 |
| MD | HGB | ALL | 0.075 / 0.061 | +0.014 [0.004, 0.025] | 0.026 | 0.20 |
| MD | HGB | LIPO | 0.363 / 0.221 | +0.142 [−1.253, 0.162] | 0.099 | 0.40 |
| MD | HGB | PPARA | 0.055 / −0.000 | +0.055 [−0.372, 0.063] | 0.32 | 0.70 |

- **Primary tests with Holm p < 0.05: 0 of 12.**
- **Secondary (BH):**
  - Cholesterol synthesis, full − additive: ΔR² +0.106 to +0.213, q 0.020–0.026.
  - Other genes: q 0.020–0.061.

## 5. Held-out: HGB pair-only (FAST_eff × CIRC_se) vs additive (`results/pair_only_vs_additive_heldout.csv`)

| Diet | Group | R² pair-only / additive / full | ΔR² pair-only − additive [CI] | Share of full − additive | p | Holm p (6 tests) |
|---|---|---|---|---|---|---|
| CD | ALL | 0.071 / 0.064 / 0.079 | +0.006 [−0.001, 0.015] | 0.45 | 0.15 | 0.44 |
| CD | LIPO | 0.292 / 0.163 / 0.294 | +0.129 [0.102, 0.173] | 0.99 | 0.0090 | 0.054 |
| CD | PPARA | 0.252 / 0.213 / 0.345 | +0.039 [−0.047, 0.099] | 0.30 | 0.31 | 0.52 |
| MD | ALL | 0.066 / 0.061 / 0.075 | +0.005 [0.002, 0.008] | 0.35 | 0.025 | 0.13 |
| MD | LIPO | 0.304 / 0.221 / 0.363 | +0.083 [−0.395, 0.124] | 0.58 | 0.031 | 0.13 |
| MD | PPARA | 0.044 / −0.000 / 0.055 | +0.044 [−0.052, 0.046] | 0.80 | 0.26 | 0.52 |

- **Pair-only tests with Holm p < 0.05: 0 of 6.**
- **Secondary (BH):** cholesterol synthesis +0.111 (CD, q 0.0096) and +0.024 (MD, q 0.0096); other genes q 0.54 (CD) and 0.11 (MD).

## 6. Result against the pre-set criteria
- **Reproducible training-fold interaction:** one pair, HGB only: FAST_eff × CIRC_se (fasting reference effect × circadian reference standard error).
- **Effect × effect pairs** (interactions between two reference signatures): none reproducible for either model.
- **Held-out benefit:** none of the 18 primary held-out tests (12 full-vs-additive + 6 pair-only) has Holm p < 0.05. The smallest is 0.054 (lipogenesis, CD, pair-only).

## 7. Limits
- Exploratory: these splits have been reused since Phase 4.
- The lipogenesis and PPARα groups span 5 and 9 modules, so the module-cluster tests have 4 and 8 df.
- The CatBoost additive comparator (depth 1) also differs in tree capacity, not only in the absence of interactions.
- Only 10 null replicates per split, so per-split exceedance corresponds to p < 1/11.
- SHAP interaction values have no separate null.
- No new hyperparameters were tuned; all held-out genes were used only for prediction.
