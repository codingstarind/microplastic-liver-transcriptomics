# Roadmap Phase 4: gene-level multi-axis prediction with nested grouped CV (2026-11-22)

## Main sources (`research/roadmap_phase4/`)
- `research_question_and_prediction_targets.md`
- `evaluation_protocol.md` (LOCKED; `results/protocol_lock.json`)
- `model_search_spaces.csv`
- `literature_model_justification.md`
- `phase4_decision_report.md`
- `results/taskA/*.csv`, `results/taskB/*.csv`
- `figures/taskA_*.png`, `figures/taskB_classification_demo.png`
- Scripts: `p4_data.py`, `p4_models.py`, `p4_run.py`, `p4_analyze.py`, `p4_classify.py`, `p4_depth_thinning.R`, `p4_alt_universe_refs.R`

## 1. Question
Do independent MASH, fasting, PPARα and time-of-day reference effects predict held-out genes' PE-associated log2FC better than:
- a null model;
- single-axis models;
- linear baselines?

Do flexible ML models (trees, forests, boosting) beat a linear model?

## 2. Why this was necessary
- Phase 3 single-axis correlations were weak (|ρ| ≤ 0.25), and the MASH axis was fragile to single mice.
- The open question was whether a combined, SE-weighted, out-of-sample model adds information.
- The Phase 3 spec memo (`roadmap_phase3/reports/phase4_specification_memo.md`) fixed targets, weights, CV grouping and decision rules **before** any fit.

## 3. Design
| Element | Choice |
|---|---|
| Rows | 11,758 genes (main universe, complete features). Genes are **not** independent units: every gene's effect comes from the same 9–10 mice. |
| Targets | GSE294180 DESeq2 MLE log2FC: PE_in_CD and PE_in_MD (separately). Diet_in_Veh = positive-control target. |
| Weights | w = 1/(s² + τ²_train). τ² by DerSimonian–Laird method of moments on training rows only. |
| Features | **F8** = log2FC + SE of MASH1 (GSE316023), FAST (GSE107787, ZT-averaged), PPARA (GSE282998 WT 20v0), CIRC_T1 (GSE117134 ZT6 − ZT2). **F4** = the four effects only. |
| Groups | 20 Ward co-expression modules (GSE117134, 54 samples, z-scored log2 RPKM). Sensitivity: chromosome. |
| Outer CV | 5 folds × 3 repeats (seeds 101–103) = 15 splits. Identical held-out genes for every model. |
| Inner CV | GridSearchCV + GroupKFold(4), scoring = neg-MSE |
| Models | null; ols_MASH; ols_FAST; ols_PPARA; ols_CIRC; ols_all; ridge; elastic net; CART tree; RF; HistGradientBoosting; CatBoost. AdaBoost and neural networks deliberately **not** run. |
| Metrics | RMSE (primary); MAE; R²_oos = 1 − SSE/SSE(training mean); Spearman; sign agreement; weighted R² |
| Paired comparisons | ΔRMSE with the Nadeau–Bengio corrected resampled-t interval (14 df). "Better" requires CI < 0 **and** a win in ≥ 12/15 splits. |
| Robustness | 8 variants (V1 MASH2, V2 FAST_clean, V3 PPARA 10v0, V4 CIRC_T2, V5 strict universe, V6 alternative muscle universe, V7 chromosome groups, V8 effects-only); 19 leave-one-mouse-out (LOO) targets |
| Negative controls | target permuted across genes; binomial depth-thinned vehicle target (`p4_depth_thinning.R`) |
| Leakage tests | `tests/test_phase4.py` (no gene or group in both training and test; inner folds ⊂ outer training; no forbidden features) |

## 4. Results
Sources: `results/taskA/summary_by_run_model.csv`, `paired_comparisons.csv`, `robustness_counts.csv`, `negative_controls_r2_oos.csv`, `coefficient_stability.csv`.

**Held-out R²_oos (mean of 15 splits):**

| Model | PE_in_CD | PE_in_MD | Diet_in_Veh |
|---|---|---|---|
| ols_MASH | 0.031 | −0.002 | 0.206 |
| ols_FAST | 0.022 | 0.045 | 0.012 |
| ols_PPARA | 0.000 | −0.001 | 0.001 |
| ols_CIRC | 0.022 | 0.005 | 0.040 |
| **ols_all (F8)** | **0.066** | **0.070** | 0.233 |
| ridge | 0.066 | 0.070 | 0.233 |
| elastic net | 0.063 | 0.068 | 0.232 |
| tree | 0.038 | 0.053 | 0.179 |
| RF | 0.073 | 0.081 | 0.247 |
| HGB | 0.078 | 0.079 | 0.285 |
| CatBoost | 0.086 | 0.079 | 0.306 |

**RMSE, null → best:**
- PE_in_CD: 0.326 → 0.312
- PE_in_MD: 0.309 → 0.296
- Diet_in_Veh: 0.490 → 0.408

**Paired ΔRMSE, ols_all vs comparators (NB-corrected 95% CI):**

| Comparison | PE_in_CD | PE_in_MD |
|---|---|---|
| ols_all vs null | −0.0107 [−0.017, −0.004] | −0.0112 [−0.025, +0.003] (not distinguishable) |
| ols_all vs best single axis | −0.0058 [−0.011, −0.0003] (vs MASH) | not distinguishable (vs FAST) |

**Tree models vs ols_all:**
- RF, HGB and CatBoost all have CIs including 0 for both PE targets.
- On the positive-control target (Diet_in_Veh), boosting does beat linear: HGB −0.015 [−0.026, −0.004]; CatBoost −0.022 [−0.037, −0.006].

**Robustness:**
- ols_all beats null in 8/8 variants and 19/19 LOO targets.
- LOO R²_oos ranges: CD 0.048–0.082; MD 0.053–0.109.
- ols_MASH is fragile: CD range −0.016 to 0.051; MD below null in 12/19 targets.
- ols_FAST is stable.
- HGB had lower mean RMSE than ols_all in most variants, but its primary paired CI includes 0. The difference is ≈ 0.001–0.002 log2.

**Negative controls (R²_oos):**
- permuted target: −0.002 to 0.000;
- depth-thinned target: −0.005 to +0.001 (thinning produced 0 DEGs).

**OLS coefficients (per SD; sign stable in 15/15 splits):**

| Feature | PE_in_CD | PE_in_MD |
|---|---|---|
| MASH | +0.060 | +0.040 |
| FAST | +0.035 | +0.058 |
| CIRC T1 | −0.041 | +0.031 |
| PPARα | −0.012 | −0.021 |

**Task B (demonstration only; no claim):**
- GAN-vs-chow classifier trained on GSE316023 (16 + 16 mice), tested externally on GSE246328 (7 v 8).
- External balanced accuracy: logistic 1.00, HGB 1.00, RF 0.94, NB 0.94, trees 0.875.
- This is a ceiling effect, so the classifiers cannot be ranked.

## 5. Failures and limitations
- Complex models never beat linear on the PE targets.
- The PE_in_MD improvement over null is not distinguishable under the corrected interval.
- **Held-out genes are not held-out mice.**
- The depth-thinning control does not test yield or quality.
- The MASH axis is conditional (Phase 3 c1 FAIL).

**Deviations:**
- D1: variants, controls and LOO used 1 repeat.
- D2: Task B HGB `min_samples_leaf` changed 20 → 3.
- D3: pandas parsed the model name "null" as NA.

## 6. Transition
About 7% of held-out variance was robustly predictable. Phase 4B (exploratory) asked whether engineered representations of the references could raise this.
