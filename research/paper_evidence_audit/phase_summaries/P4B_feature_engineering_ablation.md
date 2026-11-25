# Roadmap Phase 4B: feature engineering and block ablation (EXPLORATORY; 2026-11-22)

## Main sources (`research/roadmap_phase4b/`)
- `feature_engineering_plan.md` and `comparison_protocol.md`: locked before any new feature was scored on PE, but written **after** the Phase 4 results were known.
- `decision_memo.md`
- `results/metrics_main.csv`, `all_paired_comparisons.csv`, `extra_vs_tree_baselines.csv`, `shortlist_robustness.csv`, `control_shuffled_B.csv`, `ols_coefficient_stability.csv`
- `figures/fig1_representation_x_model_r2.png`, `fig2_feature_block_ablation.png`, `fig3_stability_per_split.png`
- Scripts: `p4b_features.py` (`RepBuilder` transformer), `p4b_run.py`, `p4b_analyze.py`, `p4b_sensitivity.py`

## 1. Question
- Do engineered representations of the **reference** data improve held-out PE-effect prediction over F8?
- Is any improvement reliable, interpretable and worth the complexity?

## 2. Why
- Phase 4 found a small, robust linear signal (R² ≈ 0.07). Tree models were statistically tied with the linear one.
- **Target-free diagnostics** (`results/feature_diagnostics.txt`) suggested three ideas:
  - the PPARα reference has low reliability (0.34);
  - reference SEs act as abundance proxies (|ρ| up to 0.93);
  - pathway context or interactions might matter.

## 3. Feature blocks
All blocks are built from reference data only, inside the fitted pipeline.

| Block | Features | Count |
|---|---|---|
| A: reliability shrinkage | shr_k = eff_k · τ̂²_k / (τ̂²_k + se_k²) | 4 |
| B: pathway membership | binary indicators for 12 frozen sets: LIPOGENESIS_CORE, CHOL_SYNTHESIS_CORE, SANDERSON_PPARA_TARGETS + 9 Hallmark. Covers 1,029 / 11,758 genes (8.7%). | 12 |
| C: response shape | FASTxPPARA and MASHxFAST (z-score products); FAST_morning_dev = FAST_clean − FAST; PPARA_consistent (agonist effect kept only where the dose response is sign-consistent) | 4 |

- **Representations:** F4, F8, A4, F8+A, F8+B, F8+C, FULL (F8+A+B+C; 28 features), FULL−A, FULL−B, FULL−C.
- **Models:** OLS, ridge, elastic net, HGB, CatBoost. RF dropped: slow, and tied with the boosters in Phase 4.
- **Locked improvement rule vs F8:** NB-corrected CI < 0 **and** ≥ 12/15 splits better **and** ΔRMSE ≤ −0.005 **and** ΔR² ≥ +0.01.
- **Controls:**
  - permuted target;
  - block-B membership shuffled across genes (seed 13).

## 4. Results
Source: `decision_memo.md`. Exact per-model values are in `../model_comparisons.md`, extracted from `results/metrics_main.csv`.

**FULL-linear vs F8-linear:**
- PE_in_CD: R²_oos 0.066 → **0.104**; ΔRMSE −0.0065 [−0.012, −0.0006]; better in 15/15 splits.
- PE_in_MD: 0.070 → **0.104–0.108**; ΔRMSE −0.0053 to −0.0060.
- The locked rule was met. The gain held in 7/7 variants and 19/19 LOO targets (ΔR² +0.008 to +0.055).

**FULL-linear vs the strongest tree baseline:**
- Tied with CatBoost-F8: ΔRMSE −0.003 (CD) and −0.004 (MD), both CIs include 0.
- Also tied with CatBoost-FULL.
- Practical size: RMSE 0.316 → 0.309 log2, a typical error about 2% smaller.

**Ablations:**

| Block | Verdict | Evidence |
|---|---|---|
| A | not useful | A4 ≈ F8 (ΔR² ≤ +0.004). Creates near-duplicate columns: FAST raw −0.32 vs shrunken +0.36 per SD inside FULL. |
| B | helps linear models for CD; not individually significant | F8+B gives CD ΔR² +0.019 (CI includes 0). The shuffled-membership control removes the gain (ridge CD 0.083 real vs 0.066 shuffled). The largest coefficient is a negative LIPOGENESIS_CORE offset, so B re-learns PE's own pathway shifts. |
| C | most consistent contributor | Removing C from FULL significantly worsens every linear model for MD (ΔRMSE +0.003 to +0.004, CI > 0) and ridge for CD. Stable coefficients: FAST_morning_dev (+0.050 MD, −0.019 CD); FASTxPPARA (+0.024 CD, +0.013 MD). |

- Tree models gained only ≈ +0.01–0.02 R² from any block and never met the rule.
- **Reproducibility:** F8 rerun through the new code reproduces Phase 4 exactly for OLS, ridge, elastic net and HGB. CatBoost differs slightly on 9/30 splits (max |ΔRMSE| 0.0016).

## 5. Caveats
- **Exploratory.** The features were chosen after the Phase 4 results, and the same splits were reused, so the intervals are optimistic.
- FULL coefficients are not interpretable (diagnosed in Phase 4C).
- Block B is effectively target-informed: it encodes known PE pathway shifts.

## 6. Key message for the paper
- Engineered features raised the **linear** model to the level of the best **tree** model. They did **not** produce a better predictor than the tree models.
- The interpretable choice is F8+C (response shape).
- The decision memo recommended stopping feature tuning on these splits and moving to diagnostics and an independent study.

## 7. Transition
FULL's coefficients flipped between folds. A diagnostic audit was therefore needed before freezing any model: collinearity, residuals, target noise and leakage (Phase 4C).
