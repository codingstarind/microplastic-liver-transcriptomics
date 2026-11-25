# Phase 4B decision memo: feature engineering and model ablation

**Status: EXPLORATORY second-stage analysis.** The feature families were chosen after the Phase 4 results were known, and the Phase 4 validation splits were reused. The protocol was locked before any new feature was scored on PE (`results/protocol_lock.json`). One smoke split of FULL-OLS was seen before the full run. Nothing here is an untouched confirmation, and nothing here was preregistered before Phase 4.

## Question
Did any engineered representation materially improve held-out PE-effect prediction (genes as rows, held-out co-expression modules) over the Phase 4 linear and tree baselines? Is that improvement reliable, interpretable and worth the complexity?

## Answer in one paragraph
**Yes, modestly and reliably, for linear models; no, against the best existing tree baseline.**
- **Linear model:** adding all three engineered blocks to F8 (representation FULL) raised out-of-sample R² for OLS/ridge/elastic net:
  - PE_in_CD: from 0.066 to **0.104** (ΔRMSE −0.0065 log2, 95% CI [−0.012, −0.0006], better in 15/15 splits);
  - PE_in_MD: from 0.070 to **0.104–0.108** (ΔRMSE −0.0053 to −0.0060).
- **Robustness:** the gain met the locked rule (CI < 0, ≥ 12/15 splits, |ΔRMSE| ≥ 0.005, ΔR² ≥ 0.01). It held in **7/7 sensitivity variants and 19/19 leave-one-mouse-out targets** (ΔR² +0.008 to +0.055).
- **Against tree baselines:** FULL-linear is **statistically tied** with the Phase 4 CatBoost-F8 baseline (ΔRMSE −0.003 CD, −0.004 MD; CIs include 0) and with CatBoost-FULL.
- **In practice:** RMSE moves from about 0.316 to 0.309 log2, roughly a 2% smaller typical error. The engineered features let a *simple, transparent* model reach what boosting already reached. They do not produce a better predictor.

## What each block contributed (ablations, `feature_ablation_results.csv`, `figures/fig2`)

| Block | Verdict | Evidence |
|---|---|---|
| **A: reliability shrinkage** | **Not useful** | A4 ≈ F4/F8 (ΔR² ≤ +0.004); removing A from FULL never hurt; for MD it slightly helped. Reason: shrunken effects are almost collinear with the raw ones for MASH/FAST (reliability 0.92/0.98), which also makes the FULL coefficients unstable (raw FAST −0.32 vs shrunken +0.36 per SD). Recommend dropping it. |
| **B: pathway membership** | **Helps linear models for CD; not individually significant** | Adding B to F8: ΔR² +0.019 (CD linear), CI includes 0. The shuffled-membership control removes the gain (ridge CD 0.083 real vs 0.066 shuffled = F8), so the gain is real membership information, not a set-size artefact. **Interpretation:** the largest stable coefficient is a negative lipogenesis-core offset, i.e. the model learns pathway-level PE shifts already found by Phase 2 enrichment, not reference similarity. |
| **C: response shape** | **The most consistent contributor, especially MD** | Removing C from FULL significantly hurt every linear model for MD (ΔRMSE +0.003 to +0.004, CI > 0) and ridge for CD. Stable coefficients: FAST_morning_dev (+0.050 MD; −0.019 CD) and FASTxPPARA (+0.024 CD, +0.013 MD). |

Tree models gained much less from any block (HGB/CatBoost ΔR² ≈ +0.01–0.02, never meeting the rule). They could already approximate the interactions that block C makes explicit for linear models.

## Reliability and limits
- **Controls:** a permuted target gives R²_oos ≈ 0 (−0.002 to 0.000), so there is no leakage signal. Phase 4's depth-thinning control (R² ≈ 0) still applies. **No representation addresses the PE ≡ depth/yield confound.**
- **Reproducibility:** F8 rerun through the new code reproduces Phase 4 exactly for OLS, ridge, elastic net and HGB. CatBoost differs slightly on 9/30 splits (max |ΔRMSE| 0.0016, mean 0.00008), because the added scaler shifts its float quantisation borders; on 1 split a near-tied inner hyperparameter choice flipped.
- **Scope:** held-out genes, not held-out mice or laboratories. Genes are correlated, so the effective n is far below 11,758. The splits were reused from Phase 4, so the intervals are optimistic.
- **Interpretability:**
  - FULL contains near-duplicate columns (A vs raw effects). Individual coefficients are therefore not interpretable even though predictions are stable.
  - Block-B coefficients mostly restate PE pathway shifts.
  - **The biological reading should come from simpler models** (e.g. F8+C), not from FULL coefficients.
- **MASH axis:** still conditional (Phase 3 positive control FAIL on c1).

## Recommendation
1. **Keep the simple linear model as the descriptive model.** If the team wants the extra 3–4 points of R², use **F8 + C** (ridge): 4 interpretable response-shape features, no near-duplicate columns, and C is the one block whose removal reliably hurts. **FULL is not recommended:** its gain is modest and its coefficients are uninterpretable.
2. **Do not adopt block A**; treat **block B** as an optional pathway-offset term (it encodes known PE enrichment, not reference similarity).
3. **Do not keep tuning features on these splits.** Further gains would be overfitting the validation design. **Move on to independent cross-study evaluation** (Phase 5, Track B), where the question is whether the same reference-based description holds for other plastic studies.

## Course concepts that actually mattered
- **Feature engineering:** interaction and response-shape features (block C) helped linear models; reliability shrinkage did not.
- **Regression and regularisation:** ridge and elastic net ≈ OLS for prediction, but they handle FULL's collinearity better.
- **Model selection and nested CV:** the inner tuning and outer grouped splits kept the comparison fair.
- **Trees, bagging and boosting:** these already captured most nonlinear structure, which is why engineered features helped them little.
- **Hypothesis testing:** corrected paired intervals and the ablations separated real block effects from noise.
- **Enrichment:** pathway membership as a feature (ii) is different from post-hoc enrichment (i), which Phase 2 already did. Here B effectively re-learns (i).
- **Sample-level pathway scores (iii)** (GSVA, singscore, PROGENy) define a different, pathway-by-sample task and were deliberately **not** mixed in.
