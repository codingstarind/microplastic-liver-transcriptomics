# 02: Domain shift and error diagnosis (Stage 2)

**Script:** `b03_gap_diagnostics.py`; mouse-resampled DESeq2 refits from `b02_external_target_resampling.R`.
**Tables:** `tables/gap_*.csv`. **Figure:** `figures/fig1_gap_targets_labellings.png`.
**No model was fitted.** Everything uses frozen predictions.

## 1. Train vs outer CV vs external (same metric, same gene-effect task; `gap_train_cv_external.csv`)

| Model (trained on) | Train ρ | Outer-CV ρ (15 folds) | Hong PE ρ | GSE245069 PS ρ | Train R² | CV R² | Hong R² (vs train mean) |
|---|---|---|---|---|---|---|---|
| F8-OLS (PE_in_CD) | 0.274 | 0.264 | **0.219** | −0.069 | 0.070 | 0.066 | 0.018 |
| F8+C-Ridge (PE_in_CD) | 0.281 | 0.270 | 0.142 | −0.080 | 0.087 | 0.080 | −0.005 |
| F8-OLS (PE_in_MD) | 0.324 | 0.309 | **0.303** | 0.181 | 0.073 | 0.070 | 0.047 |
| F8+C-Ridge (PE_in_MD) | 0.389 | 0.374 | **0.417** | 0.227 | 0.096 | 0.096 | 0.076 |

- **Bias–variance:** train ≈ CV, so the models are high-bias / low-variance. They do not overfit; they are weak.
- **In rank terms, Hong transfer is about as good as internal CV** (0.22 vs 0.26; 0.30 vs 0.31). The R² collapse is mostly a target-scale/noise effect: Hong SD 0.54 vs 0.33, median SE 0.25 vs 0.145.
- The "failure" in Hong is therefore **not a generalisation gap in the usual sense**. It is a weak predictor meeting a noisier target. The real problem is specificity (§4).
- **GSE245069 (PS):** the CD-trained models reverse sign (cross-polymer heterogeneity). The MD-trained models stay positive.

## 2. Targets: noise, tails, reliability (`gap_target_distributions.csv`, `gap_reliability_and_disattenuation.csv`)

| Target | SD | Robust SD | Excess kurtosis | Median SE | Signal share (SE-based) |
|---|---|---|---|---|---|
| GSE294180 PE_in_CD | 0.33 | 0.21 | 12.4 | 0.145 | 0.36 |
| Hong PE vs CON | 0.54 | 0.37 | 17.8 | 0.249 | 0.47 |
| Hong PS vs CON | 0.44 | 0.23 | 32.6 | 0.250 | **0.20 (noise-dominated)** |
| GSE245069 PS pooled | 0.37 | 0.22 | 18.0 | 0.128 | 0.65 |

- **Robust delete-one-mouse jackknife reliability:** Hong 0.66, GSE294180 PE_in_CD 0.55.
  - The mean-variance version (0.07) is invalidated by heavy tails and is reported only for transparency.
- **Between-study agreement of the observed PE vectors:** GSE294180 PE_in_CD vs Hong PE ρ = 0.10 (PE_in_MD vs Hong: 0.26).
  - **Disattenuated for both reliabilities: about 0.17.** Even after correcting for measurement noise, the two PE experiments share little gene-level effect. This points to domain shift or genuinely different biology, not only noise.

## 3. Where the errors are (`gap_strata_hong_pe.csv`)
- **Low-expression / high-SE genes carry no transfer:**
  - Hong baseMean lowest tercile ρ = 0.06; highest SE tercile ρ = 0.10;
  - GSE294180 high-SE tercile ρ = 0.04, vs 0.32 in its low-SE tercile.
- **Calibration slope** is about 1.2 in the reliable strata, so there is no systematic scale or sign shift, and about 0.3–0.4 in the noisy strata (attenuation by noise).
- **Mean residuals are near 0 in every stratum:** no systematic bias.

## 4. The decisive diagnostic: what the external metric measures (`gap_hong_mouse_pca.csv`, `gap_hong_20_labellings.csv`)
- **Mouse PCA of the 15 Hong mice:** PC1 explains 18% of variance. Its gene loadings correlate with the F8 prediction (ρ = 0.34) and with the reference fasting effect (ρ = 0.22).
- **Across all 20 relabellings of the 6 PE/CON mice:** transfer ρ is **89% explained (R² = 0.889)** by each labelling's difference along mouse PC1–PC3.
  - The true PE labelling ranks **5th of 20** (ρ_proxy 0.18). Four arbitrary mixes of PE and control mice score higher.
- **Interpretation:** the frozen models predict the dominant axis of between-mouse variation, a fasting / feeding-state-like programme. Any grouping of mice that happens to split along it "validates".
  - The PE label in Hong partly aligns with that axis. So does every polymer against the same shared controls (Phase 5: PP ρ 0.39, PB 0.19).

## 5. External target instability (mouse is the unit; `gap_hong_target_resampling.csv`)
- 15 mouse-resampled DESeq2 refits (6 leave-one-out, 9 two-vs-two) agree with the full target at ρ 0.87–0.96.
- **Yet F8 transfer ρ ranges 0.10–0.31 (CD-trained) and 0.14–0.39 (MD-trained).** A single mouse moves the headline metric by about ±0.1.
- With 3 v 3, no gene-level metric can be precise.

## 6. Study shift (`gap_study_shift.csv`)
- **Not identifiable** (one study per level):
  - polymer and particle;
  - dose and route (gavage vs food vs water);
  - duration (8 vs 71 vs 12 wk) and age;
  - control diet (purified vs chow);
  - strain (C57BL/6 vs ICR);
  - quantifier and annotation (Salmon/GENCODE vs StringTie/RefSeq vs vendor/NCBI).
- **Measured:**
  - **Baseline metabolic state of the controls differs strongly.** Mean FAST singscore: GSE294180 CD+Veh **−0.112**, Hong CON **+0.101**, GSE245069 CON +0.044. MASH: −0.44 / −0.33 / −0.45. External controls sit in a more fasted-like state (chow, older age, unknown sacrifice time).
  - Depth is balanced externally; in GSE294180, PE is confounded with lower depth/yield.
  - Count unit (Hong 2×).
  - Heavy-tailed source target.
  - Muscle CPM similar (1,076 / 946 / 844).

## 7. Ranked causes (evidence-weighted; uncertainty retained)

| Rank | Cause | Evidence | Strength |
|---|---|---|---|
| 1 | **Non-specific shared axis:** the models capture a fasting-like axis of between-mouse variation, so external agreement is not PE-specific | 20-labelling R² 0.89, true rank 5/20; cross-polymer ρ; shared controls | strong |
| 2 | **External target noise and n = 3 v 3** | transfer ρ 0.10–0.31 across single-mouse changes; high-SE strata ρ ≈ 0.1; Hong PS signal share 0.20 | strong |
| 3 | **Weak predictor (high bias):** internal R² is only 0.07 | train ≈ CV; noise ceiling (Phase 4C) | strong |
| 4 | **Biological / design domain shift** (diet, duration, dose, polymer, baseline state) | disattenuated between-study ρ 0.17; baseline FAST difference 0.21; CD→PS sign reversal | moderate, **not separable** |
| 5 | **Outlier-dominated source target** | top-1% share 0.77 vs 0.33 | moderate (tested by I2a) |
| 6 | **Pipeline mismatch (count unit)** | 2× unit; negligible effect on ρ | verified but **minor** |
| 7 | **Model misspecification** | trees/HGB no better (Phase 4, and B2 here); flat residual trends | weak |

**Bottom line:** the external result is explained mainly by limited specificity and precision of the test, a weak source model and an unidentifiable domain shift. It is **not** explained by a fixable pipeline error. That does not prove generalisation is impossible; the available data cannot separate PE biology from mouse-state variation.
