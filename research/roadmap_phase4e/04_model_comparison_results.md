# Model comparison results (study-held-out MASH-state task)

**Task B1:** GAN vs chow mice; train on one MASH study, test once on the other. Two directions only.
- S1: train GSE316023 TA (16 mice) → test GSE246328 (15).
- S2: the reverse.

Full table: `results/external_metrics.csv`. Predictions: `results/external_predictions.csv`. Figures: `figures/fig1–fig3`.

## Baseline → feature engineering → tuning (same target, same held-out mice)

| Stage | Representation / model | S1 AUC | S1 balanced acc. | S1 Brier | S2 AUC | S2 balanced acc. | S2 Brier |
|---|---|---|---|---|---|---|---|
| null | majority | 0.50 | 0.50 | 0.250 | 0.50 | 0.50 | 0.251 |
| **baseline** | R0 signed score, training-midpoint rule | 1.00 | 1.00 | 0.000 | 1.00 | **0.81** | **0.188** |
| + feature engineering | R1 (12 pathways), logistic C = 1 | 1.00 | 1.00 | 0.010 | 1.00 | 0.94 | 0.042 |
| + feature engineering | R2 (R0 + R1), logistic C = 1 | 1.00 | 1.00 | 0.009 | 1.00 | 0.94 | 0.049 |
| + tuning | R2, L2-logistic (inner CV chose C = 10) | 1.00 | 1.00 | 0.002 | 1.00 | **1.00** | **0.028** |
| comparator | PCA (3 PCs, training-only loadings), L2-logistic | 1.00 | 0.69 | 0.187 | 1.00 | 0.75 | 0.158 |

**Other models** (tuned by inner CV):
- Gaussian NB on R1/R2: AUC 1.0 in both directions.
- Depth-1/2 trees: Gini always selected; R1 AUC 1.0, R2 0.80 / 0.63, PCA 0.69 / 0.63.
- Random forest: R1/R2 1.0; PCA 0.91 / 0.84.
- Sensitivity S1b (training on 21 °C + 29 °C mice, 32 total): same picture.

## Reading
- **Discrimination is at ceiling** for every non-PCA representation in both directions: GAN vs chow is easy. AUC cannot rank models here.

  The bootstrap CI collapses to [1, 1], which understates the uncertainty of 15–16 test mice. A perfect score on 15 mice is still only 15 mice.
- **What feature engineering and tuning improved is cross-study calibration.**
  - The R0 rule's threshold, learnt in one study, misclassified 3 of 16 S2 mice (balanced accuracy 0.81, Brier 0.19), because rank-score scales shift between platforms.
  - Pathway features plus L2-logistic gave balanced accuracy 1.0 and Brier 0.03.
  - This is the BT3041 point: the improvement is in **threshold transfer and probability quality**, not in AUC.
- **PCA learnt study-specific axes and transferred worst** (balanced accuracy 0.69–0.75), as expected when loadings come from one platform.
- **Permutation null** (200 training-label permutations, R2 logistic):
  - S1 null mean AUC 0.51, 95th percentile **1.00**, p = 0.10;
  - S2 null mean 0.52, 95th percentile 0.97, p = 0.035.

  Permuted labels that keep most GAN mice together still carry the huge diet signal. This is the same limitation as the Phase 3 null, so the permutation test has little resolving power here.
- **Predictions vs library size in the test study:** |ρ| ≤ 0.43 for every R0/R1/R2 model (`rho_pred_libsize`), and the direction differs between S1 and S2. There is no consistent depth artefact.

## Fasting demonstration (B2; `results/fasting_demo_metrics.csv`)
- Trained on the GSE107787 muscle-clean subset (10 mice); tested on GSE225695 (6 mice).
- R0 rule and logistic: AUC 1.0; 6/6 correct.
- **n = 6**, sacrifice time unknown: a direction check only, not a validated classifier.

## Application to the 19 GSE294180 mice (C; only models that passed validity)

Scores are log-odds (not calibrated probabilities); 95% t-CI.

| Score | Diet_in_Veh (positive control) | PE_in_CD | PE_in_MD | Interaction | Leave-one-out signs |
|---|---|---|---|---|---|
| R0 + L2-logistic (both-study signature) | **+4.01 [2.92, 5.10]** | +0.35 [−0.74, 1.44], Holm 0.96 | +0.86 [−0.16, 1.89], Holm 0.28 | +0.51, Holm 0.96 | 19/19 |
| R2 + L2-logistic | **+5.55 [4.78, 6.33]** | +0.83 [0.06, 1.61], Holm 0.072 | **+1.66 [0.93, 2.39], Holm 0.0006** | +0.83, Holm 0.12 | 19/19 |
| Phase 4D MASH singscore (reference) | +0.215 [0.158, 0.273] | +0.025, Holm 0.75 | +0.054, Holm 0.15 | +0.029 | 19/19 |

**What drives R2's PE_in_MD shift** (`results/r2_feature_contributions_gse294180.csv`):

| Feature | Contribution to PE_in_MD (log-odds) |
|---|---|
| HALLMARK_CHOLESTEROL_HOMEOSTASIS | +1.17 |
| CHOL_SYNTHESIS_CORE | +0.93 |
| SANDERSON_PPARA_TARGETS | +0.35 |
| MTORC1 | +0.28 |
| LIPOGENESIS_CORE | −1.37 |
| the MASH program signature itself | **+0.18** |
| inflammatory, EMT | ≈ 0 |

- The classifier learnt "low cholesterol synthesis = GAN-like" from the chow-controlled studies, and PE lowers cholesterol synthesis.
- **The significant R2 "MASH-state" shift is therefore a cholesterol/PPARα pathway shift, not a disease-specific one.** This agrees with Phase 4D, where the MASH-singscore PE shift vanished without the cholesterol-synthesis genes.
- R2 scores correlate 0.89 with the Phase 4D fasting score and 0.67 with clock_phase. Depth |ρ| ≤ 0.22.
