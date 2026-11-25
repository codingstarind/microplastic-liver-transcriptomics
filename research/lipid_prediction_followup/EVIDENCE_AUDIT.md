# Evidence audit: lipid-gene vs other-gene held-out prediction

This file holds numbers and calculations only. Definitions are in `README.md`, and every value comes from the CSVs listed below.

**Abbreviations:**
- **R²g** = R²_group = 1 − SSE_model / SSE_null.
- **dMSE** = MSE_null − MSE_model.
- **Diff** = lipid group minus reference.
- **CI** = 95% module-block bootstrap (B = 2,000, 20 modules).
- **perm p** = within-module label permutation (2,000).
- **FE p** = module fixed effects with CR1 cluster-robust SE, t on 19 df.

## 0. Reproduction and inputs
- **Reproduction:** regenerated FULL-OLS and FULL-Ridge, and reused F8-OLS, F8+C-OLS and F8+C-Ridge, against stored Phase 4B per-split metrics (75 model × diet × split cells). Max |ΔRMSE| = 1.1e-16; max |ΔR²| = 3.2e-16 (`results/reproduction_check_vs_phase4b.csv`).
- **Universe:** 11,758 genes; 20 modules.
- **Group sizes:** LIPO 12 genes (5 modules), CHOL 21 (6), PPARA 14 (9), LIPID_any 47 (13), OTHER 11,711 (20), OTHER_noB12 10,729.
- **Pairwise overlap of the three lipid sets:** 0 genes.
- **Overlap with the other 9 block-B sets** (`results/lipid_sets_overlap_with_other_blockB_sets.csv`):
  - LIPO: 6/12 genes in at least one other set (MTORC1 4, cholesterol homeostasis 2, fatty-acid metabolism 2).
  - CHOL: 19/21 (cholesterol homeostasis 17, MTORC1 11).
  - PPARA: 12/14 (fatty-acid metabolism 12).
  - LIPID_any: 37/47.

## 1. Group metrics (primary models; `TABLE_results_concise.csv`)

| Diet | Model | Group | RMSE model / null | R²g [CI] | dMSE [CI] | genes improved |
|---|---|---|---|---|---|---|
| CD | F8-OLS | LIPO | 1.461 / 1.629 | 0.196 [0.061, 0.759] | 0.519 [0.132, 0.986] | 0.917 |
| CD | F8-OLS | CHOL | 0.510 / 0.860 | 0.648 [0.596, 0.851] | 0.479 [0.168, 0.644] | 1.000 |
| CD | F8-OLS | PPARA | 0.346 / 0.379 | 0.171 [−0.025, 0.494] | 0.025 [−0.003, 0.089] | 0.643 |
| CD | F8-OLS | LIPID_any | 0.835 / 1.025 | 0.337 [0.097, 0.594] | 0.354 [0.081, 0.529] | 0.872 |
| CD | F8-OLS | OTHER | 0.316 / 0.325 | 0.052 [0.033, 0.074] | 0.0055 [0.0039, 0.0069] | 0.549 |
| CD | F8+C-Ridge | LIPO | 1.366 / 1.629 | 0.297 [0.180, 0.754] | 0.789 [0.183, 1.192] | 0.833 |
| CD | F8+C-Ridge | CHOL | 0.442 / 0.860 | 0.737 [0.724, 0.861] | 0.545 [0.178, 0.786] | 1.000 |
| CD | F8+C-Ridge | PPARA | 0.318 / 0.379 | 0.300 [−0.073, 0.683] | 0.043 [−0.008, 0.135] | 0.786 |
| CD | F8+C-Ridge | LIPID_any | 0.770 / 1.025 | 0.436 [0.214, 0.685] | 0.458 [0.105, 0.655] | 0.894 |
| CD | F8+C-Ridge | OTHER | 0.314 / 0.325 | 0.062 [0.038, 0.087] | 0.0066 [0.0045, 0.0083] | 0.544 |
| CD | FULL-OLS | LIPO | 0.956 / 1.629 | 0.656 [−2.905, 0.834] | 1.741 [−0.553, 2.309] | 0.833 |
| CD | FULL-OLS | CHOL | 0.286 / 0.860 | 0.889 [0.845, 0.953] | 0.658 [0.182, 0.943] | 0.905 |
| CD | FULL-OLS | PPARA | 0.251 / 0.379 | 0.562 [0.077, 0.873] | 0.081 [0.006, 0.192] | 0.714 |
| CD | FULL-OLS | LIPID_any | 0.537 / 1.025 | 0.726 [0.356, 0.854] | 0.763 [0.071, 1.228] | 0.830 |
| CD | FULL-OLS | OTHER | 0.312 / 0.325 | 0.073 [0.051, 0.099] | 0.0077 [0.0060, 0.0094] | 0.547 |
| MD | F8-OLS | LIPO | 0.988 / 1.199 | 0.321 [0.017, 0.703] | 0.461 [0.005, 0.623] | 0.917 |
| MD | F8-OLS | CHOL | 0.617 / 0.976 | 0.600 [0.449, 0.704] | 0.572 [0.129, 0.645] | 1.000 |
| MD | F8-OLS | PPARA | 0.406 / 0.506 | 0.358 [0.337, 0.931] | 0.092 [0.003, 0.200] | 0.714 |
| MD | F8-OLS | LIPID_any | 0.684 / 0.932 | 0.461 [0.304, 0.595] | 0.400 [0.134, 0.529] | 0.894 |
| MD | F8-OLS | OTHER | 0.301 / 0.310 | 0.057 [0.013, 0.102] | 0.0054 [0.0012, 0.0092] | 0.574 |
| MD | F8+C-Ridge | LIPO | 1.099 / 1.199 | 0.160 [0.061, 0.697] | 0.229 [0.038, 0.441] | 0.750 |
| MD | F8+C-Ridge | CHOL | 0.673 / 0.976 | 0.524 [0.410, 0.619] | 0.499 [0.128, 0.562] | 1.000 |
| MD | F8+C-Ridge | PPARA | 0.487 / 0.506 | 0.074 [−0.002, 0.703] | 0.019 [−0.000, 0.038] | 0.714 |
| MD | F8+C-Ridge | LIPID_any | 0.763 / 0.932 | 0.331 [0.114, 0.512] | 0.287 [0.049, 0.430] | 0.851 |
| MD | F8+C-Ridge | OTHER | 0.297 / 0.310 | 0.081 [0.045, 0.124] | 0.0078 [0.0047, 0.0108] | 0.589 |
| MD | FULL-OLS | LIPO | 0.837 / 1.199 | 0.512 [−5.041, 0.745] | 0.736 [−0.390, 1.115] | 0.833 |
| MD | FULL-OLS | CHOL | 0.379 / 0.976 | 0.849 [0.705, 0.911] | 0.809 [0.168, 0.911] | 0.905 |
| MD | FULL-OLS | PPARA | 0.485 / 0.506 | 0.083 [−0.831, 0.607] | 0.021 [−0.005, 0.046] | 0.571 |
| MD | FULL-OLS | LIPID_any | 0.560 / 0.932 | 0.639 [0.233, 0.813] | 0.556 [0.083, 0.774] | 0.787 |
| MD | FULL-OLS | OTHER | 0.297 / 0.310 | 0.078 [0.046, 0.121] | 0.0075 [0.0046, 0.0103] | 0.594 |

**Mean |log2FC|** (`results/group_metrics_all.csv`):

| Group | CD | MD |
|---|---|---|
| LIPO | 1.456 | 1.040 |
| CHOL | 0.773 | 0.896 |
| PPARA | 0.319 | 0.268 |
| LIPID_any | 0.812 | 0.746 |
| OTHER | 0.213 | 0.192 |

**R²g range across the 3 outer repeats, F8 and F8+C** (`results/group_metrics_per_repeat.csv`):

| Model | Diet | LIPO | CHOL | PPARA | OTHER |
|---|---|---|---|---|---|
| F8-OLS | CD | 0.173–0.224 | 0.625–0.681 | 0.134–0.202 | 0.051–0.053 |
| F8-OLS | MD | 0.307–0.342 | 0.592–0.613 | 0.345–0.370 | 0.056–0.058 |
| F8+C-Ridge | CD | 0.278–0.310 | 0.671–0.778 | 0.206–0.367 | 0.059–0.064 |
| F8+C-Ridge | MD | 0.114–0.183 | 0.449–0.578 | 0.021–0.178 | 0.080–0.082 |

## 2. Lipid group vs OTHER, F8 and F8+C (`results/group_vs_reference_comparisons.csv`)

Columns are Diff R²g [CI] (perm p), then Diff dMSE (module-FE) [CI] (FE p).

**CD**

| Model | Group | Diff R²g [CI] (perm p) | Diff dMSE, module FE [CI] (FE p) |
|---|---|---|---|
| F8-OLS | LIPO | +0.144 [0.015, 0.704] (0.448) | +0.513 [0.068, 0.958] (0.026) |
| F8-OLS | CHOL | +0.596 [0.544, 0.809] (0.003) | +0.470 [0.357, 0.583] (<0.0001) |
| F8-OLS | PPARA | +0.119 [−0.082, 0.444] (0.543) | +0.019 [−0.027, 0.066] (0.395) |
| F8-OLS | LIPID_any | +0.285 [0.052, 0.538] (0.039) | +0.345 [0.122, 0.568] (0.004) |
| F8+C-Ridge | LIPO | +0.235 [0.121, 0.686] (0.270) | +0.780 [0.421, 1.139] (0.0002) |
| F8+C-Ridge | CHOL | +0.674 [0.644, 0.804] (0.0015) | +0.533 [0.393, 0.672] (<0.0001) |
| F8+C-Ridge | PPARA | +0.237 [−0.147, 0.637] (0.294) | +0.035 [−0.035, 0.106] (0.309) |
| F8+C-Ridge | LIPID_any | +0.373 [0.154, 0.624] (0.032) | +0.446 [0.202, 0.690] (0.001) |

**MD**

| Model | Group | Diff R²g [CI] (perm p) | Diff dMSE, module FE [CI] (FE p) |
|---|---|---|---|
| F8-OLS | LIPO | +0.264 [−0.024, 0.649] (0.317) | +0.449 [0.248, 0.649] (0.0002) |
| F8-OLS | CHOL | +0.543 [0.393, 0.661] (0.025) | +0.557 [0.400, 0.714] (<0.0001) |
| F8-OLS | PPARA | +0.301 [0.249, 0.867] (0.209) | +0.082 [−0.047, 0.212] (0.199) |
| F8-OLS | LIPID_any | +0.404 [0.238, 0.547] (0.023) | +0.387 [0.182, 0.591] (0.0009) |
| F8+C-Ridge | LIPO | +0.078 [−0.030, 0.611] (0.711) | +0.219 [0.047, 0.391] (0.016) |
| F8+C-Ridge | CHOL | +0.443 [0.327, 0.554] (0.059) | +0.482 [0.352, 0.612] (<0.0001) |
| F8+C-Ridge | PPARA | −0.007 [−0.078, 0.625] (0.974) | +0.007 [−0.014, 0.028] (0.498) |
| F8+C-Ridge | LIPID_any | +0.249 [0.033, 0.424] (0.134) | +0.272 [0.053, 0.490] (0.017) |

**Other comparisons:**
- **Within-module permutation p for Diff dMSE:** 0.0005–0.001 for LIPO, CHOL and LIPID_any in all 4 model × diet cells, except F8+C-Ridge MD LIPO (0.0095). For PPARA it is 0.0225 (F8-OLS MD), 0.147, 0.283 and 0.423.
- **Joint regression** with all three indicators plus module fixed effects (`results/joint_overlap_regression_moduleFE.csv`): coefficients are within 0.001 of the single-indicator module-FE values (max difference 0.0008); the three sets share no genes.
- **Reference OTHER_noB12** (none of the 12 block-B sets): Diff R²g values are within 0.0072 of the OTHER values (all 5 models) (same file, `reference == "OTHER_noB12"`).

## 3. Lipid group vs OTHER genes matched on null error (k = 10; `results/group_vs_nullerror_matched_other.csv`)

**Matching balance:** median null squared error, group vs matched (`results/nullerror_matching_balance.csv`).

| Group | CD | MD |
|---|---|---|
| LIPO | 2.52 vs 1.90 | 0.914 vs 0.912 |
| CHOL | 0.546 vs 0.546 | 0.723 vs 0.684 |
| PPARA | 0.0716 vs 0.0716 | 0.0165 vs 0.0165 |
| LIPID_any | 0.426 vs 0.427 | 0.519 vs 0.518 |

**Results:**

| Diet | Model | Group | R²g group / matched | Diff R²g [CI] | Diff dMSE, module FE [CI] (FE p) |
|---|---|---|---|---|---|
| CD | F8-OLS | LIPO | 0.196 / 0.073 | +0.123 [−0.008, 0.671] | +0.406 [−0.008, 0.820] (0.054) |
| CD | F8-OLS | CHOL | 0.648 / 0.073 | +0.575 [0.503, 0.797] | +0.416 [0.224, 0.609] (0.0002) |
| CD | F8-OLS | PPARA | 0.171 / 0.048 | +0.123 [−0.112, 0.455] | +0.029 [−0.037, 0.095] (0.375) |
| CD | F8-OLS | LIPID_any | 0.337 / 0.074 | +0.263 [0.034, 0.518] | +0.229 [0.011, 0.447] (0.041) |
| CD | F8+C-Ridge | LIPO | 0.297 / 0.123 | +0.174 [0.029, 0.656] | +0.569 [0.119, 1.019] (0.016) |
| CD | F8+C-Ridge | CHOL | 0.737 / 0.093 | +0.644 [0.600, 0.787] | +0.466 [0.228, 0.704] (0.0006) |
| CD | F8+C-Ridge | PPARA | 0.300 / 0.066 | +0.234 [−0.169, 0.636] | +0.047 [−0.041, 0.136] (0.276) |
| CD | F8+C-Ridge | LIPID_any | 0.436 / 0.112 | +0.324 [0.106, 0.589] | +0.282 [0.006, 0.558] (0.045) |
| MD | F8-OLS | LIPO | 0.321 / 0.008 | +0.312 [0.012, 0.690] | +0.357 [0.063, 0.651] (0.020) |
| MD | F8-OLS | CHOL | 0.600 / 0.041 | +0.559 [0.394, 0.675] | +0.401 [0.208, 0.594] (0.0003) |
| MD | F8-OLS | PPARA | 0.358 / 0.021 | +0.336 [0.242, 0.976] | +0.058 [−0.088, 0.204] (0.418) |
| MD | F8-OLS | LIPID_any | 0.461 / 0.028 | +0.433 [0.257, 0.591] | +0.283 [0.100, 0.465] (0.004) |
| MD | F8+C-Ridge | LIPO | 0.160 / 0.034 | +0.126 [0.007, 0.663] | +0.203 [0.008, 0.397] (0.042) |
| MD | F8+C-Ridge | CHOL | 0.524 / 0.050 | +0.474 [0.349, 0.585] | +0.317 [0.186, 0.449] (0.0001) |
| MD | F8+C-Ridge | PPARA | 0.074 / 0.012 | +0.062 [−0.032, 0.705] | −0.003 [−0.040, 0.035] (0.892) |
| MD | F8+C-Ridge | LIPID_any | 0.331 / 0.043 | +0.288 [0.062, 0.476] | +0.164 [−0.007, 0.335] (0.060) |

**FULL-OLS vs matched OTHER, Diff R²g [CI]:**

| Group | CD | MD |
|---|---|---|
| LIPO | +0.562 [−2.984, 0.755] | +0.480 [−5.092, 0.712] |
| CHOL | +0.795 [0.742, 0.876] | +0.800 [0.636, 0.888] |
| PPARA | +0.500 [−0.014, 0.832] | +0.074 [−0.855, 0.618] |
| LIPID_any | +0.630 [0.261, 0.765] | +0.599 [0.183, 0.779] |

## 4. Model vs model, gene-paired (`results/model_pair_gain_lipid_vs_other.csv`)

**Gain** = MSE_B − MSE_A per gene. Each row gives the mean gain in the lipid group [CI], the mean gain in OTHER, and Diff (lipid − OTHER) [CI] (FE p).

**F8+C-Ridge vs F8-OLS**

| Diet | Group | Mean gain, lipid [CI] | Mean gain, OTHER | Diff [CI] (FE p) |
|---|---|---|---|---|
| CD | LIPO | +0.270 [0.038, 0.398] | +0.0011 | +0.269 [0.037, 0.397] (0.003) |
| CD | CHOL | +0.066 [0.002, 0.157] | +0.0011 | +0.065 [0.001, 0.157] (0.0005) |
| CD | PPARA | +0.019 [−0.016, 0.063] | +0.0011 | +0.018 [−0.018, 0.063] (0.409) |
| CD | LIPID_any | +0.104 [0.021, 0.193] | +0.0011 | +0.103 [0.020, 0.192] (0.018) |
| MD | LIPO | −0.231 [−0.325, 0.075] | +0.0023 | −0.234 [−0.328, 0.072] (0.003) |
| MD | CHOL | −0.073 [−0.085, −0.001] | +0.0023 | −0.075 [−0.088, −0.003] (<0.0001) |
| MD | PPARA | −0.073 [−0.161, 0.000] | +0.0023 | −0.075 [−0.163, −0.003] (0.166) |
| MD | LIPID_any | −0.113 [−0.186, −0.040] | +0.0023 | −0.116 [−0.187, −0.043] (0.004) |

**FULL-OLS vs F8-OLS, Diff [CI]**

| Group | CD | MD |
|---|---|---|
| LIPO | +1.220 [−0.700, 1.770] | +0.273 [−0.383, 0.493] |
| CHOL | +0.176 [0.006, 0.297] | +0.235 [0.029, 0.281] |
| PPARA | +0.054 [−0.015, 0.135] | −0.073 [−0.155, −0.006] |
| LIPID_any | +0.406 [−0.015, 0.844] | +0.153 [−0.059, 0.269] |

**FULL-OLS vs F8+C-Ridge, Diff [CI]**

| Group | CD | MD |
|---|---|---|
| LIPO | +0.951 [−0.737, 1.383] | +0.507 [−0.434, 0.767] |
| CHOL | +0.112 [0.001, 0.160] | +0.310 [0.041, 0.359] |
| PPARA | +0.037 [−0.004, 0.076] | +0.003 [−0.006, 0.009] |
| LIPID_any | +0.304 [−0.033, 0.656] | +0.269 [0.020, 0.424] |

## 5. Calculation details
- **Per-gene error:**
  - mse_g = mean over r = 1..3 of (y_g − ŷ_g,r)²;
  - mse0_g = mean over r of (y_g − ȳ_train,r(g))², where ȳ_train is the DL-weighted training mean of the split that held out g in repeat r;
  - d_g = mse0_g − mse_g.
- **Group metrics:** R²g = 1 − Σ mse_g / Σ mse0_g, and dMSE = mean d_g.
- **Bootstrap:** module counts are drawn as Multinomial(20, 1/20) and the group sums (n, Σmse, Σmse0, Σd) are reweighted. Resamples that contain no group genes are dropped (NaN-quantile).
- **Permutation:** within each module, group / reference labels are permuted among genes in group ∪ reference; p = (count of |null| ≥ |obs| + 1) / 2,001.
- **Module FE:** d_g = β·1[g ∈ group] + Σ_m α_m + ε, using genes in group ∪ reference. The CR1 variance is (G/(G−1))·((n−1)/(n−k))·(X′X)⁻¹ Σ_m X_m′e_m e_m′X_m (X′X)⁻¹ with G = 20, and the test uses t on 19 df.
- **Seed:** 20261125 for the bootstrap, permutations and matching (the matching is deterministic).
