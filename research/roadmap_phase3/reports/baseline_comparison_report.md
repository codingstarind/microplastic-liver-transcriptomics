# Baseline comparisons: PE effects vs independent references

**Script:** `scripts/f01_compare_and_baselines.py`. **Table:** `results/comparisons/baselines.csv` (every target × reference pair). **Figure:** `figures/baseline_comparisons.png`.

- **Targets:** existing GSE294180 DESeq2 effects (unshrunk MLE), not re-estimated: **PE_in_CD**, **PE_in_MD**, and **Diet_in_Veh** (the positive control).
- **Gene set:** the fixed 11,888-gene universe (sensitivity: 8,360).
- **DEG rule:** padj < 0.05 and |log2FC| ≥ 0.38.
- **Uncertainty:** the range over the 19 existing leave-one-mouse-out refits of the target ("LOO").

**These are baselines, not findings about PE.** Every PE effect remains aligned with sequencing depth and yield in GSE294180.

## Baseline 2: effect correlation (Spearman ρ, observed)

| Reference | Diet_in_Veh (positive control) | PE_in_CD [LOO range] | PE_in_MD [LOO range] |
|---|---|---|---|
| MASH1 (GSE316023) | **0.355** [0.30, 0.37] | 0.211 [**0.03**, 0.24] | 0.157 [**0.02**, 0.24] |
| MASH2 (GSE246328) | 0.323 [0.24, 0.36] | 0.027 [−0.18, 0.09] | 0.041 [−0.13, 0.15] |
| FAST (GSE107787) | 0.001 | 0.135 [0.10, 0.20] | **0.254** [0.22, 0.27] |
| FAST2 (GSE225695, check) | −0.064 | 0.082 | 0.082 |
| PPARα WT 20v0 | 0.043 | −0.048 [−0.07, −0.01] | −0.086 [−0.12, −0.05] |
| PPARα-dependent genes only (n = 1,893) | 0.160 | 0.016 | −0.080 |
| CIRC T1 (ZT6−ZT2) | −0.191 | −0.109 [−0.14, −0.04] | 0.106 [0.03, 0.16] |
| CIRC T2 (ZT12−ZT0) | −0.151 | −0.069 | 0.181 [0.13, 0.22] |
| *Historical same-data split-reference (PE_in_CD vs Diet_in_Veh)* | — | *0.33 (0.29–0.40)* | — |

**Observed points:**
- **The PE–MASH1 correlations depend on single mice.**
  - PE_in_CD vs MASH1 falls from 0.21 to **0.03** without CD_MP_2.
  - PE_in_MD vs MASH1 falls to **0.02** without MD_MP_4.
  - PE vs the second MASH study (MASH2) is ≈ 0.
  - So the MASH similarity of the PE effect is **not reproducible across the two MASH references and not robust to single mice**.
- **PE vs fasting is small but stable** under LOO: 0.10–0.20 (CD) and 0.22–0.27 (MD).
- **PE vs the PPARα agonist is slightly negative** in both diets.
- The independent-reference PE_in_CD–MASH correlation (0.21) is lower than the historical same-data split-reference value (0.33), which shared mice with its target.

## Baseline 1: DEG overlap (observed)

| Reference | Target | DEGs (target / ref) | Overlap | Jaccard | Concordant / discordant |
|---|---|---|---|---|---|
| MASH1 | Diet_in_Veh | 1,363 / 4,777 | 747 | 0.139 | 631 / 116 |
| MASH1 | PE_in_CD | 430 / 4,777 | 201 | 0.040 | 134 / 67 |
| MASH2 | PE_in_CD | 430 / 4,297 | 209 | 0.046 | **84 / 125** |
| FAST | PE_in_CD | 430 / 4,009 | 290 | 0.070 | 211 / 79 |
| FAST | PE_in_MD | 190 / 4,009 | 122 | 0.030 | 112 / 10 |
| PPARα 20v0 | PE_in_CD | 427 / 3,588 | 198 | 0.052 | 95 / 103 |
| CIRC T2 | PE_in_CD | 430 / 1,369 | 109 | 0.064 | 19 / 90 |

- **Almost every overlap is "significant"** by the hypergeometric test (p < 0.01), **including overlaps where most shared genes change in opposite directions** (MASH2: 125 of 209 discordant; CIRC T2: 90 of 109 discordant).
- This shows the main weakness of DEG overlap. Responsive, high-variance genes appear in many DEG lists, so overlap enrichment says "both conditions perturb some of the same genes", not "the states are similar".
- Overlap also depends on power: the PE lists (n = 4–5 mice) are much shorter than the reference lists (n = 8–18).
- Jaccard values are all ≤ 0.14, even for the positive control.

## Comparison of the two baselines
- **Baseline 1 (DEG overlap)** rewards shared responsiveness, ignores direction unless forced, and depends on thresholds and power.
- **Baseline 2 (whole-universe correlation)** uses direction and every gene, but is dominated by the many small effects. It is fragile when one target mouse carries much of the signal.
- Neither has a valid gene-level null, because genes are correlated. The LOO range is the most honest uncertainty available here.
- **No method is ranked by its largest correlation.** These values are the benchmark that a Phase 4 model would have to improve on, with the same targets, universe and uncertainty treatment.
