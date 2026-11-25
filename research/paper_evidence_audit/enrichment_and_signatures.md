# Enrichment, pathway and biological-signature evidence (P0 to Phase 5B)

**Enrichment and signature similarity describe expression patterns.** They are **not** demonstrated mechanisms, and they do not measure pathway activity, protein, or flux. Every PE result from GSE294180 is aligned with lower sequencing depth and yield (PE ≡ depth).

---

## 1. Gene-set resources

| Resource | Version / definition | Used in | Source file |
|---|---|---|---|
| MSigDB mouse Hallmark (50 sets) | MSigDB 2026.1.Mm via msigdbr | P0 (primary family), P3, P4B (9 sets as block B), P4E (R1), P5 (panel) | `research/roadmap_phase3/results/gene_sets.csv` |
| PPAR sets | SANDERSON_PPARA_TARGETS, WP_PPAR_SIGNALING_PATHWAY, REACTOME regulation of lipid metabolism by PPARα | P0 (separate family), P3 | same |
| Custom sets (defined before any Phase 3 fit) | LIPOGENESIS_CORE (12 genes), CHOL_SYNTHESIS_CORE (21 genes) | P3, P4B, P4E, P5 | same |
| Fixed 12-set panel ("SETS12") | LIPOGENESIS_CORE, CHOL_SYNTHESIS_CORE, SANDERSON_PPARA_TARGETS + 9 Hallmark sets (fatty-acid metabolism, cholesterol homeostasis, inflammatory response, EMT, xenobiotic, bile acid, peroxisome, OXPHOS, MTORC1) | P4B block B, P4E R1, P5 pathway concordance | `roadmap_phase5/scripts/p5_common.py` (`SETS12`) |
| Frozen signatures (Phase 4D) | MASH_main_n50, FAST_main_n50, PPARA_main_n50 (50 up + 50 down each) | P4D, P5, P5B | `research/roadmap_phase4d/signatures/*_main_n50.csv` + `frozen_signature_manifest.csv` |
| Reactome | – | P0 exploratory only | – |

## 2. Methods: which question each one answers

| Method | Question | Unit | Correction | Where used |
|---|---|---|---|---|
| limma **fry** (self-contained rotation, directional) | Are this set's genes, taken together, changed? | mouse | BH within family | P0 primary; P3 |
| limma **camera** (competitive, inter-gene correlation adjusted; inter.gene.cor = 0.01 in P5) | Is this set more changed than other genes? | gene, with correlation adjustment | BH within contrast | P0 supporting; P3; P5 |
| fgsea (preranked on the Wald statistic) | Rank enrichment | gene (anti-conservative) | BH | P0 supporting only |
| Mean-z per mouse + two-way linear model | Effect size of the set score with a CI | mouse | – | P0 (fig16) |
| Hypergeometric DEG overlap | Do DEG lists share genes? | gene | descriptive | P3 Baseline 1 (shown to be uninformative) |
| singscore on frozen signatures | Per-mouse signature score → mouse-level contrasts | mouse | Holm (6 primary tests) | P4D, P5 |
| Classifier feature decomposition | Which pathway features move a frozen classifier score | mouse | descriptive | P4E, P5 |

> **Gene-level enrichment vs per-mouse signature scoring.**
> - Enrichment (fry/camera/fgsea) tests whether a set's genes changed **between groups**.
> - Per-mouse singscore gives each mouse a score. The 2 × 2 model then tests group differences with the **mouse as the replicate** (n = 19), with LOO and exact permutation.
> - The two agree qualitatively here, but they are different estimands.

---

## 3. Key pathway findings in GSE294180

### P0 (`reports/phase2/pathway_report.md`; `reports/final_report_draft.md` §4.4)

| Programme | PE under CD | PE under MD | Diet | Interaction | Robustness (`reports/phase2/robustness_report.md`) |
|---|---|---|---|---|---|
| Cholesterol homeostasis | **Down** (fry p 0.004, FDR 0.11; camera FDR < 0.001; NES −2.3) | **Down** (fry p 0.005, FDR 0.09; NES −2.6) | Down (FDR < 0.001) | none (fry p 0.84; mean-z −0.003 [−0.41, 0.40]) | 19/19 LOO; **lost under the yield covariate** |
| PPARα targets (Sanderson) | **Up (fry FDR 0.041)** | weaker (fry p 0.35) | Up | camera FDR 0.003; fry p 0.21 (exploratory) | 19/19 LOO same direction; **lost under depth/yield** |
| Lipogenesis genes (Fasn, Acly, Acaca) | down (−0.7 to −1.8 log2) | down | – | – | 19/19 LOO; 2,000/2,000 bootstrap signs |
| Inflammation (4 sets), xenobiotic, EMT/fibrosis | n.s. | n.s. | mixed | none by fry | – |
| Pathway-level interaction (53 sets) | – | – | – | **none at FDR < 0.05 (fry)** | 18/19 LOO |

### Phase 3 reference programmes (median log2FC)
Source: `research/roadmap_phase3/reports/reference_comparison_report.md`.

| Programme | Fasting (GSE107787) | GW7647 (PPARα agonist) | MASH1 (GAN vs chow) | Our Diet_in_Veh |
|---|---|---|---|---|
| Lipogenesis core | **−2.79** | **+1.15** | **+0.44 (fry ↑, FDR 0.011)** | **−1.43 (fry ↓, FDR 6e-6)** |
| Cholesterol synthesis core | −2.76 | +0.67 | −1.99 | ↓ (agrees) |
| PPARα targets | +0.82 | +1.50 | +0.14 | ↑ (agrees) |

**Interpretive keys:**
- "PPARα targets ↑ with lipogenesis and cholesterol synthesis ↓" is a **fasting-like** combination. It is **not** what a PPARα agonist produces: the agonist raises lipogenesis and cholesterol synthesis.
- Lipogenesis is opposite between the MASH-vs-**chow** references and our MASH-vs-**purified** diet, so lipogenesis cannot be used to judge "MASH-likeness".

### Phase 4C: residuals by pathway
Source: `tables/subgroup_errors.csv`.
- LIPOGENESIS_CORE mean residual: −1.24 (CD), −0.74 (MD).
- CHOL_SYNTHESIS_CORE: −0.42 / −0.53.
- The PE decrease in these programmes is **stronger than any independent reference predicts**.
- Block B (Phase 4B) gained performance by re-learning this offset.

### Phase 4D: per-mouse signatures
Source: `research/roadmap_phase4d/results/contrast_estimates_primary.csv`. Estimates with 95% t CI; n = 19.

| Score | PE_in_CD | PE_in_MD | Interaction | Diet_in_Veh |
|---|---|---|---|---|
| **FAST** | **+0.065 [0.022, 0.108]**, Holm 0.036, exact perm p 0.040 | **+0.056 [0.015, 0.096]**, Holm 0.056, perm 0.024 | −0.010 | +0.130 |
| MASH | +0.025 [−0.033, 0.082] | +0.054 [−0.0004, 0.108] | +0.029 | **+0.215 [0.158, 0.273]** (positive control PASS) |
| PPARα (secondary) | +0.031 | +0.006 | – | – |

**Sensitivities:**
- **The FAST shift is robust** to input (CPM/VST), signature size (n25: +0.094 / +0.077) and the muscle rule.
- **The MASH PE shift is fragile.** It vanishes without cholesterol-synthesis genes (+0.003 / +0.021) or without MASH–FAST shared genes.
- The Holm-adjusted exact-permutation p-values are 0.20 / 0.14 (correction recorded in Phase 4E).

**Technical:**
- clock_phase ρ: MASH 0.64, FAST 0.55;
- CD+Veh has the lowest FAST score (≈ −0.11);
- depth/yield/muscle |ρ| ≤ 0.30.

### Phase 4E: classifier decomposition
Source: `results/r2_feature_contributions_gse294180.csv`.
- The R2 "MASH-state" PE_in_MD shift (+1.66 log-odds, Holm 0.0006) is carried by cholesterol homeostasis (+1.17), cholesterol synthesis (+0.93), PPARα targets (+0.35) and MTORC1 (+0.28). Lipogenesis is −1.37; the MASH programme feature contributes **+0.18**.
- R2 scores correlate 0.89 with the 4D FAST score.

---

## 4. Cross-study pathway evidence (Phase 5, fixed 12-set panel)
Source: `research/roadmap_phase5/pathway_concordance.csv`, panel `primary12`. Values are median log2FC with camera direction and p.

| Set | GSE294180 PE_in_CD | REF fasting | REF MASH1 | **Hong PE** | Hong PS | Hong PP | Hong PB | GSE245069 PS pooled | GSE277772 PS* |
|---|---|---|---|---|---|---|---|---|---|
| LIPOGENESIS_CORE | −1.59 | −2.79 | +0.44 | **−0.82 (Down, p 0.006)** | −0.71 (p < 0.001) | −0.52 (p < 0.001) | −0.46 (p 0.12) | +0.06 (Up, p 0.023) | +0.15 |
| CHOL_SYNTHESIS_CORE | −0.74 | −2.76 | −1.99 | **−0.42 (Down, p 0.001)** | −0.36 (p 0.001) | −0.09 (p 0.10) | +0.16 (p 0.82) | −0.01 | +0.18 |
| SANDERSON_PPARA_TARGETS | +0.23 | +0.82 | +0.14 | **−0.10 (Down, p 0.26)** | −0.24 (Down, p 0.001) | −0.10 | −0.03 | **+0.41 (Up, p < 0.001)** | −0.14 |
| HALLMARK_CHOLESTEROL_HOMEOSTASIS | −0.21 | −0.66 | −0.32 | −0.22 (p < 0.001) | −0.15 | −0.04 | −0.12 | −0.13 (p 0.007) | −0.05 |
| HALLMARK_INFLAMMATORY_RESPONSE | +0.04 | −0.09 | +0.53 | +0.10 (n.s.) | +0.07 | +0.20 (p 0.047) | +0.01 | +0.03 | +0.04 |
| HALLMARK_EMT | +0.02 | −0.15 | +0.57 | −0.11 (n.s.) | −0.01 | +0.09 | −0.21 | −0.20 | −0.06 |
| HALLMARK_OXPHOS | −0.06 | −0.11 | −0.17 | −0.28 (p < 0.001) | −0.10 | −0.25 | −0.16 | +0.00 | +0.00 |

\* GSE277772: normalised input, exploratory only.

**Observed patterns** (data, not mechanism):
- **Lipogenesis and cholesterol synthesis are lower** in Hong PE, the same direction as GSE294180 PE_in_CD and the fasting reference. **Hong PS and PP show the same**, against the same 3 shared controls.
- **The PPARα-target direction is not consistent:**
  - up in GSE294180 PE_in_CD (+0.23) and GSE245069 PS (+0.41, camera p < 0.001; PPARα singscore +0.053, exact p 0.001);
  - **down** in Hong PE (−0.10) and Hong PS (−0.24).
- **Inflammation and EMT (MASH-like programmes):** not up in Hong PE or GSE245069. There is no MASH-type pathway shift in any external study.
- **Hong shows broad down-shifts** across many Hallmark sets (OXPHOS, peroxisome, xenobiotic, MTORC1) for all four polymers vs the same controls. This is consistent with an atypical control group (Phase 5: the controls have the lowest FAST mean) or a generic effect. **Not resolvable.**
- **Panel-level concordance** (Spearman over 12 sets, Phase 5 decision report): Hong PE vs GSE294180 PE_in_CD ρ = 0.66; vs fasting 0.55. These are 12 correlated sets, so the numbers are descriptive.

### Mouse-level frozen scores, externally
Source: `research/roadmap_phase5/sample_level_external_results.csv`.

| Score | Study / contrast | Estimate [95% CI] | Exact p (min) | LOO | GSE294180 PE_in_CD |
|---|---|---|---|---|---|
| FAST | Hong PE vs CON (3 v 3) | **+0.060 [−0.004, 0.124]**, g 2.5 | 0.10 (0.10) | 6/6 | +0.065 |
| FAST | Hong PP / PB / PS vs CON | +0.088 / +0.038 / +0.024 | – | – | – |
| FAST | Hong non-PE pooled vs CON | +0.050 | 0.04 | – | – |
| FAST | GSE245069 PS pooled (12 v 4) | +0.036 | 0.24 | – | – |
| PPARα singscore | GSE245069 PS pooled | **+0.053** | **0.001** | 16/16 | +0.031 |
| MASH singscore | Hong PE | −0.021 (opposite sign) | – | – | +0.025 |
| R2 logistic | Hong PE | +0.94 (69% cholesterol/fasting-like) | 0.30 | – | +0.83 |

## 5. Robustness classification

| Finding | Status | Why |
|---|---|---|
| Cholesterol programme lower with PE in GSE294180 (both diets) | **Robust within the dataset** (LOO, methods); **not separable from depth/yield** | P0 robustness matrix |
| Lipogenesis genes lower with PE in GSE294180 | Robust within the dataset; **direction conflicts with Djouina/Okamura**; agrees with Hong (all polymers) | P0 replication panel; P5 |
| PPARα targets up with PE under CD | **Fragile:** lost under depth/yield; opposite in Hong; up in GSE245069 PS | P0, P5 |
| FAST-like per-mouse shift | **Most reproducible:** P4D both diets; same direction and size in Hong PE | P4D, P5 |
| FAST-like shift is PE-specific | **Not supported:** present for PS/PP/PB vs shared controls | P5 |
| MASH-like PE shift | **Not supported:** fragile in P4D; opposite in Hong; classifier shifts are cholesterol-driven | P4D, P4E, P5 |
| Pathway-level diet × PE interaction | **None detected** (fry) | P0 |
| Gene-level interactions | Mostly **muscle artefact** | P0 |

## 6. Sensitivity to filtering, muscle and technical confounders

**Muscle exclusion:**
- The 309-gene rule removes liver-relevant genes (Pdk4, Plin4, Col1a1/2, Lpl).
- The exploratory 88-gene rule leaves the positive control unchanged (0.355 → 0.358).
- The P4D FAST/MASH results are identical under the alternative rule.

**Universe:**
- The strict 8,360-gene universe lowers the MASH positive control to ρ 0.252.
- In Phase 4 the strict universe **raises** PE_in_MD R² (ols_all 0.180 vs 0.070; `taskA_variant_stability.png`).

**Depth/yield:** covariates remove the P0 PPARα and cholesterol pathway signals, but these covariates are almost perfectly collinear with PE, so this neither confirms nor refutes a biological effect.

**Clock:** the CD+Veh clock state correlates with the signature scores (ρ 0.55–0.67). PE_in_CD estimates partly reflect the unusual reference group; PE_in_MD estimates do not depend on CD+Veh.
