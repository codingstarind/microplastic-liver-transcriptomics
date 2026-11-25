# Roadmap Phase 5: external validation of frozen models on independent microplastic studies (2026-11-23)

## Main sources (`research/roadmap_phase5/`)
- `locked_plan.md` (LOCKED before any external effect was computed; `results/plan_lock.json`)
- `phase5_decision_report.md`, `validation_diagnostics.md`
- `dataset_provenance.csv`, `external_contrasts.csv`
- `gene_level_external_results.csv`, `sample_level_external_results.csv`, `per_sample_scores.csv`, `pathway_concordance.csv`
- `frozen/frozen_manifest.json` (SHA-256 of every frozen artefact)
- `figures/` (fig1–fig4)
- Scripts: `s00_lock_plan.py`, `s01_prepare_external.py`, `s02_freeze_models.py`, `s03_external_de.R`, `s04_gene_level.py`, `s05_sample_level.py`, `s06_pathways_diagnostics.py`, `s07_figures.py`, `p5_common.py`
- Tests: `tests/test_phase5.py` (10 tests)

## 1. Questions
- **Q1 (gene level):** do the frozen gene-level models predict an independent study's plastic-vs-control log2FC vector?
- **Q2 (mouse level):** do the frozen 4D signatures and 4E classifiers move in the same direction in independent plastic-exposed mice? Is any shift fasting/cholesterol-like or MASH-like?

## 2. Data (`dataset_provenance.csv`)

| Unit | Mice | Input | Role |
|---|---|---|---|
| **HONG** (GSE261893 PE + GSE261725 PS + GSE286116 PP + GSE262399 PB) | **15** (24 GSMs); 3 controls byte-identical across the 4 series | integer counts (StringTie prepDE) | primary PE test |
| GSE245069 (Li et al. 2025): PS 80 nm / 0.5 µm / 5 µm, 1 mg/L in water, 12 wk, ICR | 16 (4 shared water controls) | estimated counts | secondary, cross-polymer |
| GSE277772: PS, oral, 8 wk | 6 | normalised values only | exploratory |

**Comparability with GSE294180:** Hong PE was given in food for 71 weeks, with size and dose not stated, on grain chow. GSE294180 used 10–150 µm PE by gavage for 8 weeks on a purified diet.

## 3. Methods
**Frozen artefacts:**
- F8-OLS (PE_in_CD primary; PE_in_MD secondary);
- F8+C-Ridge (α = 1000);
- FAST-only and MASH-only OLS;
- the 4D singscores;
- the 4E R0/R2 L2-logistic models, refitted deterministically and reproducing the saved scores to 4e-16.

**External effects:** DESeq2 Wald MLE with the same prefilter rule. Hong was fitted as one 5-group model so the shared controls are modelled once.

**Gene-level metrics:**
- weighted R² vs zero, vs the training mean and vs the study's own mean;
- Spearman with a module-block bootstrap CI;
- sign agreement;
- calibration slope;
- noise ceiling;
- **mouse-label permutation** over all relabellings (20 for 3 v 3), using a log2CPM mean-difference proxy.

**Mouse-level metrics:**
- Welch CI and Hedges g;
- exact permutation;
- leave-one-mouse-out (LOO) signs;
- correlation with technical variables;
- R2 decomposition;
- a coverage control that re-scores GSE294180 on each external gene background.

**Pathways:** camera on the 12-set panel, and concordance with GSE294180, the fasting reference and the MASH reference.

## 4. Results

**Q1, Hong PE vs CON** (`gene_level_external_results.csv`, PE_in_CD-trained models):

| Model | weighted R² vs own mean | Spearman [gene-block CI] | perm p (20 labellings) | Locked verdict |
|---|---|---|---|---|
| **F8-OLS** | **0.041** | 0.22 [0.14, 0.28] | **0.25** | **not passed** |
| F8+C-Ridge | 0.004 | 0.14 [0.08, 0.20] | 0.30 | same |
| FAST-only OLS | 0.045 | 0.29 [0.23, 0.34] | 0.05 (minimum) | PASS-direction |
| Observed GSE294180 vector | −0.053 | 0.10 | 0.25 | — |

- **Specificity:** the same predictions match PP vs CON better (ρ 0.39). PB is similar.
- **GSE245069 pooled PS:** F8-OLS ρ = −0.07, R² −0.14 → **FAIL**.
- **Signal share:** Hong PE 0.47; PS pooled 0.65.

**Q2, FAST singscore** (`sample_level_external_results.csv`):
- **Hong PE − CON = +0.060 [−0.004, 0.124]**; Hedges g 2.5; exact p = 0.10 (the minimum possible for 3 v 3); LOO 6/6.
- GSE294180 PE_in_CD was +0.065.
- Locked verdict: **DIRECTIONALLY CONSISTENT**.
- **Not PE-specific:**
  - PP +0.088, PB +0.038, PS +0.024;
  - pooled non-PE polymers vs CON +0.050 (exact p 0.04);
  - the shared controls have the lowest FAST mean of all 5 groups;
  - control-free PE minus the other polymers is +0.010.

**Other scores:**
- **MASH singscore:** Hong PE −0.021; R0 −0.56 log-odds. Both are opposite in sign to GSE294180.
- **R2:** Hong PE +0.94 (exact p 0.30), 69% cholesterol/fasting-like.
- **GSE245069 PS pooled:**
  - FAST +0.036 (p 0.24);
  - **PPARα singscore +0.053 (exact p 0.001; 16/16 LOO; strongest at 80 nm)**;
  - R2 +1.93, 70% cholesterol/PPARα-like;
  - MASH ≈ 0.
- **GSE277772 (exploratory):** FAST −0.013; MASH −0.027 (opposite in sign).

**Pathways (Hong PE):**
- Lipogenesis down (camera −0.82) and cholesterol synthesis down (−0.42), the same direction as GSE294180.
- Across the 12-set panel, ρ = 0.66 with GSE294180 PE_in_CD and 0.55 with fasting.
- Hong PS/PP/PB also lower lipogenesis (−0.46 to −0.71).

**Controls:**
- External depth is balanced (exposed/control ratio 0.95–1.0).
- Coverage control: restricting the gene background changes the GSE294180 estimates by ≤ 9%.

## 5. Corrections recorded in this phase
- Phase 4E's "every non-PCA reached AUC 1.0" is false (trees 0.63–0.94).
- The R0 "Brier" is a hard-label misclassification rate.

## 6. Decision
- **Gene-level external-validation claims: NO-GO.**
- **Narrow GO for write-up:** a fasting/cholesterol-type shift that is directionally consistent in Hong PE but not PE-specific.
- No PE-specific, causal or MASH claim.

## 7. Transition
The team asked why transfer failed and whether it could be fixed, which led to Phase 5B.
- Was the failure a pipeline error?
- What is the true gap?
- Do principled fixes help?
