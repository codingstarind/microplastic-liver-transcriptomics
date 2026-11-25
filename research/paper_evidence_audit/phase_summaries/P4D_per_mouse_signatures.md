# Roadmap Phase 4D: per-mouse signature scores in the 2 × 2 design (2026-11-23)

## Main sources (`research/roadmap_phase4d/`)
- `analysis_plan.md` (LOCKED; `results/plan_lock.json`)
- `representation_go_nogo.md`, `source_inventory.md`, `decision_report.md`, `phase5_readiness.md`
- `signatures/*_main_n50.csv` and `frozen_signature_manifest.csv`
- `results/contrast_estimates_primary.csv`, `contrast_estimates_sensitivity.csv`, `permutation_exact.csv`, `loo_contrasts.csv`, `loo_summary.csv`, `null_gene_sets.csv`, `qc_associations.csv`, `per_mouse_scores.csv`
- `figures/` (4 PNGs)
- Scripts: `scripts/` (3 files; `s01_build_signatures.py` builds the signatures)

## 1. Question
Do frozen, independently derived MASH and fasting signatures shift **per mouse** with PE in GSE294180, with the mouse as the unit of analysis?

## 2. Why
- Phases 4–4C were gene-level. Their rows are not independent, and they cannot address inference about mice.
- A per-mouse score creates an estimand with the correct experimental unit.

## 3. Methods

**Signatures** (from Phase 3 references only; frozen before any mouse was scored):
- **MASH_main_n50:** 50 up + 50 down genes. Each must agree in **both** MASH studies (padj < 0.05, |log2FC| ≥ 1).
- **FAST_main_n50:** each gene must also agree in the muscle-clean ZT0/ZT4 subset.
- **PPARA_main_n50:** secondary.

**Scoring:**
- singscore (Foroutan 2018), re-implemented in Python from the published formula.
- Input: TPM-proxy = count / vM35 median transcript length.
- Background: 13,777 filtered genes.

**Model:**
- OLS per score: score ~ MD × PE.
- Contrasts: PE_in_CD, PE_in_MD, Interaction, Diet_in_Veh.
- 95% t intervals (df 15).
- Holm correction over the 6 primary tests ({MASH, FAST} × 3 PE contrasts).

**Controls:**
- exact within-diet permutation: C(9,5) = 126 (CD) and C(10,5) = 252 (MD) relabellings;
- 19 leave-one-mouse-out (LOO) refits;
- 1,000 abundance-matched random signed null sets;
- QC correlations with depth, yield, muscle and clock_phase.

**Sensitivities:**
- CPM and VST inputs;
- n25 signatures;
- alternative muscle rule;
- MASH without CHOL_SYNTHESIS_CORE genes;
- 14 MASH/FAST shared genes removed.

## 4. Results (`results/contrast_estimates_primary.csv`)

| Score | Contrast | Estimate [95% CI] | p | Holm | Exact perm p | LOO sign stable |
|---|---|---|---|---|---|---|
| MASH | Diet_in_Veh (positive control) | **+0.215 [0.158, 0.273]** | < 1e-4 | single test | — | 19/19 |
| MASH | PE_in_CD | +0.025 [−0.033, 0.082] | 0.37 | 1.00 | 0.14 | 19/19 |
| MASH | PE_in_MD | +0.054 [−0.0004, 0.108] | 0.052 | 0.21 | 0.15 | 19/19 |
| MASH | Interaction | +0.029 [−0.050, 0.109] | 0.44 | 1.00 | 0.49 | 19/19 |
| **FAST** | **PE_in_CD** | **+0.065 [0.022, 0.108]** | 0.006 | **0.036** | 0.040 | 19/19 |
| **FAST** | **PE_in_MD** | **+0.056 [0.015, 0.096]** | 0.011 | 0.056 | 0.024 | 19/19 |
| FAST | Interaction | −0.010 [−0.069, 0.050] | 0.74 | 1.00 | 0.82 | 17/19 |
| FAST | Diet_in_Veh | +0.130 [0.087, 0.173] | < 1e-4 | descriptive | — | 19/19 |
| PPARα | PE_in_CD / PE_in_MD | +0.031 / +0.006 | 0.12 / 0.72 | 0.36 / 0.72 | 0.19 / 0.66 | 19/19 / 17/19 |

> **Correction recorded in Phase 4E (`00_task_and_estimands.md`).** The exact permutation p-values (0.040 / 0.024) were not Holm-adjusted in Phase 4D. Adjusted over the same 6-test family they become **0.20 / 0.14**.

**Null gene sets.** Observed shifts lie outside the 95% null range, e.g. FAST PE_in_CD +0.065 vs null −0.011 to +0.011. This shows specificity across genes. It is not replication, because the mouse is still the unit.

**FAST is robust:**

| Variant | PE_in_CD | PE_in_MD |
|---|---|---|
| Original (n50) | +0.065 | +0.056 |
| n25 | +0.094 | +0.077 |
| CPM / VST input | +0.067 | +0.060 |
| Shared genes removed | +0.054 | +0.031 (CI includes 0) |

**MASH is fragile:**
- Without CHOL_SYNTHESIS_CORE genes, PE_in_CD is +0.003 and PE_in_MD +0.021 (both CIs include 0).
- Without the MASH/FAST shared genes: +0.010 / +0.030.

**QC:**
- Depth separates PE completely (log2 depth 23.69–23.98 vs 24.02–24.19).
- Scores vs depth, yield and muscle: |ρ| ≤ 0.30.
- Scores vs clock_phase: ρ 0.64 (MASH) and 0.55 (FAST).
- CD+Veh has the lowest FAST scores.

## 5. Interpretation boundary
- **What the data show:** a fasting-like signature shift *associated* with PE in both diets.
- **What they cannot show:** attribution to PE (PE ≡ depth), or a clean CD-diet estimate (the CD+Veh reference group has an unusual clock state).
- **The MASH-like shift** is carried by cholesterol/fasting-shared genes.
- GSVA was not run (not installed).

## 6. Transition
The next question was a MASH-state model validated **across studies**. Would such a model read PE mice as "MASH-like"? This is Phase 4E.
