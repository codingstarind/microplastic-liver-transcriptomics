# Phase 4 specification memo (for approval; nothing here has been fitted)

**Date:** 2026-11-22. **Purpose:** fix every Phase 4 modelling choice **before** any multi-axis model is fitted to PE effects. No threshold below was chosen by looking at PE fit. All Phase 3 baseline numbers are already reported, and none of them is used to set anything here.

## 1. Fixed targets (outcomes)
- **y_g:** GSE294180 DESeq2 unshrunk MLE log2FC with lfcSE s_g, from `results/deseq2/provisional_deposited/` (`~ diet * microplastic`, deposited counts, no offset).
  - **Primary targets:** `PE_in_CD`, `PE_in_MD`, fitted separately.
  - **Positive-control target:** `Diet_in_Veh`. It is run first and must behave as in §7 before any PE result is read.
- **Mouse-level uncertainty:** the 19 existing leave-one-mouse-out refits (`results/deseq2/loo/`), used unchanged.
- **Genes:** the main universe (11,888; locked 309 muscle exclusion). Sensitivity universes: 8,360 (strict); 12,010 (exploratory muscle-restricted, flagged post-hoc).

## 2. Reference-axis effect vectors (independent mice only)

| Axis | Vector | Status |
|---|---|---|
| MASH | MASH1 GSE316023 `log2FC`; MASH2 as a **replacement** sensitivity (never both in one model) | **Exploratory** (positive control FAIL on c1) |
| Fasting | FAST GSE107787 ZT-adjusted `log2FC` (ZT-averaged) | Stable; time-dependent |
| PPARα | PPARA WT 20v0 `log2FC` (agonist-associated) | Robust; targeted assay |
| Time of day | **CIRC T1 = ZT6 − ZT2** only; T2 as a replacement sensitivity (ρ 0.86 with T1) | Phase information; female |

**Scale harmonisation:** each reference vector is z-scored over the universe genes (mean 0, SD 1). Coefficients then mean "change in PE log2FC per SD of reference effect". **Sensitivity:** rank-based inverse-normal scores, which are robust to the technology-specific tails of RPKM and TempO-Seq. Magnitudes are never compared across technologies.

## 3. Candidate models (complete list; no search)
For each target T ∈ {PE_in_CD, PE_in_MD, Diet_in_Veh}:

| Model | Terms |
|---|---|
| M0 | intercept only |
| M1 | MASH |
| M2 | FAST |
| M3 | PPARA |
| M4 | CIRC_T1 |
| M5 | MASH + FAST + PPARA |
| M6 | M5 + CIRC_T1 |

There is no interaction, no variable selection and no penalisation tuning on PE.

## 4. Weighting
- WLS with w_g = 1 / (s_g² + τ²).
- τ² is the between-gene excess variance, estimated **once per target** from M0 by the method of moments (DerSimonian–Laird form). It is then **held fixed for all models of that target**, so every model is compared on identical weights.
- This stops very precise, high-expression genes from dominating, while still down-weighting noisy genes.
- **Sensitivity:** unweighted OLS.

## 5. Noisy reference coefficients (errors in variables)
- Reference effects are estimates, so regressing on them attenuates the coefficients.
- **Primary:** report the naive WLS coefficients **and** reliability-corrected coefficients:
  - for a single-axis model, β / λ_k, with λ_k = 1 − mean(SE_k²) / var(x_k) computed on the universe;
  - for a multi-axis model, a classical measurement-error correction using the diagonal error covariance.
- The corrected values are reported as **sensitivity**, not as the headline.
- No reference gene filtering by significance: that would be selection.

## 6. Dependence among genes and validation design
- Genes are not independent (co-regulation), so the nominal regression SEs and AIC/BIC are invalid. Use:
  - **Module-blocked cross-validation:**
    - Genes are partitioned into **20 co-expression modules**, using hierarchical clustering (average linkage, 1 − Pearson) of log2(RPKM + 1) for the 54 GSE117134 liver samples (all genotypes and ZTs). This expression source is independent of the PE target and of the MASH, fasting and PPARα references.
    - The modules are a **partition**: each gene belongs to exactly one, so overlapping gene-set membership cannot leak a gene into both training and test folds.
    - Leave-one-module-out; module construction is seeded (20261122) and fixed before any fit.
  - **Known residual leakage:** co-regulation between modules, and the CIRC axis coming from the same dataset as the modules. Both are stated in the results.
  - **Cluster-robust (module) sandwich SEs** for the coefficients.
  - **Mouse-level sensitivity:** refit every model on each of the 19 LOO target vectors and report the ranges.
  - **Why not gene-set-blocked CV with Hallmark sets:** the sets overlap (a gene can be in several), which leaks information between folds. If used, it is a secondary analysis with each gene assigned to one set by a fixed priority order (core sets first, then alphabetical), with the overlap rate reported.

## 7. Evaluation metrics and decision rules (fixed now)
- **Primary metric:** held-out weighted R² of the predicted y_g (module-blocked CV), and its change relative to M0.
- **Secondary:** held-out weighted Spearman; held-out sign accuracy on target DEGs (padj < 0.05, |log2FC| ≥ 0.38).
- **Comparable baselines on the same held-out genes:**
  - **Baseline 2** = single-axis models M1–M4, which are the regression form of the Phase 3 correlations;
  - **Baseline 1** = a DEG-overlap predictor: predicted sign = the reference sign if the gene is a reference DEG, otherwise 0. It is scored by sign accuracy only, because it makes no magnitude prediction.
- **"Multi-axis is better"** only if M5 (or M6) beats the best single-axis model in held-out weighted R² in **both** the main and strict universes **and** in ≥ 16/19 LOO target refits.
- **Joint coefficients** are interpreted only if all VIFs ≤ 5 (observed Phase 3 VIFs ≈ 1.0–1.1 for MASH, FAST and PPARA).
- **Positive control first:**
  - For Diet_in_Veh, the MASH coefficient in M1 and M5 must be positive, with a module-robust 95% CI excluding 0, and M1 (MASH) must have the highest single-axis held-out R².
  - If not, MASH-axis results for PE are reported as uninterpretable.

## 8. Limits on any claim
- Every PE target comes from **one experiment in which PE is perfectly aligned with lower depth and yield**. Phase 4 can describe which reference pattern the *PE-associated effect vector* resembles. It **cannot** attribute that pattern to PE rather than to library properties.
- **No external generalisation claim is possible:** a single PE dataset, with n = 4–5 per group.
- A depth control is advisable before interpretation: run the same models on a **binomially thinned** vehicle-vs-vehicle contrast (vehicle libraries thinned to PE depth). This was proposed in Phase 1 but is not yet built. It needs approval and implementation before or alongside Phase 4.
- **Fit is expected to be low.** Single-axis |ρ| ≤ 0.25 in Phase 3. A small R² is the expected outcome, not a failure of the method.

## 9. Outputs Phase 4 would produce
- Per target × model: coefficients (naive and reliability-corrected), module-robust CIs, VIFs, CV metrics, LOO ranges.
- One comparison table against the Phase 3 baselines.
- No new reference estimation.
