# 03: Resampling and nested-CV protocol (LOCKED before any Phase 5B model was scored)

**Status of evidence:** Phase 5 external outcomes (Hong PE, GSE245069 PS) have already been seen. Every Phase 5B comparison on those studies is therefore **post-external exploratory re-development**, not fresh validation.
- **No untouched independent study exists.** Phase 2 lists no other PE- or PS-exposed mouse-liver RNA-seq with counts.
- GSE277772 was inspected in Phase 5, and its counts would need SRA re-quantification.
- A design for a genuinely new validation is in `phase6_readiness_memo.md`.

## A. Gene-effect models (unit = gene for fitting; unit = mouse for inference)
1. **Internal estimate (reused, not rerun):**
   - Phase 4 locked nested grouped CV: 20 Ward co-expression modules from GSE117134; 5 outer folds × 3 repeats (seeds 101–103); inner GroupKFold(4).
   - Phase 4C fold-level metrics are reused for F8-OLS, F8+C and F8+B+C.
   - New candidates use **exactly the same outer folds and inner scheme** (`p4_data.outer_splits`).
2. **Common internal outcome for every candidate:** held-out genes' GSE294180 PE_in_CD DESeq2 MLE log2FC (the Phase 4 target).
   - **Primary internal metric:** Spearman. It is scale-free, so candidates trained on different target scales stay comparable.
   - **Secondary:** unweighted R² vs the training mean, reported only for candidates on the log2FC scale.
3. **Mouse-level bootstrap (B = 1,000, seed 20261122, paired across candidates):**
   - Resample mice **within each arm** in GSE294180 (4/5/5/5) and in Hong (3 per group).
   - Re-estimate source and external targets with fast estimators:
     - "arith" = log2 ratio of group-mean CPM (+1), the MLE surrogate;
     - "robust" = mean of log2(CPM + 1).
   - Refit each candidate with its **selected** hyperparameters (no re-tuning inside the bootstrap; stated limitation).
   - Recompute the Hong transfer Spearman.
   - Report percentile intervals and the paired P(Δρ > 0) vs F8-OLS.
   - **Limitation:** with 3 mice per arm, within-arm resampling repeats mice often and understates variance by about (n−1)/n. The intervals are descriptive.
4. **Exhaustive relabelling (Hong PE vs CON):**
   - All C(6,3) = 20 labellings, using the arith proxy.
   - One-sided p = rank of the true labelling / 20 (minimum 0.05; conventional two-sided minimum 0.10).
   - Exchangeability is plausible within one lab and one run, but it does **not** remove the shared-control or husbandry confounding.
5. **Module-block bootstrap of genes** (2,000 resamples of the 20 modules) for descriptive gene-level CIs only. It is never used as a test of a PE effect.
6. **No naive gene-wise CIs.**

## B. Mouse-level MASH classification
- Keep the Phase 4E study-held-out design: train GSE316023 → test GSE246328 and the reverse. That gives **2 independent direction tests**.
  - Inner CV, signatures and scaler were fitted inside training studies (verified in Phase 4E). No refit.
- **New in 5B:** a stratified mouse bootstrap (2,000 resamples within class) of the saved external predictions, for:
  - AUC;
  - balanced accuracy at 0.5;
  - Brier (probabilistic models only).
- The R0 rule's 0/1 output is reported as a **misclassification rate**, never as a Brier comparison.
- The fraction of bootstrap replicates with AUC = 1 shows how imprecise "perfect" is with 15–16 mice.

## C. Multi-study microplastic model redevelopment: NOT done
- **Independent labs with plastic-exposed liver counts:** GSE294180 (PE), Hong (PE/PS/PP/PB, one shared control) and GSE245069 (PS). GSE277772 is normalised-only.
- That is **2 labs with PE** and 3 with any plastic. Leave-one-lab-out would give 2–3 folds of different polymers, doses and durations.
- Leaving one Hong polymer out while keeping its identical control mice in training is invalid.
- **No credible internal–external CV exists.** We therefore keep the source-trained model plus an external diagnostic (Steyerberg & Harrell 2016; Riley et al. 2016; Bernau et al. 2014).

## D. Hyperparameters and seeds
- **Grids:** in `04_targeted_interventions.md`.
- **Selection:** inner GroupKFold(4) on training genes, scoring = negative MSE (Spearman for I3, which is fit on normal scores).
- **Final models for external use:** refit on all genes with the modal inner choice across the 15 outer folds.
- **Seeds:** 20261122 for everything, plus the Phase 4 outer seeds 101–103.
- **Workers:** 6 loky workers (about 4 GB free RAM), BLAS/OMP threads = 1, single-threaded estimators.
