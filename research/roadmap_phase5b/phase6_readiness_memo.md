# Phase 6 readiness memo (for team review; Phase 6 NOT started)

## What can be claimed (with the exact qualifier)
1. **Internal (held-out genes, GSE294180):** reference-derived features predict about 7% (R²) / ρ ≈ 0.26 of PE gene-effect variance. Train ≈ CV, so there is no overfitting. This is a weak, stable predictor. It is **not** a statement about new mice.
2. **External rank transfer roughly matches internal CV in Hong PE** (ρ 0.22 vs 0.26 for the CD-trained model). The external "failure" was largely the R² metric meeting a noisier, wider target.
   - Phase 5B shows this agreement is **not PE-specific**. The same predictions correlate with any split of Hong mice along their main variation axis (89% of relabelling ρ explained by mouse PCs; the true labelling ranks 4th–5th of 20) and with the other polymers against the shared controls.
3. **Mouse level (Phases 4D/5):** a reproducible fasting/cholesterol-type shift with PE in GSE294180. The same direction appears in Hong PE (3 v 3, exact p = 0.10), but non-specifically across polymers.
4. **Technical audit:**
   - no pipeline error drove the Phase 5 result;
   - one verified count-unit mismatch (Hong about 2×) was corrected, with negligible effect on transfer.
5. **MASH-state classifier (Phase 4E), study-held-out in 2 directions:** discrimination at ceiling. The tuned R2 logistic gives balanced accuracy 1.00 in both directions. The R0 rule gives 0.81 in S2, and its error is a misclassification rate, not a Brier score. With 15–16 test mice, "perfect" is imprecise, and two directions are not a population of cohorts.

## What cannot be claimed
- That frozen gene-level models **validated** on independent mice. The external Spearman is not PE-specific and is not supported by the mouse-label null.
- That any Phase 5B modification **generalises better**. All 5B gains were measured on already-seen external data (post-external exploratory). None met all locked criteria (`05_…md`).
- That PE **causes** the fasting-like shift, or **induces MASH**.
- That generalisation is **impossible**. The data cannot separate PE biology from mouse-state variation, diet, duration and dose. "No evidence" is not "proved absent".
- That batch correction or resampling fixed the GSE294180 PE ≡ depth/yield confound, or the Hong design differences (dose and size unknown, 71 wk, chow, shared controls). **They did not and cannot.**

## Recommended framing for the final report
> A weak but stable gene-level predictor of PE effects, and a fasting/cholesterol-type mouse-level shift, were internally reproducible. When transferred to an independent PE experiment, both moved in the expected direction. However, the external experiment (3 v 3, shared controls, different diet, duration and dose) cannot distinguish a PE-specific effect from between-mouse metabolic-state variation, so external validation is inconclusive rather than negative.

## Design for a genuinely untouched validation study (none exists now)
- **Animals and arms:**
  - C57BL/6J males, randomised at 8 weeks;
  - **PE 10–150 µm, the same dose and gavage route as GSE294180**, vs vehicle;
  - **a matched second polymer (PS) with its OWN vehicle control**, so specificity is testable without shared controls;
  - both a purified control diet **and** grain chow (2 × 2 × 2 if feasible);
  - n ≥ 8 per arm. With 8 v 8, C(16, 8) = 12,870 labellings, so exact p resolution is about 1e-4; the observed FAST effect size g ≈ 2 gives > 95% power at α = 0.05.
- **Husbandry and confounders:**
  - fixed sacrifice ZT and recorded fasting duration (the fasting axis is the main nuisance);
  - library prep and sequencing balanced across arms, with depth recorded;
  - the same quantification as training (Salmon / GENCODE vM35, fragment counts).
- **Pre-registered frozen artefacts:**
  - F8-OLS (Phase 5 frozen);
  - FAST/MASH singscores (Phase 4D);
  - R0/R2 logistic (Phase 4E).
- **Primary endpoints:**
  - **(i)** FAST singscore PE − vehicle > 0 (exact p < 0.05);
  - **(ii)** gene-level ρ(F8-OLS, PE log2FC) above the 95th percentile of the full mouse-relabelling null;
  - **(iii)** PE effect > PS effect, each against its own control.
- **Cheaper partial option:** re-quantify GSE277772 (PS, 3 v 3, about 15–20 GB FASTQ) for count-based use. It was only rank-inspected in Phase 5, but it is PS, tiny and already viewed, so it would be a **sensitivity** analysis, not a fresh test. This needs team approval for the download.

## Go / no-go
- **GO** to a reporting Phase 6 that uses the framing above, with the Phase 5B diagnostics as a methods result.
- **NO-GO** for any external-validation, PE-specificity, causal or MASH claim.
