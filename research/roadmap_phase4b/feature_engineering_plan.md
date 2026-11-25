# Phase 4B feature-engineering plan (EXPLORATORY second-stage analysis)

**Status:** written 2026-11-22, **after** the Phase 4 results were known and **before** any new feature was scored against a PE target. This is **not** a preregistration of the original analysis. The feature families below were chosen with knowledge of the Phase 4 outcomes. Locked with `comparison_protocol.md` (hashes in `results/protocol_lock.json`).

Only **target-free** diagnostics (no PE quantity) were inspected when choosing them (`results/feature_diagnostics.txt`):
- reference SE vs reference abundance: MASH ρ = −0.70, PPARα −0.93, fasting +0.46, CIRC +0.39;
- reference reliability = 1 − mean(SE²)/var(effect): MASH 0.92, FAST 0.98, **PPARα 0.34**, CIRC 0.77;
- pathway membership coverage and module spread;
- correlation of candidate products with existing features.

## Fixed comparators (Phase 4, unchanged)
- **F4:** MASH_eff, FAST_eff, PPARA_eff, CIRC_eff (reference log2FC).
- **F8:** F4 plus the four reference SEs.

## Block A: reliability-aware (empirical-Bayes shrunken) reference effects (4 features)
- **Formula, per axis k:** shr_k = eff_k · τ̂²_k / (τ̂²_k + se_k²), with τ̂²_k = max(0, var_train(eff_k) − mean_train(se_k²)).
  - This is the normal–normal posterior mean with prior mean 0.
  - τ̂²_k is estimated **on training genes only**, inside the pipeline.
- **Why:**
  - Single-experiment effects are noisy (PROGENy: single-experiment footprints perform no better than random; Sweeney 2017: effect-size reliability matters).
  - The PPARα reference has reliability 0.34: low-abundance probes give huge, unreliable log2FCs.
  - Shrinkage replaces the raw SE columns, which are largely abundance proxies (|ρ| up to 0.93), with one reliability-weighted effect per axis.
- **Rejected alternative:** Wald z (eff/se). It grows without bound for precise, highly expressed genes (|z| up to 40) and is not comparable across assays with different n and noise models.
- **Limitation:** a single global prior per axis; no cross-assay calibration.

## Block B: prior biological programme membership (12 binary features)
- Membership indicators from the frozen Phase 3 gene-set file (`roadmap_phase3/results/gene_sets.csv`, MSigDB 2026.1.Mm plus 2 custom sets defined before any Phase 3 fit):
  LIPOGENESIS_CORE, CHOL_SYNTHESIS_CORE, SANDERSON_PPARA_TARGETS, HALLMARK_FATTY_ACID_METABOLISM, HALLMARK_CHOLESTEROL_HOMEOSTASIS, HALLMARK_INFLAMMATORY_RESPONSE, HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, HALLMARK_XENOBIOTIC_METABOLISM, HALLMARK_BILE_ACID_METABOLISM, HALLMARK_PEROXISOME, HALLMARK_OXIDATIVE_PHOSPHORYLATION, HALLMARK_MTORC1_SIGNALING.
- **Why these 12:** they are the liver-metabolism, inflammation and fibrosis programmes named in the project aims and used in Phases 2–3. They were chosen by biological role, **not** by PE significance. Phase 2/3 PE pathway results exist, but the list is the pre-existing "core + metabolic" list, not a PE-ranked selection.
- **Coverage:** 1,029 of 11,758 genes (8.7%) are in at least one set (mean 1.24 sets per member gene).
  - Each Hallmark set spans 18–20 of the 20 held-out modules (largest share in one module 0.11–0.32).
  - The small custom sets span 5–9 modules.
  - Membership therefore cannot act as a module identifier.
- **Overlap:** handled by keeping the sets as separate binary columns (no single-assignment needed for features). Overlap is reported.
- **Why:** pathway context may modulate how a gene's PE effect relates to the reference effects (Vacca 2024; Lee 2008).
- **Leakage check:** no PE-derived annotation. A **membership-shuffle control** (rows of the 12 columns permuted across genes, seed 13) is run.
- **Expected limitation:** low coverage, so little information. Most of any gain would have to come from interactions learnt by trees.

## Block C: response-shape features (4 features, from reference data only)
1. **FASTxPPARA** = z(FAST_eff) · z(PPARA_eff): a gene co-regulated by fasting and PPARα agonism (shared PPARα-mediated fasting response).
2. **MASHxFAST** = z(MASH_eff) · z(FAST_eff): a gene responding in both the disease and the feeding state.
   - z uses training-gene means and SDs.
3. **FAST_morning_dev** = FAST_clean_eff − FAST_eff: how much the fasting effect at ZT0/ZT4 (Phase 3 clean subset) deviates from the ZT-averaged effect. This captures time-specificity; the raw clean effect is 0.92-correlated with FAST_eff and would duplicate it.
4. **PPARA_consistent** = PPARA_eff · 1[sign(5v0) = sign(10v0) = sign(20v0)]: an agonist effect kept only where the dose response is directionally consistent.

**Duplication check (target-free):** max |ρ| with F8 features is 0.06 (FASTxPPARA), 0.18 (MASHxFAST) and 0.11 (consistency indicator). FAST_morning_dev was constructed to avoid the 0.92 duplication.

## Representations evaluated (no others)

| Code | Features | Count |
|---|---|---|
| F4 | fixed | 4 |
| F8 | fixed | 8 |
| A4 | Block A only (shrunken effects instead of raw effects + SE) | 4 |
| F8+A | | 12 |
| F8+B | | 20 |
| F8+C | | 12 |
| FULL | F8 + A + B + C | 28 |
| FULL−A | | 24 |
| FULL−B | | 16 |
| FULL−C | | 24 |

## Not done, and why
- **Sample-level scores** (GSVA, singscore, PROGENy activities): these produce SAMPLE × pathway matrices and do not define a per-gene feature without inventing a construct.
  - A pathway-level task would need a new target (for example per-pathway PE effects with pathways as rows: about 50 rows, overlapping, very low n). It would be evaluated separately. **Not run in Phase 4B.**
- **PE-derived lists** (Li 2025 hubs, PE DEG lists): leakage. Not used.
