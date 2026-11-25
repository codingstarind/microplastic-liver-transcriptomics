# Phase 5 readiness (planning only)

## Frozen artefacts for independent microplastic testing

| Task | Frozen model / score | Status |
|---|---|---|
| Gene level | F8-OLS (primary), F8+C-Ridge (secondary) (Phase 4C) | frozen; validated on held-out genes only |
| Mouse level, signature | Phase 4D MASH_main_n50 / FAST_main_n50 singscore | frozen; the MASH signature used both MASH studies (fine for testing *new* microplastic studies) |
| Mouse level, learned | R0 + L2-logistic and R2 + L2-logistic trained on both MASH studies (`e03`; C = 10) | **valid as a MASH-state discriminator** (study-held-out). R2's response to PE is pathway-driven, so it must be reported with the feature decomposition |
| Fasting classifier | — | **not frozen** (validation n = 6) |

## Independent microplastic studies (Phase 2)
- **Primary:** Hong-lab PE (GSE261893; 3 controls shared by four series, so the four count as one study).
- **Secondary:** GSE245069 (PS, one shared control).
- **Not usable without re-quantification:** GSE277772.

## Valid Phase 5 evaluation
1. Score each external study's mice with the frozen signatures and models (no refitting).
2. Test plastic vs control within study; report estimates with CIs (n = 3 per arm, so the minimum exact two-sided p is 0.1).
3. Gene level: predict each study's plastic log2FC vector with F8-OLS / F8+C; R², Spearman ρ and sign agreement vs null.
4. For every learned-score shift, report the feature decomposition (as in `r2_feature_contributions`).

## Failure criteria (fixed now)
- Shift sign opposite to GSE294180 or a CI centred on 0.
- Gene-level R² ≤ 0.
- A "plastic" effect that is equally present for every polymer in the Hong lab (no specificity).
