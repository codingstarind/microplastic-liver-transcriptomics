# Phase 4E decision memo (one page)

## Did rebuilding the models on the new representation improve the same sample-level task on independent data?

**Partly, and not in the way AUC would suggest.**
- On the study-held-out GAN-vs-chow task (2 studies, 15–16 test mice each), every non-PCA representation reaches AUC 1.0 in both directions. Discrimination is at ceiling, so nothing can be "improved" there.
- What changed is **cross-study calibration**. The baseline signed score with a training-study threshold misclassified 3 of 16 mice in one direction:

  | Model | Balanced accuracy | Brier |
  |---|---|---|
  | Baseline (R0 rule) | 0.81 | 0.19 |
  | Curated pathway scores + tuned L2-logistic | 1.00 | 0.03 |
  | PCA comparator | 0.69–0.75 | 0.16–0.19 |

- These are two studies only, with uncertainty that bootstrap intervals understate at ceiling. The training-label permutation p-values (0.10, 0.035) are weak because of the huge diet signal.

## Is any model reliable enough to apply to GSE294180?
**Yes, as a MASH-state discriminator** (it passed the predeclared rule in both directions), **but its PE readout is not MASH-specific.**
- R0 (the program score): PE shifts are not significant (PE_in_MD +0.86 log-odds [−0.16, 1.89], Holm 0.28).
- R2 (signature + pathways): PE_in_MD +1.66 [0.93, 2.39], Holm 0.0006. This is driven by cholesterol-synthesis/homeostasis and PPARα-target pathway features, which PE lowers or raises; the MASH program signature contributes +0.18.
- This reproduces the Phase 4D finding by an independent route. The PE-associated change looks like a cholesterol/fasting-type metabolic shift, not a disease-type (inflammation/fibrosis) shift.

## What to freeze for independent microplastic testing
- **Gene level:** F8-OLS (primary), F8+C-Ridge.
- **Mouse level:** the Phase 4D FAST and MASH singscores, plus R0 + L2-logistic. R2 + L2-logistic only with its feature decomposition reported.
- **Fasting classifier:** not frozen (external n = 6).

## What we cannot fix
- PE ≡ lower depth/yield; the CD+Veh clock state; n = 4–5 per cell.
- Chow vs purified control.
- Only two independent MASH studies.
- No PE classifier was or should be built.
