# Three different tasks: do not mix their numbers

| | (A) Gene-level regression (Phases 4–4C) | (B) Sample-level state prediction (this phase) | (C) Sample-level PE response (Phase 4D + Step 5 here) |
|---|---|---|---|
| Row | gene (~11,758) | mouse | mouse (19) |
| Outcome | estimated PE log2FC (noisy; the same 9–10 mice for every gene) | known, independently assigned label: GAN diet vs chow (or fasted vs fed) | none predicted; a frozen score is the *response variable* in score ~ diet × PE |
| Estimand | how well reference effects predict held-out genes' PE effects *within one experiment* | out-of-study accuracy of a MASH-state (or fasting-state) predictor | group differences / interaction of a frozen score among the 19 mice |
| Validation | held-out gene modules (cannot generalise to new mice) | held-out **study** (train on one MASH study, test on the other) | none (inference); leave-one-out, permutation, null gene sets |
| Independent units | 9–10 mice per PE contrast; genes are not independent | **2 MASH studies** (16 + 15 mice); fasting: 1 training study (36 livers, 10 clean) + 1 tiny test (6) | 19 mice in 4 cells (4/5/5/5); PE ≡ lower depth |
| Realistic capacity | 4–24 features, linear or shallow trees | ≤ 1–13 features, heavily regularised; 2 studies give 2 estimates, not a precise generalisation error | 1 score per model; ≤ 4 coefficients |

Changing the representation does **not** increase the number of mice. Gene-level R², sample-level AUC and factorial p-values are not comparable on one leaderboard.

**Correction to the summary in the request:** in Phase 4D the t-based Holm-adjusted p for FAST PE_in_CD was 0.036 (PE_in_MD 0.056). The exact *permutation* p-values (0.040 / 0.024) were not Holm-adjusted in Phase 4D. Adjusted over the same 6-test family they become 0.20 / 0.14, so the permutation evidence does **not** survive multiplicity. Both statements stand; nothing in Phase 4D is changed.
