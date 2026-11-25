# Phase 4D analysis plan (locked before scoring any mouse; SHA-256 in results/plan_lock.json)

**Estimand:** per-mouse signature scores (singscore total score) analysed with the original 2×2 design. **This is not gene-level prediction.**

## Signatures (frozen; `frozen_signature_manifest.csv`)
- **Primary:** MASH_main_n50, FAST_main_n50.
- **Secondary:** PPARA_main_n50 (agonist-associated; PPARα dependence provisional).
- Circadian is not scored (it is a time-pattern problem).

## Scoring
- singscore per mouse (Foroutan 2018):
  - genes ranked ascending within each mouse over the background;
  - up-set mean rank and down-set *reversed* mean rank, each normalised to the theoretical minimum and maximum and centred at 0;
  - total = up + down;
  - signature genes absent from the background are dropped (reported).
- **Background:** 13,777 GSE294180 filtered genes with a GENCODE length.
- **Primary input:** TPM-proxy. **Sensitivities:** log-CPM, VST.

## Model
- OLS per score: S = b0 + bD·MD + bP·PE + bDP·MD·PE. Treatment coding (CD, Veh reference). The mouse is the unit.
- **Contrasts:**

  | Contrast | Expression |
  |---|---|
  | PE_in_CD | bP |
  | PE_in_MD | bP + bDP |
  | Interaction | bDP |
  | Diet_in_Veh | bD (positive control) |

- 95% t-intervals (df 15). Effects in score units, plus standardised effect = estimate / residual SD.
- **Multiplicity:**
  - **primary family:** {MASH, FAST} × {PE_in_CD, PE_in_MD, Interaction} = 6 tests, Holm;
  - **positive control:** MASH Diet_in_Veh, a single test; **PASS** if the estimate > 0 and its 95% CI excludes 0;
  - **secondary family:** PPARα × 3, Holm;
  - fasting Diet_in_Veh: descriptive.
- An interaction p > 0.05 is not evidence of no interaction.

## Controls (all fixed now)
1. **LOO:** 19 refits; report the range of each contrast and any sign change.
2. **Exact within-diet label permutation:**
   - PE relabelled within CD (C(9,5) = 126) for PE_in_CD;
   - within MD (C(10,5) = 252) for PE_in_MD;
   - jointly (126 × 252 = 31,752) for the interaction;
   - two-sided, |statistic| ≥ |observed|.

   **Caveat:** PE is aligned with library depth and yield, so exchangeability under H0 is not guaranteed.
3. **Null gene sets:** 1,000 random signed sets per signature (seed 20261122), each gene replaced by a random non-signature gene from the same TPM-proxy mean-abundance decile; same up/down sizes. Report the percentile of each observed contrast within the null distribution. **Not biological replication.**
4. **QC:** scores vs depth, quality_proxy (yield), muscle and clock_phase, by group (plots plus Spearman correlations; descriptive, no covariate "correction").
5. **Sensitivities:**
   - input (CPM, VST);
   - size (n25);
   - muscle rule (alt 88-gene universe, n50);
   - MASH without CHOL_SYNTHESIS_CORE genes;
   - FAST and MASH with their 14 shared genes removed.

   GSVA is not run (not installed).

## Stopping rule
- If the MASH positive control fails, the MASH score is reported as failed.
- Null or unstable PE results are reported as they are.
- No gene, threshold or score is changed after scoring.
