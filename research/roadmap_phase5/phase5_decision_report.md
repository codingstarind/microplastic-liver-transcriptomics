# Phase 5 decision report: external validation of frozen models

**Plan:** locked before any external effect was computed (`locked_plan.md`, `results/plan_lock.json`).
- No model was refitted, retuned or reoriented on external data.
- The frozen R0/R2 models reproduce the Phase 4E GSE294180 scores to 4e-16 (`frozen/frozen_manifest.json`).
- Tests: 10/10 pass. Phase 4D lock verified.

## Inputs actually available (`dataset_provenance.csv`)

| Unit | Mice | Input | Role |
|---|---|---|---|
| **Hong lab:** GSE261893 PE + GSE261725 PS + GSE286116 PP + GSE262399 PB | **15** (24 GSMs). The 3 controls are byte-identical in all 4 series, with identical SRA read counts | integer counts | primary PE test |
| **GSE245069** (PS 80 nm / 0.5 µm / 5 µm, ICR) | 16 (4 water controls shared) | estimated counts | cross-polymer |
| **GSE277772** (PS) | 6 | normalised values only, no counts on GEO | exploratory, not validation |

- **AIBD GSE245069 artefacts:** not found anywhere on this machine, so that analysis was done fresh here.
- **Comparability:**
  - Hong PE is in food for 71 weeks, size/dose not stated, grain chow.
  - GSE294180 PE is 10–150 µm by gavage for 8 weeks, purified control diet.
  - So Hong is "PE vs control on a control diet", not a replica of CD+PE vs CD+Veh.

## (A) Numerically reproduced (internal)
- **Frozen 4E R0/R2 models:** reproduce the saved GSE294180 log-odds exactly.
- **Frozen 4D singscores:** reproduce the Phase 4D PE_in_CD estimates exactly (< 1e-9).
- **F8-OLS:** full-gene coefficients have the same signs as the Phase 4C fold medians.
- **Corrections to the Phase 4E wording:**
  - It is **false** that every non-PCA model reached AUC 1.0. For example, the R0 tree / random forest reached 0.81 in S2, and the R2 tree 0.80 / 0.63. The ceiling applies to the shortlisted logistic / pathway models.
  - The R0 rule's Brier score is a hard-label misclassification rate. It is **not** a like-for-like probabilistic calibration comparison with the logistic models.

## (B) External evidence

### Q1: gene level, Hong PE vs CON (`gene_level_external_results.csv`)

| Frozen model (trained on PE_in_CD) | weighted R² vs own mean | Spearman [gene-block CI] | mouse-label permutation p (20 labellings) | locked verdict |
|---|---|---|---|---|
| **F8-OLS (primary)** | **0.041** | 0.22 [0.14, 0.28] | **0.25** | **not passed:** positive, but not supported by permutation |
| F8+C-Ridge | 0.004 | 0.14 [0.08, 0.20] | 0.30 | same |
| FAST-only OLS | 0.045 | 0.29 [0.23, 0.34] | 0.05 (minimum attainable) | PASS-direction |
| observed GSE294180 PE vector | −0.053 | 0.10 | 0.25 | R² < 0 |

- **Why the CI and the permutation disagree:** the module-block bootstrap treats genes as the units. The permutation shows that random relabellings of the 6 Hong mice give Spearman up to 0.34 (95th percentile) with the same predictions. Most of the apparent gene-level agreement is therefore generic between-mouse (largely fasting-axis) variation, not a plastic effect.
- **Specificity check:** the same predictions match **PP vs CON better** (F8 ρ = 0.39) and PB similarly. These contrasts share the controls.
- **GSE245069 pooled PS:** F8-OLS(CD) ρ = −0.07, R² −0.14. This is a **FAIL**; only the FAST-only axis is positive (0.19, p = 0.13).
- **Noise ceiling:** signal share is 0.47 (Hong PE) and 0.65 (PS pooled), so the targets are interpretable but noisy.

### Q2: mouse level (`sample_level_external_results.csv`, `per_sample_scores.csv`, fig1)
- **Primary (locked): Hong PE − CON, FAST singscore = +0.060** [−0.004, 0.124]. Hedges g 2.5, exact p = 0.10 (the minimum for 3 v 3), leave-one-out 6/6. GSE294180 PE_in_CD was +0.065.
  - → **DIRECTIONALLY CONSISTENT** under the locked rule.
- **The shift is not PE-specific:**
  - PP +0.088, PB +0.038, PS +0.024. Pooled non-PE polymers vs CON +0.050 (exact p 0.04).
  - The 3 shared controls have the **lowest** FAST mean of all 5 groups.
  - Control-free PE minus the other polymers is only +0.010.
  - We cannot tell a generic particle/food effect from a control-group artefact (feeding state or sacrifice timing of 3 mice).
  - The pre-fixed failure criterion "effect equally present for every polymer" applies: **no PE-specific claim**.
- **MASH singscore:** Hong PE −0.021 and R0-logistic −0.56 log-odds. Both are **opposite** to GSE294180, so there is no MASH-type shift.
- **R2-logistic:** Hong PE +0.94 (exact p 0.30). The decomposition is 69% cholesterol/fasting-like features, and the disease-like share is net negative. This is the same pattern as GSE294180 (fig3).
- **GSE245069 PS pooled:**
  - FAST +0.036 (p 0.24, n.s.);
  - PPARα singscore +0.053 (exact p 0.001, 16/16 LOO; strongest at 80 nm);
  - R2 +1.93, 70% cholesterol/PPARα-like;
  - MASH ≈ 0.
- **GSE277772 (exploratory):** FAST −0.013 and MASH −0.027, in the opposite direction.
- **Controls:**
  - Depth is balanced in all external studies (exposed/control ratio 0.95–1.0). Unlike GSE294180, there is no PE ≡ depth confound.
  - Coverage control: reduced gene backgrounds change the GSE294180 estimates by ≤ 9%.

### Pathways (`pathway_concordance.csv`, fig4)
- **Hong PE:** lipogenesis (−0.82) and cholesterol synthesis (−0.42) are down, the same direction as GSE294180 PE_in_CD. Across the 12-set panel, Spearman is 0.66 with GSE294180 PE_in_CD and 0.55 with fasting.
- **But Hong PS, PP and PB also lower lipogenesis** (−0.46 to −0.71). Cross-polymer gene-effect correlations are 0.57–0.77, against about 0.5 expected from the shared control alone.

## (C) Exploratory observations
- Across two independent labs and four polymers, plastic-exposed livers tend to show a **lipid-programme / fasting-like / PPARα-type shift, not a MASH or inflammation shift**. This agrees with Phases 4D and 4E.
- Each piece rests on 3–4 mice per arm and shared controls.
- GSE245069 80 nm shows the clearest PPARα/fasting-like shift.

## (D) Limitations
- **External PE evidence is one lab, 3 v 3, and shares controls with three other polymers:**
  - the minimum exact p is 0.10;
  - specificity cannot be established;
  - the controls look atypical on the FAST axis.
- **Domain shift:** duration (71 vs 8 weeks), route, size/dose (unknown), grain chow vs purified diet, strain (ICR), and quantification pipeline.
- **Gene-level:** gene-wise CIs are anti-conservative. Only the mouse-label permutation respects the true unit (mice), and it does not support F8.
- **Unknowns:** fasting/clock state at sacrifice is unknown in every external study; muscle carry-over is similar to GSE294180, with no imbalance.
- **GSE277772:** normalised values only; never treated as validation.

## Go / no-go for Phase 6 (rule fixed in `locked_plan.md` §8)
- **Gene-level external-validation claims: NO-GO.** Q1 did not pass in Hong (permutation p 0.25), and it failed in GSE245069.
- **Narrow GO** for a write-up phase whose main result is the one below. It must include the R2 decomposition and the non-specificity finding.
- **No GO** for any PE-specific, causal or disease (MASH) claim.

## Recommended main result
> Frozen, pre-registered scores built from independent reference data show that PE exposure in GSE294180 is associated with a fasting/cholesterol-type hepatic metabolic shift, not a MASH-type shift.
>
> In an independent PE study (Hong lab), the frozen fasting score moved in the same direction with almost the same magnitude (+0.060 vs +0.065; 3 v 3, exact p = 0.10). The MASH score did not.
>
> The same shift appeared for PS, PP and PB against the same shared controls, so it is not PE-specific and may partly reflect the control group.
>
> Gene-level prediction models did not validate externally once the mouse-label permutation null was applied.

## BT3041 teaching story (technique → why it was used)
1. **Preprocessing and identifier harmonisation** (vM35 mapping, prefilters, TPM-proxy vs log-CPM). Cross-study scores are meaningless without a fixed gene universe and comparable inputs.
2. **Count models** (DESeq2 Wald MLE; limma-voom / trend sensitivity). The external target must be the same estimand as the training target (MLE log2FC + SE).
3. **Feature engineering, before vs after:**
   - Gene level: F8 → F8+C added +0.014 to +0.026 internal R², but **less** external agreement in Hong PE (ρ 0.14 vs 0.22).
   - Sample level: R0 threshold rule → R1/R2 logistic + tuning fixed cross-study threshold transfer (S2 balanced accuracy 0.81 → 1.00), but R2's output is dominated by lipid pathways.
4. **Model selection and tuning frozen before testing** (locked plan + hashes). This prevents selecting on the test set.
5. **Validation at the right unit:**
   - held-out genes (Phase 4), then held-out study (Phase 4E), then independent experiment (here);
   - permuting **mice**, not genes, exposed an over-optimistic gene-bootstrap CI.
6. **Hypothesis testing with small n:** exact permutation (3 v 3 → p ≥ 0.10), leave-one-out, and effect sizes with CIs instead of p-value chasing.
7. **Enrichment:** camera (inter-gene correlation) on a pre-fixed panel; concordance with reference axes (fasting vs MASH) to interpret the direction.
8. **Diagnostics:**
   - shared-control correlation (expected 0.5);
   - control-group rank;
   - depth/muscle balance;
   - coverage control;
   - noise ceiling.
9. **Honest interpretation:** report the failure (Q1), the non-specificity (Q2), and what cannot be separated.

**The team decides** whether this narrow GO is acceptable and how to word it. Nothing was committed, and Phase 6 was not started.
