# 04: Targeted interventions (≤ 3 families; LOCKED with 03 before scoring)

Each family answers a problem **measured in Stages 1–2**. Nothing here is chosen to maximise a Hong score. Because Phase 5 outcomes are known, any Hong or GSE245069 gain is **exploratory**.

## Comparators (not interventions)

| ID | Model | Training target | Hyperparameters |
|---|---|---|---|
| **B0** | F8-OLS (frozen, Phase 5) | GSE294180 PE_in_CD MLE, w = 1/(s² + τ²) | none |
| **B1** | F8+C-Ridge (frozen) | same | α = 1000 |
| **B2** | HGB on F8 (the one ensemble comparator; CatBoost is not installed in this environment) | same, weighted | inner grid: learning_rate {0.03, 0.1} × min_samples_leaf {50, 200}; max_leaf_nodes 8, max_iter 300 (the Phase 4 modal choice is inside the grid) |

## I1: Exact pipeline harmonisation (verified mismatch: count unit)
- **Problem (Stage 1):** Hong counts are about 2× fragment-equivalent (1.64 counts per read pair, vs 0.78 in GSE294180 and 0.73 in GSE245069). StringTie prepDE derives counts from coverage of both mates.
- **Fix, fixed a priori (not tuned):**
  - multiply Hong counts by **0.5** and round;
  - rerun the identical Phase 5 DESeq2 pipeline (5-group model, same prefilter rule);
  - sensitivity factor 0.775/1.641 = 0.472, which matches GSE294180's counts per read pair.
- **No model change.** B0/B1 are re-evaluated on the harmonised Hong target.
- **Expected:** SEs and weights change, DEG counts drop, and point estimates barely move. Rank transfer should be nearly unchanged.
- **Not done:** batch correction or cross-study normalisation. Study is confounded with treatment, so it would erase biology.

## I2: Reliability-aware training target (measured problem: outlier-dominated source target, few mice)
- **Problem:** the GSE294180 PE_in_CD MLE is heavy-tailed. The top 1% of genes carry 77% of the squared effect (Hong 33%), and r(MLE, mean-of-logs) = 0.63 (Hong 0.96). These are single-mouse spikes (muscle/blood carry-over). Only 4 v 5 mice inform it.
- **Variants (bounded 3-cell ablation; F8 features; model = Ridge with inner-CV α ∈ {0.1, 1, 10, 100, 1000}, which includes ≈ OLS):**
  - **I2a robust estimator:** source target = mean of log2(CPM + 1), PE − Veh within CD, with a limma-style moderated SE (d0 = 4). This is the same mathematical estimator later used for the external robust proxy.
  - **I2b more mice:** source target = PE main effect averaged over diets, (PE_in_CD + PE_in_MD)/2, from the saved 19-mouse DESeq2 fit (MLE contrast 0,0,1,0.5).
  - **I2c both:** the mean of the two robust simple effects.
- **Weights:** 1/(s² + τ²_DL), computed from each target's own SEs.
- **Disclosure:** Phase 5 showed that PE_in_MD-trained models transferred better to Hong. I2b is partly motivated by that **seen** result, so it is post-hoc. I2a is motivated only by Stage 1/2 evidence on the training data.
- **Literature:** precision weighting and random-effects τ² (DerSimonian & Laird 1986); outlier-robust estimation (Huber 1964); effect-size reproducibility improves with more samples and studies (Sweeney et al. 2017).

## I3: Transfer-robust loss and representation (measured problem: heavy tails and scale shift)
- **Problem:** the target SDs differ (0.33 vs 0.54), and R² is scale-sensitive while rank transfer is not. Heavy tails distort least squares.
- **Model:** sklearn HuberRegressor on F8 (standardised), with target = rank-based inverse-normal scores of the source PE_in_CD MLE.
  - Inner grid: epsilon {1.1, 1.35, 1.75} × alpha {1e-4, 1e-2}; inner scoring = Spearman.
- Predictions are on a normal-score scale, so **only Spearman/sign metrics are compared**; no R² superiority claim.
- **Literature:** Huber 1964; rank-based signed scores that are independent of other samples (Foroutan et al. 2018, singscore) motivate scale-free comparisons across platforms.

## Primary comparison and acceptance criteria (fixed now)
- **Primary outcome (unchanged from Phase 5):** Spearman(prediction, Hong PE − CON DESeq2 MLE log2FC) over the shared genes.
- A candidate is an **exploratory external improvement** only if **all** of the following hold:
  - **A1:** Δρ vs B0 ≥ **+0.05**.
  - **A2:** paired mouse bootstrap P(Δρ > 0) ≥ **0.90**.
  - **A3:** internal outer-CV Spearman (held-out genes, PE_in_CD target) is not more than 0.02 below B0.
- Add **"with PE-specificity"** only if **A4** also holds: the true labelling ranks **1st of 20**, **or** Spearman with the control-free gene vector PE − mean(PS, PP, PB) is > 0 with bootstrap P > 0.90.
- **Secondary (reported, not decisive):**
  - GSE245069 pooled-PS Spearman (cross-polymer);
  - weighted R² vs own mean (log2FC-scale candidates only);
  - the I1-harmonised target.
- **Verdict labels:**
  - "technical mismatch fixed";
  - "true positive transfer improvement (exploratory)";
  - "internal improvement only";
  - "inconclusive external power";
  - "no biological transfer detected".
