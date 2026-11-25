# Phase 4: research question, prediction targets and feasibility gate

**Written 2026-11-22, before any Phase 4 model was fitted** (locked together with `evaluation_protocol.md` and `model_search_spaces.csv`; hashes in `results/protocol_lock.json`).

## Status of inputs carried from Phase 3
- **MASH reference (GSE316023):** **CONDITIONAL / exploratory.** The locked positive-control criterion c1 failed (exact permutation p = 0.036 > 0.01; Phase 3 review addendum). It is not repaired here, and the exact permutation audit is already complete. Every result that uses the MASH axis carries this status.
- **Other references:**
  - fasting (ZT-averaged, time-dependent);
  - PPARα (WT GW7647 20v0; agonist-associated, PPARα dependence provisional);
  - time of day (T1 = ZT6 − ZT2; female).
- **GSE294180:** PE is completely aligned with sequencing depth and yield. **No model in this phase can separate PE from that alignment.**

## A. Primary task: continuous gene-effect prediction (FEASIBLE; carried out)

| Element | Definition |
|---|---|
| **Question** | Do independent MASH, fasting, PPARα and time-of-day reference effects predict the held-out GSE294180 PE-associated gene effects better than null, single-axis and linear baselines? |
| **Rows** | Genes in the fixed Phase 2 main universe (11,888; locked 309-gene muscle exclusion) with a non-missing target and all primary features. **Rows are not independent biological replicates.** |
| **Target y** | GSE294180 DESeq2 MLE log2FC: **PE_in_CD** and **PE_in_MD**, modelled separately (existing Phase 2 results; not re-estimated). |
| **Training weights** | w = 1/(s² + τ²), with s = target lfcSE and τ² estimated (DerSimonian–Laird) **on training rows only**. Not a feature. |
| **Features** (primary set F8) | log2FC and SE of MASH1 (GSE316023), FAST (GSE107787, ZT-averaged), PPARA (GSE282998 WT 20v0) and CIRC_T1 (GSE117134), i.e. 4 effects + 4 uncertainties. **Excluded:** gene IDs, target significance, target baseMean, gene-set labels, any PE-derived quantity. |
| **Independent units** | There are none at the level of the claim. Held-out genes test **within-experiment generalisation across gene modules**, not held-out mice or laboratories. Effective sample size is far below the number of genes, because gene effects are correlated. |
| **Pooled CD+MD model** | **Rejected as primary.** The two targets share gene rows, and their errors are correlated through shared normalisation and the common GSE294180 dispersion fit. A pooled model would double-count genes. Each target is reported separately. |

## B. Optional mouse-level classification (FEASIBILITY: DEMONSTRATION ONLY)

**Question:** can a biological-state classifier (GAN/MASH diet vs control) trained on independent mice transfer to other mice?

| Check | Finding |
|---|---|
| Comparable MASH-vs-control groups | GSE316023 wk23 at 21 °C (8 v 8) and 29 °C (8 v 8): both classes appear at both temperatures, so class is **not** confounded with temperature. GSE246328 24 w (7 v 8) in a different lab. |
| Extra groups to increase n | GAN wk12/16/32 (GSE316023) and 8 w/16 w (GSE246328) have **no matched controls**. Adding them would confound class with diet duration and age. GSE294629 re-deposits GSE246328 (same mice). **No valid way to increase n.** |
| Class vs study | Both classes are present in both studies; training on one study and testing on the other does not make class equal to study. |
| Separability | GAN vs chow explains ~50% of expression variance (PC1 r = 0.97 with diet). Near-perfect accuracy is expected from **any** classifier, so differences between classifiers will be uninformative (ceiling). |

**Decision:**
- Run a **small demonstration**: train on GSE316023 (32 mice), test externally on GSE246328 (15 mice).
- Apply descriptively to GSE294180 **vehicle** mice only (CD 4, MD 5; purified control diet, so the control is not comparable).
- Accuracy is reported as **inconclusive for model comparison**. It is not a research claim.

**Prohibited and not done:**
- a PE-vs-vehicle classifier on the 19 mice;
- a three-state MASH/fasting/PPARα classifier (each class from a different laboratory would learn study identity).

## What success would and would not mean
- **Would mean:** within this one experiment, the pattern of PE-associated gene effects is partly predictable from independently measured biological-state effects. The robust per-axis contributions describe which state the PE vector resembles.
- **Would not mean:**
  - that PE causes that state;
  - that the result generalises to other mice or laboratories;
  - that a technical (depth/yield) explanation is excluded;
  - that coefficients or importances identify biology (errors in variables, cross-technology scales).
