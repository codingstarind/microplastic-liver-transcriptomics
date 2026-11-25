# Research progression: from the factorial re-analysis to the Phase 5B transfer audit

This document separates two things:
- **chronology:** what was done, in order, and what was known at that time;
- **scientific reasoning:** why each step followed from the evidence then available.

Later corrections are flagged as **[later correction]** at the point they were discovered. They are never presented as known earlier.

Per-phase detail is in `phase_summaries/`. Numbers come from the files cited there.

> **Phase 5C** (search for new datasets; in progress) is **excluded** from this paper.

---

## Timeline at a glance

| Step | Date | Question | Main data | Key output | Why we moved on |
|---|---|---|---|---|---|
| **P0** pre-roadmap factorial re-analysis (course "Phase 1–2" branches) | to 2026-11-02 | Does PE resemble MASH? Does the PE effect depend on diet? Is it robust? | GSE294180 (19 mice) | Strong diet effect; PE lowers lipogenic and cholesterol genes in both diets; no robust interaction (mostly muscle artefact); PE ≡ depth | The "PE resembles MASH" test used a **same-data** diet axis; fasting/time of day is an alternative explanation |
| **R1** problem framing | 2026-11-21 | What should the restart test? | Literature (14 plastic + 17 methods papers), GEO screen | Track A: independent reference axes. Track B: cross-study directions. Track C (PE classifier) rejected | Track A needs verified independent references |
| **R2** reference datasets | 2026-11-22 | Which independent MASH, fasting, time-of-day and PPARα datasets are usable, and can they be mapped without bias? | 6 references + 3 plastic studies | Fixed universe (11,888), 309-gene muscle list, verified designs; Hong controls shared | Signatures could now be estimated |
| **R3** reference signatures + baselines | 2026-11-22 | Are the signatures reliable? Does our diet effect match independent MASH (positive control)? What do baselines say? | 6 references | MASH positive control **FAIL** (c1) → exploratory; fasting stable but ZT-dependent; axes nearly orthogonal; PE–MASH ρ fragile, PE–FAST small but stable | Single axes are weak (|ρ| ≤ 0.25): does a multi-axis, held-out model add anything? |
| **R4** gene-level prediction | 2026-11-22 | Do the references jointly predict held-out PE gene effects? Do trees beat linear? | 11,758 genes × F8 | R²_oos ≈ 0.07 (robust vs null); trees tied with linear on PE | Can representations improve the linear model? |
| **R4B** feature engineering (exploratory) | 2026-11-22 | Do engineered blocks A/B/C help? | same splits | FULL raises linear R² 0.066 → 0.104 but ties CatBoost; C most useful; A useless; B re-learns PE | FULL coefficients flip sign → need diagnostics before freezing |
| **R4C** diagnostic audit | 2026-11-22 | Why are fits low? Are the FULL gains trustworthy? Which models should be frozen? | saved predictions | FULL VIF 233; target noise ceiling (signal share 0.36–0.39); single-mouse influence; freeze F8-OLS (+ F8+C) | Gene-level CV cannot address mice → need a per-mouse estimand |
| **R4D** per-mouse signatures | 2026-11-23 | Do frozen MASH/FAST signatures shift per mouse with PE? | 19 mice | FAST +0.065 (CD), +0.056 (MD); MASH positive control passes, but the PE MASH shift is fragile (cholesterol genes) | Is a MASH-state readout valid across studies? |
| **R4E** study-held-out classifier | 2026-11-23 | Does a MASH-state classifier transfer across MASH studies, and what does it say about PE mice? | 2 MASH studies + 19 mice | AUC at ceiling; engineering improves threshold transfer (S2 balanced accuracy 0.81 → 1.00); the R2 PE shift is cholesterol/PPARα-driven | All internal routes exhausted → test on independent plastic studies |
| **R5** external validation | 2026-11-23 | Do frozen models transfer to independent mice? | Hong (15), GSE245069 (16), GSE277772 (6) | Gene level: Hong ρ 0.22 but permutation p 0.25 (not passed); GSE245069 fail. Mouse level: Hong FAST +0.060, directionally consistent but not PE-specific | Why did it fail? A pipeline error, or something fixable? |
| **R5B** transfer audit and rebuild | 2026-11-23 | What failed, why, and do principled fixes help? | same (post-external) | No pipeline error (count unit 2×, negligible effect); transfer tracks the mouse-state axis (R² 0.89); no intervention met the locked criteria; verdict "inconclusive external power" | Stop; team review; design of an untouched study |

---

## Narrative with reasoning

### P0: Factorial re-analysis of the source study
**Known at the time:** only GSE294180 and the Jung 2026 paper.

**What was done:**
- The original paper's pairwise DEG lists were replaced by a DESeq2 `~ diet * PE` model, so the PE effect within each diet and the interaction could be estimated.
- QC exposed three problems:
  - muscle carry-over (PC1 = muscle genes);
  - complete PE ≡ depth separation;
  - a distinct clock state of CD+Veh.

**What it found:**
- Diet dominates (2,262 DEGs).
- PE lowers lipogenic and cholesterol programmes in both diets.
- No robust interaction: 87 → 16 genes with the muscle covariate, and 1 in edgeR/limma.

**The key limitation that drove everything after:**
- "PE resembles MASH" was tested against the diet effect **from the same mice**. The split-reference ρ was 0.33, and only under CD.
- The pattern "lipogenesis/cholesterol ↓, PPARα targets ↑" is also what **fasting or a different sacrifice time** looks like, and the CD+Veh group was in an unusual clock state.
- A same-data reference could not distinguish these explanations.

### Roadmap Phase 1: Reframing the question
**Purpose:** decide what to test and with which data.

**Inputs:**
- a literature verification (14 plastic–liver papers, 17 methods papers, each with its access level);
- a GEO feasibility screen.

**Decision:**
- **Track A:** compare PE effects with **independent** reference states (MASH from other labs on the same GAN diet product, fasting, time of day, PPARα activation).
- **Track B:** cross-study directions.
- **No PE-vs-vehicle classifier:** n is too small, and it would learn depth.

**Literature finding:** the direction conflict in lipogenic genes. **[later correction, recorded in this phase]** Jung's Table S7 does contain an interaction analysis, but it was never reported, and it is dominated by muscle genes.

### Roadmap Phase 2: Building the independent reference library
**Purpose:** verify candidate references before any effect was estimated. No PE result was used in this phase.

**What was done:**
- Designs were read from GEO SOFT records and the papers.
- Identifiers were mapped to GENCODE vM35 with MGI synonym rescue (e.g. Arntl → Bmal1).
- Technical risks were documented:
  - muscle bleed aligned with feeding state (GSE107787) and with dose group (GSE282998);
  - depth clustered by dose.
- Fixed in advance:
  - the 11,888-gene universe;
  - a 309-gene muscle exclusion rule derived from independent tissue data;
  - the eligible-sample list.

**Corrections recorded:**
- GSE246328's "8-week GAN" group is actually 34 wk of pre-feeding.
- GSE294629 is the same set of mice as GSE246328.
- GSE277772 deposits normalised values, not counts.

**Hong lab:** the four series share 3 byte-identical control libraries, so they form **one study**.

### Roadmap Phase 3: Estimating signatures and the first honest test
**Estimation:**
- MASH (DESeq2), fasting (limma-trend `~ ZT + state`), time of day (harmonic regression) and PPARα (DESeq2 dose and genotype models).
- Analysis plan locked first.

**Positive control:** our Diet_in_Veh vs the independent MASH1.
- ρ 0.355 and sign agreement 0.79, but **below the locked 99th-percentile null** (c1).
- → **FAIL**, so the MASH axis was kept as exploratory.
- **[later correction, post-hoc addendum]** Exact enumeration of 12,870 allocations gives p = 0.036 and q99 0.390. The proposed balanced-null replacement was withdrawn.

**Key biology of the references:**
- Fasting and PPARα agonism are **opposite** for lipogenesis and cholesterol synthesis (fasting −2.79, agonist +1.15).
- So "PPARα targets ↑ with lipogenesis ↓" is fasting-like, not agonist-like.
- Lipogenesis is opposite between MASH-vs-chow and our MASH-vs-purified diet.

**Baselines:**
- PE–MASH1 ρ 0.21, but it collapses to 0.03 without one mouse.
- PE–FAST ρ 0.14 (CD) / 0.25 (MD), stable.
- DEG overlap is "significant" even when most shared genes disagree in direction.

**Consequence:** single axes are weak and the axes are nearly orthogonal (VIF ≈ 1). That motivated a joint, weighted, held-out regression, specified **before fitting** (`phase4_specification_memo.md`).

### Roadmap Phase 4: Gene-level prediction
**Design:**
- Rows are genes; the target is the GSE294180 PE log2FC.
- Features: the effects and SEs of the 4 references (F8).
- Precision weights 1/(s² + τ²).
- Nested CV with **whole co-expression modules held out**.
- Baselines: null and single-axis models.
- Flexible comparators: tree, RF, HGB, CatBoost.
- Paired comparisons with the Nadeau–Bengio correction.

**Results:**
- Multi-axis OLS explains about 7% of held-out variance (CD 0.066, MD 0.070).
- It is robust across 8 variants and 19 LOO targets.
- Trees are **not** better on the PE targets, though they are better on the diet positive control.

**Limits:**
- These are held-out **genes**, not mice.
- Depth thinning controls depth only.

### Roadmap Phase 4B: Feature engineering (explicitly exploratory)
**Question:** can better representations of the references help? Three blocks were tested:
- **A:** reliability shrinkage;
- **B:** pathway membership;
- **C:** response-shape interactions and time-specificity.

The protocol was locked before scoring, but after the Phase 4 results were known.

**Results:**
- FULL lifted linear R² to about 0.10 and met the locked rule vs F8.
- It was only **tied** with CatBoost-F8.
- Ablations:
  - **C** is the useful block;
  - **A** adds nothing;
  - **B**'s gain is GSE294180's own lipogenesis offset, so it is target-informed.

**Why we moved on:** FULL coefficients were unstable across folds.

### Roadmap Phase 4C: Diagnose before freezing
**Findings:**
- **FULL collinearity:** raw vs shrunken effects r = 0.998; VIF 233.
- **Noise ceiling:** signal share is about 0.36–0.39, so R² of 0.07–0.11 is near what the target allows.
- **Heteroskedastic errors:** errors are concentrated in noisy genes.
- **Influence:** single mice are influential.
- **Pathway residuals:** lipogenesis and cholesterol genes are under-predicted. The PE effect there is stronger than any reference predicts.

**Frozen for external use:** F8-OLS (primary) and F8+C-Ridge (secondary).

**The unanswered question:** do the models transfer to **new mice**? Gene-level CV cannot answer that.

### Roadmap Phase 4D: The right unit (the mouse)
**Design:**
- Frozen signatures (built only from references) were scored per mouse with singscore.
- They were analysed with the 2 × 2 model, using LOO, exact within-diet permutation and abundance-matched null gene sets.

**Results:**
- **Fasting-like shift with PE in both diets:** +0.065 / +0.056.
- **MASH positive control passes**, but the PE MASH shift vanishes without cholesterol-synthesis genes.

**Caveats:**
- The CD+Veh clock state (ρ 0.55–0.64 with clock_phase).
- **[later correction in 4E]** Holm-adjusted permutation p-values are 0.20 / 0.14.

### Roadmap Phase 4E: A cross-study MASH-state classifier
**Leakage problem:** the 4D MASH signature used both MASH studies, so it could not be validated across studies.

**Design:**
- Classifiers were trained on one MASH study and tested on the other (two directions).
- Representations were rank-based scores: R0 programme score; R1 12 pathways; R2 = R0 + R1; PCA.

**Results:**
- **AUC is at its ceiling.** What feature engineering and tuning improved was **threshold transfer and calibration**: S2 balanced accuracy went from 0.81 to 1.00, and Brier from 0.19 to 0.03.
- **PCA transferred worst.**
- **Applied to GSE294180:** R2 shows a significant PE_in_MD "MASH-state" shift, but it decomposes into cholesterol/PPARα features. The MASH programme itself contributes +0.18 log-odds.
- **[later correction in P5]** The claim "AUC 1.0 for every non-PCA representation" was false (trees 0.63–0.81). The R0 "Brier" is a misclassification rate.

### Roadmap Phase 5: The first independent-mouse test
**Locked plan:**
- frozen artefacts with hashes;
- Hong treated as one study with shared controls;
- GSE245069 secondary;
- GSE277772 exploratory.

**Gene level:**
- Hong PE ρ 0.22, but the relabelling p = 0.25 → **not passed**.
- PP matches better than PE (0.39).
- GSE245069 fails.

**Mouse level:**
- Hong PE FAST +0.060 (same direction and size as GSE294180) → directionally consistent.
- The same shift appears for PS/PP/PB → **not PE-specific**.
- The MASH score is opposite in sign.
- The GSE245069 PPARα score is up.

**Decision:** gene-level validation NO-GO; narrow GO for write-up.

### Roadmap Phase 5B: Understanding the failure without fooling ourselves
**Status:** every Hong number is post-external. The protocol and at most 3 intervention families were locked before scoring.

**Results:**
- **No pipeline error.** The one mismatch is Hong's 2× count unit, with negligible effect on transfer.
- **Rank transfer roughly equals internal CV** (0.22 vs 0.26), so the "failure" is specificity and precision:
  - transfer ρ is 89% explained by mouse PCs across the 20 labellings;
  - the mouse bootstrap gives −0.03 to 0.30;
  - the disattenuated between-study ρ is about 0.17;
  - control baselines differ.
- **Interventions:** the robust estimator, Huber regression and HGB gave no gain. Diet-averaging gave +0.07, but it is post-hoc and failed the locked criteria.
- **Verdict:** "inconclusive external power", not "generalisation impossible".
- **Memo:** a design for an untouched validation study; Phase 6 not started.

---

## The arc in one paragraph

We started by asking whether PE "looks like MASH" in one public experiment, and found the comparison was circular and confounded. We rebuilt the comparison around independent reference states, and the strongest consistent signal turned out to be **fasting- and cholesterol-type, not MASH-type**. We then:
- quantified it with held-out-gene prediction (about 7%, robust, linear as good as boosting);
- refined it with feature engineering (linear models catch up with boosting, at a cost to interpretability);
- diagnosed why fits are low (target noise, single mice);
- moved to the correct unit, the mouse (a fasting-like shift in both diets);
- validated a MASH-state readout across studies (its PE signal is cholesterol-driven);
- tested independent plastic studies (the direction holds, but specificity and precision do not).

The final audit showed this is a power and design problem, not a coding error.
