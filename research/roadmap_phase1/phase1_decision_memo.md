# Phase 1 decision memo

**For:** BT3041 project team · **Date:** 2026-11-21 · **Status:** recommendation for discussion. The team owns the questions, hypotheses and final decisions.

## 1. What problem should we investigate?

**Recommended:** *When the disease reference comes from independent mice, does the liver response to polyethylene (PE) look like MASH, or like another known state such as PPARα activation or fasting/time-of-day? And does that answer change on a MASH diet?*

This keeps our two original questions but makes RQ1 testable without circular reasoning.

**Possible second strand:** *Do the main liver metabolic programmes (fat synthesis, cholesterol synthesis, PPARα targets, inflammation) change in the same direction with plastic exposure across independent labs?*

## 2. Why is it worth studying?

- **The headline claims rest on comparisons that cannot prove them.** The field says "microplastics worsen fatty liver disease" and Jung 2026 says PE effects resemble MASH. Both rest on significant-vs-non-significant comparisons and on disease signatures built from the same mice.
- **The pattern has another obvious explanation.** The PE signature in our data is lower fat and cholesterol synthesis with higher PPARα targets. That is also what fasting or a different time of sacrifice looks like. Our control-diet vehicle mice do have a distinct clock-gene state.
- **The literature contradicts itself on direction.**
  - Fat-synthesis genes go **down** in Jung (PE), Makdissi and Chiu (both polystyrene).
  - They go **up** in Djouina (PE), Okamura and Xu P (both polystyrene).
  - PPARα is induced in some studies and suppressed in one that fasted its mice overnight.

## 3. What have researchers already done?

- **Study size and reporting.** Most studies use 3–6 mice per group. They report DEG lists at very different thresholds; two called DEGs without multiple-testing correction.
- **Design versus analysis.**
  - Nine studies use a diet × plastic design.
  - None reports and interprets a diet × plastic interaction for molecular data.
  - **Correction to our earlier review:** Jung did compute one (Table S7), but never mentioned it. It is dominated by skeletal-muscle genes, the same artefact we found.
- **What our own Phase 2 already showed.**
  - A strong diet effect.
  - Fasn, Acly and Acaca reliably lower with PE in both diets.
  - Cholesterol programme lower with PE in both diets.
  - No robust diet × PE interaction.
  - These PE results cannot be separated from the fact that every PE library was sequenced more shallowly.
- **Only one other PE liver transcriptome is public:** GSE261893. It has 3 vs 3 mice, a 71-week chow diet and no linked paper. Its controls are byte-identical to those of three other series from the same lab.

## 4. What could we improve?

| Track | Improvement | Status |
|---|---|---|
| **A (primary)** | Build independent reference axes: MASH from GSE316023, which uses **the same MASH diet as Jung (D09100310)**, plus GSE246328; and PPARα from GSE282998. Then use regression and model selection to ask which axis best explains the PE response in each diet. Report interaction estimates with confidence intervals. | Candidate |
| **B (complementary)** | Re-analyse the ~4 independent public plastic datasets with one pipeline. Compare programme-level directions, not DEG lists or hub genes, and check each dataset's sequencing depth against treatment. | Candidate |
| **C (machine learning)** | Exposure classifiers are **rejected**: they would learn sequencing depth or the lab, and n is far too small. The only sensible ML is a MASH-state score trained on external data, which is really an implementation of Track A. | Merge into A |

## 5. What would we compare against?

- **Track A:**
  - **Baselines:** Jung's DEG-overlap approach, and our own same-data split-reference correlation (ρ ≈ 0.33).
  - **Positive control:** our diet contrast must align with the external MASH axis. If it doesn't, the axis is unusable.
  - **Technical controls:** computationally downsample vehicle libraries to PE sequencing depth; refit with the muscle-contamination covariate.
- **Track B:**
  - **Baseline:** gene-level DEG overlap and hub-gene agreement between studies.
  - **Reference for "how good can agreement get":** agreement between the two independent MASH-diet studies.

## 6. Which data would show whether our approach is better?

- **Already public and screened:**
  - GSE294180 (ours);
  - GSE316023 and GSE246328 (independent MASH references);
  - GSE282998 (PPARα reference);
  - GSE245069, GSE277772, GSE261893 and siblings (independent plastic studies).
- **Still missing:**
  - **A fasting or time-of-day liver reference with counts.** Not yet searched; the fasting comparison depends on it.
  - **Confirmation that GSE277772 has raw counts.** Our two checks disagreed.
- **Approach A "is better" if:**
  - it passes the positive control;
  - it explains held-out pathways' PE effects better than the same-data baseline;
  - its conclusion survives the depth and muscle controls.
- **Approach B "is better" if:** programme-level agreement between labs is clearly higher than gene-level agreement.

## 7. Limits we must keep stating

1. PE is perfectly aligned with lower sequencing depth and yield. **No new analysis removes this.**
2. Control-diet vehicle mice have a distinct clock state.
3. Skeletal-muscle contamination in some samples.
4. Deposited counts lack length offsets (assumptions A1–A3).
5. Cage is the same as treatment group.
6. External studies differ in polymer, strain, diet, route and duration, so they can never be pooled into one dataset.

## 8. Decisions the team needs to make

1. **Scope versus deadline.** The final-report draft says it is due **2026-11-25**. Is this restart for that submission, or for later work? Track A is feasible in days; Track B needs more.
2. **Adopt Track A as primary?** If yes, which reference axes (MASH only, + PPARα, + fasting if a dataset is found)?
3. **Include Track B?** If yes, in full or only as a small direction-of-effect check plus the depth QC?
4. **Formally drop exposure classification** (Track C) as a model aim?
5. **Reword the novelty claim:** "first *reported and interpreted* diet × plastic interaction, with contamination and depth checks", not "first interaction test".
6. **Resolve the Ppara discrepancy before any Ppara claim.** On the same counts, Jung's tables give +0.53 log2 (padj 0.011) for PE under the control diet; our model gives +0.16 (FDR 0.62).
7. **Housekeeping:** commit `research/` and `presentation/` (currently only on one laptop); obtain the Chiu 2025 and Xu 2024 full texts through the library; decide whether to email Jung's authors about sacrifice times, food intake and the Xenium slide map.

*Supporting files in `research/roadmap_phase1/`: `existing_work_inventory.md`, `paper_methods_matrix.csv`, `dataset_feasibility_screen.csv`, `candidate_model_tracks.md`, `references.md`.*
