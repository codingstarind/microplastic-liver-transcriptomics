# Roadmap Phase 1: problem framing, literature and dataset screen (2026-11-21)

## Main sources
- `research/roadmap_phase1/phase1_decision_memo.md`
- `existing_work_inventory.md`
- `candidate_model_tracks.md`
- `references.md`
- `paper_methods_matrix.csv`
- `dataset_feasibility_screen.csv`
- `working_notes/` (agent notes; the team has not checked them)

## 1. Question
What problem should the research restart investigate, and with which data? No new analysis was run in this phase.

## 2. Why it was necessary (what P0 had shown)
- **The similarity claim was circular.** P0's "PE resembles MASH" test (split-reference ρ ≈ 0.33) used a diet axis built from the same 19 mice.
- **There is an alternative explanation.** The PE pattern (lipogenesis and cholesterol synthesis down, PPARα targets up) is also produced by fasting or time of day. CD+Veh mice have a distinct clock state.
- **The literature is inconsistent on direction.**
  - Lipogenic genes go **down** in Jung (PE), Makdissi (PS) and Chiu (PS).
  - They go **up** in Djouina (PE), Okamura (PS) and Xu P (PS).
- **Interactions are computed but not reported.** Nine published studies use a diet × plastic design, but none reports and interprets a molecular interaction. Jung's Table S7 contains one (76 genes, mostly muscle), but it is never mentioned in the text. This was corrected in the memo.

## 3. Decision (recommended to the team)
- **Track A (primary):** test PE effects against **independent** reference axes (MASH, PPARα, fasting, time of day) using regression and model selection.
- **Track B (complementary):** compare programme-level directions across independent plastic studies.
- **Track C (exposure classifiers): rejected.** n is far too small, and a classifier would learn sequencing depth or lab. The only acceptable ML is a MASH-state score trained externally.

## 4. Data screened
`dataset_feasibility_screen.csv`:
- GSE294180;
- MASH: GSE316023, GSE246328;
- PPARα: GSE282998;
- plastic studies: GSE245069, GSE277772, and GSE261893 with its siblings;
- **missing at this point:** a fasting or time-of-day reference.

## 5. Literature basis
`references.md`: 14 microplastic–liver papers and 17 methods papers, each with its access level (FT / SI / ABS). Key papers:
- Jung 2026 (source of our dataset; FT + SI);
- Djouina 2023, Li 2025 (GSE245069), Chiu 2025 (GSE277772; ABS only);
- Makdissi 2026; Wu & Smyth 2012 (camera); Wu 2010 (ROAST/fry);
- Bernau 2014; Varma & Simon 2006; Ambroise & McLachlan 2002; Leek 2010; Zhang 2020 (ComBat-seq).

## 6. Limitations stated
1. PE ≡ lower depth and yield.
2. CD+Veh clock state.
3. Muscle carry-over.
4. Missing length offset.
5. Cage ≡ group.
6. External studies can never be pooled into one dataset.

## 7. Open questions the memo left for the team
- Unexplained Ppara gap: Jung +0.53 vs ours +0.16.
- Deadline 2026-11-25.
- Whether to commit `research/`.

## 8. Transition
Track A needs independent reference datasets that are verified, gene-mapped and QC-checked. That is the purpose of roadmap Phase 2.
