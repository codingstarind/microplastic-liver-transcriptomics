# Research phases and follow-up analyses

These folders are listed in the order the work was done. Each folder has its own README with the exact commands, inputs, outputs and decision documents. Earlier folders are read-only inputs to later ones, and no earlier output was overwritten. Folders whose names start with `roadmap_` are the planned phases; the others are focused follow-ups.

| Order | Folder | Dates (2026) | Content | Used in report |
|---|---|---|---|---|
| 1 | [`literature/`](literature/) | 11-20 | Literature review, research-gap audit, dataset audit, paper notes | §1 |
| 2 | [`roadmap_phase2/`](roadmap_phase2/) | 11-22 | Reference dataset discovery, download, QC, gene universe | §2.1 |
| 3 | [`roadmap_phase3/`](roadmap_phase3/) | 11-22 | Reference signatures (MASH, FAST, CLOCK, PPARα), positive control, baselines | §2.1 |
| 4 | [`roadmap_phase4/`](roadmap_phase4/) | 11-22 | Gene-level regression (F4/F8; OLS, Ridge, Elastic Net, CatBoost, HGB) on 20 co-expression modules; locked protocol | §3.4–3.5 |
| 5 | [`roadmap_phase4b/`](roadmap_phase4b/) | 11-22 | Feature blocks A/B/C, FULL (28 features), ablation | §3.4, Table 4 |
| 6 | [`roadmap_phase4c_diagnostics/`](roadmap_phase4c_diagnostics/) | 11-22 | Collinearity, leakage and error diagnostics; out-of-fold predictions | §3.4 |
| 7 | [`roadmap_phase4d/`](roadmap_phase4d/) | 11-22 to 11-23 | Per-mouse signature scores (singscore) | §3.2 |
| 8 | [`roadmap_phase4e/`](roadmap_phase4e/) | 11-23 | Cross-study MASH-state classifier; frozen scores on Main PE | §3.3, Table 3 |
| 9 | [`roadmap_phase5/`](roadmap_phase5/) | 11-23 | External validation of frozen models (Hong-PE and others); locked plan | §3.5, §4.2 |
| 10 | [`roadmap_phase5b/`](roadmap_phase5b/) | 11-23 | Transfer audit and diagnosis (post-external, exploratory) | §5.3 |
| 11 | [`enrichment_audit/`](enrichment_audit/) | 11-24 | Gene-mapping audit and fry/camera sensitivity for 50 Hallmark + 3 PPAR sets | §3.1 |
| 12 | [`pathway_followup/`](pathway_followup/) | 11-25 | Eight pathway scores and tests; classifier pathway ablation | §3.2, Table 2, §4.2 |
| 13 | [`lipid_prediction_followup/`](lipid_prediction_followup/) | 11-25 | Lipid and cholesterol gene groups vs other genes | Table 4 |
| 14 | [`nonlinear_chol_followup/`](nonlinear_chol_followup/) | 11-25 | CatBoost and HGB vs linear models for 21 cholesterol-synthesis genes | Table 5 |
| 15 | [`interaction_followup/`](interaction_followup/) | 11-25 | Interaction diagnostics with permutation controls | §4.2 |
| 16 | [`section3_1_figures/`](section3_1_figures/) | 11-25 | Report Figures 1 and 2 | Fig. 1–2 |

**Interpreters.**
- Python steps use `.venv`.
- Steps that fit CatBoost use the interpreter that has catboost 1.2.10 installed (`C:\Python314\python.exe` in our runs).
- R steps run in the WSL conda environment `mpliver-r`.

Versions: see the Software section of the [main README](../README.md#software).

**Large inputs.** Reference-study matrices are downloaded by `roadmap_phase2/scripts/` into `D:/mpliver_work/phase2_refs/`, which can be overridden with the environment variable `MPLIVER_P2_DATA`. They are not committed.

**Kept locally, not in the repository.** Four side folders that the report does not use are kept on the team's machine only: the research-strategy notes (`roadmap_phase1`), a published-gene replication check (`cross_study_replication`), the evidence index used while writing (`paper_evidence_audit`) and a paused external study (`roadmap_phase5c`, PRJNA916285, nothing scored). So are the `*.zip`/`*.7z` archive copies.
