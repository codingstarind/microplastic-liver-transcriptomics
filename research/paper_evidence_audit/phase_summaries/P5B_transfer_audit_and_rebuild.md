# Roadmap Phase 5B: cross-study transfer audit, diagnosis and targeted rebuild (2026-11-23)

## Main sources (`research/roadmap_phase5b/`)
- `README.md`
- `01_pipeline_equivalence_audit.md`
- `02_domain_shift_and_error_diagnosis.md`
- `03_resampling_and_nested_cv_protocol.md` (LOCKED)
- `04_targeted_interventions.md` (LOCKED)
- `05_results_and_generalization_verdict.md`
- `phase6_readiness_memo.md`
- `pipeline_crosswalk.csv`
- `tables/`: `equivalence_checks`, `count_unit_per_sample`, `gap_*`, `intervention_verdicts`, `cv_fold_metrics`, `bootstrap_summary`, `labelling_ranks`, `mash_classifier_bootstrap`, `hyperparameter_choices`, `i1_unit_harmonisation_effect`, `determinism_check`, `candidate_coefficients_per_SD`
- `results/protocol_lock.json`, `results/hong_resampling/`, `results/intervention_targets/`, `results/models/`
- `figures/fig1–fig4`
- Scripts: `b00_common.py` through `b07_figures.py`
- Tests: `tests/test_phase5b.py` (10 tests)

## 1. Question
What exactly failed to generalise in Phase 5, and why? Was it a pipeline error? Do principled, pre-locked fixes help?

**Every Hong and GSE245069 number in this phase is post-external exploratory.** Those outcomes were already seen in Phase 5.

## 2. Stage 1: pipeline equivalence (`01_…md`)
- **Checked and found correct (16 automated checks):** sign, IDs, feature order, leakage, estimator and gene mapping. Mapping coverage was 99.99% (Hong) and 99.74% (GSE245069).
- **M1, count unit:** Hong counts are about 2× the fragment-equivalent unit.
  - Counts per read pair: GSE294180 0.775, GSE245069 0.734, Hong 1.641.
  - Correcting with × 0.5 changes the DEG count from 1,133 to 951 but leaves ρ essentially unchanged (0.219 → 0.220).
- **M2, outlier-dominated training target:** the top 1% of genes carry 77% of the squared effect (Hong 33%). r(MLE, mean-of-logs) = 0.63; excess kurtosis 12.4.
- **M3, external-target uncertainty was not propagated in Phase 5:** transfer ρ ranges 0.10–0.31 across 15 mouse-resampled DESeq2 refits.

## 3. Stage 2: gap diagnosis (`02_…md`)

**Train vs CV vs external** (`gap_train_cv_external.csv`):

| Model | Train ρ | CV ρ | Hong ρ | GSE245069 ρ |
|---|---|---|---|---|
| F8-OLS (CD-trained) | 0.274 | 0.264 | 0.219 | −0.069 |
| F8+C-Ridge (CD-trained) | 0.281 | 0.270 | 0.142 | −0.080 |
| F8-OLS (MD-trained) | 0.324 | 0.309 | 0.303 | 0.181 |
| F8+C-Ridge (MD-trained) | 0.389 | 0.374 | 0.417 | 0.227 |

Train ≈ CV, so the models are high-bias and do not overfit.

**Target noise:**
- Hong SD 0.54 vs 0.33 (GSE294180).
- Robust jackknife reliability: Hong 0.66; GSE294180 0.55.
- Between-study PE vectors agree at ρ 0.10, about **0.17** after disattenuation.

**Where the errors are:** in high-SE and low-expression genes (ρ ≈ 0.06–0.10).

**Mouse PCA of Hong:**
- PC1 accounts for 18% of variance.
- PC1–PC3 explain **89%** of transfer ρ across the 20 labellings.
- The true labelling ranks **5th of 20**.

**Baseline shift:** control FAST singscore is −0.112 (GSE294180 CD+Veh), +0.101 (Hong) and +0.044 (GSE245069).

**Ranked causes:**
1. a non-specific shared mouse-state axis;
2. target noise with 3 v 3 mice;
3. a weak predictor;
4. domain shift that cannot be separated;
5. the outlier-dominated source target;
6. the count unit (minor);
7. model misspecification (weak evidence).

## 4. Stage 3: locked resampling protocol (`03_…md`)
- Reuse the Phase 4 nested grouped CV.
- **Mouse-level paired bootstrap:** B = 1,000; resample mice within arm in both GSE294180 and Hong; refit each candidate with its selected hyperparameters.
- Exhaustive relabelling of the 20 Hong labellings.
- Module-block gene bootstrap as a descriptive check only.
- **MASH classifier bootstrap:** 2,000 stratified resamples of the saved 4E external predictions.
- No multi-study or leave-one-study-out redevelopment. There are only 2 PE labs, and the Hong polymers share controls.

## 5. Stage 4: interventions (locked before scoring; `04_…md`)
- **Comparators:** B0 F8-OLS; B1 F8+C-Ridge; B2 HGB-F8 (CatBoost is not installed in `.venv`).
- **I1:** Hong count unit × 0.5 (sensitivity factor 0.472).
- **I2a:** robust mean-of-logs source target with a moderated SE (d0 = 4).
- **I2b:** diet-averaged PE MLE target from all 19 mice. **Post-hoc:** motivated by the Phase 5 observation that the MD-trained model transferred better.
- **I2c:** robust estimator + diet-averaging.
- **I3:** Huber regression on rank-based normal scores.
- **Acceptance criteria:**
  - A1: Δρ ≥ +0.05;
  - A2: bootstrap P(Δρ > 0) ≥ 0.90;
  - A3: internal CV ρ no more than 0.02 below B0;
  - A4 (specificity): true labelling ranked 1st of 20, or control-free ρ > 0 with P > 0.90.

## 6. Stage 5: results (`05_…md`, `tables/intervention_verdicts.csv`)

| Candidate | Train ρ | CV ρ (SD) | CV R² | Hong ρ | Δ vs B0 | Mouse-bootstrap 95% CI | P(Δ > 0) | Rank /20 | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| B0 F8-OLS | 0.274 | 0.264 (0.040) | 0.066 | 0.219 | — | **−0.03 to 0.30** | — | 4 | comparator |
| B1 F8+C-Ridge | 0.281 | 0.270 (0.040) | 0.081 | 0.142 | −0.076 | −0.06 to 0.29 | 0.14 | 5 | comparator |
| B2 HGB-F8 | **0.378** | 0.276 (0.046) | 0.079 | 0.206 | −0.012 | −0.04 to 0.30 | 0.37 | 4 | variance gap |
| I2a robust | 0.255 | 0.246 | 0.063 | 0.217 | −0.002 | — | 0.56 | 3 | no gain |
| I2b diet-averaged | 0.242 | 0.234 | 0.052 | **0.293** | **+0.074** | 0.01 to 0.35 | 0.87 | 4 | fails A2, A3 |
| I2c robust + averaged | 0.224 | 0.217 | 0.048 | 0.281 | +0.062 | 0.02 to 0.35 | 0.88 | 4 | fails A2, A3 |
| I3 Huber | 0.276 | 0.265 | n/a | 0.220 | +0.001 | — | 0.21 | 4 | no gain |

- **No candidate met all locked criteria.**
- **Mechanism of the I2b gain:** the CIRC coefficient shrinks from −0.042 to about −0.004. This is consistent with the CD+Veh clock anomaly, but it is post-hoc.
- **Ridge:** α = 1000 again, at the grid maximum, with predictions equal to OLS.

**MASH classifier bootstrap** (`mash_classifier_bootstrap.csv`):
- AUC = 1 in 100% of replicates for the shortlisted models.
- S2 balanced accuracy: R0 0.81 [0.63, 0.94]; R1 0.94 [0.81, 1.00]; R2 1.00.
- S2 Brier: R2 0.028 [0.008, 0.056].

**Determinism:** a rerun reproduced every output exactly (maximum difference 0).

## 7. Verdicts
| Question | Verdict |
|---|---|
| Pipeline mismatch? | Technical mismatch fixed (count unit); negligible effect on transfer. |
| Did the fixes help? | Internal improvement only / inconclusive external power. |
| Gene-level generalisation? | **Inconclusive external power.** |
| PE-specific biological transfer? | Not detected. This is **not** proof of absence. |

**Deviation:** the protocol specified 6 workers and the run used 10. The results are seed-identical, and determinism was verified.

## 8. Transition
The memo states that no untouched validation study exists, and it gives a design for one. Phase 6 (reporting) was not started. Phase 5C (new datasets) is ongoing and is **excluded from the current paper**.
