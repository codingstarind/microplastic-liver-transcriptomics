# Figure and table recommendations

**Basis:** every figure was inspected visually. Raw notes are in `scripts/figure_notes_raw.md`; the full inventory is in `figure_catalog.csv` (60 rows covering 76 distinct figures plus 182 repeated DESeq2 diagnostic plots).

The proposals below **reuse existing plots and scripts**. Where a plot needs a change, the change is named here. No analysis was altered.

## Main figures (proposed)

### Fig. 1: Study design and analysis workflow (NEW schematic; no data)
- **Content:**
  - the 2 × 2 GSE294180 design (CD/MD × Veh/PE; n = 4/5/5/5);
  - independent reference datasets (MASH ×2, fasting, time of day, PPARα);
  - the three estimands: gene-level prediction, mouse-level signatures, study-held-out classifier;
  - external plastic studies (Hong, GSE245069, GSE277772);
  - which steps were locked and frozen.
- **Why needed:** the project has three different units of analysis, and readers confuse them. A schematic prevents mixing their numbers.
- **Build:** as a vector diagram (draw.io or similar). Take every n from `research_progression.md`.

### Fig. 2: Data quality and confounding in GSE294180
| Panel | Content | Source figure |
|---|---|---|
| A | library size and quality | F001 + F006 |
| B | PCA before/after muscle-gene removal | F010 + F013 |
| C | PC1 loadings (muscle genes) | F011 |
| D | muscle marker screen | F009 |
| E (optional) | interaction genes track muscle | F018 |

- **Message:** PE ≡ lower depth and quality; muscle carry-over; diet is the dominant biology.
- **Rebuild:** recombine the existing PNGs as one vector figure from the existing scripts (`scripts/04–07`). Use one colour scheme.

### Fig. 3: Independent reference axes and the positive control
| Panel | Content | Source figure |
|---|---|---|
| A | reference similarity heatmap | F033 |
| B | positive control | F030; add the exact-permutation null (12,870 allocations, `roadmap_phase3/results/posthoc_audit/exact_permutation.json`) |
| C | fasting reference ZT dependence | F029, panel 3 |

- **Caption must say:** "locked criterion c1 failed; MASH axis treated as exploratory".

### Fig. 4: Gene-level prediction and feature engineering (the main ML figure)
| Panel | Content | Source |
|---|---|---|
| A | model-comparison matrix: R²_oos and Spearman for 12 models × 3 targets, with NB CIs | `tables/T01` + `T02` (replaces F034/F037 visually) |
| B | representation × model R² heatmap | F039 |
| C | block ablation | F040 |
| D | collinearity, VIF and condition number | F042 |

- **Message:**
  - multi-axis linear ≈ best tree models;
  - FULL helps linear models but only ties CatBoost;
  - FULL coefficients are unidentifiable.
- **Rebuild A** from the tables. This is a new plot, but it uses existing numbers only.

### Fig. 5: Per-mouse signatures and the MASH-state classifier
| Panel | Content | Source figure |
|---|---|---|
| A | FAST/MASH singscores per mouse with contrast CIs | F045 + F047 |
| B | study-held-out classifier calibration (before/after) | F050 |
| C | R2 classifier decomposition | F055 |

- **Message:**
  - a fasting-like shift with PE in both diets;
  - AUC is at its ceiling, and engineering improved threshold transfer;
  - the "MASH-state" PE shift is cholesterol/PPARα-driven.
- **Retitle F049** if it is used.

### Fig. 6: External studies and why transfer is not PE-specific
| Panel | Content | Source |
|---|---|---|
| A | external per-mouse FAST scores (Hong) | F053, top-left panels |
| B | gene-level transfer per polymer with the relabelling null | F054, right panel; replace the gene-block CI with the 5B mouse-bootstrap interval |
| C | 20 relabellings vs PC1 | F057, panel 3 |
| D | pathway concordance heatmap | F056 |

- **Message:**
  - the directionally consistent fasting/cholesterol-type shift also appears for other polymers;
  - gene-level agreement tracks between-mouse variation.

## Main tables (proposed)
| Table | Content | Source |
|---|---|---|
| **T1** | Datasets: role, n, design, platform, independence, caveats | P2 dataset selection report; P5 `dataset_provenance.csv` |
| **T2** | Gene-level model comparison | `model_comparisons.md` §1 |
| **T3** | Feature engineering vs F8 and vs the strongest tree baseline | §2 |
| **T4** | Per-mouse signature contrasts, GSE294180 and external | `enrichment_and_signatures.md` §3–4 |
| **T5** | Study-held-out MASH classifier, before/after, with bootstrap | `model_comparisons.md` §4 |
| **T6** | External transfer and post-external interventions (explicitly labelled exploratory) | `model_comparisons.md` §5 |

## Supplementary (suggested set)
- **S1:** EDA panels F002–F005, F007, F008, F012, F014, F015.
- **S2:** DESeq2 diagnostics (F022).
- **S3:** pathway forest (F016) and split-reference concordance (F017, retitled).
- **S4:** reference QC 6-panel (F024) and reference diagnostics (F026–F028, F031, F032).
- **S5:** Phase 4 LOO and variant stability (F035, F036); Task B demo (F038).
- **S6:** Phase 4B per-split stability (F041); Phase 4C residuals and influence (F043, F044).
- **S7:** Phase 4D technical covariates (F046) and the sensitivity forest (rebuilt F048).
- **S8:** Phase 4E permutation null (F051) and GSE294180 classifier scores (F052).
- **S9:** Phase 5B domain shift (F058), interventions (F059, labelled exploratory) and coefficients (F060).
- **S10:** cross-study gene panel forest (F021).

## Exclude
- F019 and F020 (8-gene heatmap/PCA; descriptive selection).
- F023 (182 sensitivity-fit diagnostics; cite the folder).
- F037 (superseded by the model-comparison matrix).

## Missing visualisations worth building (existing numbers only)
1. **Study-design / workflow schematic** (Fig. 1).
2. **Model-comparison matrix** (Fig. 4A) from `T01`/`T02`: models × targets, R²_oos with NB CI, plus a "better than F8-OLS?" marker.
3. **Cross-study validation summary:** one forest plot of the frozen FAST effect across GSE294180 PE_in_CD, PE_in_MD, Hong PE/PS/PP/PB and GSE245069. Data: `roadmap_phase5/sample_level_external_results.csv`.
4. **Evidence-strength ladder:** a table-like graphic of internal → held-out gene → held-out study → independent mice, with the result at each level. Optional; text-heavy.
5. **Rebuilt Phase 4D sensitivity forest** (F048) with CIs from `roadmap_phase4d/results/contrast_estimates_sensitivity.csv`.

## Fixes needed before reuse
1. **F049 title is false.** Retitle, e.g. "AUC is at ceiling for most models; trees vary".
2. **F017 title** implies an independent MASH axis. Use "same-study diet effect".
3. **F054 error bars** are gene-block CIs (anti-conservative). Show the mouse bootstrap (−0.03 to 0.30) or the permutation null only.
4. **F045 / F046 colour schemes differ.** Unify the 4-group palette across all figures.
5. **Remove "PROVISIONAL"** from P0 figure titles once the team accepts assumptions A1–A3.
6. **Export all main figures as vector** (PDF/SVG). Most roadmap figures exist only as PNG.
