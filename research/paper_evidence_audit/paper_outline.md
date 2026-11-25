# Paper outline: section skeleton with phase mapping (no prose drafted)

**Working title (team to decide):** *Re-examining polyethylene-associated liver transcriptomics with independent reference states, held-out validation and external microplastic studies.*

**Scope:** P0 + Roadmap Phases 1–5B. Phase 5C is excluded.

## Abstract
Write last. Components:
- background;
- aim;
- data;
- approach (three estimands);
- the main numbers (one per estimand);
- the external result;
- limitations (PE ≡ depth; 3 v 3).

## 1. Introduction
| Paragraph | Content | Source phase | Key citations (`citation_map.md`) |
|---|---|---|---|
| 1.1 | Microplastics and liver; the source study's claims (PE hepatotoxicity, MASH aggravation, PPARα–Anxa2) | P1 | jung2026, horvatits2022* |
| 1.2 | Conflicting directions across plastic studies; diet × plastic designs without interpreted interactions | P1 | djouina2023, okamura2023, makdissi2026, chiu2025, ahn2025 |
| 1.3 | Methodological problem: same-data disease signatures (circularity) and the fasting/time-of-day alternative explanation | P0, P1 | – |
| 1.4 | Aims: (i) factorial re-analysis with QC; (ii) comparison against independent reference states; (iii) held-out gene and mouse-level evaluation; (iv) external validation in independent plastic studies. State the course context | P1 | – |

## 2. Materials and Methods
| Subsection | Content | Phase | Evidence files |
|---|---|---|---|
| 2.1 Data | GSE294180 design; reference datasets; external plastic studies (Table 1) | P0, R2, R5 | `phase_summaries/P2_*`, `roadmap_phase5/dataset_provenance.csv` |
| 2.2 Preprocessing and QC | filters; normalisation; EDA/PCA; muscle screen; depth/yield; quantification pilot (A1–A3) | P0, R2 | `methods_inventory.csv` M01–M10 |
| 2.3 Factorial DE and pathway tests | DESeq2 ~ diet*PE; contrasts; BH; fry/camera/fgsea; robustness | P0 | M11–M13, M18–M21 |
| 2.4 Independent reference signatures | MASH, fasting, time of day, PPARα models; positive-control criteria (locked) | R3 | M14, M15, M22 |
| 2.5 Gene-level prediction | targets; weights; features F8; module-grouped nested CV; models; metrics; NB-corrected comparisons; variants, LOO, negative controls | R4 | M29–M32, M37, M38, M45 |
| 2.6 Feature engineering and diagnostics | blocks A/B/C; ablation rule (exploratory); collinearity/VIF; residual and target-quality diagnostics; freezing rule | R4B, R4C | M27, M28, M35, M46 |
| 2.7 Per-mouse signature scoring | signature construction; singscore; 2 × 2 OLS; Holm; exact permutation; null gene sets | R4D | M23, M24, M41 |
| 2.8 Study-held-out MASH-state classifier | representations R0/R1/R2/PCA; models; inner CV; metrics; validity rule; application | R4E | M34, M36, M39, M42 |
| 2.9 External validation | locked plan; frozen artefacts; external DE; gene- and mouse-level metrics; relabelling null | R5 | M16, M40 |
| 2.10 Transfer audit | equivalence checks; count-unit; gap diagnostics; mouse bootstrap; interventions (post-external) | R5B | M17, M43, M47–M49 |
| 2.11 Reproducibility and AI use | locks/hashes; seeds; software versions; tests; AI usage declaration | all | `reports/ai_usage_declaration.md` |

## 3. Results
| Subsection | Message (evidence, not wording) | Figure / table | Phase |
|---|---|---|---|
| 3.1 Data quality and confounding | PE ≡ depth; muscle carry-over defines PC1; diet dominates biology | Fig. 2, F001/F006/F009–F013 | P0 |
| 3.2 Factorial effects | 858 / 472 / 2,262 / 87 DEGs; lipogenic and cholesterol genes ↓ with PE in both diets; interactions mostly muscle; no pathway interaction | F016, F018 | P0 |
| 3.3 Reference states | References reproducible (MASH1 vs MASH2 ρ 0.56) and near-orthogonal; fasting ≠ PPARα for lipogenesis; positive control FAIL on the locked criterion | Fig. 3, F025, F029, F030, F033 | R3 |
| 3.4 Baselines | PE–MASH fragile; PE–FAST small and stable; DEG overlap uninformative | F025, Table | R3 |
| 3.5 Gene-level prediction | R² ≈ 0.07 robust; trees not better on PE (better on diet) | Fig. 4A, T2 | R4 |
| 3.6 Feature engineering | FULL linear 0.104 vs F8 0.066; tied with CatBoost; block C useful, A useless, B target-informed; FULL coefficients unidentifiable (VIF 233); noise ceiling | Fig. 4B–D, T3 | R4B, R4C |
| 3.7 Per-mouse signatures | FAST +0.065 / +0.056; MASH PE shift fragile | Fig. 5A, T4 | R4D |
| 3.8 Cross-study MASH-state classifier | AUC ceiling; threshold transfer improved (0.81 → 1.00); PE readout cholesterol-driven | Fig. 5B–C, T5 | R4E |
| 3.9 External studies | Gene level not passed (p 0.25; GSE245069 fail); FAST directionally consistent in Hong (+0.060) but not PE-specific; pathway directions | Fig. 6, T6 | R5 |
| 3.10 Why transfer was inconclusive | no pipeline error; rank transfer ≈ CV; agreement tracks mouse PCs; mouse bootstrap −0.03 to 0.30; interventions fail locked criteria (post-external) | Fig. 6C, F057–F059 | R5B |

## 4. Discussion
| Paragraph | Content | Evidence |
|---|---|---|
| 4.1 Principal findings | A fasting/cholesterol-type rather than MASH-type pattern, reproducible internally and directionally present in one independent PE study, but not PE-specific | `enrichment_and_signatures.md` §5 |
| 4.2 Alternative explanations | library depth/yield; feeding/sacrifice time (CD+Veh clock); cell composition; shared controls | `limitations_and_corrections.md` A1–A3, A9 |
| 4.3 Comparison with the source study and literature | lipogenic directions; Ppara discrepancy; interaction artefacts; "compounded" MASH effect not supported | P0, P1, P3 |
| 4.4 Methodological lessons | circular references; DEG overlap; unit of analysis (genes vs mice); held-out genes ≠ new mice; AUC ceilings vs calibration; feature engineering vs interpretability; post-hoc vs frozen | `model_comparisons.md` §0 |
| 4.5 Limitations | A1–A13 | `limitations_and_corrections.md` |
| 4.6 What would resolve it | untouched validation design (≥ 8 per arm, own controls per polymer, fixed ZT, pre-registered frozen endpoints) | `roadmap_phase5b/phase6_readiness_memo.md` |

## 5. Conclusions
Three to four bullet points, each tied to one estimand and one limitation.

## Back matter
- Data and code availability.
- AI usage declaration.
- Author contributions ([TEAM]).
- Supplementary figures S1–S10 (`figure_recommendations.md`) and supplementary tables (`tables/T*.csv`).

## Phase → section map
| Phase | Methods | Results | Discussion |
|---|---|---|---|
| P0 | 2.2, 2.3 | 3.1, 3.2 | 4.2, 4.3 |
| R1 | 1.2–1.4 (motivation) | – | 4.3 |
| R2 | 2.1, 2.2 | (Table 1) | 4.5 |
| R3 | 2.4 | 3.3, 3.4 | 4.1, 4.4 |
| R4 | 2.5 | 3.5 | 4.4 |
| R4B | 2.6 | 3.6 | 4.4 |
| R4C | 2.6 | 3.6 | 4.4, 4.5 |
| R4D | 2.7 | 3.7 | 4.1 |
| R4E | 2.8 | 3.8 | 4.1, 4.4 |
| R5 | 2.9 | 3.9 | 4.1, 4.2 |
| R5B | 2.10 | 3.10 | 4.4, 4.6 |

\* `horvatits2022`: verify the polymer list before citing.
