# Roadmap Phase 3: independent reference signatures, MASH positive control and conventional baselines (2026-11-22)

## Main sources (`research/roadmap_phase3/`)
- `analysis_plan.md` (locked; SHA-256 in `results/plan_lock.json`)
- `reports/phase3_decision_report.md`
- `reports/phase3_review_addendum.md` (post-hoc audit)
- `reports/phase4_specification_memo.md`
- `reports/{mash,fasting,circadian,ppara}_reference_report.md`
- `reports/{reference_comparison,baseline_comparison,count_offset_validation}_report.md`
- `reference_effects/*.csv.gz`
- `results/**`
- `figures/*.png`
- `tests/test_phase3.py`

## 1. Questions
1. Can reliable signatures for MASH, fasting, time of day and PPARα activation be estimated from **independent** mice?
2. Does our own diet effect (Diet_in_Veh) agree with an independent MASH axis (the positive control)?
3. What do conventional baselines (DEG overlap, effect correlation) say about PE?

## 2. Why
P0's "PE resembles MASH" test was circular: same mice and a shared reference. Phase 2 supplied verified independent datasets, a fixed universe of 11,888 genes, and a 309-gene muscle list.

## 3. Methods implemented (scripts in `roadmap_phase3/scripts/`)
| Reference | Model | Script |
|---|---|---|
| MASH1 GSE316023 (8 v 8) | DESeq2 `~ diet` on rounded estimated counts; no length offset (primary); FPKM-implied `avgTxLength` offset as sensitivity; TN 29 °C pair as sensitivity | `a01_offset_audit.py`, `b01_mash_references.R` |
| MASH2 GSE246328 (7 v 8) | DESeq2 `~ diet`; leave-one-out (LOO) | `b01` |
| FAST GSE107787 (18 v 18) | limma-trend on log2(RPKM + 1), `~ ZT + state`; ZT × state interaction check; clean subset ZT0 + ZT4 muscle-low (4 v 6) | `d01_other_references.R` |
| FAST2 GSE225695 (3 v 3) | direction check | `d01` |
| CIRC GSE117134 (female, 6 ZT × 3) | limma-trend harmonic `~ cos + sin`; rhythmicity F-test; amplitude/phase; signed contrasts T1 = ZT6 − ZT2 and T2 = ZT12 − ZT0 | `d01` |
| PPARA GSE282998 | DESeq2 `~ dose` (factor) on TempO-Seq probe counts: 5/10/20 v 0, trend, muscle-adjusted; genotype × dose interaction (WT − KO) | `d01` |
| Ppara discrepancy | DESeq2 configurations on GSE294180 Ppara | `e01_ppara_discrepancy.R` |
| Positive control | Spearman, sign agreement, hypergeometric DEG overlap, fry pathways, relabelling null (1,000 draws), LOO, offset, strict universe | `c01_positive_control.py` |
| Comparisons and baselines | pairwise Spearman matrix, VIF, DEG overlap (hypergeometric, Jaccard), effect correlation with LOO ranges | `f01_compare_and_baselines.py` |
| Post-hoc audit | exact enumeration of C(16, 8) = 12,870 allocations; pre-registered narrower muscle rule (88 genes) | `g01_exact_permutation_audit.py`, `g02_*` |

**Multiple testing:**
- genes: BH within each (dataset, contrast);
- pathways: BH within each (contrast, family);
- overlaps and correlations: descriptive.

## 4. Key results
### Reference signatures
| Signature | Result |
|---|---|
| MASH1 | 6,260 genes padj < 0.05 |
| MASH2 | 5,252 genes padj < 0.05 |
| MASH1 vs MASH2 | ρ = 0.557; 90% sign agreement among 3,306 shared significant genes; 24/24 pathways agree |
| FAST | 8,552 padj < 0.05; 4,296 genes with ZT × state FDR < 0.05; 45% of fasting DEGs are ZT-dependent; clean-subset ρ 0.919; GSE225695 ρ 0.489 |
| CIRC | 2,914 rhythmic genes (female); median phase difference vs male data 1.4 h; T1 vs T2 ρ = 0.86 |
| PPARA WT 20v0 | 3,604 padj < 0.05; ρ vs dose trend 0.978; vs muscle-adjusted 0.946 |
| PPARα dependence | 1,893 genes formally dependent (the stocks differ, so this is provisional except for canonical targets). Canonical targets WT vs KO: Cyp4a14 +5.7 vs +1.7; Acot1 +5.8 vs −0.3 |

### Programme medians (log2)
| Programme | Fasting | GW7647 |
|---|---|---|
| Lipogenesis | −2.79 | +1.15 |
| Cholesterol synthesis | −2.76 | +0.67 |
| PPARα targets | +0.82 | +1.50 |

"PPARα targets up with lipogenesis down" is therefore a **fasting-like** combination, **not** an agonist-like one.

### Collinearity
- VIF: MASH1 1.03, FAST 1.07, PPARA 1.02, CIRC T1 4.07, T2 4.04.
- MASH vs FAST ρ 0.06; FAST vs PPARA 0.04.

### Positive control (Diet_in_Veh vs MASH1)
| Measure | Value |
|---|---|
| ρ (11,886 genes) | 0.355 |
| Sign agreement (1,247 genes) | 0.787 |
| DEG overlap | 747 (p = 4 × 10⁻³¹) |
| Pathways | 2/3 core agree; **lipogenesis opposite** (GAN vs chow ↑, our MD vs purified CD ↓) |
| Monte Carlo null q99 | 0.382 → criterion c1 **not met** → **formal FAIL** |
| Exact enumeration | p = 0.036; q99 0.390; balanced 4 + 4 maximum 0.288 |
| Strict universe | ρ = 0.252 |
| MASH2 vs Diet | 0.323 |

**Verdict:** the MASH axis may be used as **EXPLORATORY** only.

### Baseline 2: effect correlation (Spearman, [LOO range])
| Reference | PE_in_CD | PE_in_MD |
|---|---|---|
| MASH1 | 0.211 [**0.03**, 0.24] | 0.157 [**0.02**, 0.24] |
| MASH2 | 0.027 | 0.041 |
| FAST | 0.135 [0.10, 0.20] | 0.254 [0.22, 0.27] |
| PPARA 20v0 | −0.048 | −0.086 |
| CIRC T1 | −0.109 | 0.106 |

The PE–MASH correlations depend on single mice: without CD_MP_2 or MD_MP_4 they fall to 0.03 and 0.02.

### Baseline 1: DEG overlap
Almost every overlap is "significant" by the hypergeometric test, even when most shared genes change in opposite directions:
- PE_in_CD vs MASH2: 84 concordant / 125 discordant;
- vs CIRC T2: 19 / 90.

Overlap is therefore **not informative** about state similarity.

### Ppara discrepancy
Every DESeq2 configuration gives +0.16 to +0.17 for PE_in_CD, against Jung's +0.53. DESeq2 settings do not explain the gap.

## 5. Failures and limitations
- **MASH positive control: formal FAIL.**
- The MASH references use a chow control, so lipogenesis disagrees.
- **Fasting:** RPKM only; time-dependent; a 24 h fast is extreme.
- **Circadian:** female; RPKM.
- **PPARα:** targeted FFPE assay; dose-linked batch; different breeding stocks for WT and KO.
- **The muscle list removes liver-relevant genes:** Pdk4, Plin4, Col1a1/2, Lpl.
- **MASH1 library order is aligned with diet.**
- The GSE294180 limitations are unchanged.

## 6. Transition
Phase 3 delivered independent and nearly orthogonal reference vectors (VIF ≈ 1). Single-axis baselines were weak (|ρ| ≤ 0.25). The question became whether a **multi-axis, SE-weighted regression with held-out-gene validation** explains PE effects better than single axes. `phase4_specification_memo.md` fixed this before fitting:
- models M0–M6;
- DL-τ² weights;
- 20 co-expression modules from GSE117134 for blocked CV;
- a decision rule requiring a win in both universes **and** in ≥ 16/19 LOO refits.
