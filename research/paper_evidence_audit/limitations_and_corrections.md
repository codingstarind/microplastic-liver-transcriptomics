# Limitations, corrections, failed experiments and unresolved issues

Each item is split into **technical observation**, **numerical evidence** (with source) and **biological interpretation**. Interpretation is the team's to finalise.

---

## A. Limitations

### A1. PE is confounded with sequencing depth and yield in GSE294180
- **Observation:** every PE library is shallower than every vehicle library. All 19 libraries come from one flowcell.
- **Evidence:**
  - read pairs 13.5–16.5 M (PE) vs 17.0–19.1 M (vehicle) (`reports/final_report_draft.md` §4.1);
  - log2 depth 23.69–23.98 vs 24.02–24.19 (P4D `decision_report.md`);
  - depth/yield covariates have R² ≈ 0.72–0.76 on PE; adjusting cuts PE_in_CD DEGs from 858 to 116 / 66 and removes the PPARα and cholesterol pathway signals (P1 `existing_work_inventory.md` §3);
  - chrM fraction is higher in PE libraries (`figures/eda/fig06`).
- **Interpretation:** no model in any phase can separate PE from library properties.
  - Binomial thinning (P4) tests **depth only**: it produced 0 DEGs and R² ≈ 0, and does not simulate yield or quality.
  - External studies have balanced depth (P5: ratio 0.95–1.0), which is why they matter.

### A2. Skeletal-muscle contamination and filtering sensitivity
- **Observation:** muscle CPM is MD_1 5,190, CD_MP_1 539, MD_2 316, MD_3 161; PC1 (30%) is muscle.
- **Evidence:**
  - 23/87 interaction genes are muscle-tracking. 16 interaction genes remain with the muscle covariate, and 1 in edgeR or limma-voom.
  - The 309-gene exclusion list removes Pdk4, Plin4, Col1a1/2, Lpl, Cpt1b, Cav1 and Vldlr (P3 decision report §8).
  - The alternative 88-gene rule leaves the positive control unchanged (0.355 → 0.358) and the P4D scores identical.
- **Same artefact elsewhere:** GSE107787 FED livers (aligned with feeding state) and GSE282998 (aligned with dose group); the published Jung Table S7 interaction list.
- **Interpretation:** gene-level interactions are not evidence of diet-dependent PE biology.

### A3. Unusual circadian / feeding state of the CD+Veh mice
- **Evidence:**
  - pre-defined clock-phase score: CD+Veh −2.38 vs +0.08 to +1.01 in the other groups;
  - P4D scores vs clock_phase: ρ 0.64 (MASH), 0.55 (FAST);
  - P4E R2 vs clock_phase ρ 0.67;
  - CD+Veh FAST score −0.112, the lowest of any group across studies (P5B `gap_study_shift.csv`).
- **Interpretation:** PE_in_CD estimates partly reflect an atypical reference group. Sacrifice time is unknown. The PE_in_MD estimates do not depend on CD+Veh, and they show the same FAST direction.

### A4. Count units and length normalisation
- **Evidence:**
  - The deposited GSE294180 counts lack a length offset (assumptions A1–A3; the 5-sample pilot was "not identifiable", with totals 4.5–7.2% higher).
  - GSE316023 FPKM-implied lengths correlate only 0.497 with annotation, so the offset was rejected as primary. 883 offset-sensitive genes; no significant sign reversals.
  - Hong counts are about 2× the fragment-equivalent unit (1.641 vs 0.775 counts per read pair). Correcting changed DEGs 1,133 → 951 and transfer ρ 0.219 → 0.220 (P5B).
  - The P4D TPM-proxy uses median transcript length, which is constant per gene.
- **Interpretation:** rank-based comparisons are robust to these; magnitude comparisons across technologies are not attempted.

### A5. Small biological sample sizes
| Dataset | n per group | Consequence |
|---|---|---|
| GSE294180 | 4/5/5/5 | interaction CI width 0.80 log2 |
| Hong | 3 v 3 | minimum exact two-sided p = 0.10 |
| GSE245069 | 4 water controls | – |
| MASH references | 7–8 | – |
| Fasting check (GSE225695) | 3 | – |

Single mice are influential:
- dropping CD_MP_2 changes the target by RMS 0.19 log2 and F8-OLS R² from 0.066 to 0.048;
- PE–MASH1 ρ falls from 0.21 to 0.03 without CD_MP_2.

### A6. Uncertainty in the independent reference signatures
- **MASH1** positive control formal FAIL: c1, exact p = 0.036 > 0.01; strict-universe ρ 0.252. The MASH references use a chow (not purified) control, so lipogenesis is opposite.
- **Fasting:** RPKM only, ZT-averaged (45% of DEGs ZT-dependent), extreme 24 h fast.
- **Circadian:** female; RPKM.
- **PPARα:**
  - targeted FFPE assay;
  - dose-linked depth and muscle batch;
  - WT and KO are different breeding stocks, so the 1,893 "PPARα-dependent" genes are provisional;
  - reliability 0.34 (P4B diagnostics).
- **Errors in variables** attenuate every reference coefficient (Fuller 1987).

### A7. Dependence between genes measured in the same mice
- **Evidence:** random-gene CV R² 0.069 vs module-grouped 0.066 (P4C). Gene-block bootstrap CIs (e.g. Hong ρ 0.22 [0.14, 0.28]) are far narrower than mouse-level intervals (−0.03 to 0.30).
- **Interpretation:** the effective n is far below 11,758. Only mouse-level resampling and permutation respect the experimental unit.

### A8. Cross-study differences in diet, dose, route, polymer and duration
| Study | Polymer | Size | Route / dose | Duration | Diet | Strain |
|---|---|---|---|---|---|---|
| GSE294180 | PE | 10–150 µm | gavage, 2 mg/day | 8 wk | purified | – |
| Hong | PE | NOT STATED | in food, dose NOT STATED | 71 wk | grain chow | – |
| GSE245069 | PS | 80 nm / 0.5 µm / 5 µm | water | 12 wk | – | ICR |
| GSE277772 | PS | – | – | – | – | (normalised values only) |

- Disattenuated between-study PE ρ ≈ 0.17.
- Control FAST baselines differ: −0.112 / +0.101 / +0.044.
- **Not separable** (one study per level).

### A9. Reused controls in Hong's polymer experiments
- **Evidence:** CON1–3 count columns are byte-identical across GSE261893, GSE261725, GSE286116 and GSE262399, and the control SRA read counts are identical. So 12 GSMs = 3 mice.
- **Interpretation:** the four polymer contrasts are correlated (expected r ≈ 0.5 under the null), and specificity cannot be established. The shared controls have the lowest FAST mean of all groups.

### A10. Held-out-gene validation is not new-mouse validation
Phases 4–4C validate generalisation across gene modules **within one experiment**. Only Phases 4E (held-out study, MASH state) and 5 (independent plastic studies) test new animals.

### A11. Models that looked good internally did not transfer
- F8+C-Ridge had the highest internal CV ρ (0.270; R² 0.080–0.096) but the lowest Hong ρ (0.142).
- HGB had train ρ 0.378 vs CV 0.276 (variance gap) and Hong 0.206.
- F8-OLS (CD-trained) had GSE245069 ρ −0.069.

### A12. Post-hoc vs frozen analyses
| Status | Analyses |
|---|---|
| **Frozen / locked** | P3 analysis plan; P4 protocol; P4D plan + signatures; P4E plan; P5 locked plan + frozen manifest; P5B protocol (03/04) |
| **Post-hoc / exploratory** | P3 review addendum (exact permutation, alternative muscle rule); all of P4B (features chosen after P4); P4C fix experiment (reused folds); P5B interventions (Hong already seen); I2b diet-averaging (motivated by a seen result) |

### A13. Resampling with only three mice per arm
- Within-arm bootstrap duplicates mice: B0 control-free ρ is 0.002 observed vs a bootstrap median of 0.049.
- Intervals are descriptive.
- Exhaustive relabelling allows a minimum one-sided p of 0.05 (20 labellings).
- The permutation nulls in P3 and P4E have low resolving power when the reference effect dominates.

---

## B. Corrections: later reports corrected earlier ones (use the corrected value)

| # | Earlier statement | Corrected statement | Where corrected |
|---|---|---|---|
| C1 | "No study fitted a diet × plastic interaction" | Jung's Table S7 contains one (76 genes, mostly muscle) but never reports or interprets it | P1 inventory §4 |
| C2 | Jung Ppara induction "+0.53" vs ours | Jung within-CD +0.53 (padj 0.011), within-MD +0.34 (n.s.); ours +0.16 to +0.17 in every DESeq2 configuration. **Unexplained** | P1 §4; P3 `results/ppara_discrepancy/` |
| C3 | GSE316023 length "median CV 1%" | Across the 16 contrast samples, 603 genes differ > 1.2× between GAN and chow | P2 signoff §2 |
| C4 | GSE246328 "8-week GAN group" | Pre-fed 34 wk GAN; the 24 w arm ≈ 58 wk GAN; only the 24 w pair is valid | P2 dataset selection §1.2 |
| C5 | GSE277772 raw counts (two checks disagreed) | Normalised values, not counts | P2 external status |
| C6 | Positive-control balanced 4+4 maximum 0.277 (sampled) | 0.288 (exact enumeration); proposed null replacement **withdrawn**; c1 stays FAIL | P3 review addendum |
| C7 | P4D "exact permutation p 0.040 / 0.024" | Not Holm-adjusted; adjusted over the 6-test family → 0.20 / 0.14 | P4E `00_task_and_estimands.md` |
| C8 | P4E "every non-PCA representation reached AUC 1.0" (text and fig1 title) | False: trees on R2 0.80 / 0.63; R0 tree S2 0.81 | P5 decision report (A) |
| C9 | P4E R0 rule "Brier 0.188" | A hard-rule misclassification rate, not a probabilistic Brier | P5 decision report; P5B |
| C10 | P5 Hong ρ 0.22 "[0.14, 0.28]" | That interval is a gene-block bootstrap (anti-conservative). The mouse-level bootstrap gives −0.03 to 0.30 | P5B `05_*.md` §3 |
| C11 | P5B rank of the true Hong labelling | "5th of 20" (log2CPM proxy, `02_*.md`) vs "4th of 20" (arith proxy, intervention table). Report both with their proxy | this audit |
| C12 | P5B mean-variance jackknife reliability 0.065 | Invalid under heavy tails; the robust jackknife gives Hong 0.66, GSE294180 0.55 | P5B `02_*.md` §2 |
| C13 | P5B sign check of GSE294180 target vs simple mean-of-logs (r = 0.59) | Due to an MLE vs mean-of-logs mismatch; the corrected check (Spearman vs log ratio of group-mean CPM) passes (0.77). This led to the heavy-tail finding | P5B `01_*.md` |
| C14 | P3 offset-audit pass mark ρ ≥ 0.50 | Not in the locked plan; the decision rests on undocumented provenance | P3 README deviations |
| C15 | P5B protocol said 6 workers | The run used 10; outputs are seed-identical and determinism was verified | P5B README |

---

## C. Failed or negative experiments (keep in the paper)
1. **MASH positive control:** formal FAIL on c1 (P3).
2. **DEG-overlap baseline:** "significant" even with mostly discordant directions (P3).
3. **Complex models (RF/HGB/CatBoost):** never beat linear on PE targets (P4).
4. **Block A (reliability shrinkage):** not useful, and it caused collinearity (P4B/4C).
5. **FULL representation:** improved R² but its coefficients are uninterpretable (P4C).
6. **PCA representation:** worst cross-study transfer (P4E).
7. **Phase 5 gene-level external validation:** not passed (Hong permutation p 0.25); failed in GSE245069 (ρ −0.07).
8. **MASH singscore externally:** opposite sign in Hong.
9. **PE specificity:** the FAST shift is present for PS/PP/PB vs the shared controls.
10. **Phase 5B interventions:** none met the locked criteria.
11. **Fasting classifier:** external n = 6, not frozen (P4E).
12. **HGB classifier in P4 Task B:** degenerate at the default min_samples_leaf = 20 (deviation D2).

## D. Unresolved technical issues
- The Ppara discrepancy with Jung.
- Assumptions A1–A3 (deposited counts, no length offset). The full re-quantification of 14 runs was never run.
- Sacrifice time and food intake in GSE294180 are unknown; the authors were not contacted.
- Xenium slide confounding (barcode prefix split by treatment) is open.
- Dose and particle size in Hong are not stated, and there is no linked paper.
- Not installed: GSVA (P4D/4E) and CatBoost in the Phase 5B environment.
- The P0 counts were never independently checked by the team (the final-report draft status box says so).
