# Phase 3 analysis plan (pre-specified)

**Written:** 2026-11-22, **before** any Phase 3 model was fitted. The SHA-256 of this file is recorded in `results/plan_lock.json` before the first model run. Any later change is logged in §11 (Deviations) with a reason; the original text is not edited.

**Fixed inputs (Phase 2, not re-selected):**
- `research/roadmap_phase2/metadata/phase3_gene_universe.csv`: main universe 11,888 genes; sensitivity universe 8,360 genes.
- `metadata/muscle_exclusion_genes.csv` (309 genes; already removed from both universes).
- `metadata/phase3_eligible_samples.csv` and `reference_contrast_plan.csv`.
- Prepared matrices in `D:/mpliver_work/phase2_refs/<GSE>/prepared/`.

**Software:** R 4.4.3 (WSL env `mpliver-r`): DESeq2 1.46.0, limma 3.62.1, edgeR 4.4.0. Python 3.14 (pandas, numpy, scipy, matplotlib) for comparisons and figures. Seed: **20261122**.

**Target effects (GSE294180), not re-estimated:**
- Taken from `results/deseq2/provisional_deposited/results_{Diet_in_Veh,PE_in_CD,PE_in_MD}.csv.gz`: unshrunk MLE `log2FC`, `lfcSE`, `padj`.
- Model: DESeq2 `~ diet * microplastic`, 13,777 filtered genes, deposited counts, no length offset.
- Leave-one-mouse-out refits: existing `results/deseq2/loo/`.

## 1. Common rules
- **Gene sets for reporting and comparisons:** the main universe intersected with the genes each dataset can estimate. Reference effects are reported for **every** universe gene present in the dataset, never only for significant genes.
- **Multiple testing:** BH within each (dataset, contrast) over the universe genes reported. Independent filtering is **off** so that every universe gene has a padj. Pathway tests: BH within each (dataset, contrast, gene-set family).
- **DEG definition (all baselines and the positive control):** padj < 0.05 **and** |log2FC| ≥ 0.38 (1.3-fold; the Phase 2 threshold). Sensitivity: padj < 0.05 only.
- **Effect scale:** log2 fold change for counts; difference of log2(RPKM+1) means for RPKM. Effects from different technologies are compared **only by rank and sign**, never by magnitude.
- **No pooling** of expression matrices across studies.

## 2. Stage A: GSE316023 length-offset decision rule (fixed before any fit)

**Audit checks** (`scripts/a01_offset_audit.py`):
- (i) the count and FPKM matrices have identical samples, genes and zero pattern;
- (ii) implied length L = count·10⁹ / (FPKM · N) is in the bp range and correlates with GENCODE transcript length;
- (iii) between-sample variation is not merely a column constant (identifiable up to sample scaling);
- (iv) the share of universe gene–sample cells where 4-decimal FPKM rounding gives > 1% length error;
- (v) DESeq2 normalisation factors built from L are finite, positive and row-centred (geometric mean 1);
- (vi) rounding of counts changes < 0.01% of any library.

**Decision:**
- **The offset model is primary only if** (i)–(vi) all pass **and** the FPKM provenance is documented by the authors.
- If the provenance is undocumented but (i)–(vi) pass, the **no-offset model is primary** and the offset model is a **co-reported approximate sensitivity model**.

  Reason: an undocumented normalisation cannot be verified as the per-sample effective length that tximport expects. The target contrast (GSE294180) also has no length offset.
- Genes whose call (DEG status or sign at padj < 0.05) differs between the two models are flagged `offset_sensitive`.

## 3. Stage B: MASH references (DESeq2)
- **GSE316023:**
  - samples: `MASH_primary` (GAN wk23 TA ×8 vs chow wk23 TA ×8);
  - `~ diet`, reference level `chow`, coefficient `diet_GAN_vs_chow`;
  - fitted on universe genes present (size factors from those genes);
  - `cooksCutoff` default; outlier replacement as DESeq2 default (≥ 7 replicates) recorded.
- **Outputs per gene:** log2FC (MLE), lfcSE, 95% Wald CI, stat, p, padj, baseMean, ashr-shrunk log2FC (for display only), detected flag (≥ 10 counts in ≥ 4 of the 16 samples), Cook's max.
- **Diagnostics:** dispersion plot, MA plot, Cook's distance per sample.
- **Sensitivities:**
  - offset alternative;
  - `MASH_sens_TN` contrast (8 v 8, 29 °C);
  - leave-one-mouse-out (16 refits; Spearman of LOO log2FC vs full);
  - sensitivity universe.
- **GSE246328:** `MASH_secondary` (Vehicle 24 w ×7 vs chow 24 w ×8), `~ diet`, the same outputs and diagnostics, and LOO.

## 4. Stage C: MASH positive control (pre-specified criteria)

**Inputs:**
- target = GSE294180 `Diet_in_Veh` (MLE log2FC, padj);
- reference = GSE316023 `diet_GAN_vs_chow` from the **primary** Stage A model;
- genes = main universe ∩ non-NA in both.

**Metrics:**
- **M1:** Spearman ρ of log2FC.
- **M2:** sign agreement among genes with padj < 0.05 in **both**.
- **M3:** DEG overlap with a fixed-universe hypergeometric test (one-sided), Jaccard, and direction-concordant count.
- **M4:** pathway agreement in 6 **core sets**: LIPOGENESIS_CORE, CHOL_SYNTHESIS_CORE, SANDERSON_PPARA_TARGETS, HALLMARK_FATTY_ACID_METABOLISM, HALLMARK_INFLAMMATORY_RESPONSE, HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.
  - Method: limma `fry` (self-contained, directional) on log-normalised expression with the dataset's own design; `camera` (competitive) reported alongside.
  - Agreement = same direction in sets with fry FDR < 0.05 in both.
- **M5 uncertainty:**
  - **sample-relabelling null:** 1,000 random balanced relabellings of the reference mice (8 v 8), recomputing the reference effect as the difference of mean log2 normalised counts and its ρ with the target. Label permutation keeps gene–gene correlation, unlike gene permutation;
  - **LOO ranges** in reference (16) and target (19 existing);
  - offset/no-offset;
  - sensitivity universe.

**Acceptance criteria:**

| Outcome | Condition |
|---|---|
| **PASS** | ρ ≥ 0.30 **and** ρ above the 99th percentile of the relabelling null; **and** M2 sign agreement ≥ 0.70; **and** among core sets significant (fry FDR < 0.05) in both, ≥ 2/3 agree in direction (at least 2 such sets required); **and** ρ ≥ 0.20 in every LOO refit and under both offset choices |
| **CONDITIONAL** | the first two conditions hold but the pathway or LOO/offset condition fails |
| **FAIL** | otherwise; diagnosed, not tuned |

**Why these thresholds:**
- ρ ≥ 0.30 is a modest bar for two diet studies that differ in control diet (chow vs purified) and duration (23 vs 8 wk). Cross-model mouse liver comparisons typically show moderate genome-wide correlation, so a higher bar would reject a usable axis for known design reasons.
- 0.70 sign agreement among jointly significant genes is far above the 0.5 chance level and robust to scale differences.
- The pathway rule checks that agreement holds in the programmes we care about, not only genome-wide.

## 5. Stage D: second MASH reference
GSE246328 is compared with GSE316023 and with Diet_in_Veh using M1–M4 and the relabelling null (7 v 8). This is reported as reproducibility, with no pass/fail. It is the benchmark for "how similar are two independent MASH references".

## 6. Stage E: fasting (GSE107787, RPKM)
- **Main model:** limma-trend on log2(RPKM + 1), universe genes present; `~ ZT (factor) + state`, reference `FED`, coefficient `stateFAST`.
  - **Check:** mean–variance trend plot (the trend must be smooth and the residual SD must not be dominated by the RPKM floor).
  - **Interaction:** `~ ZT * state`, F-test of the 5 interaction terms; report the number of genes with interaction FDR < 0.05. If > 10% of fasting-DEGs show interaction FDR < 0.05, the fasting effect is reported as time-dependent and the single signed effect is labelled an *average over ZT*.
- **Clean subset:** `in_clean_sensitivity` (ZT0+ZT4, muscle-low: 4 FAST vs 6 FED), `~ ZT + state`.
- **Stability rule:** the fasting reference is labelled **CONDITIONAL** if Spearman(main, clean) < 0.60 **or** sign agreement < 0.80 among main-model DEGs.
- **GSE225695:** DESeq2 `~ state` (3 v 3), reference `Fed`. **Direction check only:** ρ with the main model, plus sign agreement on main DEGs.
- **Limitations stated:** RPKM loses depth information; no NB model; a 24 h fast is an extreme state.

## 7. Stage F: circadian (GSE117134 WT, female, RPKM)
- limma-trend on log2(RPKM + 1); `~ cos(2πZT/24) + sin(2πZT/24)`.
  - Rhythmicity: F-test on both terms, BH.
  - Amplitude = √(a² + b²) (log2 units, half peak-to-trough); peak phase = atan2(b, a)·24/2π (h).
  - Phase SE by the delta method.
- **Sensitivity:** categorical `~ factor(ZT)` F-test; agreement of rhythmic calls.
- **Signed time contrasts** (predicted log2 difference from the harmonic fit, with SE via limma contrasts):
  - **T1 = ZT6 − ZT2** (4-h drift inside a typical morning necropsy window);
  - **T2 = ZT12 − ZT0** (light-to-dark transition).

  These are the only circadian quantities used in signed comparisons. Rhythmic amplitude is never used as a signed effect.
- **Male sensitivity (descriptive):** GSE107787 FED livers, the same harmonic model; phase agreement for genes rhythmic in both (circular correlation). This is not independent of the fasting axis and is muscle-contaminated at ZT8–20.

## 8. Stage G: PPARα (GSE282998, TempO-Seq probe counts)
- **WT B6129SF2/J:** DESeq2 `~ dose` (factor 0/5/10/20, reference 0), contrasts 5v0, 10v0, 20v0. Trend: a separate model `~ dose_mgkg` (numeric; log2FC per 10 mg/kg).
- **Muscle sensitivity:** `~ muscle_high + dose`.
- **Technical:** size factors and depth by dose reported.
- **Primary signed PPARα reference = WT 20v0.**
- **Consistency flag:** the sign of 20v0 agrees with 5v0 and 10v0 (point estimates).
- **Genotype:**
  - one model on WT B6129SF2/J + Ppara-null samples, `~ genotype + genotype:dose`;
  - interaction contrast = (WT 20v0) − (KO 20v0), with SE.
- **Classes:**
  - `PPARa_dependent` = WT 20v0 padj < 0.05 **and** interaction padj < 0.05 with the same sign as the WT effect;
  - `agonist_associated_dependence_unresolved` = WT padj < 0.05 without a significant interaction;
  - `not_responsive` otherwise.

  Different breeding stocks (B6129SF2/J vs B6;129S4-Ppara^tm1Gonz/J) are a stated caveat for the interaction. J:ARC(S) is not used in the main reference.
- **Effects are probe-panel, FFPE, targeted.** They are compared with RNA-seq by rank and sign only.

## 9. Stage H: reference comparison
- **Signed references:** MASH1, MASH2, FAST (main), FAST2 (GSE225695), PPARA_20v0 (WT), PPARA_dependent subset, CIRC_T1, CIRC_T2.
- **Pairwise:** Spearman ρ on the pairwise-overlapping main universe; sign agreement among genes significant in both; fry/camera direction per core pathway.
- **Collinearity on genes present in all signed primary references:** correlation matrix, and the variance inflation of each reference regressed on the others (R² → VIF). Pairs with |ρ| ≥ 0.7 are flagged as possibly not separable.

## 10. Stage I: baselines (after the references are documented)
- **Targets:** PE_in_CD, PE_in_MD (and Diet_in_Veh as the positive control).
- **Baseline 1 (DEG overlap):**
  - per reference: Jaccard, one-sided hypergeometric (background = universe genes tested in both);
  - direction-concordant / discordant counts.
- **Baseline 2 (effect correlation):**
  - Spearman ρ over the shared universe;
  - sign agreement among reference-significant genes;
  - target LOO range (existing 19 refits);
  - a reference-label relabelling null for count references.
- **Historical comparator:** split-reference ρ ≈ 0.33 (same mice; Phase 2).
- **Circadian:** T1/T2 only.
- **No method ranking by the largest correlation.**

**Ppara discrepancy (added at the team's request):**
- Reproduce Jung's pairwise approach (DESeq2 on the CD samples only, default settings).
- Compare with our joint-model estimate.
- Isolate the effect of: the sample subset (size factors, dispersion), prefiltering, shrinkage and outlier replacement.

## 11. Deviations
(Filled in during execution.)
