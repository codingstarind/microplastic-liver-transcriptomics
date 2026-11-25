# Phase 3 readiness report

**Date:** 2026-11-22 · **Branch:** `research/phase2-datasets` (not merged) · **Status:** Phase 2 data package complete except for the items under "Not done". Nothing in this phase estimated a signature, fitted a regression or interpreted biology.

## 1. What datasets are ready

"Ready" means:
- downloaded and checksummed;
- the design was verified against GEO and (where stated) the paper;
- sample metadata match the matrix columns (tested);
- genes are mapped to the GSE294180 annotation with every decision recorded;
- QC has been run.

| Axis | Dataset | Planned contrast (biological n) | Values | Ready? |
|---|---|---|---|---|
| MASH 1° | GSE316023 | GAN 23 wk vs chow 23 wk, 21 °C (8 v 8); TN pair (8 v 8) as sensitivity | integer (rounded estimated counts) | **Ready** |
| MASH 2° | GSE246328 | GAN vehicle 24 w (58 wk GAN) vs chow 24 w (7 v 8) | integer counts | **Ready** |
| Fasting 1° | GSE107787 | 24 h fasted vs fed, matched at 6 ZTs (18 v 18) | **RPKM** | **Ready for log-scale linear models only** |
| Fasting 2° | GSE225695 | fasted vs fed (3 v 3), sacrifice time unknown | integer (rounded) | Ready as a direction check only |
| Circadian 1° | GSE117134 | WT female, 6 ZT × 3 | **RPKM** | **Ready for log-scale harmonic models only** |
| Circadian (male sens.) | GSE107787 FED | 6 ZT × 3 | RPKM | Ready; not independent of the fasting control arm |
| PPARα 1° | GSE282998 | B6129SF2/J GW7647 5/10/20 vs vehicle (12 v 4); Ppara-null arm (12 v 4) | integer probe counts | **Ready** (targeted assay) |
| PPARα 2° | GSE8295 | Wy14643 vs control, 129S1 | GCRMA array | **Not prepared** (audited only) |

Files:
- prepared matrices: `D:/mpliver_work/phase2_refs/<GSE>/prepared/`;
- metadata and mapping: `research/roadmap_phase2/datasets/<GSE>/`;
- eligible samples: `metadata/reference_contrast_plan.csv`.

## 2. Data preparation that remains (not done)

1. **AIBD GSE245069 analysis files:** not found locally. The team needs to supply them.
2. **GSE8295:** not prepared. It would need the GPL probe annotation and an array-specific pipeline; optional.
3. **Count-level data for the fasting and circadian references:** not available on GEO. Options:
   - (a) use RPKM with limma-trend (planned);
   - (b) re-quantify the 36 GSE107787 liver runs (about 1.2 G SE100 reads, roughly 60–80 GB FASTQ) with the existing Salmon/M35 pipeline on D:.

   Option (b) is **not started** and needs team approval.
4. **GSE277772 raw counts:** need SRA re-quantification (6 runs). Not started.
5. **A fixed muscle-gene exclusion list and a fixed gene universe:** these are decisions for the start of Phase 3 (see the gene mapping report, §6). They are not data steps.

## 3. Which reference signatures can be estimated in Phase 3

| Signature | Estimable? | Method (planned; see `reference_contrast_plan.csv`) |
|---|---|---|
| MASH (GAN vs chow) | Yes, from 2 independent studies | DESeq2 on counts → per-gene log2FC + SE |
| Fasting vs fed, time-of-day adjusted | Yes | limma-trend on log2(RPKM+1), `~ ZT + state`; muscle genes excluded |
| Time of day | Yes (diurnal, female), plus a male sensitivity check | harmonic regression on log2(RPKM+1) |
| PPARα agonist (WT), with KO specificity | Yes | DESeq2 on probe counts per genotype; 10 vs 0 as a batch-matched sensitivity |

## 4. Valid comparisons and misleading ones

**Valid (with stated limits):**
- Projecting the GSE294180 PE effects (PE_in_CD, PE_in_MD) and the Diet_in_Veh effect onto **independently estimated** reference effects, on a **pre-fixed common gene universe**.
- **Positive control:** our Diet_in_Veh vs the GSE316023 MASH axis (and vs GSE246328).
- MASH-axis agreement between GSE316023 and GSE246328, as a benchmark for how similar two independent references of the same state are.
- Fasting vs PPARα axis collinearity, computed on the shared universe before any PE regression.

**Misleading; do not do:**
- Using GSE294629 together with GSE246328 (same animals), or counting the 4 Hong-lab series as 4 studies.
- Rounding RPKM, FPKM or TPM and running DESeq2/edgeR (GSE107787, GSE117134, GSE245069's FPKM/TPM, GSE277772's normalised values).
- Pooling raw matrices across studies, or batch-correcting across studies with treatment nested in study.
- Calling a GSE225695 fasted-vs-fed signature a pure fasting effect (sacrifice time unknown).
- Using GSE246328's 8 w or 16 w GAN groups against its 24 w chow group (age-confounded), or describing them as "8 weeks of GAN".
- Mixing GSE282998 genotypes, strains or drugs in one contrast without a genotype or strain term. Its absolute fold-changes should not be compared with RNA-seq fold-changes; use directions and ranks.
- Treating similarity to any reference as evidence of a shared mechanism.
- Treating any reference as a fix for the PE ≡ sequencing-depth confounding in GSE294180. **None of these datasets addresses it.**
- Including skeletal-muscle genes in any similarity score. Muscle carry-over is present in GSE294180, GSE107787 (aligned with feeding) and GSE282998 (aligned with dose group).

## 5. Baselines to reproduce next (Phase 3), and the data each needs

| Approach | Needs | Available? |
|---|---|---|
| **Baseline 1:** published-style DEG overlap (PE DEG list ∩ reference DEG list at fixed thresholds, with a hypergeometric test on a fixed universe) | PE DEG lists (Phase 2 results, existing); a reference DEG list per axis; a fixed universe | Yes, after the references are estimated |
| **Baseline 2:** directional signature correlation (Spearman of log2FCs on the fixed universe), including our existing split-reference ρ ≈ 0.33 as the same-data comparator | PE log2FCs (existing); reference log2FCs | Yes |
| **Candidate model:** regression of PE gene effects on MASH, fasting, circadian and PPARα reference effects, SE-weighted, with model selection | PE effects + SEs (existing); reference effects + SEs (Phase 3); collinearity diagnostics; gene-set-blocked cross-validation | Inputs available. The model form is to be decided in Phase 4. |

## 6. Validation criteria to carry into Phases 3 and 4 (pre-specify before seeing PE results)

1. **Positive control:** Diet_in_Veh must align with the external MASH axis (direction agreement and a correlation clearly above a gene-permutation null on the fixed universe). If it fails, the MASH axis is not used.
2. **Technical controls:**
   - binomial thinning of vehicle libraries to PE depth, which must not reproduce the PE projection pattern on its own;
   - a refit excluding muscle-tracking genes and with a muscle covariate;
   - the same muscle handling in GSE107787 and GSE282998.
3. **Reference independence:** no GSE294180 sample in any reference; no reference sample used in two axes that are compared as independent. The GSE107787 FED arm is shared by the fasting axis and the male circadian sensitivity axis, so these two are not independent of each other.
4. **Gene-universe consistency:** one universe, fixed in advance, used for every method. Never "genes significant in both".
5. **Collinearity:** report the fasting–PPARα–circadian–MASH axis correlations and variance-inflation factors. If fasting and PPARα are too collinear, report their joint contribution instead of separate coefficients.
6. **Model evaluation:** identical observations and targets across methods; prespecified metrics; held-out gene sets (blocked by pathway); no selection by the largest apparent correlation.

## 7. Unresolved technical limitations

- GSE294180: PE ≡ lower depth and yield; CD+Veh clock state; muscle carry-over; unknown sacrifice time. None of these can be resolved with external data.
- Chow-vs-purified control diet in both MASH references.
- Female circadian reference; RPKM-only fasting and circadian references.
- Targeted FFPE assay for PPARα; dose-group batch pattern in GSE282998.
- Treatment-aligned muscle signal in GSE107787 FED livers.

## 8. Recommended Phase 3 starting point

1. Fix the gene universe and the muscle exclusion list (a short written decision).
2. Estimate the two MASH axes and run the **positive control** first (Diet_in_Veh vs GSE316023). This is the go/no-go for Track A.
3. Estimate the fasting, circadian and PPARα axes; report their collinearity.
4. Reproduce Baselines 1 and 2 on the fixed universe.

**Are we ready to start baseline modelling?** Yes, for the MASH and PPARα axes (count data). Yes with a stated method restriction (log-scale linear models) for the fasting and circadian axes. The team should first confirm the deadline and decide the items in §2.3 and in the gene mapping report §6.
