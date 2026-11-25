# Dataset selection report (Phase 2)

**Date:** 2026-11-22 · **Status:** recommendation for team review. "Verified" means checked directly against the GEO SOFT records, the deposited files and (where stated) the paper's methods. Fields the sources do not state are marked **NOT STATED**, not filled in.

## Summary of selections

| Axis | Selected | Role | Contrast (biological replicates) | Data format |
|---|---|---|---|---|
| MASH | **GSE316023** | Primary | GAN D09100310 23 wk vs LabDiet 5001 chow 23 wk, 21 °C (8 vs 8) | Non-integer estimated counts (+ FPKM) |
| MASH | **GSE246328** | Secondary | GAN vehicle 24 w (**58 wk GAN in total**) vs chow vehicle 24 w (7 vs 8) | Integer STAR counts |
| Fasting | **GSE107787** | Primary | 24 h fasted vs fed, matched at 6 Zeitgeber times (18 vs 18) | **RPKM only** |
| Fasting | GSE225695 | Secondary (check only) | 24 h fasted vs fed (3 vs 3); sacrifice time NOT STATED | StringTie-derived counts |
| Circadian | **GSE117134** | Primary | WT female liver, ZT0–ZT20 every 4 h (6 × 3) | **RPKM only** |
| Circadian | GSE107787 FED arm | Male sensitivity | FED liver, 6 ZT × 3 | RPKM only; shares mice with the fasting control arm |
| PPARα | **GSE282998** | Primary | B6129SF2/J WT: GW7647 5/10/20 mg/kg/day vs vehicle (4 per dose vs 4), plus Ppara-null specificity arm | TempO-Seq targeted probe counts |
| PPARα | GSE8295 | Audited, **not prepared** | Wy14643 vs control in 129S1 WT (4 vs 4) | Affymetrix GCRMA |
| — | GSE294629 | **Rejected** | Contains 54 of the 55 GSE246328 libraries | Salmon counts |

---

## 1. MASH references

### 1.1 GSE316023 (primary): verified design

Source: GEO SOFT (77 GSMs) and Nat Commun 2026;17:6964 (PMC13392390, Methods read in browser).

| Item | Verified value |
|---|---|
| Diet | GAN **D09100310** (Research Diets), the same catalogue product as the GSE294180 MASH arm. Composition is not given in this paper; the Gubra papers (GSE246328/GSE294629) give 40 kcal-% fat (46% saturated FA by weight), 22% fructose, 10% sucrose, 2% cholesterol, 4.49 kcal/g. |
| Control diet | **LabDiet 5001 normal chow** (grain-based). **Not** the purified control D09100304 used in GSE294180. |
| Strain / sex / age | C57BL/6J (Jackson), male, 7–8 wk at diet start ("Age: 7 weeks upon start of GAN diet", Fig. 2 legend) |
| Duration | GAN sampled at 12, 16, 23 and 32 wk; **chow controls only at wk 23** |
| Housing | TA 21 ± 1 °C ("A" in sample names) or TN 29 ± 1 °C ("B"). The Table 1 footnote gives the ranges as 20–22 and 28–30 °C. Lights on 06:00. |
| Sacrifice | Cardiac puncture under isoflurane **after a 5 h fast**. Clock time NOT STATED. |
| Chow vs GAN at wk 23 sampled at comparable times? | Both groups are labelled wk23 with the same protocol. Same-day or same-session sampling is **NOT STATED**. |
| Pooling | No biological pooling. "Pooled into a single lane" refers to library multiplexing only. 77 GSMs map to 77 individual libraries (two GSMs have 2 SRA runs each, which are technical sequencing units). |
| Assay | Bulk poly-A RNA-seq (TruSeq stranded mRNA), NovaSeq 6000 PE, mm10. Counting software NOT STATED. |
| Files | `COUNT_…csv.gz` (non-integer estimated counts) and `FPKM_…csv.gz`; raw reads in SRA (PRJNA1400485) |

**Complete sample membership of the selected contrast** (from GEO, not from the Phase 1 summary):
- **GAN wk23 TA (case):** SD_wk23_A_R1–R8 = GSM9443427–GSM9443434
- **Chow wk23 TA (control):** NC_wk23_A_R1–R8 = GSM9443443–GSM9443450
- **TN sensitivity pair:** GSM9443435–442 (GAN) vs GSM9443451–458 (chow), 8 vs 8

**Why it is selected:** independent mice; the same GAN product; male C57BL/6J; n = 8 per group; housing temperature recorded; matched 5 h fast; deep libraries (median 50 M read pairs); no technical depth difference between GAN and chow at wk23 TA (ratio 0.93, Mann–Whitney p = 0.13, no separation).

**Differences from GSE294180 that must stay visible:**

| | GSE294180 (ours) | GSE316023 TA wk23 |
|---|---|---|
| MASH diet | D09100310 | D09100310 (same) |
| Control diet | **D09100304 purified** | **LabDiet 5001 grain chow** |
| Diet duration | 8 wk | **23 wk** |
| Age at sampling | NOT STATED in our audit (diet started in young adults) | ≈ 30–31 wk |
| Strain / sex | C57BL/6, male | C57BL/6J (JAX), male |
| Housing | NOT STATED | 21 °C |
| Pre-sacrifice fast | NOT STATED | 5 h |
| Platform | NextSeq 2000, Salmon/GENCODE M35, integer counts | NovaSeq 6000, mm10 (tool NOT STATED), estimated counts |
| Gene coverage | 13,777 filtered genes | 12,586 of these present; 12,259 detected in the contrast samples |

The chow-vs-purified difference is the largest threat to this reference. Grain chow differs from purified diets in fibre, phytoestrogens and fat source, so a GAN-vs-chow signature contains a "purified vs grain diet" component that our CD-vs-MD contrast does not. The positive control in Phase 3 (our Diet_in_Veh effect against this axis) tests whether the axis is still usable. It cannot remove this difference.

### 1.2 GSE246328 (secondary): verified design

Source: GEO SOFT (55 GSMs), BMC Gastroenterol 2024;24:59 Methods (read in browser), GSE294629 SOFT.

**Correction to the Phase 1 screen:** the "8-week GAN group" is **not** 8 weeks of GAN diet. Per the paper: "Animals were fed the GAN diet (vehicle and chow reversal groups), or chow (chow vehicle) for 34 weeks". This was followed by liver biopsy, 4 weeks of recovery, and then the 8, 16 or 24 week intervention period ("fed chow or the GAN diet for up to 58 weeks").

Totals below assume the 34-week pre-feeding ends at intervention start. The paper's phrasing allows that the 4-week post-biopsy recovery may add to it, which would make each total up to 4 weeks longer.

| GEO group | n | Total diet history |
|---|---|---|
| Vehicle 8 w | 9 | ≈ 42 wk GAN |
| Vehicle 16 w | 8 | ≈ 50 wk GAN |
| Vehicle 24 w | 7 | ≈ 58 wk GAN |
| Chow-reversal 8/16/24 w | 7/8/8 | 34 wk GAN, then 8/16/24 wk chow |
| Chow vehicle 24 w | 8 | ≈ 58 wk chow (Altromin 1324) |

- Strain C57BL/6J from Janvier (C57BL/6JRj per GSE294629), 5–6 wk at start, 21 °C, lights on 03:00. All arms were dosed s.c. with vehicle once daily. GAN mice were biopsied and single-housed. Whether chow mice were biopsied is NOT STATED. Terminal blood was taken after a 4 h fast.
- **Sex:** not stated in GSE246328 or its paper. GSE294629 (the same libraries) states male.
- **Only the 24 w pair is a valid contrast.** The 8 w and 16 w GAN groups have no age-matched chow control.
- **Overlap with GSE294629 confirmed:** 54 of 55 GSE246328 library IDs (e.g. `16174`) appear as GSE294629 FASTQ names. GSE294629 adds semaglutide arms. These are the **same animals**, so only GSE246328 is used. Neither may serve as validation of the other.

**Why secondary:** an independent lab, the same diet product, integer counts and no depth difference between arms (ratio 1.02). But it represents **advanced MASH after 58 weeks** with biopsy and single housing, which is far from our 8-week exposure. Use it as a second MASH axis and as the "how well do two MASH references agree" benchmark, not as the primary axis.

---

## 2. Fasting and circadian references

The axes are kept separate. A fasting axis is valid only if fasted and fed mice were sampled at the same clock times.

### 2.1 GSE107787 (primary fasting; male circadian sensitivity)

Source: GEO SOFT (72 GSMs) and Kinouchi 2018 Cell Rep STAR Methods (PMC6433478).

- 8-wk male C57BL/6J (JAX), normal chow, 12:12 light–dark, single-housed for 1 wk before harvest.
- **Fasting protocol:** "The fasted mice were placed in new cages without food 24-hr before sacrifice". Every fasted mouse was therefore fasted exactly 24 h, and fasting onset shifts with harvest time.
- Liver and gastrocnemius harvested at **ZT0, 4, 8, 12, 16 and 20 in both FAST and FED** (n = 3 per ZT per state). **Feeding state and time of day are fully crossed and balanced,** so the fasting effect can be estimated with time of day as a separate term (and vice versa).
- Liver samples: 36 (FAST 18, FED 18). Muscle samples are not used.
- Assay: TruSeq stranded mRNA, HiSeq 2500 SE100; 72 samples multiplexed in 3 pools of 24; ELAND/RefSeq in-house pipeline. **Deposited data are RPKM only** (no counts on GEO).

**Limitations:**
- RPKM only, so no negative-binomial count model. Phase 3 must use log-scale linear models (limma-trend) or re-quantify from SRA (36 liver runs, about 1.2 G reads in total).
- **Muscle-gene signal in FED livers at ZT8–ZT20** (median summed muscle RPKM about 50 vs about 0.8 in fasted livers; see QC report). This is a treatment-aligned technical artefact, most likely cross-sample bleed in mixed liver and muscle pools. Muscle-marker genes must be excluded, with a muscle-score sensitivity analysis.
- A 24 h fast is an extreme state. Similarity to it says "fasting-like direction", not "comparable magnitude".

**Recommendation: primary fasting reference.** It is the only verified dataset in which a fasting effect is separable from time of day.

### 2.2 GSE225695 (subseries of GSE225697): secondary fasting check

Male C57BL/6J (SLAC), 10 wk, fed / 24 h fasted / 2 h refed (n = 3). It deposits FPKM plus StringTie-derived counts. **The clock time of sacrifice is NOT STATED** for any group, so fasting may be confounded with time of day. With n = 3 it serves only as a count-based direction check of GSE107787, not as a reference axis.

### 2.3 GSE117134 (primary circadian)

Source: GEO SOFT. WT, Bmal1-KO and liver-reconstituted mice, **female**, 8–12 wk, ad libitum chow, livers at ZT0–ZT20 every 4 h (n = 3 each). **Only the 18 WT livers are used.** RPKM only, from the same core and pipeline as GSE107787. These are independent mice from the fasting reference, which is why the time-of-day axis is taken from here.

**Limitations:**
- Female (our data are male). The male GSE107787 FED arm is the sensitivity check, but it shares samples with the fasting control arm.
- Light–dark (diurnal) sampling, so clock and feeding rhythms are combined. A "circadian" axis here means "time of day under ad libitum feeding".
- RPKM only.

**What the circadian axis can and cannot do:** GSE294180's sacrifice time is **unknown** (Phase 1). A time-of-day axis can show whether PE-associated changes resemble differences between times of day. It **cannot** establish that our groups were sampled at different times.

---

## 3. PPARα reference

### 3.1 GSE282998 (primary): verified design

Source: GEO SOFT (192 GSMs) and Heintz et al. 2025 Toxicol Sci (PMC12198672, Methods read in browser).

| Item | Verified value |
|---|---|
| GW7647 arm | 48 samples: 3 strains × 4 doses (0, 5, 10, 20 mg/kg/day) × 4 male mice |
| Strains | **B6129SF2/J** (WT, the matched background for the knockout), **B6;129S4-Ppara^tm1Gonz/J** (Ppara-null), **J:ARC(S)** (Swiss outbred, CD-1 equivalent) |
| Route / vehicle / duration | Oral gavage in **1% carboxymethylcellulose**; 5 days; euthanised **24 h after the last dose** |
| Age / housing | 9–12 wk; 20–26 °C; ~12 h light cycle; fasting and clock time NOT STATED (APAP arm only was fasted) |
| Controls | Separate 0 mg/kg vehicle group per strain within the GW7647 arm (n = 4); vehicle groups are not shared with other chemicals' arms |
| Assay | **TempO-Seq** whole-transcriptome mouse panel on **FFPE** liver: 30,146 probes for 21,398 genes; integer probe counts; NovaSeq X SE |
| GEO metadata error | For the knockout-strain samples, the GEO "genotype" field says "WT" in the GW7647 arm. The "strain" field (`PPARa KO`) and the paper are used; recorded in `datasets/GSE282998/sample_metadata.csv`. |

**Selected contrast:** B6129SF2/J GW7647 (5/10/20, 12 mice) vs vehicle (4 mice), within genotype and within arm.
- **Specificity arm:** the same comparison in Ppara-null mice. Genes that respond in WT but not in the null mice define the PPARα-dependent subset.
- **J:ARC(S):** strain-robustness sensitivity only.
- **No mixing:** genotypes, strains and drugs are not mixed. Rosiglitazone, HFPO-DA and APAP arms were not downloaded.

**Technical findings:**
- **Muscle-gene signal** in 2 of 4 mice in every 0 and 10 mg/kg group (summed CPM about 200–270), and in none of the 5 or 20 mg/kg mice, in all three strains. This pattern suggests a processing batch shared by the 0 and 10 mg/kg groups. The 10-vs-0 comparison is batch-matched for this feature; 5 and 20 vs 0 are not.
- Sequencing depth also varies between dose groups (medians about 15 M vs 25–29 M).

**Gene coverage:** 11,966 of our 13,777 filtered genes are present; 10,060 are detected in the primary-contrast samples.

**Recommendation: suitable as the primary PPARα reference,** with three caveats:
1. It is a targeted FFPE assay. Compare directions and ranks, not absolute fold-change scales, against bulk poly-A RNA-seq.
2. It is a strong pharmacological agonist dose. Similarity means "PPARα-target-like direction", not equivalent PPARα activity, and it says nothing about PPARα protein activity in our mice.
3. The dose-group batch pattern must be handled.

### 3.2 GSE8295 (secondary candidate, not prepared)

Male 129S1/SvImJ WT vs 129S4/SvJae Ppara-null, Wy14643 0.1% in food for 5 days (n = 4), Affymetrix, GCRMA.

It is not prepared:
- it uses a different technology (log intensities);
- it uses a different agonist;
- WT and null are different 129 substrains;
- probe-to-gene annotation would need the platform file.

Only the normalised series matrix was downloaded (checksummed). It could be prepared later as an independent agonist check if the team wants one.

---

## 4. Not found locally

- **The AIBD report and its GSE245069 analysis** (16 PS-exposed livers; Welch and permutation tests, sensitivity analyses) were **not found** on this machine. The repository, Desktop, Downloads and course folders were searched by name and content for `GSE245069`, `AIBD`, `Welch` and `polystyrene`. The only local GSE245069 mentions are in Phase 1 notes and literature files. The one PDF in Downloads (`BT3041_Microplastics_Codebase_Research_Report.pdf`) mentions GSE245069 once, as a candidate. **This does not mean the analysis was never performed;** its files are just not on this machine.
- **Course submission instructions:** none found locally (course folder `Desktop/BT3041` holds tutorials only; lecture PDFs in Downloads contain no deadline). The only deadline seen is the final-report draft header ("Final report due 2026-11-25"). **Please confirm the actual deadline.**
