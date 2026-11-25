# Literature review: microplastics, liver transcriptomics and metabolic dysfunction

**Prepared for:** BT3041 project team (IIT Madras), Phase 2 planning
**Date:** 20 November 2026
**Scope:** 19 papers (16 read as full text, 3 as abstract only), published 2022–2026.

> **How to read this document.** It is organised by *scientific question*, not by paper. Every
> factual claim carries a traceable source: a `[shortid]` pointing to a per-paper record in
> [`papers/`](papers/), which holds the verbatim quotes and their locations (figure, table, section).
> The comparison matrix is [`paper_matrix.csv`](paper_matrix.csv); repository checks are in
> [`dataset_audit.md`](dataset_audit.md); open questions are in [`research_gaps.md`](research_gaps.md).
>
> **Labels used throughout.** *VALIDATED* = supported by an orthogonal non-omics assay (qPCR,
> western, histology, enzyme activity, knockout, rescue). *OMICS-ONLY* = from sequencing/enrichment
> alone. *COMPUTATIONAL-PREDICTION* = from docking, network topology or TF-activity inference.
> *NOT VERIFIED* = we could not read the source for that field.
>
> **Nothing here is a hypothesis for the team.** Where the literature is ambiguous, this document
> says so rather than resolving it.

---

## 1. What was searched, and what was included

Europe PMC was searched on 20 November 2026 across title/abstract for combinations of
(microplastic | nanoplastic | polystyrene | polyethylene | PET | PVC) × (liver | hepatic |
hepatocyte) × (transcriptome | RNA-seq) × (mouse | rat), plus targeted searches for
microplastic × high-fat/MASLD/MASH, for null or contrasting results, for human liver detection,
and for single-cell/spatial studies. Nineteen papers were selected to cover:

| Category | Papers |
|---|---|
| The project's own source paper | `jung2026` |
| Other **polyethylene** oral exposure, liver | `ees2023_pe_liver_fibrosis`, `stoten2024_pe_tlr2` |
| Plastic × diet (or × genotype) factorial designs | `ehp2023_ps_hfd`, `ehp2024_ps_hfd_obesity`, `ijms2024_ps_lep`, `liverint2025_np_masld`, `ijbs2025_psnp_lysosome`, `advsci2026_psnp_circadian`, `jnanobiotech2026_ncor1_ppara`, `liverint2026_lowdose_alleviate`, `stoten2024_scrna_nafld` |
| Chow-only liver transcriptomes (comparators) | `sci_rep_2025_ps_hubgenes`, `biology_2026_ps_rorg`, `toxsci2024_obesogenic_ps`, `advsci2026_pet_gutliver` |
| Mechanistic / multi-omic | `natmetab2026_size_macrophage` |
| Purely computational | `plosone2026_pet_comp` |
| Human tissue | `horvatits2022_human_cirrhosis` |

Three papers (`stoten2024_pe_tlr2`, `stoten2024_scrna_nafld`, `liverint2026_lowdose_alleviate`)
are behind subscriptions; only their abstracts could be read. They are used for direction and
context only, never for numbers.

---

## 2. Question 1 — Does microplastic exposure change the liver transcriptome at all, and how large is the effect?

**Established.** Oral exposure to plastic particles produces measurable hepatic transcriptional
change in every in-vivo study that looked, across polymers and sizes. But the effects are
consistently **modest next to diet**, and the reported DEG counts are not comparable between
papers because the thresholds differ by orders of magnitude.

| Study | Polymer / size | Model | DEG count (as reported) | Threshold used |
|---|---|---|---|---|
| `jung2026` | PE 10–150 µm | PE × MASH diet | 503 CD-only, 182 MD-only, 97 shared | \|log2FC\|≥0.38, FDR<0.05 |
| `sci_rep_2025_ps_hubgenes` | PS 80 nm / 0.5 / 5 µm | chow | 333 / 602 / 758 | \|log2FC\|>1, adj. p<0.05 |
| `biology_2026_ps_rorg` | PS 80 nm | chow | 689 (389 up, 300 down) | \|log2FC\|>1, **raw p<0.05** |
| `toxsci2024_obesogenic_ps` | PS 0.5 / 5 µm | chow | 283 / 49 | q<0.05, no fold-change filter |
| `ijbs2025_psnp_lysosome` | PS 50 nm | PS × HFD | 361 / 628 / 1501 | \|log2FC\|>0.5, **raw p<0.05** |
| `jnanobiotech2026_ncor1_ppara` | PS 50 nm | PS × HFD | 236 (95 up, 141 down) | \|log2FC\|>1, adj. p<0.05 |
| `advsci2026_pet_gutliver` | PET 10–200 µm | chow | 28 at 8 wk; 83 at 29 wk | **not reported at all** |

**Three cautions the team should carry forward.**

1. **Two of these DEG lists are uncorrected.** `biology_2026_ps_rorg` and `ijbs2025_psnp_lysosome`
   both called genes at raw p < 0.05 with no FDR control, and neither names its DE tool. Their
   counts are not comparable to FDR-controlled counts, including the team's own.
2. **One paper reports DEG counts with no stated method.** `advsci2026_pet_gutliver` gives counts
   but no DE tool, fold-change cut-off, FDR threshold or read depth, so those numbers cannot be
   reproduced from the paper.
3. **The diet effect dwarfs the particle effect wherever both are measurable.** `ijbs2025_psnp_lysosome`
   reports 1501 DEGs between HFD+PS and CD+PS, versus 361 and 628 for the particle contrasts —
   and conspicuously **omits the HFD-vs-CD contrast entirely**, so its readers cannot judge the
   relative size of the two effects. The team's own factorial design does not have this blind spot.

**Reporting quality is poor across the board.** Sequencing depth is reported in only one of the
nine transcriptome studies (`advsci2026_psnp_circadian`: median 43.78 M clean reads, alignment
rates given, single batch, no correction applied). Per-group RNA-seq *n* is missing from
`jung2026`, `biology_2026_ps_rorg` and `ijbs2025_psnp_lysosome`.

---

## 3. Question 2 — Do polyethylene findings replicate, and does Jung et al. agree with independent PE work?

This matters most for RQ1, because the project's dataset is PE-specific and PE is
under-studied relative to polystyrene.

**Only two independent PE oral-exposure mouse liver studies were found.** One
(`stoten2024_pe_tlr2`) is abstract-only. The other, `ees2023_pe_liver_fibrosis` (Djouina et al.
2023), is the only substantive independent comparison available — and it used **qPCR panels, not
a transcriptome**.

### Where Jung and Djouina agree

- **`Ppara` is induced by PE.** Djouina reports "a two-fold increase in expression of the
  peroxisome proliferator-activated receptor alpha (Ppara) gene ... after exposure to both types
  of PE microbeads" (6 weeks), and again at 9 weeks for the larger beads. Jung reports Ppara
  induced in PE groups. Two independent labs, two routes, two sexes, same direction. This is the
  single most reproducible PE-specific transcriptional observation in the literature.
- **Liver injury markers rise.** Both report increased serum transaminases with PE.

### Where they disagree — and it is a real disagreement in direction

**De novo lipogenesis moves in opposite directions.**

| | `jung2026` (PE 10–150 µm, 2 mg/day gavage, male, 8 wk, CD and MASH diet) | `ees2023_pe_liver_fibrosis` (PE 36 and 116 µm, 100 µg/g food, **female**, 6–9 wk, chow) |
|---|---|---|
| *Fasn* | **DOWN** (shared PE DEG) | **UP** in all PE groups |
| *Acly* | **DOWN** (shared PE DEG) | **UP** in all PE groups |
| *Acaca* | not reported | **UP** in all PE groups, both time points |
| Liver TG | UP under MASH diet | UP (combined beads) |
| *Ppara* | UP | UP (~2-fold) |

Both studies report increased hepatic lipid content while disagreeing about the lipogenic
transcripts underlying it. **This should not be called a contradiction without stating the
differences**, which are substantial: sex (male vs female — and Djouina themselves note *Cyp2a4*
is ~1000-fold female-biased), route (daily gavage of a fixed 2 mg vs dietary admixture at
100 µg/g food), particle size distribution, and the presence of a MASH diet arm in one and not
the other. Djouina used no multiple-testing correction across a large qPCR panel, which also
weakens individual gene calls.

**A third relevant observation:** `natmetab2026_size_macrophage` reports, from single-nucleus
RNA-seq, that **de novo lipogenesis is suppressed while fatty-acid uptake rises** — the same
direction as Jung, in polystyrene, with no diet challenge. So Jung's lipogenic down-regulation is
not isolated; it has independent support from a different polymer and platform.

### A mechanistic disagreement worth the team's attention

Djouina's authors conclude that PE particles of 36–116 µm probably **never reached the liver**:
*"it seems unlikely that PE MPs reached the liver in our experimental conditions. Therefore, the
effects are more likely to be mediated indirectly, in particular through the effects of MPs on
dysbiosis and intestinal immunity."* Jung, by contrast, reports PE-like particles detected in
MD+PE liver tissue by micro-Raman and O-PTIR.

These are not necessarily incompatible (different size ranges, different detection methods), but
they imply **different mechanisms** — direct hepatic particle action versus a gut-liver axis
effect — and therefore different expectations for what a liver transcriptome should show. Three
other papers support the indirect route: `ehp2024_ps_hfd_obesity` abolished the liver phenotype
with antibiotics and transferred it with faecal transplant; `advsci2026_pet_gutliver` found a
full NAFLD-like liver phenotype while **pyrolysis-GC/MS detected no PET in liver**; and
`natmetab2026_size_macrophage` found no particles inside hepatocytes at all, with Kupffer cells
as the reservoir.

---

## 4. Question 3 — Does a metabolically stressed liver respond differently to plastic? (RQ2)

**This is the central methodological finding of the review, and it is unambiguous.**

> **Not one study in this set formally tested a plastic × diet interaction.**

Nine papers use a 2×2 (or 2×2-like) design. Every one of them analyses it with pairwise
comparisons, and most with one-way ANOVA across the four groups:

| Study | Design | Analysis actually used | Interaction term? |
|---|---|---|---|
| `jung2026` | diet × PE | one-way ANOVA + Tukey; pairwise DESeq2 | **No** |
| `ehp2023_ps_hfd` | diet × PS | one-way ANOVA + Holm–Šidák | **No** |
| `ehp2024_ps_hfd_obesity` | diet × PS (+ABX, +FMT) | one-way ANOVA + Tukey | **No** (the only two-way is group × *time*) |
| `liverint2025_np_masld` | diet × NP (unbalanced 6/12/6/12) | one-way ANOVA + Tukey, all contrasts vs Chow+Veh | **No** |
| `ijbs2025_psnp_lysosome` | diet × PS | one-way ANOVA; three separate pairwise DEG contrasts | **No** (and no HFD-vs-CD contrast) |
| `ijms2024_ps_lep` | genotype × PS | t-tests / one-way ANOVA | **No** |
| `jnanobiotech2026_ncor1_ppara` | diet × PS | two-way ANOVA **fitted** | Model exists; **no interaction F or p ever printed**; omics cover only 2 of 4 groups |
| `advsci2026_psnp_circadian` | diet × PS (× age) | two-way ANOVA **fitted** for phenotypes | Model exists; **no interaction F or p printed**; transcriptomics pairwise only |
| `natmetab2026_size_macrophage` | particles × HFD arm | ANOVA/Kruskal–Wallis + FDR | **No** |

The two papers that come closest (`jnanobiotech2026_ncor1_ppara`,
`advsci2026_psnp_circadian`) do fit two-way models for phenotypes, but neither prints the
interaction statistic, and both analyse their transcriptomes as pairwise contrasts only.

**Why this matters.** Nearly every "exacerbation" claim in this literature rests on the pattern
*significant in the high-fat arm, not significant in the control arm* — the difference-of-significance
fallacy. The clearest illustration is `ehp2023_ps_hfd`, which concludes that microplastics caused
metabolic disturbance "only in mice fed a high-fat diet", yet its own data show microplastic
effects **under the normal diet too**: raised iPGTT glucose and reduced goblet cells/Muc2 mucus
thickness in both diets, plus non-significant trends in serum lipids and Oil Red O that the
authors themselves attribute to sample size. `ijbs2025_psnp_lysosome` states in its abstract that
there were "no effects ... on a normal chow diet" while its own results report increased
macrophages, neutrophils, collagen deposition and 361 DEGs in the chow+particle group.

**So: the team's pre-specified factorial interaction model addresses something the published
literature has not done.** That is a defensible, verifiable novelty claim — but it is a claim
about *statistical approach*, not about discovering a new biological effect, and it should be
worded that way.

### What the diet-interaction literature does agree on

Despite never testing it formally, the studies converge on one qualitative pattern:

- **Plastic alone, on a control diet, does little to the liver.** `liverint2025_np_masld`:
  chow+NP mice "showed normal liver morphology". `jnanobiotech2026_ncor1_ppara`: PSNPs
  "without phenotypically affecting ND-fed mice". `ehp2024_ps_hfd_obesity`: no body-weight effect
  under normal diet. `liverint2026_lowdose_alleviate` (abstract): "no obvious pathological changes
  were observed in NP group". Even the contrarian paper agrees on this arm.
- **On a high-fat background, a liver phenotype appears** in most studies — but see §6 for the
  exception that reverses direction.

### And one study where the effect reverses sign by metabolic background

`ijms2024_ps_lep` crossed polystyrene with *genetic* obesity (leptin knockout) rather than diet.
Hepatic steatosis and *PPARg*/*C-EBPa* went **up in wild-type** mice and **down in leptin-knockout**
mice; serum lipids rose only in wild-type. This is a direct demonstration that metabolic
background can flip the sign of a microplastic effect — but, again, the genotype dependence was
never tested as an interaction, only as separate within-genotype comparisons.

---

## 5. Question 4 — PPARα: induced or suppressed?

This is the sharpest biological disagreement in the set, and it bears directly on Jung's
mechanistic claim.

| Study | Polymer / size / route | Direction |
|---|---|---|
| `jung2026` | PE 10–150 µm, gavage, 8 wk | *Ppara* **UP** (vs CD+Veh); hub gene of the MD+PE WGCNA module |
| `ees2023_pe_liver_fibrosis` | PE 36/116 µm, diet, 6–9 wk | *Ppara* **UP** ~2-fold |
| `sci_rep_2025_ps_hubgenes` | PS 80 nm–5 µm, water, 12 wk | *Ppara* **UP**; *Acox1*, *Cpt1a*, *Mcad* up |
| `jnanobiotech2026_ncor1_ppara` | PS 50 nm, gavage 150 mg/kg/day, 90 d | PPARα signalling **SUPPRESSED**: *Cpt1a*, *Acox1*, *Acsl3*, *Cyp4a31* down at mRNA; ACSL3/CPT1A down at protein |
| `advsci2026_psnp_circadian` | PS 100 nm, water, 10 wk | PPARα not measured; FAO capacity (CPT1 activity) **DOWN** under combined exposure |
| `ijbs2025_psnp_lysosome` | PS 50 nm, 8–20 wk | "PPAR signalling" enriched, but it is **PPARγ** that is measured and up; PPARα never reported |

`jnanobiotech2026_ncor1_ppara` is explicit that it contradicts the prevailing direction:
*"the observed suppression (rather than activation) of this pathway"*.

**How the team should treat this.** It is a genuine direction conflict, but the conditions are not
comparable, and three differences are documented and material:

1. **Polymer and size.** 50 nm polystyrene nanoparticles versus 10–150 µm polyethylene — different
   uptake routes entirely (the nanoparticle study shows hepatic delivery; the PE studies may not
   reach the liver at all).
2. **Dose and route.** 150 mg/kg/day by gavage for 90 days versus a fixed 2 mg/day for 8 weeks.
3. **An overnight fast before sacrifice** in `jnanobiotech2026_ncor1_ppara`. Fasting is itself a
   strong physiological activator of hepatic PPARα and its FAO targets, so **all four groups in
   that study sit in a PPARα-primed state** — which changes what a measured suppression means.

Note also that **neither** paper measured what it most needed to: `jnanobiotech2026_ncor1_ppara`
never measured *Ppara* itself (only targets and the corepressor), and `jung2026` reports *Ppara*
induction relative to CD+Veh, a contrast that mixes diet and PE.

---

## 6. Question 5 — Contrasting and negative results

The team asked specifically for these. They are a real and under-cited part of this literature.

**A direction reversal.** `liverint2026_lowdose_alleviate` (Liver International 2026, Wistar rats,
90 days) reports that **low-dose** polystyrene nanoplastics *reduced* high-fat-diet-induced lipid
accumulation, inflammation and fibrosis, with fatty-acid oxidation up, AMPK up, mTOR down and
lipophagy proposed as the mechanism. That is the exact inverse of `ijbs2025_psnp_lysosome`, which
reports mTOR **up** and autophagic flux **blocked** under a similar exposure. *Caveat we cannot
resolve:* the paper is abstract-only, and **the numeric dose is never stated in the abstract** —
so the dose at which the direction flips cannot be established. It is also the only rat study here.

**A size-dependent reversal within one study.** `natmetab2026_size_macrophage` found that under a
high-fat diet, 500 nm microplastic **increased** hepatic lipid while 50 nm nanoplastic
**decreased** it, in the same experiment. The same paper found nanoplastic *reduced* weight gain
and fat mass, and *increased* oxygen consumption.

**Well-documented nulls.** These deserve as much weight as the positives:

- `toxsci2024_obesogenic_ps` (chow, 13 weeks): liver histology unchanged, hepatic triglycerides
  unchanged, ALT/AST unchanged, hepatic GSH/GSSG unchanged, ileum histology normal, plasma
  adipokines unchanged — alongside 283 DEGs and raised hepatic cholesterol.
- `sci_rep_2025_ps_hubgenes`: oxidative stress explicitly **not** changed — "the levels of GPX,
  MDA, and SOD were not increased by exposure to PS-NMPs of any of the three diameters" — which
  the authors note is "contrary to previous studies".
- `jnanobiotech2026_ncor1_ppara`: no liver injury, no inflammation, no fibrosis — a clean
  steatosis-only model.
- `natmetab2026_size_macrophage`: no hepatocellular injury, no fibrosis, and explicitly **no
  inflammation** (no myeloid expansion, no serum cytokine rise).
- `ehp2023_ps_hfd`: under normal diet, the entire liver qPCR panel was null (0.86–1.22× control).

**Two papers that disagree with each other head-on.** `ehp2023_ps_hfd` (0.45–0.53 µm PS-COOH,
4 weeks) reports **no body-weight change** and **lower** absolute liver weight in the
high-fat+plastic group. `ehp2024_ps_hfd_obesity` (5 µm PS, 8 weeks) reports **>17% greater**
body-weight gain and **higher** liver weight. The 2024 paper cites the 2023 one as the
contrasting result and attributes the difference to particle size and properties.

---

## 7. Question 6 — Inflammation, oxidative stress, mitochondria and fibrosis

The pathway-level picture is far less consistent than review articles suggest.

**Inflammation.** Reported up in `ehp2023_ps_hfd` (hepatic *Tnfa*, *Il6*, *Il1b* 5.5–6.2×, high-fat
arm only), `liverint2025_np_masld`, `advsci2026_pet_gutliver`, `sci_rep_2025_ps_hubgenes` and
`ijbs2025_psnp_lysosome` (the last on **both** diets). But explicitly **absent** in
`natmetab2026_size_macrophage` and `jnanobiotech2026_ncor1_ppara`, and in `jung2026` the lobular
inflammation score was only "an increasing trend". `ees2023_pe_liver_fibrosis` found inflammation
up in combined-bead groups while *Il1b* went **down** in single-bead groups.

**Oxidative stress.** Up in `biology_2026_ps_rorg` (ROS +46%, MDA +35%, SOD/CAT/GSH down) and
`liverint2025_np_masld`. **Explicitly unchanged** in `sci_rep_2025_ps_hubgenes` and
`toxsci2024_obesogenic_ps` (the latter despite NRF2/glutathione pathway *enrichment*, which the
authors attribute to glutathione S-transferases). **Not measured at all** in `jung2026`,
`ehp2023_ps_hfd`, `ehp2024_ps_hfd_obesity`, `jnanobiotech2026_ncor1_ppara` or
`advsci2026_pet_gutliver`. `ees2023_pe_liver_fibrosis` found MDA *down* at 6 weeks and *up* at
9 weeks — direction depends on duration.

**Mitochondria.** The best-supported mitochondrial evidence is `advsci2026_psnp_circadian`
(reduced cristae density, membrane potential, ATP/ADP, CPT1 activity, and respiration by
Seahorse) and `biology_2026_ps_rorg` (ATP and Complexes I/III/V down). Most other papers infer
mitochondrial dysfunction from transcripts alone; `advsci2026_pet_gutliver` claims it in the
abstract with no locatable supporting panel.

**Fibrosis.** Up in `sci_rep_2025_ps_hubgenes`, `liverint2025_np_masld`, `advsci2026_pet_gutliver`
(29 weeks only) and `ees2023_pe_liver_fibrosis` (clearly only on a CCl₄ background). **Absent** in
`jung2026` ("despite the absence of overt fibrosis"), `natmetab2026_size_macrophage` and
`jnanobiotech2026_ncor1_ppara`. `ijbs2025_psnp_lysosome` reports collagen deposition up by
Picro-Sirius Red while its own transcriptome shows fibrosis genes unchanged — and does not
reconcile the two.

---

## 8. Question 7 — Technical factors: what the literature does and does not control

The project's Phase 2 found several PE-associated results sensitive to technical factors. It is
worth knowing how the published field handles the same issues. **Mostly, it does not.**

| Control | How many of the 16 full-text studies |
|---|---|
| Sampling time / circadian phase controlled | **1** (`advsci2026_psnp_circadian`: all collections in a ZT4–ZT6 window, 7 days light-dark acclimation) |
| Sequencing depth and alignment rates reported | **1** (`advsci2026_psnp_circadian`) |
| Full particle characterisation (measured size + zeta + PDI) | **3** (`advsci2026_psnp_circadian`, `jnanobiotech2026_ncor1_ppara`, `natmetab2026_size_macrophage` — the last reports values only in an unread supplement) |
| Endotoxin testing of particles | **0** |
| Blinding stated for at least histology | **5** |
| Individual housing (removes cage confounding) | **2** (`ehp2023_ps_hfd`, `ijbs2025_psnp_lysosome`) |
| Matched purified control diet | **2** (`ehp2024_ps_hfd_obesity` D12450B vs D12492; `jung2026` D09100304 vs D09100310) |

**Cage confounding is pervasive and almost never modelled.** `jung2026` housed mice in groups of
4–5 per treatment, so cage is co-extensive with group. `toxsci2024_obesogenic_ps` allocated by
cage (2 cages per group). `ehp2024_ps_hfd_obesity` used 3 mice/cage for a microbiota-driven
phenotype without any cage term. `advsci2026_psnp_circadian` declares the mouse the statistical
unit while delivering exposure through shared drinking water in group housing.

**Circadian phase is almost entirely ignored** — notable because several of the genes at issue
(*Fasn*, *Scd1*, *Elovl6*) are strongly diurnal, and because `biology_2026_ps_rorg` builds its
entire mechanism on RORγ, a clock-linked nuclear receptor, without reporting sampling time.
`jnanobiotech2026_ncor1_ppara` fasted mice overnight before sacrifice, which independently
activates the PPARα programme it then reports as suppressed.

**`advsci2026_psnp_circadian` is the methodological benchmark of this set** — and even its authors
refuse to read a single-time-point clock measurement as rhythm disruption: "interpretation of this
single-time-point measurement requires caution because circadian function is determined by dynamic
parameters, including phase, amplitude, and periodicity."

---

## 9. Question 8 — What is validated versus predicted

The team asked not to confuse computational predictions with experimentally validated findings.
The distinction turns out to be load-bearing in this literature.

**Strongest experimental validation:**
- `natmetab2026_size_macrophage` — ex-vivo co-culture showing Kupffer cells from exposed mice are
  **sufficient** to induce lipid droplets in naive hepatocytes; in-vivo phagocytosis and bacterial
  clearance assays; genetic fate-mapping; diphtheria-toxin Kupffer-cell depletion.
- `ehp2024_ps_hfd_obesity` — antibiotic depletion abolished the phenotype; faecal transplant
  transferred it.
- `liverint2025_np_masld` — SHH-neutralising antibody attenuated stellate-cell activation
  (in vitro).
- `advsci2026_psnp_circadian` — pharmacological rescue (AICAR, SR9009) both in vitro **and in vivo**.

**Mechanisms that are inferred, not demonstrated:**
- `jung2026`'s PPARα → *Anxa2* axis: the ChIP-qPCR is real, but the **Ppara-null evidence comes
  from an external dataset** (GSE8295, different strain, no PE exposure), and the agonist,
  antagonist and siRNA work is all in primary hepatocytes. **There is no in-vivo PE × Ppara-loss
  experiment.**
- `biology_2026_ps_rorg`'s RORγ causality: no knockout, knockdown, agonist or antagonist rescue.
- `jnanobiotech2026_ncor1_ppara`: no genetic manipulation at all; the NCoR1–PPARα interaction
  central to the model is **never measured** (no co-immunoprecipitation, no ChIP); the single
  causal experiment is rapamycin in a cell line, in a supplementary figure.
- `sci_rep_2025_ps_hubgenes`'s hub genes: validated at mRNA level only, with no functional test.
- `plosone2026_pet_comp`: **entirely computational** — and it models PET as a soluble monomer-like
  fragment, not a particle, which makes its conclusions about microplastic exposure a non-sequitur
  however sound the docking is.

---

## 10. Where the course's statistical toolkit genuinely fits

The team asked where BT3041 methods could be scientifically useful — without forcing machine
learning into the project. Honest assessment:

**Strong fit, and justified by the gap in §4:**
- **Factorial modelling with an interaction term** (regression/ANOVA). This is the single clearest
  methodological gap in the literature. It is exactly what the course covers, and it needs no
  machine learning.
- **Hypothesis testing and multiple-testing correction.** Two papers in this set called DEGs at
  raw p < 0.05; one reports no thresholds at all. Pre-specifying and reporting a single correction
  scheme is a genuine contribution to rigour.

**Reasonable fit:**
- **PCA and clustering** for QC and confounder detection — already used in Phase 1, and the
  literature shows almost nobody does this (no paper in this set reports sample-level QC,
  library-size balance or outlier analysis for its transcriptome).
- **Cross-validation / model selection** if, and only if, the team builds a predictive model;
  there is no obligation to.

**Weak fit — proceed with caution:**
- **Classification of exposure status from expression.** With *n* = 4–5 per group and technical
  factors collinear with treatment (the team's own Phase 2 finding), a classifier would most
  likely learn library depth, not biology. `plosone2026_pet_comp` is a cautionary example of how
  impressive-looking computational output can rest on nothing.

---

## 11. Summary table: what is settled, what is preliminary, what is open

| Claim | Status | Key sources |
|---|---|---|
| Oral plastic exposure changes the liver transcriptome | **Established** (direction varies) | all in-vivo transcriptome papers |
| Plastic alone on a control diet produces little liver pathology | **Established** across polymers | `liverint2025_np_masld`, `jnanobiotech2026_ncor1_ppara`, `ehp2024_ps_hfd_obesity`, `liverint2026_lowdose_alleviate` |
| *Ppara* is induced by **polyethylene** specifically | **Replicated twice**, independently | `jung2026`, `ees2023_pe_liver_fibrosis` |
| PPARα is *suppressed* by polystyrene nanoplastics | **Preliminary**, one strong study, conditions differ | `jnanobiotech2026_ncor1_ppara` |
| Plastic "exacerbates" diet-induced liver disease | **Widely claimed, never formally tested** | all nine 2×2 studies |
| Lipogenic genes go down with PE while lipid content rises | **Preliminary**; supported by one independent polymer/platform, contradicted by the one independent PE study | `jung2026`, `natmetab2026_size_macrophage` vs `ees2023_pe_liver_fibrosis` |
| Effects are gut-mediated rather than direct hepatic | **Preliminary but converging** | `ees2023_pe_liver_fibrosis`, `ehp2024_ps_hfd_obesity`, `advsci2026_pet_gutliver`, `natmetab2026_size_macrophage` |
| Microplastic exposure causes hepatic oxidative stress | **Contested** — explicit nulls exist | `biology_2026_ps_rorg` vs `sci_rep_2025_ps_hubgenes`, `toxsci2024_obesogenic_ps` |
| Low doses may be protective rather than harmful | **Open**, one abstract-only report, dose unstated | `liverint2026_lowdose_alleviate` |
| Microplastics are present in human liver | **Established for cirrhotic liver** in one small case series; causality open | `horvatits2022_human_cirrhosis` |

---

## 12. Explicitly unverified items

Stated here so nothing in this review is mistaken for a checked fact:

- Three papers were read as **abstracts only**: `stoten2024_pe_tlr2`, `stoten2024_scrna_nafld`,
  `liverint2026_lowdose_alleviate`. Their design fields, sample sizes, doses and statistics are unknown.
- **Supplementary files were not read** for most papers. Where a number lives only in a supplement,
  the per-paper record says so.
- For `horvatits2022_human_cirrhosis`, **which six polymers were detected is not verified here** —
  relevant because Jung motivates its PE focus by asserting PE is the most abundant polymer in
  human liver, a claim that should be traced to its own source.
- Whether the 19 Xenium samples in `jung2026` are the same 19 mice as the bulk libraries
  **is not stated in the paper**.

---

## References

Full bibliographic details, DOIs and links are in [`paper_matrix.csv`](paper_matrix.csv)
(columns `title`, `authors`, `journal`, `doi`, `doi_link`, `pmid`, `pmcid`) and at the top of each
per-paper record in [`papers/`](papers/). Metadata were retrieved from the Europe PMC REST API on
20 November 2026 and are archived in [`papers/_citations.json`](papers/_citations.json).
