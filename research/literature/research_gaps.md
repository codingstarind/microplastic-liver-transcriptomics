# Research gaps

**Date:** 20 November 2026. Companion to [`literature_review.md`](literature_review.md),
[`paper_matrix.csv`](paper_matrix.csv) and [`dataset_audit.md`](dataset_audit.md).

> **What this document is.** A sorted list of what the literature has and has not settled, in
> three tiers of confidence. **It does not propose hypotheses and it does not recommend a
> dataset.** Those are the team's decisions. Items are phrased so that each one can be checked
> against the per-paper records in [`papers/`](papers/) before anyone acts on it.

---

## Tier A — Established gaps
*Things the published literature demonstrably has not done. We checked every paper in the set;
the absence is verifiable, not an impression.*

### A1. No study has formally tested a plastic × diet interaction

Nine papers use a 2×2 design (or 2×2-like, with genotype or age as the second factor). **None
reports an interaction coefficient, F-statistic or p-value.** Seven analyse the four groups by
one-way ANOVA with pairwise post-hoc tests; two (`jnanobiotech2026_ncor1_ppara`,
`advsci2026_psnp_circadian`) fit a two-way model for phenotypes but never print the interaction
term, and both analyse their transcriptomes as pairwise contrasts only.

Consequence: essentially every "exacerbation" or "synergy" claim in this field rests on the
difference-of-significance pattern (*significant on high-fat diet, not significant on control
diet*). Where the data are visible, that pattern is sometimes contradicted by the papers' own
results — `ehp2023_ps_hfd` reports microplastic effects under the normal diet for glucose
tolerance and intestinal mucus; `ijbs2025_psnp_lysosome` claims "no effects ... on a normal chow
diet" in its abstract while reporting 361 chow-arm DEGs plus increased macrophages, neutrophils
and collagen.

**Status:** this is the single clearest, most defensible gap identified. It is a gap in
*statistical approach*, not a missing biological observation — and it should be described that
way to avoid overclaiming.

### A2. Sampling time is almost never controlled, and circadian confounding is unexamined

One study out of sixteen full-text papers reports a controlled collection window
(`advsci2026_psnp_circadian`: all collections at ZT4–ZT6 after 7 days of light–dark
acclimation). No other paper in the set reports time of sacrifice at all.

This matters concretely, not abstractly:
- Several genes central to this literature (*Fasn*, *Scd1*, *Elovl6*) are strongly diurnal.
- `biology_2026_ps_rorg` builds its entire mechanism on RORγ — a clock-linked nuclear receptor —
  without reporting sampling time.
- `jnanobiotech2026_ncor1_ppara` fasted mice overnight before sacrifice, which independently
  activates the hepatic PPARα programme it then reports as suppressed.

No published study has examined whether microplastic-associated hepatic signatures are confounded
by, or interact with, circadian phase. The team has already identified a clock-phase difference in
its own control group, so this is a gap the project is unusually well placed to speak to.

### A3. Technical covariates of RNA-seq are essentially unreported

Across nine transcriptome studies:
- **Sequencing depth is reported once** (`advsci2026_psnp_circadian`).
- **Per-group RNA-seq *n* is missing** from `jung2026`, `biology_2026_ps_rorg` and
  `ijbs2025_psnp_lysosome`.
- **DE tool, thresholds and correction are entirely absent** from `advsci2026_pet_gutliver`,
  whose DEG counts therefore cannot be reproduced.
- **Two studies called DEGs at uncorrected p < 0.05** (`biology_2026_ps_rorg`,
  `ijbs2025_psnp_lysosome`).
- **No study reports sample-level QC**, library-size balance across groups, or outlier analysis
  for its transcriptome.

No published microplastic liver transcriptome has been assessed for the kind of
depth/quality-versus-treatment collinearity the team found in GSE294180.

### A4. Cage is confounded with treatment and never modelled

`jung2026` housed mice in groups of 4–5 per treatment (cage ≡ group). `toxsci2024_obesogenic_ps`
allocated by cage (2 cages per group). `ehp2024_ps_hfd_obesity` used 3 mice per cage for an
explicitly microbiota-mediated phenotype with no cage term. `advsci2026_psnp_circadian` declares
the mouse the statistical unit while delivering exposure through shared drinking water in group
housing. Only `ehp2023_ps_hfd` and `ijbs2025_psnp_lysosome` housed animals individually.

### A5. No endotoxin testing anywhere

**Zero of nineteen papers report endotoxin or LPS testing of their particles.** This is most
striking in `ehp2024_ps_hfd_obesity`, where **serum LPS is a primary outcome** and the mechanism
is metabolic endotoxaemia — yet the particles themselves were never screened.

### A6. Polyethylene is severely under-studied relative to polystyrene

Of nineteen papers, **two** study polyethylene by oral exposure in rodent liver
(`ees2023_pe_liver_fibrosis`, `stoten2024_pe_tlr2` — the latter abstract-only), plus the project's
source paper. Everything else is polystyrene (ten papers), PET (two), or not applicable. PE is
simultaneously the polymer Jung et al. motivate as most abundant in human liver and the one with
the thinnest independent evidence base.

---

## Tier B — Possible limitations requiring verification
*Things that look like problems but that we could not fully check. Each carries the specific
action that would settle it. Do not treat any of these as established.*

### B1. Is Xenium slide confounded with PE exposure in the Jung spatial data?

In the deposited `cell_metadata.csv` (930,743 cells), barcode prefix `G12_` contains **only the
two vehicle groups** and `G34_` **only the two PE groups**. The paper states the 19 samples sat on
**two Xenium slides**. If `G12`/`G34` are those slides, slide is completely confounded with PE —
structurally the same problem as the depth/quality collinearity the team found in the bulk data.

**Not confirmed.** The 419 MB `.h5` was not fully retrieved and the record has no README.

**To settle it:** download `MP_PE_Mouse_Liver_04162026.h5` and inspect its metadata slots for a
sample or slide field; or ask the corresponding authors for the sample-to-slide mapping.

### B2. Are the 19 Xenium samples the same 19 mice as the bulk libraries?

The paper says only "19 samples from two Xenium slides" and refers to "matched mouse samples".
It never states the correspondence. Until confirmed, Fig. 4A (pseudo-bulk versus bulk correlation)
should not be read as a within-animal validation.

**To settle it:** ask the authors, or check whether the `.h5` carries sample identifiers that map
to GSM titles.

### B3. What exactly are the two re-used PPARα ChIP-seq runs?

`SRX5028180` and `SRX5142528` are labelled `PPARa_ChIPseq_NCOR_WT_1` and
`PPARa_ChIPseq_GPS2_WT_3`, and both belong to study **SRP140448**, a GPS2 liver-knockout
cistrome study. Jung uses them to argue PPARα binds the *Anxa2* enhancer/promoter, names no peak
caller and no input control, and lists these accessions only in Methods — not in its Data
Availability statement.

**To settle it:** read the source paper (Liang et al., *Nat Commun* 10:1684, 2019) and the GSM
records for genotype, treatment and matched input.

### B4. Does the lipogenesis direction conflict survive the condition differences?

`jung2026` reports *Fasn* and *Acly* **down** with PE; `ees2023_pe_liver_fibrosis` reports *Fasn*,
*Acly* and *Acaca* **up** with PE. Both report increased hepatic lipid content. The studies differ
in sex (male vs female), route (gavage of a fixed 2 mg/day vs 100 µg/g in food), particle size
distribution, and the presence of a MASH-diet arm. Djouina applied no multiple-testing correction
across a large qPCR panel.

`natmetab2026_size_macrophage` independently reports lipogenesis **suppressed** with uptake up —
agreeing with Jung, in a different polymer and platform.

**To settle it:** nothing in the literature can; it would need an experiment matching sex and
route. What the team *can* do is state the conflict precisely and not resolve it by assertion.

### B5. Is PPARα induced or suppressed, and does fasting state explain the discrepancy?

Three studies report *Ppara* up (`jung2026`, `ees2023_pe_liver_fibrosis`,
`sci_rep_2025_ps_hubgenes`); one reports PPARα signalling suppressed
(`jnanobiotech2026_ncor1_ppara`). The suppression study fasted all mice overnight before
sacrifice — a strong PPARα activator — and **never measured *Ppara* itself**, only targets and the
corepressor NCoR1. Jung reports induction relative to CD+Veh, **a contrast that mixes diet and PE**.

**To settle it:** check whether a within-diet PE contrast in GSE294180 shows *Ppara* induction at
all (the team's Phase 2 comparison suggests it does not), and record fasting state as a covariate
when comparing across studies.

### B6. Is the hepatic phenotype direct or gut-mediated?

Four independent lines point away from direct hepatic particle action:
`ees2023_pe_liver_fibrosis` (authors judge particles probably never reached the liver);
`ehp2024_ps_hfd_obesity` (antibiotics abolished the phenotype, faecal transplant transferred it);
`advsci2026_pet_gutliver` (full NAFLD-like phenotype with **no PET detectable in liver** by
pyrolysis-GC/MS); `natmetab2026_size_macrophage` (no particles in hepatocytes; Kupffer cells are
the reservoir, and Kupffer cells alone are *sufficient* to induce hepatocyte lipid droplets).

Against this, `jung2026` reports PE-like particles in MD+PE liver by micro-Raman and O-PTIR — but
with no procedural blanks, no CD-group tissue examined, no particle quantification, and the Raman
performed on slides already processed for Xenium.

**To settle it:** not resolvable from a liver transcriptome alone. Relevant because it changes what
a bulk liver signature *means* — a direct hepatocellular response versus a downstream consequence
of gut and immune changes.

### B7. Does the "low-dose protective" result hold, and at what dose?

`liverint2026_lowdose_alleviate` reports steatosis, inflammation and fibrosis all **reduced** by
low-dose polystyrene nanoplastics on a high-fat diet, with AMPK up, mTOR down and lipophagy
proposed — the inverse of `ijbs2025_psnp_lysosome`'s mTOR-up/autophagy-blocked mechanism.

**Unverifiable as it stands:** abstract-only, **no numeric dose given**, particle size not stated,
group structure ambiguous (Methods name two groups, Results imply a third), and it is the only rat
study in the set.

**To settle it:** obtain the full text (Wiley, subscription) and extract dose, size, route, *n* and
group structure.

### B8. Was polyethylene among the polymers detected in human liver?

`jung2026` motivates its PE focus by asserting PE is the polymer most abundant in human liver.
`horvatits2022_human_cirrhosis` detected **six polymers, 4–30 µm**, in cirrhotic liver, but we did
not verify which six. Its abbreviation list includes PS, PET, PVC, PMMA, POM and PP.

**To settle it:** read the Results/Table of the EBioMedicine paper, and trace Jung's claim to the
source it actually cites.

---

## Tier C — Unanswered questions
*Open scientific questions the field has not addressed. Listed as questions, deliberately not as
hypotheses.*

**C1.** Does a formally modelled diet × plastic interaction exist at gene or pathway level, or do
the apparent "exacerbations" dissolve into two additive main effects? No published study provides
the answer, in any polymer.

**C2.** Are microplastic-associated hepatic expression changes separable from circadian phase and
feeding rhythm? High-fat diets themselves shift murine feeding rhythms, which shift the liver
clock — so diet, clock and treatment are entangled in every study that did not control sampling
time.

**C3.** Does the direction of the lipid-metabolism response depend on polymer, particle size,
dose, route, sex or metabolic background — and in what order of importance? Individual reversals
are documented for **size** (`natmetab2026_size_macrophage`: 500 nm increased, 50 nm decreased
hepatic lipid under HFD), for **metabolic background** (`ijms2024_ps_lep`: steatosis up in
wild-type, down in leptin-knockout), and possibly for **dose**
(`liverint2026_lowdose_alleviate`). No systematic study crosses these factors.

**C4.** How much of a bulk liver "plastic effect" is a change in cell-type composition rather than
per-cell expression? `natmetab2026_size_macrophage` shows particles sit in Kupffer cells, not
hepatocytes, and that Kupffer-derived signals alone suffice to make hepatocytes accumulate lipid;
`stoten2024_scrna_nafld` reports macrophage subpopulation shifts; `jung2026` reports a midlobular
hepatocyte increase and Kupffer/stellate/lymphocyte decrease. Bulk data cannot distinguish the two
by construction.

**C5.** Is there a reproducible PE-specific hepatic signature at all, independent of laboratory?
With only two independent PE oral studies (one abstract-only, one qPCR-only), this cannot
currently be answered.

**C6.** Do microplastic-associated mouse liver genes correspond to anything in human MASLD liver?
Three public human cohorts exist (GSE37031, GSE63067, GSE89632) and only one computational paper
has used them in this context — with PET modelled as a soluble monomer rather than a particle.

**C7.** What is the minimum particle dose that produces a detectable hepatic transcriptional
effect? Only `natmetab2026_size_macrophage` ran a dose titration (100 / 10 / 1 mg/kg), finding
glycogen loss at 100 and 10 but not 1 mg/kg. Every other study used a single dose level.

---

## What this review does **not** establish

Stated plainly so the team is not misled by its own evidence base:

- **It does not identify a "best" dataset.** `dataset_audit.md` lists what exists and what each
  would permit; the choice is the team's.
- **It does not establish that any finding in the team's Phase 2 is correct or incorrect.** The
  literature supplies context and directional priors, not validation of a specific re-analysis.
- **It does not prove Jung et al. are wrong about anything.** Several of their claims rest on
  evidence types this review flags as inferential (external datasets, in-vitro pharmacology,
  pairwise contrasts). That is a statement about *evidence strength*, not about correctness.
- **It cannot settle the direction conflicts** in §B4 and §B5 — the experiments that would settle
  them have not been done.
- **Three papers were read as abstracts only**, so the field may contain relevant detail this
  review does not capture.
- **Supplementary files were not read** for most papers. Where a load-bearing number lives only in
  a supplement, the per-paper record says so explicitly.

---

## Suggested next decisions for the team (not conclusions)

These are the choices the evidence base now makes possible. Each is a decision, not a
recommendation:

1. **How to word the novelty claim.** The interaction-testing gap (A1) is real and verifiable, but
   it is methodological. Decide whether the report claims "first formal interaction test in this
   literature" (defensible, narrow) or something broader (not supported).
2. **Whether to pursue the Xenium slide-confounding question (B1).** It would either strengthen
   the team's technical-confounding argument considerably or close it off. One download settles it.
3. **Whether to use GSE245069 as an independent comparison.** It is the only independent liver
   transcriptome with public raw counts, but it is a different polymer, strain and design.
4. **Whether to obtain the three subscription papers**, particularly
   `liverint2026_lowdose_alleviate` (B7), which is the only direction-discrepant report found.
5. **How to report the PPARα conflict (B5)** — as a genuine conflict with documented condition
   differences, or as a condition-dependent effect. Both framings are defensible; asserting either
   without the caveats is not.
