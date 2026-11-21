# Dataset audit: what public data actually exist, and what each would permit

**Date of checks:** 20 November 2026. Every record below was queried directly (GEO SOFT/FTP,
NCBI E-utilities, Zenodo API, PRIDE API, Metabolomics Workbench REST, Dryad API, GitHub API,
DOI resolution). **"Verified" means we retrieved the record itself**, not that the paper
mentioned it.

> **Headline.** Of 19 papers, only **four** deposit liver molecular data that could be
> re-analysed, and only **three** of those are transcriptomes. Eight papers generated omics and
> deposited nothing. One cites an accession that is **not public**.

---

## 1. Summary table

| Accession | Paper | Type | Status | Reanalysable? |
|---|---|---|---|---|
| **GSE294180** | `jung2026` | Mouse liver bulk RNA-seq | ✅ Public | **Yes** — the project's own dataset |
| **Zenodo 19603085** | `jung2026` | Xenium spatial, 930,743 cells | ✅ Public, CC-BY | **Partly** — see §3 |
| **GSE245069** | `sci_rep_2025_ps_hubgenes` | Mouse liver bulk RNA-seq | ✅ Public | **Yes** — raw counts included |
| **GSE291616** | `natmetab2026_size_macrophage` | Mouse liver snRNA-seq | ✅ Public | **Yes**, but n = 2/group |
| **PXD072158** | `natmetab2026_size_macrophage` | Kupffer-cell proteomics | ✅ Public (PRIDE) | Yes |
| **ST003776** | `natmetab2026_size_macrophage` | Liver lipidomics, 28 samples | ✅ Public (Metab. Workbench) | Yes |
| **bonndata FK2/M42B8O** | `natmetab2026_size_macrophage` | Histology images | ✅ DOI resolves | Images only |
| **PRJNA990701** | `ehp2024_ps_hfd_obesity` | 16S rRNA, mouse feces | ✅ Public, 20 runs | Microbiome only |
| **Dryad h44j0zpsj** | `toxsci2024_obesogenic_ps` | "Polystyrene-induced changes in adipose transcriptome" | ⚠️ Public but **adipose only** | **No liver data** |
| **GSE282998** | re-used by `jung2026` | PPARα/γ agonist TempO-Seq, 192 samples | ✅ Public | Yes |
| **GSE8295** | re-used by `jung2026` | WT vs *Ppara*-null microarray, 16 samples | ✅ Public | Yes |
| **SRX5028180 / SRX5142528** | re-used by `jung2026` | PPARα ChIP-seq | ✅ Public | Yes — **but see §5** |
| **GSE37031 / GSE63067 / GSE89632** | `plosone2026_pet_comp` | Human NAFLD liver arrays | ✅ All public (15 / 18 / 63 samples) | Yes |
| **GitHub + Zenodo 18925269** | `plosone2026_pet_comp` | Analysis code | ✅ Both resolve | Yes |
| **PRJNA1305997** | `biology_2026_ps_rorg` | Liver RNA-seq (claimed) | ❌ **NOT PUBLIC** | **No** |
| *(none)* | 8 further papers | — | ❌ No deposit | **No** |

---

## 2. GSE294180 — the project's own dataset (`jung2026`)

**Verified.** Series public; BioProject PRJNA1248512 with **19 SRA experiments**, matching the
19 GSM records. Supplementary file: `GSE294180_MP_raw_counts_gene_level.csv.gz`.

Already audited in depth by the team in
[`reports/dataset_understanding.md`](../../reports/dataset_understanding.md) — not repeated here.

**One item this literature audit adds:** the paper **never states the per-group RNA-seq *n***
(Fig. 1B says n = 5/group for phenotypes; Methods say "n = 4 or 5 per group"). The GEO record is
the authoritative source for the sequencing design (4 CD+Veh, 5 each for the others), and the team
has already recorded this discrepancy.

---

## 3. Zenodo 19603085 — the Xenium spatial data (`jung2026`)

**Verified open, CC-BY 4.0.** Record title: *"Spatial Ppara-Anxa2 Axis in Polyethyelen-Driven Liver
Damage"* (the typo is in the record). Three files:

| File | Size |
|---|---|
| `MP_PE_Mouse_Liver_04162026.h5` | 418.6 MB |
| `cell_metadata.csv` | 58.2 MB |
| `umap_embeddings.csv` | 48.2 MB |

**We downloaded and parsed `cell_metadata.csv` in full.** Contents:

- **930,743 cells.** Columns: `barcode`, `nCount_Xenium`, `nFeature_Xenium`, `Condition`,
  `SCT_snn_res.0.8`, `Main`, `Broad`.
- Cell counts by condition: Chow_Veh 210,258 · HFHF_Veh 204,674 · Chow_MP 290,882 · HFHF_MP 224,929.
- Broad cell types: Hepatocytes 629,882 · EC 113,112 · Myeloid 92,958 · Stellate 42,841 ·
  Cholangiocytes 26,963 · Hepato_Cholangio 12,974 · Lymphocyte 11,264 · Mesothelial 749.

### ⚠️ Two findings the team should weigh

**(a) There is no mouse identifier.** The metadata carry a `Condition` label only. Barcodes are
prefixed `G12_` or `G34_`, and those prefixes partition the data **exactly by treatment**:

| Barcode prefix | Conditions present | Cells |
|---|---|---|
| `G12_` | Chow_Veh, HFHF_Veh — **vehicle groups only** | 414,932 |
| `G34_` | Chow_MP, HFHF_MP — **PE groups only** | 515,811 |

The paper's Methods state: *"A total of 19 samples from two Xenium slides were integrated."* If
`G12`/`G34` correspond to those two slides, then **Xenium slide is completely confounded with PE
exposure** — the same class of technical confounding the team found in the bulk data
(library depth and quality collinear with PE).

**This is not yet confirmed.** We could not complete the 419 MB `.h5` download to check for a
per-sample or slide field inside the object, and the record carries no README. Two ways to settle
it: (i) download the `.h5` and inspect its metadata slots; (ii) email the corresponding authors
(taegyu-oh@ou.edu, aditya-joshi@tamu.edu) and ask for the sample-to-slide mapping.

**(b) Without a mouse ID, the deposited metadata cannot support mouse-level statistics.** The
paper's own diversity and cell-proportion analyses use mouse as a random effect and explicitly
state they do so "to avoid pseudoreplication" — meaning the authors held sample identity that the
deposit does not expose. Any re-analysis of these files at cell level would be pseudo-replicated.

**What it would permit as deposited:** descriptive comparison of cell-type composition between
conditions; verification of the reported midlobular-hepatocyte increase and Kupffer/stellate
decrease as *counts*; UMAP reproduction. **Not** valid per-mouse inference.

---

## 4. GSE245069 — independent PS liver transcriptome (`sci_rep_2025_ps_hubgenes`)

**Verified public** (released 23 Jan 2025). BioProject PRJNA1026899, **16 SRA runs**, platform
GPL28330 (DNBSEQ-T7).

- **16 liver samples, 4 per group:** H2O, 80 nm, 0.5 µm, 5 µm. Male ICR mice, 1 mg/L in drinking
  water, 12 weeks.
- Supplementary `GSE245069_gene_expression.txt.gz` (2.1 MB, 20,926 rows, 63 columns) contains
  **`read_count_*` columns for all 16 samples**, plus FPKM and TPM and per-group averages.
- A second file `GSE245069_rna_expression_.txt.gz` has transcript-level FPKM/TPM (no counts).

**This is the most reusable independent dataset found.** It would permit: a count-based
re-analysis with the team's own pipeline and thresholds; a direct test of whether *Ppara*
induction replicates in an independent polymer and lab; and a size-trend analysis the original
paper did not perform (it compared each size to control pairwise only).

**Caveats to carry:** outbred ICR mice (not C57BL/6), no diet factor, no vehicle/surfactant
control, and the original DE used limma with |log2FC| > 1.

---

## 5. Data re-used by `jung2026` — all verified, with one provenance flag

| Accession | What it actually is | Verified contents |
|---|---|---|
| **GSE282998** | TempO-Seq, FFPE liver, WT and **PPARα-knockout** mice given HFPO-DA, GW7647 (PPARα agonist), rosiglitazone or APAP | Public, **192 samples**, 4 gene-count CSVs, BioProject PRJNA1191425 |
| **GSE8295** | Affymetrix 430 2.0, **129S1/SvImJ WT vs *Ppara*-null**, Wy14643 0.1% in food, 5 days, n = 4/group | Public, 16 samples, RAW tar |
| **SRX5028180** | `GSM3478065: PPARa_ChIPseq_NCOR_WT_1` | Public; run SRR8208984 |
| **SRX5142528** | `GSM3517892: PPARa_ChIPseq_GPS2_WT_3` | Public; run SRR8330968 |

### ⚠️ Provenance flag on the ChIP-seq

Both SRX runs belong to study **SRP140448**, whose title is *"Transcriptome comparison between WT
and GPS2 liver knockout (LKO) livers and GPS2 NCOR PPARa cistrome and epigenome analysis in
livers."* `jung2026` cites them as reference 58 (a GPS2-knockout paper) and uses them to argue
that PPARα binds the *Anxa2* enhancer/promoter.

The two runs are labelled `_WT_` so they are plausibly the wild-type arms — but the team should
confirm genotype, treatment and input controls before citing them as evidence of
wild-type hepatic PPARα occupancy. `jung2026` names **no peak caller and no input control** in
its readable Methods; the figure shows signal tracks in IGV.

Also note: **these two accessions appear only in Jung's Methods, not in its Data Availability
statement.**

---

## 6. `natmetab2026_size_macrophage` — the best-documented deposit in the set

All four stated repositories verified:

| Accession | Verified contents |
|---|---|
| **GSE291616** | Public (released 9 Jun 2026). **6 samples: 2 PBS, 2 × 50 nm, 2 × 500 nm**, single-nucleus RNA-seq of C57BL/6J liver after 12 weeks of weekly gavage. Files: barcodes/features/matrix per sample (10x format), ~562 MB RAW tar. BioProject PRJNA1234520, 6 SRA runs. |
| **PXD072158** | Public at PRIDE. Title: *"Kupffer cells exposed to 12 weeks of microplastic (500nm), nanoplastic (50nm) or PBS control."* Mouse liver; published 13 Aug 2026. |
| **ST003776** | Public at Metabolomics Workbench. *"Investigation of Hepatic Lipid Alterations Following Micro- and Nanoplastic-Ingestion"*, Mus musculus, **28 samples (9 PBS, 9 × 50 nm, 10 × 500 nm)**, University of Bonn, CC BY-NC-ND. |
| **bonndata 10.60507/FK2/M42B8O** | DOI resolves to the Bonn repository (histology images). |

**What it would permit:** an independent test of whether hepatocyte de-novo-lipogenesis
suppression (which agrees in direction with `jung2026`) reproduces — but **n = 2 per group** makes
the snRNA-seq unsuitable for formal differential testing on its own. The lipidomics (n = 9–10) is
much better powered. No code repository.

---

## 7. ⚠️ PRJNA1305997 (`biology_2026_ps_rorg`) — the accession is not public

The paper's Data Availability Statement reads, verbatim: *"RNA-seq data are available in the SRA
under BioProject accession number PRJNA1305997."*

**Checked 20 November 2026:**
- NCBI E-utilities BioProject search for `PRJNA1305997` → **0 results** ("No items found").
- NCBI SRA search linked to that BioProject → **0 experiments**.
- The BioProject web page returns **"is not public"**.

So the only stated accession for that paper's 689-DEG transcriptome **cannot be retrieved**. The
paper also shows RORγ ChIP-seq signal tracks with no methods and no accession at all. Treat this
study as **not reanalysable**, and note that its DEG list was called at raw p < 0.05 with no
FDR control and no named DE tool.

---

## 8. ⚠️ Dryad h44j0zpsj (`toxsci2024_obesogenic_ps`) — real, but not the liver data

The paper's Data Availability points to `https://doi.org/10.5061/dryad.h44j0zpsj`.

**Verified.** The record exists (CC0, published 1 Feb 2024, 386 KB) and is titled
**"Polystyrene-induced changes in adipose transcriptome"**. Its complete contents:

```
README.md                  1,147 bytes
TO-DELPvsDMSO_down.xlsx  173,844 bytes
TO-DELPvsDMSO_up.xlsx    211,020 bytes
```

That is **adipose** gene lists from the delphinidin-versus-vehicle comparison. The paper's
**liver mRNA-seq (n = 10/group, 283 and 49 DEGs) is not deposited anywhere** — no GEO, no SRA, no
raw counts. Combined with the fact that the paper never names its DE tool, platform or depth, its
liver transcriptome results cannot be independently reproduced.

---

## 9. `plosone2026_pet_comp` — fully reproducible inputs, unvalidated biology

An instructive contrast: the weakest biology in the set has the best code availability.

- **GSE37031** — public, 15 samples (8 NASH, 7 control), GPL14877.
- **GSE63067** — public, 18 samples, steatosis and NASH, GPL570.
- **GSE89632** — public, 63 samples (20 SS, 19 NASH, 24 healthy donors), GPL14951.
- **GitHub `pipipeng288/Code`** — exists, created 9 March 2026, five folders: `CIBERSORT`,
  `Differential expression analysis`, `GO`, `KEGG`, `WGCNA`.
- **Zenodo 10.5281/zenodo.18925269** — resolves to record **18925270**, containing
  `pipipeng288/Code-v1.0.zip` (29 KB).

The analysis is reproducible; the **conclusions** are not experimentally supported (see the review,
§9).

---

## 10. Papers with **no** public data

Eight studies generated molecular data and deposited none. All state some version of "available
from the corresponding author on reasonable request":

| Paper | What was generated but not deposited |
|---|---|
| `ijbs2025_psnp_lysosome` | Full mouse liver bulk RNA-seq (4 groups, 2 time points) |
| `advsci2026_psnp_circadian` | **40 RNA-seq libraries** (8 groups × 5) |
| `jnanobiotech2026_ncor1_ppara` | Liver RNA-seq (n = 4 × 2 groups) **and** untargeted metabolomics (n = 6 × 2) |
| `advsci2026_pet_gutliver` | Liver bulk RNA-seq (n = 4/group × 2 time points) **and** 16S |
| `liverint2025_np_masld` | No omics generated (histology/qPCR/western only) |
| `ees2023_pe_liver_fibrosis` | No omics generated (qPCR panels only) |
| `ijms2024_ps_lep` | No omics generated |
| `ehp2023_ps_hfd` | 16S rRNA sequencing — **no data availability statement at all** |
| `horvatits2022_human_cirrhosis` | Analytical data only |

`advsci2026_psnp_circadian` is the most regrettable: it is the methodologically strongest study in
the set (ZT-controlled sampling, full particle characterisation, reported depth and alignment
rates, two-way ANOVA) and its 40 libraries would be the ideal external comparison — but they are
not public.

---

## 11. What each usable dataset would permit, in plain terms

| Dataset | A reanalysis could ask |
|---|---|
| **GSE294180** (own) | The factorial interaction question the source paper did not test. |
| **GSE245069** | Does *Ppara* induction replicate in an independent lab, polymer and mouse strain? Is there a particle-size trend? Counts are available, so the team's own pipeline and thresholds apply. |
| **GSE291616** | Does hepatocyte lipogenesis suppression reproduce at single-nucleus resolution? (Exploratory only — n = 2/group.) |
| **ST003776** | Is there a lipid-species signature accompanying the transcriptional changes? (n = 9–10, decent power.) |
| **GSE282998 / GSE8295** | What does a *bona fide* PPARα agonist or *Ppara* loss do to the genes the team sees moving? A pharmacological/genetic reference axis for interpreting PE effects. |
| **Zenodo 19603085** | Cell-composition comparison only, pending resolution of the slide-confounding and mouse-ID questions (§3). |
| **GSE37031 / 63067 / 89632** | Human MASLD reference — useful for asking whether mouse PE-associated genes are MASLD-associated in human liver. |

**No recommendation is made here about which dataset to use.** That is the team's decision.

---

## 12. Method notes

Queries used: GEO `acc.cgi` SOFT (`targ=self`/`targ=gsm`, `form=text`), GEO FTP `suppl/` and
`filelist.txt`, NCBI E-utilities `esearch`/`esummary` (databases `sra`, `bioproject`),
`zenodo.org/api/records/{id}`, `ebi.ac.uk/pride/ws/archive/v3/projects/{id}`,
`metabolomicsworkbench.org/rest/study/study_id/{id}/summary` and `/factors`,
`datadryad.org/api/v2`, `api.github.com/repos`, and DOI resolution via `doi.org`.
Europe PMC `resultType=core` supplied bibliographic metadata
([`papers/_citations.json`](papers/_citations.json)).

The Xenium `cell_metadata.csv` was downloaded in full (58.2 MB) and parsed with Python; the
`.h5` object (418.6 MB) was **not** fully retrieved, which is why §3(a) remains open.
