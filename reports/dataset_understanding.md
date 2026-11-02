# Dataset understanding: from mouse to matrix (GSE294180)

Sources checked:

- the GEO series record and all 19 GSM sample records (SOFT file, downloaded 2 Nov 2026);
- the full text of Jung *et al.*, *Sci. Adv.* 12(25): eaec8681 (2026), read through Europe PMC (PMC13274612, CC BY licence).

Each fact below gives its source. Where GEO and the paper disagree, both are reported.

## 1. Biological motivation

- **Microplastics (MPs)** are plastic particles of 1 µm–5 mm. Polymers including polyethylene (PE), polystyrene (PS), polypropylene and PET have been reported in human blood, placenta and testes. The authors state that PE is the polymer most abundant in human liver.
- **Prior work** on MPs and liver mostly used PS beads and bulk RNA-seq. It suggested that MPs cause inflammation and disturb lipid metabolism, which could feed into **MASH** (metabolic dysfunction-associated steatohepatitis, formerly NASH).
- **Nuclear receptors** (NRs) such as PPARα control hepatic fatty-acid oxidation, inflammation and detoxification. PPARα agonists are being tested in clinical trials for MASH.
- **The authors' aim:** test whether PE causes liver injury, whether it worsens diet-induced steatohepatitis, and where in the liver lobule and in which cells this happens. They used bulk RNA-seq for whole-tissue gene expression and 10x Xenium spatial transcriptomics for cell-level, zone-resolved expression. Their mechanistic conclusion is a PPARα → *Anxa2* regulatory link.
- **Why PE:** in an in-vitro screen of seven polymers on primary hepatocytes (24 h, 100 µg/mL), PE reduced viability the most (paper Fig. S1A).

## 2. Experimental design

**A 2 × 2 factorial design.** The experimenters set two factors:

| Factor | Levels | What exactly |
|---|---|---|
| Diet | **CD**: control diet | Research Diets **D09100304**: 10 kcal% fat, sucrose-matched. A *purified* low-fat diet, although GEO calls it "Chow Diet". |
| | **MD**: MASH diet | Research Diets **D09100310**: 40 kcal% fat, 20 kcal% fructose, 2% cholesterol |
| Oral gavage | **Veh**: vehicle | 0.1% Tween 20 in molecular-grade water, daily |
| | **PE**: polyethylene | **2 mg PE per mouse per day**, suspended in the vehicle, daily |

**Microplastic.** Polyethylene microspheres, **10–150 µm** in diameter, from Cospheric (catalogue CPMS-0.96 in the in-vitro methods). The size distribution within that range is not reported. The dose is a fixed mass per mouse, not per kg body weight. For a 25–30 g mouse it is roughly 70–80 mg/kg/day; that figure is our estimate, not a reported value. The paper detected "PE-like particles" in liver tissue by micro-Raman and O-PTIR spectroscopy.

**Duration.** Diet and gavage ran for **8 weeks**. The text implies they started together but does not say so explicitly.

**Animals** (paper Methods):

- male C57BL/6 mice from the Jackson Laboratory, 8 weeks old at the start;
- 22 ± 1 °C, 25 ± 5% humidity, 12 h light/dark cycle, food and water *ad libitum*;
- "Groups of mice (n = 4 or 5 per group) were housed together". Cage is therefore probably the same as treatment group, and cage effects cannot be separated from treatment effects.

**Controls.**

- CD+Veh is the reference for both factors.
- MD+Veh is the control for PE within the MASH diet.
- The vehicle gavage controls for the handling and stress of daily gavage and for Tween 20.
- There is no sham (non-gavaged) group and no control particle (e.g. an inert non-plastic bead of the same size).

**Sex.** GEO records only "genotype: WT". Sex comes from the paper and is confirmed by our data: Xist ≈ 0 and Y-chromosome genes are expressed in all 19 samples (see the audit report).

## 3. Independent biological samples

| Group | GEO samples (GSM8899917–35) | Libraries |
|---|---|---|
| CD+Veh | 4 (Liver_CD_1–4) | S39–S42 |
| CD+PE | 5 (Liver_CD_MP_1–5) | S48–S52 |
| MD+Veh | 5 (Liver_MD_1–5) | S43–S47 |
| MD+PE | 5 (Liver_MD_MP_1–5) | S53–S57 |

- **19 mice, one liver sample each.** Each library is a distinct mouse, so these are biological replicates; there are no technical replicates.
- The paper's Fig. 1B says *n = 5 within each of the four groups* for body weight, but GEO has only 4 CD+Veh samples. Either one CD+Veh mouse was not sequenced or the legend is approximate. The paper does not say.
- The GEO extraction protocol says RNA was extracted from "four randomly selected (using Randomice tool)" livers. That fits CD+Veh (4) but not the other groups (5). This is unresolved; see the audit report.
- The paper says 19 samples were also profiled on two Xenium slides. These may be the same 19 mice, but that is not confirmed.

## 4. What was measured, and how

**Physiology** (paper Fig. 1, not deposited in GEO):

- body weight (weekly) and liver-to-body weight ratio;
- blood glucose, serum ALT, liver triglycerides;
- H&E histology with blinded NAFLD activity score (steatosis, lobular inflammation, ballooning).

The paper reports that PE raised ALT on both diets, raised liver triglycerides on the MASH diet, and increased steatosis scores. PE also increased body weight on CD (significant from week 7) but not on MD.

**Bulk RNA-seq** (the data we use):

```
Mouse liver (snap-frozen) ─► ~50 mg tissue ─► TRIzol total RNA ─► Bioanalyzer QC
   ─► NEBNext poly(A) mRNA isolation (500 ng input per GEO; "up to 1 µg" per paper)
   ─► IDT xGen Broad Range stranded library prep, indexed (S39–S57)
   ─► pooled, Illumina NextSeq 2000, 2 × 150 bp, target ~15 M read pairs/sample
   ─► fastp v0.23.4: adapter trimming, poly-G / poly-X removal
   ─► Salmon v1.10.3, selective alignment (--validateMapping), auto library type (-l A),
       --seqBias, --useVBOpt, against the GENCODE GRCm39 M35 transcriptome
       → per-transcript estimated read counts (quant.sf)
   ─► tximport (R) with a transcript-to-gene table → gene-level estimated counts
   ─► [the deposited CSV] 55,836 gene symbols × 19 samples, integer values
   ─► authors: DESeq2 (size factors, Wald tests), |log2FC| ≥ 0.38, FDR < 0.05; fgsea GSEA
```

Steps that GEO and the paper **do not** document:

- the tximport `countsFromAbundance` setting;
- whether and where counts were rounded;
- that the transcript-to-gene table maps to **gene symbols**, not Ensembl IDs. We inferred this from the file; it means genes sharing a symbol were probably summed;
- Salmon mapping rates, decoys, GC-bias correction, RIN values;
- time of day and fasting state at sacrifice.

## 5. Bulk vs spatial vs single-cell transcriptomics

| | Bulk RNA-seq (GSE294180, this project) | Single-cell RNA-seq | Spatial transcriptomics (Xenium, in the paper) |
|---|---|---|---|
| Unit measured | Whole tissue piece: an average over millions of cells | Individual dissociated cells | Individual cells *in place* in a tissue section |
| Genes | Whole transcriptome (~20k protein-coding genes detectable) | Whole transcriptome (3′ or 5′ tag) | Targeted panel: 379-gene Mouse Tissue Atlassing panel + 100 custom genes |
| Spatial context | Lost | Lost (cells dissociated) | Kept: zonation (portal ↔ central), neighbourhoods |
| What a difference means | A change in expression per cell **or** a change in cell-type composition; cannot tell which | Can attribute changes to cell types | Can attribute changes to cell types *and* liver zones |
| Data type | Counts per gene per sample | Counts per gene per cell | Transcript counts per gene per segmented cell + coordinates |
| Statistical unit | Mouse | Cell (but replicate is still the mouse) | Cell / region (replicate is still the mouse) |

Single-cell RNA-seq was **not** done in this study. The paper uses "spatial single-cell" to mean Xenium imaging. The Xenium data are on Zenodo (record 19603085), not GEO. Because bulk data average over cell types, a bulk "PE effect" could be a change in cell composition, for example more immune cells. The paper reports a shift towards midlobular hepatocytes and fewer Kupffer cells in MD+PE.

## 6. What the authors already did with these bulk data

To avoid simply repeating their work:

- **PCA** (DESeq2 `plotPCA`): "distinct transcriptomic profiles among the experimental groups" (Fig. 2A).
- **Pairwise DESeq2**, |log2FC| ≥ 0.38 and FDR < 0.05:
  - PE-specific DEGs: 503 genes under CD only, 182 under MD only, 97 shared;
  - other counts quoted: 600 PE-altered, 1,437 MD-altered and 279 MD+PE vs MD (Fig. 2B–C, tables S5–S7).
- **Named genes:**
  - down with PE: Fasn, Acly, Me1, Pnpla3, Angptl8, Lepr, Hmgcr;
  - up with PE: Plin4, mt-Nd6, Cyb5a, Cyp4a12a, Gadd45g;
  - MD+PE-specific: Ldlr, Tm7sf2, Mylip, Pask, Oit3, Nkap;
  - Ppara up in PE groups, Nr1i3 (CAR) up in MD+PE.
- **GSEA (fgsea)** with a custom "LIPID_METABOLISM_AND_TOXICITY" NR gene set: enriched for PE on both diets, strongest for MD+PE vs CD+Veh.
- **NR expression heatmap** (Fig. 2K) and a comparison with a PS-microplastic liver dataset (GSE245069): Ppara up and Esr down in both.
- **Spatial and mechanistic work** (Figs 3–6): 15 cell clusters, zonation, inflammatory "hotspots" with lower cell-type diversity, WGCNA on pseudo-bulk central hepatocytes, and Ppara–Anxa2 ChIP-seq / knockout / siRNA validation.

**Not reported by the authors for the bulk data:**

- any formal test of the **diet × PE interaction** (they used pairwise contrasts);
- sample-level QC: library sizes, contamination, outliers;
- sensitivity of the PCA to preprocessing;
- the confounding of library order and size with PE treatment;
- clustering of samples without supervision.

These gaps are where our analysis can add something.
