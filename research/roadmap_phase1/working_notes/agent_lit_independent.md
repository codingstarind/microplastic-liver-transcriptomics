# Independent methods verification: microplastics and rodent liver transcriptomics

Verified 2026-11-21 from primary sources only. The local prior notes were re-checked against the sources, not copied. Labels used below: **FT** means full text read. **ABS** means abstract only. **UNVERIFIED** means the field could not be read in any source I accessed.

---

## Part 1, Paper 1. Toxicol Sci 2024: "Obesogenic polystyrene microplastic exposures disrupt the gut-liver-adipose axis"

- **Citation:** Zhao J, ... Cave MC, O'Toole TE. Toxicol Sci 2024;198(2):210-220.
- **Links:** DOI https://doi.org/10.1093/toxsci/kfae013. PMID 38291899, PMC10964747.
- **Access: FT (PMC HTML).** I read the Methods, Results, all figure and table legends, the Discussion and the Data availability statement.
- **Supplement:** the supplementary zip (Supp Figs 1-9, Supp Tables 1-6) was **NOT read**. PMC serves it behind a proof-of-work bot check, which I did not bypass, and Europe PMC has no supplement for this article.
- **Dryad (read via API):** doi:10.5061/dryad.h44j0zpsj. The dataset is titled "Polystyrene-induced changes in adipose transcriptome". It contains only README.md, TO-DELPvsDMSO_up.xlsx and TO-DELPvsDMSO_down.xlsx, which are the adipose delphinidin-vs-DMSO DEG tables.
- **URLs read:**
  - https://pmc.ncbi.nlm.nih.gov/articles/PMC10964747/
  - https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.h44j0zpsj
  - https://datadryad.org/api/v2/versions/276616/files
  - Europe PMC REST core record

### 1. Question
Does PS in drinking water (0.5 and 5 um) disrupt the gut-liver-adipose axis in chow-fed mice? A secondary question is whether delphinidin blunts PS-induced weight gain (Abstract; Intro).

### 2. Design and n

**Main study**
- Groups: water, 0.5 um PS, 5 um PS.
- Allocation: "randomly selected cages (2 cages, n = 10 mice per treatment)". Mice were housed 5 per cage. **Treatment is therefore confounded with cage** (2 cages per arm). Source: Methods "Animals and treatment".

**Liver mRNA-seq**
- Fig 6 legend: "Cut-off values are q-value < 0.05. n = 10". Fig 7 legend: "n = 10".
- That this means 10 libraries per group, one per mouse, is implied but **UNVERIFIED**.

**Phenotype n by assay**
- Flow cytometry: n = 8-10 (Fig 1).
- Adipokines: n = 10 (Fig 2).
- Western blot: n = 4 (Fig 3).
- Adipose qPCR (Table 1): eWAT n = 6-9; PVAT n = 4-5.
- Liver histology, TG and cholesterol (Fig 5): n not stated.

**Delphinidin sub-study**
- PS-water mice received i.p. delphinidin 20 mg/kg or vehicle (20% DMSO), twice weekly for 4 weeks.
- Groups (Fig 4): U n = 4, Veh n = 8, Del n = 8.
- There was no delphinidin-only arm without PS.

**Blinding:** not reported.

### Animals and exposure
| Field | Value |
|---|---|
| Species, strain, sex, age | Male C57BL/6 (Abstract says C57BL/6J), Jackson Laboratories, 12 weeks old at purchase, 1-week acclimation |
| Polymer | Pristine PS beads, Creative Diagnostics |
| Size | "0.5 um" (0.4-0.6) and "5 um" (4.5-4.9); supplier ranges, no in-house characterisation |
| Dose | 1.0 ug/ml in drinking water, ad libitum |
| Vehicle | Filtered, autoclaved water. Stock sonicated for 30 min; water renewed every other day; bottles agitated on the other days |
| Duration | 13 weeks |
| Diet | Chow (Discussion: "fed a chow diet"); brand and composition not reported |
| Sampling | Pentobarbital euthanasia, cardiac puncture. **Fasting state and time of day: UNVERIFIED (not reported)** |
| Housing | 5 per cage, AAALAC facility. Light cycle not reported |

### 3. Datasets
- **Generated:**
  - Liver mRNA-seq (UofL Center for Genetics and Molecular Medicine).
  - eWAT RNA-seq (Novogene).
  - Bile-acid LC-MS.
- **Public:**
  - Dryad holds only the adipose delphinidin DEG spreadsheets.
  - **The liver RNA-seq raw data and counts are NOT public.** No GEO, SRA or BioProject accession appears in the paper, and none is in the Dryad record.
- **Re-used:** none.

### 4. Statistical models
- **Liver DE:** "q-value <0.05 with a log2 fold change = 0" (Methods "RNAseq analysis").
  - The DE software (DESeq2, edgeR or other) is **not named**, and the FDR method is not named. **UNVERIFIED.**
  - Platform, read depth and aligner are **not reported**.
- **Liver enrichment:** IPA core analysis per comparison. Pathways were called overrepresented if -log(p) > 1.3, i.e. unadjusted p < 0.05. No activation z-scores are reported in the text.
- **Adipose RNA-seq:** KEGG, GO and Reactome analyses by Novogene. DEG thresholds are not in the main text.
- **Phenotypes:** GraphPad Prism 9.5.1. Normality was checked with D'Agostino-Pearson and variance with Bartlett, then:
  - one-way ANOVA with Tukey, or
  - Kruskal-Wallis with Dunn, or
  - Brown-Forsythe/Welch ANOVA with Dunnett.
  - t-tests for two-group comparisons.
- **Bile acids:** one-way ANOVA with Tukey.
- **Not done:** no factorial model, no interaction term (there is only one factor), no WGCNA, PPI or hub-gene analysis, and no ML. Cage is not modelled.

### 5. Suitability of the models
- Pairwise comparisons against a single control fit a one-factor, three-arm design.
- Weaknesses:
  - Cage confounding (2 cages per arm) is ignored, so between-mouse variance is underestimated.
  - IPA pathways are called at unadjusted p.
  - The unnamed DE tool and missing FDR method make the DEG counts irreproducible.
  - With n = 10 per arm the study is adequately sized for bulk DE in principle.

### 6. Inputs and outputs
- **Inputs:** whole-liver bulk RNA, n = 10 per group. Each PS size was compared with water.
- **Outputs:** DEG lists (Supp Tables 2-5, not read), IPA canonical pathways, and phenotypes (body composition, hepatic lipids, bile acids, GSH).
- There is no predictive target.

### 7. Evaluation
- Only significance thresholds were used.
- Liver DEGs were **not validated by qPCR**.
- Adipose Akr1b7 and Fam3b were confirmed by qPCR, and Fam3b also by Western blot (Results, delphinidin section).

### 8. Validation
- **Internal only.** There is no independent cohort and no external dataset.
- The phenotypic NULLs (histology, TG, ALT/AST, GSH/GSSG)