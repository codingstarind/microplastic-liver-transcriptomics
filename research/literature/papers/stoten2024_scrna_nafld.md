# stoten2024_scrna_nafld (ABSTRACT-ONLY record; written by coordinator from Europe PMC core metadata)

## Citation
Liu W, Li M, Guo H, Wei S, Xu W, Yan Y, Shi Y, Xu Z, Chang K, Wei G, Zhao S. "Single-cell transcriptome analysis of liver immune microenvironment changes induced by microplastics in mice with non-alcoholic fatty liver." Science of the Total Environment 2024; 912:168308. DOI 10.1016/j.scitotenv.2023.168308. PMID 37977403. **No PMCID** — not available as open-access full text in Europe PMC/PMC.

## Access level read
**ABSTRACT ONLY.** Publisher full text not retrieved. Fields marked NOT VERIFIED could not be read and must not be cited from this record.

## Biological question
How the hepatic immune microenvironment changes when microplastic exposure is combined with high-fat-diet-induced NAFLD: "further investigation is required to understand the alterations in the immune microenvironment associated with non-alcoholic fatty liver disease (NAFLD) occurrence following exposure to MPs." (Abstract)

## Experimental design (groups, controls, randomisation/blinding)
**2 × 2 layout stated in the abstract:** "mice ... were given a normal chow or high-fat diet (NCD or HFD, respectively) plus free drinking of sterile water with or without MPs" — i.e. NCD, NCD+MP, HFD, HFD+MP. Randomisation, blinding, group sizes: NOT VERIFIED.
**Whether a formal diet × MP interaction was tested: NOT VERIFIED** (the abstract reports only group comparisons, e.g. "MPs plus HFD" vs control).

## Species/strain/sex/age; tissue(s)
Mice; strain, sex and age NOT VERIFIED. Tissue: liver (histology and single-cell transcriptomics).

## Biological replicates per group
NOT VERIFIED. The abstract gives a total cell count only: "a single-cell resolution transcriptomic atlas of 43,480 cells in the mice's livers of the indicated groups". **Number of mice per group — the actual statistical unit — is not stated in the abstract.** This is the key field to recover if the full text is obtained.

## Polymer; particle size; supplier; label; characterisation; endotoxin
Polymer: not named in the abstract beyond "microplastics" ("tiny particles of microplastics (MPs) at the nano-scale level"), implying sub-micron particles, most likely polystyrene — **NOT VERIFIED**. Size, supplier, characterisation, endotoxin: NOT VERIFIED.

## Dose; route; vehicle; frequency; duration
Route: oral, drinking water ("free drinking of sterile water with or without MPs"). Dose, concentration, duration: NOT VERIFIED.

## Diet / metabolic-stress model
Normal chow diet (NCD) vs high-fat diet (HFD); exact diets, % kcal fat, supplier and duration NOT VERIFIED.

## Data types
- **Single-cell RNA-seq (scRNA-seq)** of liver, 43,480 cells total; platform and chemistry NOT VERIFIED.
- Liver histology with NAFLD activity scoring ("NAFLD activity scores", "hepatic steatosis, inflammatory cell infiltrations and ballooning degeneration").
- Cell–cell communication analysis ("potential cell-to-cell cross-talk"); tool NOT VERIFIED.

## Statistics & computational methods
NOT VERIFIED (no clustering parameters, DE method, thresholds, multiple-testing handling or pseudo-replication handling stated in the abstract).

## Main findings with effect directions and statistical evidence
All from the abstract; no p-values, effect sizes or figure locations available.
- "mice in the HFD groups had remarkably greater NAFLD activity scores than those from the NCD groups" — diet effect, histology. VALIDATED (histology) per abstract.
- "administration of MPs plus HFD further worsened the histopathological changes in the mice's liver, leading to hepatic steatosis, inflammatory cell infiltrations and ballooning degeneration" — MP added to HFD worsens histology. **Whether this was tested as an interaction or as a pairwise contrast: NOT VERIFIED.**
- "MPs exacerbated the pro-inflammatory response and influenced the stemness of hepatocytes during HFD feeding" — OMICS-ONLY (scRNA-seq).
- "treatment with MPs significantly increase the infiltration of the infiltrating liver-protecting Vsig4+ macrophages ... while remarkably decreasing the angiogenic S100A6+ macrophage subpopulation" — macrophage subpopulation shifts; OMICS-ONLY (no validation method stated in the abstract).
- "mice treated with MPs plus HFD exhibited significantly increased recruitment of CD4+ cells and heightened exhaustion of CD8+ T cells than those from the control group" — OMICS-ONLY; note the comparator here is "the control group", not the diet-matched HFD group.

## Pathway-theme directions
- Lipid metabolism / lipogenesis: steatosis worse histologically; no gene-level lipid direction in the abstract.
- PPARα / FAO: not mentioned.
- Cholesterol: not mentioned.
- Inflammation: UP (pro-inflammatory response, immune infiltration, macrophage and T-cell shifts).
- Oxidative stress: not mentioned.
- Mitochondria: not mentioned.
- Fibrosis: not mentioned.
- Circadian: not mentioned.
- MP × diet interaction: claimed verbally as worsening ("further worsened", "exacerbated ... during HFD feeding"); **formal interaction testing NOT VERIFIED**.

## Author-acknowledged limitations
NOT VERIFIED (no Discussion available).

## Methodological points to verify independently
If the full text is obtained, check first: **number of mice per group** (43,480 cells may come from very few animals, and cell-level p-values across thousands of cells from a handful of mice are pseudo-replication unless a sample-level model is used); polymer identity and size; dose and duration; whether MP-only (NCD+MP) effects were reported; and whether the scRNA-seq data were deposited in GEO/SRA/GSA.

## Data availability statement (verbatim) and accessions
NOT VERIFIED — no data availability statement could be read, and no accession appears in the abstract. Whether the 43,480-cell atlas was deposited anywhere is unknown.

## Relevance to the team
The closest published design to the team's RQ2 in a different modality: microplastics on a high-fat-diet background, profiled at single-cell resolution, reporting immune-compartment shifts (macrophage subsets, T cells) rather than hepatocyte-intrinsic changes. It is directly relevant to the interpretation that bulk PE effects could reflect **cell-composition change** rather than per-cell expression change. Because only the abstract could be read — including the missing per-group n — it should be cited as a qualitative pointer only.
