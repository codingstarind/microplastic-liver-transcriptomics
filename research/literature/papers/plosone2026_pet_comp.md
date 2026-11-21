# plosone2026_pet_comp

## Citation
Integrated computational analysis identifies FABP4, PTGS2, and HPGD as key molecular targets linking PET microplastic exposure to metabolic dysfunction-associated steatotic liver disease. Zhang Y, Yu Y, Ge B, Gong D, Zheng L, Wang Y, Chen P (corresponding). PLOS One 2026; 21(7): e0354607 (published 2026-07-24; received 2026 Jan 14, accepted 2026 Jul 7). DOI 10.1371/journal.pone.0354607; PMID 42497141; PMCID PMC13399308. https://doi.org/10.1371/journal.pone.0354607 . Licence: CC BY.

## Access level read
Full text XML (Europe PMC, PMC13399308), parsed to plain text. Read: abstract, Introduction, Methods 2.1–2.12, Results 3.1–3.8, figure legends, Discussion (incl. limitations), Conclusion, Data Availability, Funding. Supplementary S1 Table (DOCX) NOT read (not in XML).

## Biological question
Entirely in silico: which predicted human protein targets of PET (as a small-molecule SMILES) overlap with MASLD-associated genes from public human liver transcriptome data, and do they dock/stay bound in molecular-dynamics simulation? Abstract Background: "molecular targets and mechanisms in MASLD are not well defined." **No PET exposure data of any kind were generated or analysed.**

## Experimental design (groups, controls, randomisation/blinding if stated)
- No wet-lab experiment, no animals, no PET-exposed samples. Ethics: "Not applicable." (2.12)
- Arm 1 (chemical target prediction, 2.1): PET structure/SMILES from PubChem; targets from ChEMBL, PharmMapper, SwissTargetPrediction, merged, human only → "229 non-redundant potential target proteins" (3.1, Fig 2B).
  Note: the SMILES given, "CC(=O)C1=CC=C(C=C1)C(=O)OCCOC", is a small monomer-like fragment, not a polymer or particle.
- Arm 2 (disease transcriptome, 2.2–2.4): GEO GSE37031, GSE63067, GSE89632, "comprising a total of 38 control and 58 NAFLD samples" (2.2); quantile normalisation + ComBat (sva); limma DE; WGCNA.
- Overlap of PET targets with MASLD genes → 19 genes (3.4, Fig 5A); STRING PPI + cytoHubba/MCODE → hubs FABP4, PTGS2, HPGD (3.5).
- CIBERSORT (LM22) immune deconvolution, MASLD vs control (2.7, 3.6).
- Docking (AutoDock Vina) + 100 ns MD (GROMACS) of PET with FABP4, HPGD, PTGS2 (2.8–2.9, 3.7–3.8).
- CTD comparison with CCl4, APAP, BPA target genes (2.10; S1 Table).
- Randomisation/blinding: not applicable / not reported.

## Species/strain/sex/age; tissue(s)
Species of the GEO datasets: human (target prediction restricted to "Homo sapiens entries", 2.1). Donor sex/age: not reported. Tissue: "liver tissue" (3.6). Platform type: the abstract calls them "three GEO microarray cohorts" while 2.7 refers to deconvolving "bulk RNA-seq profiles" (internal inconsistency, see points to verify).

## Biological replicates per group
No experimental replicates (no wet lab). Human sample counts from public data: "a total of 38 control and 58 NAFLD samples" pooled across GSE37031, GSE63067 and GSE89632 (Methods 2.2). Per-dataset breakdown is not restated in the text. Donor sex, age, batch and platform differences are handled only by quantile normalisation + ComBat.

## Polymer; particle size; supplier; label; characterisation; endotoxin
Not applicable, and **PET is not modelled as a particle**. PET enters the analysis only as a small-molecule structure/SMILES from PubChem (Methods 2.1). No particle size, shape, dose, supplier or characterisation exists anywhere in the study. This is the single most important caveat: the work models a soluble monomer-like PET fragment binding to proteins, not microplastic particles.

## Dose; route; vehicle; frequency; duration
Not applicable — no exposure of any kind. No dose-response is modelled.

## Diet / metabolic-stress model
Not applicable. The "disease" axis is human MASLD/NAFLD vs control liver from public cohorts; no diet information is used.

## Data types
Public human liver expression data (three GEO cohorts); in-silico target prediction (ChEMBL, PharmMapper, SwissTargetPrediction); PPI (STRING); immune deconvolution (CIBERSORT, LM22); molecular docking (AutoDock Vina); molecular dynamics (GROMACS, 100 ns); toxicogenomic comparison (CTD: CCl4, APAP, BPA). No new molecular data.

## Statistics & computational methods
- Harmonisation: quantile normalisation per dataset, then ComBat (sva) batch correction; PCA before/after (Methods 2.2).
- DE: limma, "Linear models were fitted with disease status (MASLD vs. control) as the primary covariate", empirical Bayes. **Threshold stated in Methods 2.3: "|log2 fold change| >= 1.5 and adjusted p-value (FDR) < 0.05, with Benjamini-Hochberg correction".**
- WGCNA: pickSoftThreshold to scale-free R² > 0.85; TOM; minModuleSize = 30, mergeCutHeight = 0.25; module eigengene–trait Pearson correlation (|R| > 0.5, P < 0.05).
- PPI: STRING confidence ≥ 0.4; cytoHubba with seven named algorithms (Degree, Closeness, Betweenness, Eigenvector, EPC, MNC, MCC) — the text calls them "Six topological algorithms" while listing seven; MCODE for sub-networks.
- Docking: AutoDock Vina; "free energies <= -5.0 kcal/mol were considered indicative of favorable interactions"; PyMOL inspection.
- MD: 100 ns all-atom, GROMACS 2022/2022.3, TIP3P water, 0.15 M NaCl.
- **No machine-learning classifier is used** (no LASSO, random forest, SVM or ROC/AUC analysis anywhere in the text).

## Main findings with effect directions and statistical evidence
Every finding below is COMPUTATIONAL-PREDICTION. There is no qPCR, western blot, cell or animal validation anywhere in the paper (searched: no "qPCR", no "western").
- PET target prediction: "229 non-redundant potential target proteins" (Results 3.1).
- MASLD DEGs: Results state "The threshold is set to |log2FC| >= 2 and p < 0.05. A total of 357 differential genes were obtained, of which 153 were up-regulated and 204 were down-regulated" (Fig 3G).
- WGCNA intersection: "the WGCNA module set included 736 genes, while 265 DEGs were identified through differential expression analysis, and their overlap yielded 95 intersecting genes" (Results 3.4, Fig 4G).
- PET-target × MASLD overlap: 19 genes, "enriched in pathways related to lipid metabolism, fatty acid degradation, glycolysis/gluconeogenesis, and chemical carcinogenesis" (Abstract; Fig 5A).
- Hub genes: "Network topology consistently highlighted FABP4, PTGS2, and HPGD as central hub genes" (Abstract; Results 3.5).
- Immune deconvolution (CIBERSORT): "increased M2 macrophages and γδ T cells, with decreased monocytes, dendritic cells, and naive B cells"; "PTGS2 and FABP4 expression showed strong correlations with innate immune cells" (Abstract; Results 3.6). **This describes MASLD vs control, not a PET effect.**
- Docking: "favorable PET binding to all three proteins (-6.3 to -6.9 kcal/mol)" (Abstract; Results 3.7).
- MD: "simulations confirmed stable complexes over 100 ns, with predominantly hydrophobic interactions" (Abstract; Results 3.8).
- Speculative mechanism in the Discussion, presented without data: "once PET microplastics accumulate in hepatocytes or Kupffer cells, they may form persistent complexes with these metabolic and inflammatory regulators, leading to long-term functional impairment".

## Pathway-theme directions
- Lipid metabolism: implicated via FABP4 (fatty-acid binding) and enrichment in lipid metabolism/fatty-acid degradation — PREDICTION only; no direction of change caused by PET is established.
- PPARα / FAO: not analysed as such (FABP4 is a PPAR-responsive gene, but no PPAR pathway test is reported).
- Cholesterol: not analysed.
- Inflammation: implicated via PTGS2 (COX-2) and HPGD (prostaglandin degradation) and the CIBERSORT shifts — these are MASLD-associated, not PET-induced.
- Oxidative stress, mitochondria, fibrosis, circadian: not analysed.
- MP × diet interaction: not applicable; no exposure factor exists.

## Author-acknowledged limitations (quotes + location)
From the Discussion ("Several limitations of our study should be acknowledged"):
1. "our analyses relied on computational prediction methods and publicly available transcriptomic datasets, which require experimental validation in cell culture and animal models."
2. "the molecular docking and dynamics simulations, while indicating favorable binding energetics, do not account for the complex in vivo environment including competitive binding, me[tabolism]…" (sentence continues beyond the extracted window).
The Introduction also concedes the targets "require rigorous experimental validation to establish their role in environmental pollutant-mediated liver disease."

## Methodological points to verify independently
1. **PET is modelled as a monomer-like small molecule** (SMILES "CC(=O)C1=CC=C(C=C1)C(=O)OCCOC"), not as a polymer or particle. Docking a soluble fragment to FABP4/PTGS2/HPGD does not model particle uptake, and the Discussion's leap to "PET microplastics accumulate in hepatocytes" is not supported by the analysis performed.
2. **DEG thresholds and counts are internally inconsistent:** Methods 2.3 says |log2FC| ≥ 1.5 with BH-adjusted p < 0.05; Results say |log2FC| ≥ 2 with p < 0.05 (unadjusted) giving 357 DEGs; Results 3.4 refers to 265 DEGs. Which analysis produced the 19 overlapping genes?
3. **Force-field conflict:** Methods 2.8 says the 100 ns MD used "CHARMM36"; Methods 2.9 says "AMBER99SB-ILDN". Both describe the same simulations.
4. "Six topological algorithms" are announced but seven are listed (cytoHubba).
5. Platform inconsistency: the abstract calls the cohorts "microarray", while Methods 2.7 refers to deconvolving "bulk RNA-seq profiles". CIBERSORT LM22 behaves differently on array vs RNA-seq input.
6. Pooling three cohorts on different platforms (GPL14877, GPL570, GPL14951) with ComBat and no cohort covariate in the limma model risks residual confounding; PCA is the only check reported.
7. No cross-validation, held-out cohort or classifier performance metric is reported, so "key targets" rest on network topology alone.
8. The CIBERSORT result characterises MASLD vs healthy liver; it carries no information about PET exposure, although the abstract presents it inside a PET-focused narrative.

## Data availability statement (verbatim) and all accessions/code links mentioned
Verbatim (Data Availability): "The gene expression datasets analyzed in this study are publicly available from the NCBI Gene Expression Omnibus (GEO) database under accession numbers GSE37031, GSE63067, and GSE89632. All original code generated for this study is openly available under the MIT License. The code can be accessed via GitHub at https://github.com/pipipeng288/Code and has been archived in Zenodo with the DOI: https://doi.org/10.5281/zenodo.18925269. The repository includes well-documented scripts for all analyses, along with a detailed README file. Any additional data or code-related inquiries can be directed to the correspond[ing author]."
Accessions/links: GEO **GSE37031**, **GSE63067**, **GSE89632**; GitHub **pipipeng288/Code**; Zenodo DOI **10.5281/zenodo.18925269**. Databases used: PubChem, ChEMBL, PharmMapper, SwissTargetPrediction, STRING, CTD (https://ctdbase.org/). Supplement: S1 Table (DOCX, not read).

## Relevance to the team
A cautionary methodological comparator rather than biological evidence. It shows the "network toxicology" genre the team could be tempted to imitate: public human MASLD cohorts plus chemical-target prediction plus docking, producing named hub genes with no exposure data and no wet-lab validation. Its code and input accessions are genuinely public (unusually good for this set), so the analysis is reproducible; but its conclusions about PET *microplastics* do not follow from particles, and its internal threshold inconsistencies are a concrete example of why the team should pre-specify and report DE thresholds exactly once.
