# sci_rep_2025_ps_hubgenes

## Citation
Li Q, Zhu K, Huang L, Niu X, Li L, Gao L, Xia Z. "Polystyrene microplastics induce liver fibrosis and lipid deposition in mice through three hub genes revealed by the RNA-seq." Scientific Reports 2025; 15(1): 2583. DOI 10.1038/s41598-025-86810-5. PMID 39833454. PMCID PMC11747203. https://europepmc.org/articles/PMC11747203 (Q. Li and K. Zhu contributed equally; corresponding: L. Gao, Z. Xia.) Received 2023 Oct 6; accepted 2025 Jan 14. Licence CC BY-NC-ND. Author list, volume/article no. from Europe PMC metadata (meta_39833454.json).

## Access level read
Full text (Europe PMC full-text, PMC11747203.txt converted from XML): Abstract, Methods, Results, figure legends, Discussion, Conclusion, Data availability. Supplementary Material 1 (Suppl. Table 1 primers, Suppl. Figs 1-4: oxidative-stress assays, volcano plots/heatmap, GO CC/MF, PPI) NOT read (not in XML text).

## Biological question
Molecular mechanism of liver fibrosis after chronic oral exposure to polystyrene nano-/microplastics (PS-NMPs) of three sizes (80 nm, 0.5 µm, 5 µm) in mice; identification of hub genes from liver RNA-seq. Intro: "the objective of this study is to investigate the impact of PS microplastics on liver injury."

## Experimental design (groups, controls, randomisation/blinding if stated)
- 4 groups: control H2O, 80 nm, 0.5 µm, 5 µm PS-NMPs (Methods "Animal model").
- Quote: "Forty mice were equally and randomly divided into four groups, including the control group (H2O) (n = 10)..." (n = 10 per group).
- Control: "The control mice (H2O) were exposed to normal water without PS-NMPs." (no vehicle/surfactant control described; particle suspension medium not stated beyond sterilized water).
- Randomisation: stated ("equally and randomly divided"); method of randomisation not reported.
- Blinding: not reported (searched Methods/histology section; histology "ten random fields were selected at × 400 magnification and analyzed by ImageJ", no blinding statement).
- No diet/metabolic-stress factor (standard housing; diet not specified). Single dose, single time point (12 weeks).
- Cage allocation / housing density: not reported.

## Species/strain/sex/age; tissue(s)
- Mouse, ICR, male, 5 weeks old at start: "Five-week-old ICR male mice were used in the present experiments." (Methods "Animal model").
- DISCREPANCY (metadata only): Europe PMC/MEDLINE MeSH indexing lists "Mice, Inbred C57BL", whereas the paper text says ICR (outbred). Paper text is authoritative; MeSH tag appears erroneous.
- Housing: SPF, 12 h light/dark, 21-23 °C, free access to diet and water. Time of sacrifice (clock time) not reported.
- Tissues: liver (histology, IHC, qPCR, WB, TEM, RNA-seq), serum (ALT/AST). Note: oxidative-stress Methods say measured "in ovaries" (see Methodological points) - in male mice; presumably a copy error for liver.

## Biological replicates per group (per assay - RNA-seq n may differ from phenotype n)
- Allocated: n = 10/group, 40 mice total (Methods "Animal model").
- RNA-seq: "Four independent biological replicate samples were sequenced for each group." (Methods "RNA-seq and data analysis") -> 4 groups x 4 = 16 libraries. Whether each sample is one mouse or pooled is not stated explicitly ("independent biological replicate samples").
- Liver/body weight, ALT, AST: "n = 5–8" (Fig. 1 legend). H&E quantification n = 8 (Fig. 1F).
- MTS (Masson) quantification: n not stated in Fig. 2 legend (2A/B). α-SMA qPCR n = 6 (Fig. 2C); α-SMA WB n = 5 (Fig. 2D); α-SMA IHC n = 8 (Fig. 2F).
- CD68 IHC n = 8; inflammatory qPCR n = 6; IL-1β and IL-6 WB n = 5 (Fig. 3 legend).
- Lipid gene qPCR n = 5–6 (Fig. 6 legend). Oil Red O quantification (Fig. 4B): n not reported. TEM: representative only.
- Hub-gene qPCR validation "by increasing the samples (n = 5–6)" (Fig. 8E–H legend).
- Oxidative stress (GPX/GSH, MDA, SOD): Suppl. Fig. 1, n NOT VERIFIED (supplement not read).
- Note: 10 allocated vs n = 5–8 analysed; attrition/exclusion criteria not reported.

## Polymer; particle size (nominal and measured); supplier; surface/fluorescent label; characterisation; endotoxin
- Polymer: polystyrene (PS).
- Nominal sizes: 80 nm, 0.5 µm, 5 µm. "PS-NMPs of three different diameters (80 nm, 0.5 µm, and 5 µm) were purchased from Baseline ChromTech Research Center (Tianjin, China)" (Methods "Materials and reagents").
- Measured size (DLS/SEM/TEM of particles): not reported. Zeta potential: not reported.
- Surface functionalisation (plain/carboxyl/amine) and fluorescent label: not reported.
- Particle tracking / tissue accumulation in liver: not measured in this study (Discussion cites refs 3,12,13 for accumulation).
- Endotoxin testing: not reported. Suspension additives (surfactant, azide) in stock: not reported.

## Dose (as reported, and units); route; vehicle; frequency; duration
- Dose: 1 mg/L in drinking water: "PS-NMPs were diluted in sterilized water to a final concentration of 1 mg/L." (Methods "Animal model").
- Route: oral, ad libitum drinking water: "PS-NMPs-exposed mice directly drank the diluted PS-NMPs water for 12 weeks."
- Vehicle: sterilized water. Frequency of water renewal: not reported. Water intake / ingested dose (mg/kg/day): not reported.
- Duration: 12 weeks. Authors' rationale: 1 mg/L "is believed to be more representative of the actual exposure levels of NMPs in the environment" (Discussion para 2).

## Diet / metabolic-stress model (exact diet, % kcal fat, supplier, control diet)
- No metabolic-stress diet. Diet composition/supplier not reported ("had free access to their diets and water", Methods). All groups on the same (unspecified) standard diet. No HFD/MASH arm.
- Introduction cites prior work that "NMPs exposure exacerbated liver fibrosis in HFD-fed mice by promoting inflammation and oxidative stress" (ref 14) - background only, not tested here.

## Data types (RNA-seq platform/depth, proteomics, metabolomics, lipidomics, scRNA, spatial, histology, serum)
- Bulk liver RNA-seq: "Transcriptome sequencing of RNA was completed by the Beijing Genomics Institution (BGI)." Platform/instrument, read length, library type, sequencing depth: not reported in main text.
- qRT-PCR (Bio-Rad CFX; normalised to β-actin; primers in Suppl. Table 1), Western blot (α-SMA, IL-1β, IL-6; β-actin loading), IHC (α-SMA, CD68), histology (H&E, Masson trichrome, Oil Red O), TEM (HT7800, Hitachi), serum ALT/AST (Jiancheng kits), oxidative-stress kits (GPX/GSH, MDA, SOD; Jiancheng).
- No proteomics, metabolomics, lipidomics, liver TG/cholesterol biochemistry, scRNA or spatial data.

## Statistics & computational methods (DE tool, thresholds, multiple testing, pathway tool, models; tests for phenotypes)
- DE: limma in R: "Differentially expressed genes (DEGs) between the H2O and NMPs groups were analyzed using linear models for the microarray data (limma) package in R." Whether voom/trend was used for count data, read alignment/quantification tool and reference genome: not reported.
- Thresholds: "A |log2 (fold change)| > 1 and an adjusted P value < 0.05 were considered the threshold." Adjustment method (BH assumed?) not named.
- Each PS group compared to H2O separately (3 pairwise contrasts); shared DEGs by Venn; UpSetR; pheatmap.
- Enrichment: DAVID (GO BP/CC/MF, KEGG, REACTOME); background/universe not stated; enrichment FDR thresholds not stated.
- PPI: STRING v11.5, minimum interaction score 0.4; Cytoscape 3.4.0 CytoHubba, "12 algorithms"; GeneMANIA neighbour network.
- Phenotypes: "mean ± S.E.M."; SPSS 23; "one-way ANOVA with multiple comparisons. Post hoc multiple comparison LSD was applied." (Statistical analysis). LSD is uncorrected for multiplicity.

## Main findings with effect directions and statistical evidence
- Liver weight and liver/body weight ratio decreased, more with larger particles; body weight, ALT, AST unchanged (Results 1; Fig. 1A–E; n = 5–8). "remarkable decreases in liver weight and liver/body weight ratio as PS-NMPs particle size increased". VALIDATED (gravimetry, serum biochemistry).
- Histological liver injury (vacuolar degeneration, edema, irregular hepatic cords, binuclear cells, inflammatory infiltration) (Fig. 1F–G; n = 8). VALIDATED (H&E, semi-quantitative).
- Fibrosis: Masson trichrome area increased for 0.5 and 5 µm vs control (Fig. 2A–B). α-SMA mRNA unchanged (Fig. 2C, n = 6); α-SMA protein increased "with the increase in PS-NMPs particle size" (WB n = 5; IHC n = 8; Fig. 2D–F). VALIDATED (histology, WB, IHC). Conclusion states α-SMA up for "all sizes".
- Inflammation: CD68 IHC increased for 80 nm and 0.5 µm (Fig. 3A, n = 8); Il6, Il1b mRNA up for 80 nm and 0.5 µm; Tnf and Cxcl2 up for 5 µm; Cxcl12 up for 0.5 µm; Areg up for 80 nm (Fig. 3B–G, n = 6); IL-1β and IL-6 protein up (Fig. 3H–I, n = 5). VALIDATED (IHC, qPCR, WB).
- Oxidative stress NOT changed: "the levels of GPX, MDA, and SOD were not increased by exposure to PS-NMPs of any of the three diameters" (Results; Suppl. Fig. 1). VALIDATED-null (biochemical kits; supplement not read).
- Lipid deposition increased (Oil Red O, Fig. 4A–B; TEM Fig. 4C); "the PS-NMPs with the smallest diameter led to more lipid deposition." VALIDATED (histology/TEM; n and statistics for 4B not reported in legend).
- RNA-seq DEGs vs H2O: "A total of 333, 602, and 758 DEGs" for 80 nm, 0.5 µm, 5 µm (Results "Bioinformatics analyses"; Suppl. Fig. 2). Up/down split not given in main text. Enrichment: lipid metabolism BP in all three; glucose homeostasis (0.5, 5 µm); unfolded protein / ER stress response (5 µm); KEGG "metabolic pathways"; REACTOME "metabolism of lipids" (Fig. 5). OMICS-ONLY.
- 87 shared DEGs (46 up, 41 down) "shared among all four RNA-seq datasets" (sic; three contrasts) (Fig. 7A); enriched in positive regulation of lipid metabolism, negative regulation of lipid biosynthetic process / lipid storage; KEGG fatty acid elongation, biosynthesis of unsaturated fatty acids; CC mitochondrial matrix. OMICS-ONLY.
- PPI (85 nodes, 51 edges) -> CytoHubba hub genes Acot3, Abcc3, Nr1i3, Fmo2 (Fig. 7E). COMPUTATIONAL-PREDICTION.
- Hub genes up in RNA-seq (Fig. 8A–D); qPCR in expanded samples (n = 5–6): Acot3, Abcc3, Nr1i3 significantly up; Fmo2 not significant (Fig. 8E–H). VALIDATED at mRNA level only (qPCR); no functional test (no knockdown/agonist/protein) linking these genes to fibrosis. Claim "are involved in PS-NMPs exposure-induced liver fibrosis" is associative.
- Lipid-gene qPCR (Fig. 6, n = 5–6): FAO genes Acox1, Ppara, Cpt1a, Mcad up (80 nm, 0.5 µm); Ppara up (5 µm); Cpt1b up (80 nm). Synthesis: Apoc4, Fas up (0.5 µm); Pparg up (0.5, 5 µm); Scd1 up (80 nm); Acc - no change described in text. Transport: Fatp1 up (80 nm, 0.5 µm); Fatp2 up (0.5 µm). VALIDATED (qPCR).

## Pathway-theme directions
- Lipid / lipogenesis: UP. Oil Red O and TEM lipid deposition up; qPCR Fas (0.5 µm), Scd1 (80 nm), Pparg (0.5, 5 µm), Apoc4 (0.5 µm) up; Acc no change described (Fig. 6B). Srebf1, Acly, Cd36 not reported.
- PPARα / FAO: UP. Ppara up (all sizes), Acox1, Cpt1a, Mcad up (80 nm, 0.5 µm), Cpt1b up (80 nm) (Fig. 6A). Acot3 (peroxisomal thioesterase) up (Fig. 8). Cyp4a not reported in main text.
- Cholesterol synthesis: not reported.
- Inflammation: UP (size-dependent). Il6, Il1b, CD68 (80 nm, 0.5 µm); Tnf, Cxcl2 (5 µm); IL-1β/IL-6 protein up (Fig. 3). NF-κB, NLRP3 not measured.
- Oxidative stress: NO CHANGE (GPX/GSH, MDA, SOD; Suppl. Fig. 1). Nrf2/Hmox1 not reported.
- Mitochondria: not directly assayed; shared DEGs enriched in CC "mitochondrial matrix" (Suppl. Fig. 4) - OMICS-ONLY.
- ER stress: unfolded protein / ER stress response BP enriched for 5 µm (Fig. 5A) - OMICS-ONLY.
- Fibrosis: UP for 0.5 and 5 µm (Masson); α-SMA protein up, mRNA unchanged (Fig. 2).
- Xenobiotic / nuclear receptors: Nr1i3 (CAR), Abcc3 up (Fig. 8).
- Circadian: not reported (no clock genes or sacrifice time mentioned).
- MP x diet interaction: not applicable (no diet factor). Only size comparisons, each vs H2O; no formal size trend test.

## Author-acknowledged limitations (quotes + location)
- No explicit limitations section. Closest statements:
  - Hub genes: "the role of ACOT3 in chronic exposure to NMPs-induced liver fibrosis still needs further study." (Discussion, ACOT3 para)
  - "its mechanism in long-term exposure to PS-NMPs needs to be deeply researched." (Discussion, ABCC3 para)
  - "Whether NR1I3 promotes liver fibrosis ... by regulating lipid metabolism remains to be further explored." (Discussion, NR1I3 para)
  - Prior studies "have typically used higher concentrations of NMPs than what animals are typically exposed to" - framed as a strength of their 1 mg/L dose (Discussion para 2).
  - Oxidative stress null is "Contrary to previous studies" (Discussion para 3).

## Methodological points to verify independently
- Particle characterisation: are measured size distribution, surface chemistry, label, endotoxin and stock additives reported anywhere (supplement)? None in main text.
- Is each RNA-seq replicate one mouse? Were the 4/group sequenced mice a random subset of the 10? Were all 16 libraries prepared/sequenced in one BGI batch?
- limma described as "for the microarray data" on RNA-seq: was voom/limma-trend applied to counts or to FPKM? Alignment/quantification pipeline and depth not reported.
- Phenotype n (5–8) < allocated n (10): reasons for exclusion/attrition? Different subsets per assay?
- Post hoc LSD without multiplicity correction; many qPCR targets x 3 sizes - check robustness of individual significant calls.
- Internal inconsistencies to check: oxidative-stress Methods state "in ovaries" (male mice); heading lists GPX but Methods list GSH; "tubulointerstitial fibrosis" (kidney term) for liver; sections "(4 mm)"; IL-6 WB antibody not listed in Methods; Discussion says 0.5/5 µm "resulted in a reduction in ... liver fibrosis" (likely means induction); "shared among all four RNA-seq datasets" for 3 contrasts; Results say liver-weight decrease with increasing size while Discussion says 80 nm caused more lipid deposition.
- Blinding of histology/IHC scoring not stated; Oil Red O n and test not given (Fig. 4B).
- Cage effects / housing per group, sacrifice time of day, fasting status: not reported.
- Hub-gene validation is mRNA-only; no functional or protein evidence for Acot3/Abcc3/Nr1i3 causality.
- Single sex (male), outbred ICR strain, single dose and time point.
- Metadata: MeSH indexes "Mice, Inbred C57BL" while text says ICR.

## Data availability statement (verbatim) and all accessions/code links mentioned
- Verbatim (Data availability): "The datasets analyzed during the current study are available in the GEO repository, and the accession number to datasets is GSE245069."
- Accessions mentioned anywhere: GEO GSE245069 (Data availability only). No SRA/BioProject, code repository or other links. Tools cited by URL: geneontology.org, kegg.jp, DAVID, STRING.
- Cross-check vs coordinator note (GSE245069: 16 samples, H2O/80 nm/0.5 µm/5 µm, 4 each, ICR male, 1 mg/L drinking water, 12 weeks): AGREES with paper text (4 groups; "Four independent biological replicate samples were sequenced for each group"; ICR male; 1 mg/L; 12 weeks). No discrepancy found between paper and GEO description as relayed. Only discrepancy found is MeSH "C57BL" vs text "ICR". GEO record itself not checked by me.

## Relevance to the team
Healthy-diet, chronic low-dose oral PS (not PE) study with public liver bulk RNA-seq (GSE245069, n = 4/group); reports PS-associated up-regulation of PPARα/FAO and lipogenic genes plus lipid deposition without oxidative stress. Provides a non-MASH, different-polymer liver transcriptome against which PE direction-of-effect can be compared; no diet interaction tested.
