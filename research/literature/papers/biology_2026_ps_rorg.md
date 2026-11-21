# biology_2026_ps_rorg

## Citation
Ding P, Muniyappan M, Zhu C, Li C, Adam SY, Wang Y, Tyasi TL, Yuan P, Hu P, Liu H, Cai D. "Polystyrene Nanoplastic Exposure Causes Reprogramming of Anti-Oxidative Genes Hmox1 and Sod3 by Inhibiting Nuclear Receptor RORγ in the Mouse Liver." Biology (Basel) 2026; 15(2): 135. DOI 10.3390/biology15020135. PMID 41594870. PMCID PMC12837483. https://europepmc.org/articles/PMC12837483 . Licence CC BY. Author list/volume/issue/article no. from Europe PMC metadata (meta_41594870.json).

## Access level read
Full text (Europe PMC XML / PMC12837483.txt): Simple Summary, Abstract, Methods 2.1-2.14, Results 3.1-3.5, all figure legends (Figs 1-6) and Tables 1-2 (parsed from XML), Discussion, Conclusions, Data Availability. Supplementary Figs S1 (GSEA) and S2 (WB images) NOT read.

## Biological question
Whether 80 nm PS nanoplastics alter hepatic transcription of the antioxidant genes Hmox1 and Sod3 via the nuclear receptor RORγ and histone modifications. Intro 1: "We aimed to examine whether PS-NPs alter their expression through mechanisms involving RORγ and associated epigenetic modifications."

## Experimental design (groups, controls, randomisation/blinding if stated)
- Two groups: vehicle (PBS) vs PS-NPs (80 nm). Quote (Methods 2.2): "twelve mice were randomly assigned to two groups (n = 6 per group): a vehicle control group, which received phosphate-buffered saline (PBS), and a PS-NP-exposed group".
- Randomisation: stated; method not reported. Blinding: not reported (searched Methods 2.2-2.14).
- Single dose, single time point (28 days). No diet/metabolic-stress factor. Cage allocation/housing: not reported.

## Species/strain/sex/age; tissue(s)
- Mouse C57BL/6, male, 6 weeks old at purchase, 7-day acclimatisation: "Six-week-old male C57BL/6 mice were purchased from Yangzhou University's Model Animal Center" (Methods 2.2).
- Housing conditions (light cycle, temperature), time of day of sacrifice, fasting status: not reported.
- Tissues: liver (histology, enzymes, ROS, mitochondrial complexes, ATP, qPCR, WB, RNA-seq, ChIP), blood (CBC) and serum (ALT/AST/ALP).

## Biological replicates per group (per assay)
- n = 6 per group throughout as stated: Fig. 1 legend "The data is presented as means ± SD (n = 6)"; same in Fig. 2, Fig. 4, Fig. 5, Fig. 6 and Table 2 legends; Abstract "(n = 6, treatment (PS-NPs) vs. vehicle group (Veh))".
- RNA-seq: number of libraries per group NOT reported (Methods 2.12 describes only enrichment analysis; no library/sequencing section). Whether RNA-seq used all 6 mice per group, a subset, or pooled RNA is not stated.
- ChIP-seq (Fig. 5A,B): n and methods for ChIP-seq are not reported (only ChIP-qPCR is described in Methods 2.13).
- Mitochondrial complexes, ROS, ATP, caspase 3/7, DNA content: n = 6 (figure legends).
- DNA content qPCR run "in triplicate" (technical replicates, Methods 2.5).

## Polymer; particle size (nominal and measured); supplier; surface/fluorescent label; characterisation; endotoxin
- Polystyrene nanoplastics, nominal 80 nm, stock 10 mg/mL: "The Baseline ChromTech Research Center (Tianjin, China) provided PS-NPs (80 nm, 10 mg/mL)" (Methods 2.1). Same supplier as the Sci Rep 2025 PS study.
- Measured size (DLS/TEM/SEM), polydispersity, zeta potential: not reported (no characterisation section; searched text for DLS/zeta/characterisation).
- Surface functionalisation / fluorescent label: not reported. Particle aggregation in PBS, sonication/dispersion protocol: not reported.
- Endotoxin testing: not reported. Stock additives (surfactant/azide): not reported.
- Tissue uptake/biodistribution of particles: not measured here (Discussion cites refs for accumulation).

## Dose (as reported, and units); route; vehicle; frequency; duration
- "For 28 days, mice in the PS-NP group were administered 200 µL of a PS-NP suspension (10 mg/mL) via oral gavage once daily, delivering a dose of 1 mg per mouse per day" (Methods 2.2).
- Route: oral gavage, once daily; vehicle: PBS 200 µL, same route; duration 28 days.
- Approx. body-weight-normalised dose (~40-50 mg/kg/day for a 20-25 g mouse) NOT stated by the authors.
- INCONSISTENCY: Discussion 4 states the study "examined the toxicological impacts of PS-NPs ingested via drinking water", contradicting the gavage protocol in Methods 2.2.

## Diet / metabolic-stress model (exact diet, % kcal fat, supplier, control diet)
- Standard chow only, both groups: "The standard diet was purchased from Jiangsu Xietong Pharmaceutical Bio-Engineering Co., Ltd." (Methods 2.1); composition in Table 1 (per kg: crude protein 200 g, crude fat 40 g, coarse fibre 50 g, ash 80 g, calcium 10.18 g, total phosphorus 12 g, lysine 13.2 g, Met+Cys 7.8 g, vitamins/minerals). % kcal fat not given (4% crude fat by weight); no HFD/MASH model; no high-fat or choline-deficient arm.
- Feed intake measured: PS-NP group -8% average daily feed intake (Results 3.1, Fig. 1B).

## Data types (RNA-seq platform/depth, proteomics, histology, serum, epigenomics)
- Bulk liver RNA-seq: platform, library prep, read length, depth, aligner and quantification NOT reported anywhere (Methods 2.12 covers only GSEA 4.1.0, Metascape, DAVID and bioinformatics.com.cn plots).
- ChIP-seq for RORγ (Fig. 5A,B signal tracks) - no ChIP-seq methods, peak caller or data source stated; possibly re-used/public data (not specified).
- ChIP-qPCR (Methods 2.13; Covaris E220 sonication; antibodies for RORγ, P300, SRC1, Pol II, Ser5-Pol II, Ser2-Pol II, H3K9ac, H3K18ac, H3K27ac, H3K4me1, H3K4me3, H3K27me3; "ChIP-re-ChIP" mentioned in Discussion but not in Methods).
- qRT-PCR (2^-ΔΔCt; reference gene not named), Western blot (GAPDH, RORγ, Sod3, Hmox1; antibody sources not given), H&E histology + damage score, caspase 3/7 assay, DNA-content qPCR, SOD/CAT/GSH/MDA/ROS kits, mitochondrial Complex I/III/V and ATP kits, CBC (Sysmex XN-9000), serum ALT/AST/ALP.
- No proteomics/metabolomics/lipidomics; no liver TG/cholesterol measurement; no fibrosis staining (no Masson/Sirius red).

## Statistics & computational methods
- Phenotypes/qPCR/ChIP-qPCR: "All data analysis was conducted using GraphPad Prism (version 9.0), following normality assessments and evaluations of variance homogeneity. The results were given as means ± standard deviation (SD). Statistical analysis was conducted using two-tailed non-parametric t tests. Significance was defined as p < 0.05" (Methods 2.14). Note the internal tension: normality/variance checks reported, yet "non-parametric t tests" (test not named, e.g. Mann-Whitney) applied; no correction for the many endpoints tested.
- DEG threshold: "689 differentially expressed genes (DEGs, |Log2 (fold change) | > 1, p < 0.05)" (Results 3.3) - raw p, NO multiple-testing correction stated; DE software (DESeq2/edgeR/limma) not named.
- Enrichment: GSEA 4.1.0; Metascape; DAVID; plots via bioinformatics.com.cn (accessed 23 May 2025). Gene-set versions, backgrounds and FDR thresholds not reported.
- Table 2 uses superscript letters ("Means with different superscripts in the same row differ significantly") for a 2-group comparison.

## Main findings with effect directions and statistical evidence
- Body weight -10%, average daily feed intake -8%, liver weight -19% in PS-NP vs Veh (Results 3.1; Fig. 1A-C, n = 6). VALIDATED (gravimetry). Note confounding by reduced food intake is not addressed.
- Liver histology: immune-cell infiltration, hepatocyte vacuolisation, displaced/pyknotic nuclei, increased sinusoidal space; damage score +22.7% (Fig. 1D,E). VALIDATED (H&E; scoring blinded? not stated).
- Apoptosis/DNA: caspase 3/7 activity +38%; DNA content -29.4% (Fig. 1F,G). VALIDATED (enzymatic assay, qPCR-based DNA assay).
- Serum liver enzymes UP: ALT +39.4%, AST +31.6%, ALP +25% (Fig. 1H-J). VALIDATED (kits). (Results text says "liver ALT, AST, and ALP levels" though assays were on serum.)
- Haematology (Table 2, n = 6): WBC 9.345 ± 0.85 -> 15.51 ± 0.72 (+17.77% per text), monocytes 5.57 ± 0.38 -> 9.32 ± 0.31 %, platelets 439.50 ± 10.31 -> 357.00 ± 5.72 (10^9/L, -18.48%); other CBC indices unchanged. VALIDATED (analyser). DISCREPANCY: text says WBCs "increased by 17.77%" but the tabulated values differ by ~66%; monocyte % change ~67% vs "30%" in text.
- Oxidative stress UP: SOD -21.42%, CAT -29.12%, GSH -23.02%, MDA +35.15%, ROS +45.87% (Results 3.2; Fig. 2A-E). VALIDATED (biochemical/fluorimetric kits). Note Fig. 2 legend lists panel C as MDA and D as GSH while Results text order differs; legend states "*** p < 0.05".
- Mitochondria DOWN: ATP -23.45%, Complex I -31.23%, Complex III -29.67%, Complex V -31.41% (Fig. 2F-I). VALIDATED (kits).
- Transcriptome: 689 DEGs (389 up, 300 down) PS-NP vs Veh (Fig. 3A). GO/KEGG enrichment in inflammation and apoptosis pathways; GSEA showed downregulated hallmarks of cell cycle, apoptosis, oxidative stress, lipid oxidation and inflammatory response (Fig. S1A-E; Figs 3B-F, 4A-E). OMICS-ONLY. Note internal tension: KEGG/GO described as enriched in inflammation/apoptosis while GSEA reports these same pathways as downregulated.
- Hmox1 and Sod3 mRNA DOWN by 45.34% and 48.62%; RORγ mRNA -19.34%; all three proteins decreased by immunoblot (Results 3.3; Fig. 4F-I, Fig. S2). VALIDATED (qRT-PCR + western blot).
- RORγ occupancy at Hmox1 and Sod3 loci reduced; co-activators P300 and SRC1 reduced; Pol II, Ser5-Pol II and Ser2-Pol II reduced (Results 3.4; Fig. 5C-H, n = 6). VALIDATED (ChIP-qPCR; ChIP-seq tracks Fig. 5A,B shown without methods). No RORγ knockout/knockdown or agonist/antagonist rescue was performed, so RORγ causality is INFERRED, not demonstrated.
- Chromatin modifiers: Hdac11 and Hdac5 mRNA DOWN; Kdm2a, Kdm3a, Kdm5c, Kdm6a, Kdm6b UP (transcriptomics + qRT-PCR validation; Fig. 6A-D). VALIDATED (qPCR) for the listed genes. (Discussion writes "Kdm5a" where Results say Kdm5c.)
- Histone marks at Hmox1/Sod3 loci: H3K9ac and H3K18ac DOWN; H3K4me3 and H3K27me3 UP; H3K27ac and H3K4me1 unchanged (Results 3.5; Fig. 6E-J). VALIDATED (ChIP-qPCR). Results text mislabels H3K4me3/H3K27me3 as "related to transcriptional activation"; Discussion's closing sentence states these changes "result in enhanced expression of these genes", the opposite of the paper's own conclusion (downregulation).

## Pathway-theme directions
- Lipid metabolism / lipogenesis: not measured directly (no liver TG, no Fasn/Srebf1/Acly/Cd36/Pparg data). A "lipid oxidation pathway" gene heatmap is shown (Fig. 4D) and GSEA reports lipid oxidation hallmarks DOWN (Fig. S1) - OMICS-ONLY; individual gene directions not stated in text.
- PPARα / FAO: not measured in this study. Discussion only cites ref [44] that "80 nm NPs lower NR PPARα activation, which in turn reduces mitochondrial fatty acid β oxidation and increases liver steatosis" - background citation, not this paper's data.
- Cholesterol synthesis: not measured; mentioned only as background (RORγ-MVA/MVK pathway, Discussion).
- Inflammation: mixed/ambiguous. Histological immune infiltration and raised WBC/monocytes (UP, VALIDATED); transcriptome inflammation pathways described as enriched but GSEA hallmark "inflammatory response" DOWN (OMICS-ONLY). No Tnf/Il1b/Il6/NF-κB/NLRP3 measurements reported.
- Oxidative stress: UP (ROS +, MDA +, SOD/CAT/GSH enzyme activity -, Hmox1 and Sod3 mRNA and protein DOWN). VALIDATED. Nrf2 itself was NOT measured (discussed only).
- Mitochondria: IMPAIRED (ATP, Complex I/III/V activities all down). VALIDATED.
- Fibrosis: not assessed (no fibrosis staining or markers).
- Circadian: not measured. RORγ is a clock-linked nuclear receptor (Intro cites its role in circadian gene expression), but no clock genes, no sampling time (ZT) reported - a notable gap given RORγ's rhythmicity.
- Apoptosis/cell cycle: apoptosis UP by caspase 3/7 (VALIDATED) while GSEA apoptosis hallmark DOWN (OMICS-ONLY) - directions conflict.
- MP x diet interaction: none; no diet factor and no interaction testing. All claims are PS-NP vs vehicle pairwise.

## Author-acknowledged limitations (quotes + location)
- No limitations section. Nearest statements:
  - "While the underlying causes are intricate and multifaceted, including inflammation, oxidative stress, disturbance of lipid metabolism, and epigenetic changes, further investigation is needed to comprehend these systems." (Discussion 4, KDM/HDAC paragraph)
  - "Further research is needed to fully elucidate these pathways and validate molecular targets" (Conclusions 5)
  - On conflicting literature: "it is reported that livers exposed to PS-NPs at a dose of 0.5 µm did not exhibit significant oxidative stress ... which, again, may result from variations in the PS-NP dosages and particle sizes." (Discussion 4) - this refers to a study like the Sci Rep 2025 PS paper.
  - "The discrepancy may therefore result from variations in PS-NP dosages and administration routes." (Discussion 4, on hepatotoxicity vs no histological damage)

## Methodological points to verify independently
- RNA-seq reporting is essentially absent: which platform, how many libraries per group, depth, aligner, quantification and DE tool? Were DEGs called on n = 6 vs 6 or on pooled samples? Were the 689 DEGs corrected for multiple testing (text gives raw p < 0.05 only)?
- ChIP-seq tracks (Fig. 5A,B) appear without any ChIP-seq method, accession or peak-calling description - is this the authors' own data (is it in PRJNA1305997?) or re-used public RORγ ChIP-seq?
- "ChIP-re-ChIP" is claimed in the Discussion but no re-ChIP experiment appears in Methods or Results.
- Statistics: "two-tailed non-parametric t tests" is not a named test; with n = 6 per group a Mann-Whitney U has limited resolution (min two-sided p = 0.002). Confirm which test produced "p < 0.001" values. No adjustment across dozens of endpoints (CBC 16 parameters, many ChIP targets).
- Percentage changes in the text vs tabulated means do not match for WBC and monocytes (Table 2 vs Results 3.1) - recompute.
- Body weight and feed intake both fell ~10%/8%: are liver-weight, enzyme and transcriptome changes confounded by reduced food intake / caloric restriction? No pair-fed control.
- Gavage of 1 mg/mouse/day (~40-50 mg/kg/day) is a high, non-environmental dose; vehicle is PBS with no dispersant control.
- Particle characterisation, aggregation state in PBS and endotoxin contamination not reported.
- Time of sacrifice / circadian phase not reported, which matters for RORγ, Hmox1 and clock-controlled transcription.
- Single sex (male), one strain, one dose, one time point; no RORγ loss-of-function or agonist rescue to establish causality.
- Reference gene for qRT-PCR, antibody catalogue numbers/validation for ChIP and WB not reported.
- Internal inconsistencies: "ingested via drinking water" (Discussion) vs oral gavage (Methods); H3K4me3/H3K27me3 described as activation marks; final Discussion sentence states the mark changes "result in enhanced expression", opposite to the reported downregulation; Kdm5c vs Kdm5a; Fig. 2 legend "*** p < 0.05"; Results say "Hmox1 and Sod3 were downregulated" but one sentence concludes "Sod3 may emerge as a key factor".

## Data availability statement (verbatim) and all accessions/code links mentioned
- Verbatim (Data Availability Statement): "RNA-seq data are available in the SRA under BioProject accession number PRJNA1305997."
- Accessions mentioned anywhere: BioProject PRJNA1305997 (SRA) only. No GEO accession, no ChIP-seq accession, no code repository.
- Supplementary materials link given: https://www.mdpi.com/article/10.3390/biology15020135/s1 (Figures S1, S2).
- Tool URLs cited: metascape.org, davidbioinformatics.nih.gov, bioinformatics.com.cn (accessed 23 May 2025); GSEA 4.1.0.
- Ethics: Yangzhou University SYXK [Su] 2023-0089, approved 26 February 2024.

## Relevance to the team
A short (28-day), high-dose gavage PS-nanoplastic study on chow-fed male C57BL/6 with liver bulk RNA-seq (SRA PRJNA1305997) reporting oxidative stress, mitochondrial impairment and transcriptional repression of Hmox1/Sod3 via reduced RORγ occupancy - a mechanism-level, chromatin-based readout absent from the team's PE x MASH dataset. It offers a contrasting polymer/route/dose and, since RORγ is clock-linked and no sampling time is reported, it is a reminder to check circadian covariates; it tests no diet interaction.
