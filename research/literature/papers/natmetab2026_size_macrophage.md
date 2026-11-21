# natmetab2026_size_macrophage

## Citation
Makdissi N, Viola MF, Steinheuer LM, Blank-Stein N, Kardinal R, Musacchio F, et al. (33 authors; N.M. and M.F.V. contributed equally). "Size-dependent plastic exposure disrupts macrophage function and tissue-specific metabolism." Nature Metabolism, 2026; volume 8, issue 10, pages 2084-2101. DOI 10.1038/s42255-026-01615-8. PMID 42811166. PMCID PMC13690586. Link: https://doi.org/10.1038/s42255-026-01615-8

## Access level read
Full text (Europe PMC full-text XML, converted to text: PMC13690586.txt / .xml). Main text, all figure legends (Figs 1-4), Extended Data Fig. 1-10 legends, full Methods, Data availability, Acknowledgements read. NOT read: Supplementary Figs 1-8 and Supplementary Tables 1-7 (not in the XML) - so particle characterisation values (Supplementary Fig. 1) and DEG counts per cell type (Supplementary Table 1) are NOT VERIFIED numerically.

## Biological question
Whether chronic, orally ingested micro- and nanoplastic particles (MNPs) accumulate in liver-resident macrophages (Kupffer cells, KCs), how this reprograms KC transcription and effector function, and whether such KC changes translate into hepatocyte/systemic metabolic dysregulation in a particle-size-dependent manner, including under a high-fat-diet (HFD) challenge. "Here we identify KCs as the primary hepatic reservoir for MNPs in young male mice" (Abstract).

## Experimental design
Polystyrene (PS) particles by oral gavage; PBS gavage = control ("Ctrl") throughout. Main paradigms:
- Acute bioavailability: single oral gavage 100 mg/kg of fluorescent PS of 50, 100, 300, 500 nm; organs 24 h later (Methods "Plastic administration model"; Extended Data Fig. 1a).
- Chronic feeding: weekly oral gavage for 12 weeks, 100 mg/kg (also 10 and 1 mg/kg dose-titration arms; "4W" = 4-week arm) - 3 groups: Ctrl (PBS), NP (50 nm), MP (500 nm) (Fig. 1a; Methods).
- Washout: 12 weeks exposure then 6-month washout (Extended Data Fig. 7a).
- HFD arm: 2 weeks control diet (CD) acclimation, 8 weeks MNP on CD, then switch to HFD with 4 more weeks of MNP, then HFD alone for 6 further weeks (Fig. 2a; Methods "Chronic plastic administration with HFD feeding").
- Gold-core MNP arm: 0.8 mg/week gold-core PS-coated MP or NP by gavage for 5 weeks (Methods; Fig. 1d legend states "0.8 mg kg-1 twice per week, alternating NPs and MPs" - note the Methods and the Fig. 1d legend disagree on dose/frequency; AMBIGUOUS).
- Genetic models: double fate-mapper Tnfrsf11aCre;Rosa26LSL-YFP;Ms4a3FlpO;Rosa26FSF-tdTomato; Cxcr4CreERT2;Rosa26LSL-tdTomato; Clec4fDTR/+ (diphtheria-toxin KC depletion) with Clec4f+/+ littermate controls; Clec4fDTR/+;Ms4a3Cre/+;Rosa26LSL-tdTomato/+.
- Ex vivo: KC-hepatocyte co-culture (KCs from Ctrl or MP mice + primary hepatocytes from untreated 5-week-old female mice).
Randomisation: not reported (looked in Methods "Mice", "Plastic administration model", "Statistics"). Blinding: only for histopathology - "a trained pathologist blinded to experimental groups performed all assessments" (Methods, Immunohistochemistry section). No statement of blinding for other readouts.

## Species/strain/sex/age; tissue(s)
Mouse, C57BL/6JRcc wild type (transgenics on C57BL/6J background crossed to C57BL/6Jrcc). "Young (8-12-week-old) male mice were used for all experiments, with the exception of experiments with Cxcr4CreERT2/+; Rosa26LSL-tdTomato/+ mice, in which both sexes were used" (Methods, Mice). HFD arm used male C57BL/6JRcc 6-8 weeks old at diet start (Methods). Hepatocytes for co-culture came from 5-week-old chow-fed FEMALE C57BL/6JRcc mice (Methods, co-culture).
Tissues: liver (main), spleen, brown adipose tissue (BAT), gonadal white adipose tissue (gWAT), brain (cortex, hippocampus), blood, bone marrow, lymph nodes, Peyer's patches, small intestine/jejunal lamina propria, pancreas, serum, faeces.
Tissues with omics: liver ONLY. snRNA-seq = liver nuclei; lipidomics = liver (10 mg tissue); proteomics = isolated liver KCs. No omics on gut, adipose, brain or muscle.

## Biological replicates per group (per assay)
n differs by assay; all values below are quoted from figure legends.
- snRNA-seq: n = 2 per condition. "snRNA-seq of livers of mice (n = 2 per condition) after chronic MNP feeding" (Fig. 3d legend). This is the ONLY n given for the snRNA-seq; Fig. 1e legend gives no n.
- Lipidomics: "n = 9 (Ctrl), n = 8 (NP) and n = 10 (MP)" (Fig. 1j legend); total-lipid panel n = 9/8/10 (ED Fig. 6g); DBI n = 9/9/10 (ED Fig. 6i).
- Proteomics of KCs: n not reported in the text read (ED Fig. 9a legend gives no n); "Three samples were excluded because of 20% lower protein identifications" (Methods, Proteomics).
- Oil Red O (chow, 12 wk): "n = 6 (Ctrl), n = 10 (NP) and n = 9 (MP)" (Fig. 1i). ORO 4 wk: n = 7/8/8 (ED Fig. 5h). ORO under HFD: "n = 5 (Ctrl), n = 7 (NP) and n = 10 (MP)" (Fig. 2i).
- PAS glycogen (100 mg/kg): n = 9/9/10 (ED Fig. 5d); 10 mg/kg: n = 4/4/4 (ED Fig. 6b); 1 mg/kg: n = 4/5/5 (ED Fig. 6e).
- Body-weight curve, chronic: "n = 35 (Ctrl), n = 32 (NP), n = 39 (MP)" (ED Fig. 4a); HFD weight gain n = 18/18/17 (Fig. 2b).
- Platelet clearance 4 wk: n = 8/7/9 (Fig. 3j); 12 wk time-course: n = 7/10/12 (Fig. 3l); Clec4fDTR arm: n = 9 (Ctrl), n = 10 (MP) (Fig. 4d).
- E. coli CFU blood: n = 7 (Ctrl), n = 8 (NP), n = 6 (MP), "pooled from two independent experiments" (Fig. 3k).
- Co-culture: "Ctrl, 30; MP, 40 [images, technical replicates] ... from a total of n = 3 (Ctrl) and n = 4 (MP) biological replicates" (Fig. 3c legend).
- Serum cytokines (Luminex): n = 5/4/4 (ED Fig. 9d). Myeloid cell counts: n = 8/7/8 (ED Fig. 9c).
Many panels state "Data are pooled from >2 independent experiments" (Fig. 1, Fig. 2, Fig. 3 legends).

## Polymer; particle size; supplier; label; characterisation; endotoxin
- Polymer: polystyrene (PS) only, carboxylated (negatively charged). No polyethylene, no PET, no PVC.
- Sizes: NP = 50 nm, MP = 500 nm for all chronic work; bioavailability screen also used 100 nm and 300 nm. Nominal sizes only; measured sizes are in Supplementary Fig. 1, NOT VERIFIED (not in the XML).
- Suppliers: "Yellow-Green Fluorescent Carboxylated PS; Polysciences, 15700-10, 17149-10, 16662-10, 2405-10; Red Fluorescent Carboxylated PS; Megasphere, CAF050NM, CAF100NM, CAF200NM, CAF500NM" (Methods). Unlabelled negatively charged PS for fate-mapper work: Polysciences 15913-10, 09836-15. Gold-core PS-coated carboxylated particles: CD Bioparticles CPGN-500 / CPGN-50. Europium MP 500 nm: Abvigen ABPST-1-0500-EuC.
- Labels: yellow-green or red fluorescence; gold core; europium; unlabelled.
- Characterisation: DLS on a Zetasizer Ultra - "Z-average hydrodynamic diameter and polydispersity index", zeta potential measured 5x per sample (Methods, "Size and charge assessment of MNPs"); values reported only in Supplementary Fig. 1 (not read).
- Endotoxin: NOT REPORTED. Searched Methods, particle sections and whole text for "endotoxin"/"LPS" - no hit. Biocorona is discussed conceptually in the introduction but not measured.
- Authors note they could NOT quantify tissue plastic burden: Europium signal was lost, and pyrolysis-GC-MS failed - "In some instances, samples from MNP-exposed mice yielded quantitative values that were comparable to or even lower than those obtained from Ctrl tissues" (Results, pyrolysis paragraph).

## Dose; route; vehicle; frequency; duration
Doses: 100 mg/kg body weight standard; titration arms at 10 mg/kg and 1 mg/kg; gold-core arm 0.8 mg/week (Methods) / "0.8 mg kg-1 twice per week" (Fig. 1d legend - inconsistent). Route: oral gavage (one i.v. arm for europium only). Vehicle: PBS (also the control gavage). Frequency: once weekly. Duration: 12 weeks (4 weeks for "4W" panels; 5 weeks for gold-core; 24 h for acute). Verbatim: "mice received PBS or MNPs once weekly by oral gavage for 12 weeks ... Doses were 1 mg kg-1 or 10 mg kg-1 of body weight for titration experiments and 100 mg kg-1 of body weight otherwise" (Methods).

## Diet / metabolic-stress model
Control diet: "a CD (13 kJ% fat; Ssniff, E15748-047)"; HFD: "an HFD (60 kJ% fat; Ssniff, E15742-347)" (Methods, Chronic plastic administration with HFD feeding). Units are kJ% fat, not kcal%. Timing: 2 wk CD acclimation -> 8 wk MNP on CD -> switch to HFD + 4 wk MNP -> 6 further weeks HFD without MNP (Fig. 2a). The non-HFD (main) experiments are chow-fed/lean. There is NO MASH/NASH diet (no fructose/cholesterol/CDAHFD) and no steatohepatitis model; authors state this explicitly as a limitation.

## Data types
- snRNA-seq (single-NUCLEUS, not bulk): liver nuclei, 10x Chromium Single Cell 3' v3.1 Dual Index, 16,000 nuclei loaded per sample, "sequenced on an Illumina NextSeq 2000 with P4 flow cells to an average depth of 30,000 reads per nucleus"; CellRanger v7.1.0, Seurat v5.0.1 (Methods, snRNA-seq). Nine cell types identified (Fig. 1e).
- Shotgun lipidomics of liver (10 mg tissue, internal standards, Thermo Q Exactive Plus, LipidXplorer) (Methods, Lipidomics).
- DIA proteomics of sorted/enriched liver KCs (Evosep One + timsTOF Pro 2, diaPASEF, DIA-NN v1.8.1, MaxLFQ) (Methods, Proteomics).
- Flow cytometry (particle uptake per cell type, myeloid counts, fate-mapping), confocal immunofluorescence, STEM/serial block-face EM, histology (ORO, PAS, H&E, Sirius Red, insulin IHC) with blinded pathology scoring.
- Serum: AST, ALT, glucagon, FGF21 ELISAs; Piccolo Lipid Plus panel (cholesterol, HDL, TG); 48-plex Luminex cytokines.
- Physiology: GTT, ITT, fasting glucose, indirect calorimetry (PhenoMaster, O2 consumption, activity, food/water intake at 22 C and 4 C), NMR body composition.
- RT-qPCR on liver (fibrosis genes, washout arm; normalised to Hprt and Tatabp, Pfaffl method).
- Pyrolysis-GC-MS and europium fluorimetry (both reported as FAILED for quantification).

## Statistics & computational methods
- Phenotypes: GraphPad Prism 10; Shapiro-Wilk for normality; one-way ANOVA + FDR or unpaired t-test (normal data); Kruskal-Wallis + FDR (non-normal); two-way ANOVA + FDR for two factors (condition x time); Wilcoxon + FDR for ranked data; "Outliers were identified using Grubbs' test (alpha = 0.05) and removed" (Methods, Statistics). Data as mean +/- s.d.
- snRNA-seq DE: "Differential expression and ORA compared MP and NP with Ctrl using a Wilcoxon test", two-sided Wilcoxon rank-sum with Benjamini-Hochberg (Fig. 3h legend). Thresholds: "adjusted P <= 0.05 with log2(fold change) >=0.75 or <=-0.75 in hepatocytes ... and adjusted P <= 0.1 with log2(fold change) >=0.25 or <=-0.25 in all other cell types" (Methods). Note Fig. 1g legend states hepatocyte DEGs used "log2(FC) above 0.75 or below -0.75 and q < 0.1" - the q threshold differs between Methods (0.05) and Fig. 1g legend (0.1); AMBIGUOUS.
- Gene filtering: "Per cell type, genes outside the 0.5th-99.5th expression percentile were excluded as outliers, and per sample, only genes expressed in >=1% of cells were retained".
- ORA: clusterProfiler v4.10.0 against Hallmark/Reactome/GO-BP via msigdbr v7.5.1, GO-BP restricted to first four hierarchy levels, gene sets 10-300 genes, FDR <= 10% with >= 3 genes.
- Module scores: Seurat AddModuleScore on curated MSigDB lipid gene sets, averaged per mouse, expressed as delta vs PBS (no significance test reported for module scores).
- TF activity: VIPER v1.36.0 with DoRothEA mouse regulons v1.14.1 on pseudobulk; Wilcoxon + BH. Explicitly: "Although no single transcription factor reached significance" (Results, VIPER paragraph on hepatocytes).
- Cell-cell communication: CellChat v1.6.0, default parameters, genes overexpressed if in >=10% cells, >=10 cells/group.
- Proteomics: Perseus v1.6.15, Welch's t-test with S0 = 0.1, FDR < 0.05, 250 randomisations; imputation from downshifted normal distribution.
- Doublets: DoubletFinder v2.0.4, assumed doublet rate 25%.
- NO mixed model / cage-level or batch covariate is described anywhere for the omics or the phenotypes.

## Main findings with effect directions, evidence and label
Distribution / uptake
- MNPs cross the gut and accumulate mainly in liver, spleen and BAT; not in brain, gWAT, blood, bone marrow, lymph nodes or small intestine after chronic feeding (Results; ED Fig. 2a-c). VALIDATED (flow cytometry + confocal IF).
- KCs are the dominant hepatic reservoir; NPs additionally in CD31+ LSECs; hepatocytes show no detectable intracellular particles (Fig. 1b-d; Supplementary Fig. 5). VALIDATED (flow cytometry, confocal IF, gold-core STEM/serial block-face EM).
- Size dependence of brain entry: NPs <=100 nm detected in microglia, MPs 300-500 nm "not detectable" (ED Fig. 1d). VALIDATED (IF), descriptive/no stats given.
- Particles persist in KCs after a 6-month washout, but the metabolic phenotype resolves (ED Fig. 7b-d). VALIDATED (IF, PAS, Sirius Red, RT-qPCR).

Liver phenotype (chow, 12 weeks, 100 mg/kg)
- INCREASED hepatic neutral lipid with MP only; NP indistinguishable from Ctrl: "Oil Red O staining showed increased neutral lipid deposition in MP-exposed livers at 12 weeks, while NP-treated mice were indistinguishable from Ctrl" (Results; Fig. 1h,i; Kruskal-Wallis + FDR). VALIDATED (histology). At 4 weeks no difference (ED Fig. 5g,h) - time-dependent.
- DECREASED glycogen (PAS) in MP but not NP at 100 and 10 mg/kg; no effect at 1 mg/kg (ED Fig. 5c,d; ED Fig. 6a-f). VALIDATED (histology). Supports a dose threshold.
- NO hepatocellular injury or fibrosis: AST/ALT unchanged and Sirius Red negative at 12 weeks (ED Fig. 5e,f). VALIDATED (serum ELISA, histology).
- Liver weight: trend to reduction in NP, not significant; no change in MP (ED Fig. 4b). VALIDATED (organ weights, Kruskal-Wallis + FDR).
- Lipidomics: MP INCREASED LPE, LPC and sphingomyelin classes; at species level "in NP, none reached significance, whereas in MP five LPC species, seven LPE species and two monogylceride species were significantly elevated" (Fig. 1j; Kruskal-Wallis + Dunn). Total lipid unchanged (ED Fig. 6g). VALIDATED (MS lipidomics; note this is an omics assay but a direct biochemical measurement, not inference).

Liver transcriptomics (snRNA-seq, n = 2/group)
- Hepatocytes, not KCs, had "by far the greatest number of differentially expressed genes" in both NP and MP (Fig. 1f; counts only in Supplementary Table 1, NOT VERIFIED). OMICS-ONLY.
- Hepatocyte ORA converged on lipid pathways: 'lipid localisation', 'lipid export from cell', 'cholesterol homoeostasis' shared; 'glycolysis' MP-specific; 'response to toxic substance' NP-specific (Fig. 1g). OMICS-ONLY (the MP glycolysis term is corroborated by PAS glycogen loss - partial histological support).
- Shared hepatocyte programme: de novo lipogenesis SUPPRESSED, fatty-acid uptake and lipoprotein-receptor expression INCREASED, in both NP and MP (Results; Fig. 1k,l). OMICS-ONLY (module scores + DEGs; no qPCR/western confirmation of individual lipogenic genes). Specific canonical genes (Fasn, Srebf1, Acly, Cd36, Ppara, Cpt1a, Acox1, Cyp4a) are NOT named individually in the readable text - gene-level detail is in Fig. 1l/Supplementary Table 1 (NOT VERIFIED).
- MP-specific: VLDL secretory machinery DOWNregulated ("comparatively preserved in NPs") and Mgll "more strongly reduced in MPs" (Results). OMICS-ONLY. Authors themselves state this is not required for the phenotype: "the co-culture experiment demonstrates that MP-KC-derived signals are sufficient to induce hepatocyte lipid droplet formation ... indicating that a VLDL secretory block is not a prerequisite" (Discussion).
- Nuclear receptors: Nr1h3 (LXRalpha) DOWN - "excluding it as a driver of the observed lipid gene programme"; bile acid synthesis cascade UP (Cyp7a1, Cyp8b1, Cyp27a1) with reduced feedback repressor Nr0b2 (Results). OMICS-ONLY; authors say other nuclear-receptor inputs "may contribute upstream but that we did not test functionally".
- Lands-cycle / ceramide enzyme changes co-occur with the lipidomic shifts, but authors flag "these transcript-lipid associations are correlative and we do not infer causation for individual species" (Results). OMICS-ONLY.
- Hepatocyte VIPER: NP engaged HNF4alpha, PPARalpha, C/EBPbeta; MP a pro-inflammatory signature (STAT1, IRF2, E2F4) with loss of antioxidant regulators MAFG/MAFK (Fig. 1m). COMPUTATIONAL-PREDICTION, and explicitly non-significant ("no single transcription factor reached significance").

KC transcriptome / proteome / function
- KC subclustering: cluster 1 enriched in Ctrl, clusters 0 and 2 in MNP-exposed (Fig. 3d). OMICS-ONLY.
- KC ORA: UP in both - 'response to xenobiotic stimulus', 'lipid localization', 'iron uptake and transport', 'chaperone-mediated autophagy'. MP-specific DOWN - 'binding and uptake of ligands by scavenger receptors', 'regulation of chemotaxis', 'Fc gamma receptor-dependent phagocytosis'; NP UP-regulated the scavenger-receptor pathway (Fig. 3f). OMICS-ONLY.
- Apolipoproteins Apoa1, Apoa2, Apoc3 UP in KCs in both conditions; 109 commonly upregulated genes (Fig. 3g,h). OMICS-ONLY at RNA level; PARTIALLY VALIDATED at protein level - "Among the proteins significantly upregulated in MP-exposed KCs were APOC3 and transthyretin (TTR) ... thereby validating key transcriptomic predictions at the protein level" (Results; proteomics, Welch's t-test, FDR < 0.05). Note Apoa1/Apoa2 protein were not stated as validated.
- KC VIPER: shared PPARalpha, C/EBPalpha, HNF4alpha, MNT, HIF-1alpha activity; NP-selective PPARgamma; MP-selective MYC, PU.1, IRF1 with loss of MAFB (Fig. 3i). COMPUTATIONAL-PREDICTION. The proposed HIF-1alpha -> apolipoprotein -> hepatocyte lipid axis is explicitly a hypothesis here: "raising the possibility that a similar HIF-1alpha-dependent mechanism drives KC apolipoprotein secretion" - HIF-1alpha was NOT manipulated in this study.
- CellChat KC->hepatocyte: increased information flow for VCAM, OCLN, AGT, SEMA4 in both conditions (ED Fig. 9b). COMPUTATIONAL-PREDICTION; the authors then tested the inflammation implication and found none (see below).
- KEY FUNCTIONAL VALIDATION 1 (paracrine sufficiency): KCs isolated from MP-exposed mice co-cultured with naive primary hepatocytes INCREASED hepatocyte lipid droplet area vs Ctrl KCs (Fig. 3b,c; two-sided Mann-Whitney; n = 3 vs 4 mice, 30/40 images). VALIDATED (ex vivo co-culture + imaging). This is the strongest causal link in the paper: KC-derived factors are SUFFICIENT for hepatocyte lipid accumulation. The specific mediator was NOT identified (no neutralisation/knockdown of APOC3, TTR or any ligand).
- KEY FUNCTIONAL VALIDATION 2 (phagocytosis): MP (not NP) IMPAIRED KC Fc-gamma-receptor-dependent clearance of opsonised platelets already at 4 weeks (Fig. 3j, Kruskal-Wallis + FDR) and sustained at 12 weeks (Fig. 3l, mixed-effects + FDR), i.e. BEFORE lipid accumulation appears. VALIDATED (in vivo anti-platelet antibody clearance assay).
- KEY FUNCTIONAL VALIDATION 3 (bacterial clearance): MP mice cleared i.v. E. coli from blood less efficiently (Fig. 3k, one-way ANOVA + FDR); liver and spleen CFU unchanged (ED Fig. 9f). VALIDATED (CFU assay).
- NO inflammation: no increase in hepatic KC/monocyte/neutrophil/DC numbers and no serum cytokine elevation (IL-1beta, IL-12p70, G-CSF, Ccl7, IL-18 shown) (ED Fig. 9c,d). VALIDATED (flow cytometry, Luminex). Authors conclude MP "reprograms KC core functions rather than triggering a classical inflammatory response".
- KC ontogeny: KCs remain ~100% yolk-sac-derived after 12 weeks of exposure; no monocyte replacement (Fig. 4b; ED Fig. 10a,b). VALIDATED (dual fate-mapping, flow cytometry).
- After DT-mediated KC depletion, MP-exposed mice showed NO recovery of phagocytic function in monocyte-derived KCs (Fig. 4d, Mann-Whitney) and significantly FEWER moKCs at 12 weeks post-DT (Fig. 4e, unpaired t-test) with fewer IBA1+ macrophages (Fig. 4f,g), while monocyte influx, pre-moKC emergence/proliferation, ontogeny and niche factors (Dll4, Csf1, Il34) were unaffected (ED Fig. 10g-j). VALIDATED (genetic depletion/fate-mapping + flow cytometry + IF). Interpretation: maintenance, not generation, of moKCs is impaired.

HFD (diet-induced obesity) arm
- MP INCREASED and NP DECREASED hepatic lipid (ORO) relative to Ctrl under HFD (Fig. 2h,i). VALIDATED (histology). This is the paper's "exacerbation" claim for MP.
- NP reduced weight gain, reduced gWAT and BAT weight; MP showed a non-significant trend (Fig. 2b,c). VALIDATED (weights).
- Fasting glucose: no effect before HFD; both NP and MP HIGHER after HFD (Fig. 2d). GTT AUC significantly increased in NP only; ITT unchanged in all groups (Fig. 2e-g). VALIDATED (physiology).
- NP INCREASED oxygen consumption at 22 C and 4 C with reduced fat mass, no change in food/water intake or locomotion (Fig. 2j-n). VALIDATED (indirect calorimetry, NMR).
- Serum total cholesterol, HDL-C and TG did not differ (Supplementary Fig. 8, legend only). NOT VERIFIED numerically.
- INTERACTION TESTING: the MP x diet "exacerbation" was NOT tested as a formal particle x diet interaction term. The HFD cohort is analysed within itself by one-way ANOVA/Kruskal-Wallis across Ctrl/NP/MP, and two-way ANOVA is used only for condition x TIME (weight curves, GTT/ITT curves). Chow and HFD cohorts are separate experiments with different timelines, so "exacerbates" rests on pairwise comparisons within the HFD cohort plus the separate chow cohort - IMPORTANT for the team.

## Pathway-theme directions
- Lipogenesis / de novo lipid synthesis: DOWN (suppressed) in hepatocytes for BOTH NP and MP (module scores + DEGs). OMICS-ONLY. Individual genes (Fasn, Srebf1, Acly, Pparg) not named in the readable text.
- Fatty-acid uptake / lipoprotein receptors: UP in hepatocytes, both sizes. OMICS-ONLY. (Cd36 not named explicitly.)
- Lipid export (VLDL machinery) and turnover (Mgll): DOWN, MP-selective. OMICS-ONLY.
- PPAR/FAO: PPARalpha appears only as VIPER-inferred activity - increased in NP hepatocytes and in KCs of both conditions; PPARgamma NP-selective in KCs. COMPUTATIONAL-PREDICTION, and non-significant in hepatocytes. No Ppara/Cpt1a/Acox1/Cyp4a transcript or protein data in the readable text; no FAO flux measured.
- Cholesterol: "cholesterol homoeostasis" enriched in hepatocyte ORA (both sizes; up/down direction of the term not stated); bile-acid synthesis genes UP (Cyp7a1, Cyp8b1, Cyp27a1) with repressor Nr0b2 DOWN; Nr1h3 (LXRalpha) DOWN. OMICS-ONLY. Serum total cholesterol and HDL-C unchanged in the HFD arm (VALIDATED, serum panel).
- Inflammation: NO classical inflammatory phenotype - no hepatic immune-cell infiltration, no serum cytokine elevation, no gWAT crown-like structures (VALIDATED by flow cytometry, Luminex, blinded histology). A pro-inflammatory signature exists only as inferred TF activity in MP (STAT1, IRF2, E2F4 in hepatocytes; MYC, PU.1, IRF1, loss of MAFB in KCs) - COMPUTATIONAL-PREDICTION. Tnf, Il1b, Il6, NF-kB and NLRP3 are not reported as regulated; serum IL-1beta and IL-18 were measured and unchanged.
- Oxidative stress: only the inferred loss of MAFG/MAFK antioxidant-regulator activity in MP hepatocytes - COMPUTATIONAL-PREDICTION. No ROS, MDA, SOD, GSH, Nrf2 or Hmox1 measurement anywhere in the paper (full-text search). This paper does not provide direct oxidative-stress evidence.
- Mitochondrial function: not assayed at organelle level in liver. Whole-body oxygen consumption INCREASED in NP under HFD (VALIDATED, indirect calorimetry); BAT lipid vacuoles smaller in both NP and MP (VALIDATED, blinded pathology scoring).
- Fibrosis: ABSENT - Sirius Red negative at 12 weeks and after washout; fibrosis-gene RT-qPCR after washout showed no pathology (VALIDATED).
- Circadian: NOT addressed. No clock genes; no Zeitgeber time or time-of-sacrifice statement in Methods (searched "circadian", "Zeitgeber", "time of day" - no hits). Housing is described only as "a 12 h light-dark cycle".
- MP x diet interaction: asserted verbally ("microplastics exacerbate hepatic lipid accumulation" under diet-induced obesity, Abstract) but NOT formally tested as an interaction term (see Main findings).

## Author-acknowledged limitations (quotes + location)
- Disease stage: "Our findings were obtained under chow-fed (lean) conditions and a relatively short HFD challenge, which models early, mild metabolic perturbation rather than advanced steatohepatitis" (Discussion).
- Unmeasured endpoints: "We also did not measure fasting insulin or non-esterified fatty acids, and therefore did not assess systemic insulin resistance or adipose tissue function directly" (Discussion).
- Sex: "our study was conducted mainly in young male mice, which develop diet-induced hepatic steatosis and metabolic dysfunction more rapidly than females. Whether females show comparable size-dependent susceptibility ... remains to be determined" (Discussion).
- Tissue plastic quantification: "reliable quantification of MNPs in tissues is currently limited by methodological constraints in detection and mass-based conversion"; "Whether the plastic is gradually excreted, or whether the reduced fluorescent signal is caused by degradation or signal leakage remains unclear" (Discussion).
- Mechanism not reduced to one ligand: "we propose that their combined perturbation, rather than any single ligand-receptor interaction, drives progressive hepatic lipid retention selectively in MP-exposed livers" (Results, end of KC section).
- Untested upstream regulators: other nuclear-receptor inputs "were altered in a manner that may contribute upstream but that we did not test functionally" (Results, hepatocyte section).
- Correlative omics: "these transcript-lipid associations are correlative and we do not infer causation for individual species" (Results).
- Competing NP explanations: "The two possibilities are not mutually exclusive, and distinguishing between them will require comprehensive in vivo tracer-based lipid flux studies" (Discussion).
- LSECs: "the contribution of LSECs to metabolic dysfunction in this context remains largely unexplored" (Discussion).
- Statistical candour on TF inference: "Although no single transcription factor reached significance" (Results, hepatocyte VIPER paragraph).
- Removed data: "behavioural experiments, which had to be removed during the revision of the paper" (Author contributions).
- Descriptive imaging: "Image is not representative of the distribution of MNPs across the whole liver but is merely descriptive" (Fig. 1d legend).

## Methodological points to verify independently
(Questions to check, not criticisms.)
1. snRNA-seq is n = 2 mice per condition (Fig. 3d legend) - 6 libraries total. Was differential expression run per NUCLEUS (Wilcoxon across nuclei) rather than per mouse? VIPER was explicitly pseudobulked by sample, DE apparently was not; nuclei-as-replicates inflates significance at n = 2.
2. Were the 6 snRNA-seq libraries processed and sequenced in a single batch, or is condition collinear with 10x lane / sequencing run? Not stated.
3. The hepatocyte DEG threshold differs between Methods ("adjusted P <= 0.05") and the Fig. 1g legend ("q < 0.1"). Which applies to the reported ORA?
4. Gold-core dosing conflicts: Methods "0.8 mg per week ... for 5 weeks" vs Fig. 1d legend "0.8 mg kg-1 twice per week, alternating NPs and MPs". Which, and mg or mg/kg?
5. Grubbs outlier removal (alpha = 0.05) was applied to phenotype data - how many points were removed per panel, and was this pre-specified? Source Data files would show it.
6. Cage/litter as the unit of analysis is not mentioned. Were Ctrl, NP and MP mice co-housed or housed by treatment group? Treatment-wise housing would make cage collinear with group, exactly the issue the team has.
7. Time of sacrifice / circadian phase is not reported. PAS glycogen, lipogenic and bile-acid transcripts are strongly diurnal - were all groups collected in the same window?
8. Measured particle size, PDI, zeta potential (Supplementary Fig. 1) and endotoxin content could not be read. Commercial carboxylated PS suspensions often contain surfactant/azide - were particles washed, and was a bead-free vehicle control matched for preservatives?
9. Is PBS gavage a sufficient control, or would a surfactant-only control be needed to separate particle effects from formulation effects?
10. n varies between panels from the same cohort (ORO n = 6/10/9 vs PAS n = 9/9/10) - which animals contributed to which assay and were any excluded?
11. In the HFD arm, MNP dosing STOPPED 6 weeks before sacrifice and the cohort differs in age and timeline from the chow cohort; was any model ever fitted with particle x diet as factors?
12. Proteomics n per group is not stated in the readable text, and 3 samples were excluded post hoc on a QC criterion - check Supplementary Table 5.
13. Co-culture hepatocytes came from 5-week-old FEMALE mice while all in vivo work is male - does the sex mismatch affect the paracrine readout?
14. "Data pooled from >2 independent experiments" - was experiment/batch ever included as a factor in the statistical models?
15. KC preparations for proteomics and co-culture used CD11b magnetic enrichment, not sorted Clec4f+/Vsig4+ KCs - what purity, and could contaminating myeloid cells drive the co-culture effect?
16. Module scores (Fig. 1k) are presented as deltas with no significance testing - are any of these differences statistically supported?

## Data availability statement (verbatim) and all accessions/code links
Verbatim (section "Data availability"): "The data included in the paper can be accessed as follows: for lipidomics, at Metabolomics Workbench, study ID ST003776; liver snRNA-seq, GSE291616; histology, bonndata-Image data related to publication 'Size-dependent plastic exposure disrupts macrophage function and tissue-specific metabolism'54; Proteomics, PRIDE PXD072158. Source data are provided with this paper."
Accessions / links mentioned anywhere (reported as stated; existence not verified here):
- Metabolomics Workbench study ST003776 (lipidomics) - Data availability.
- GEO GSE291616 (liver snRNA-seq) - Data availability. NOTE: this is NOT the team's GSE294180.
- PRIDE PXD072158 (KC proteomics) - Data availability AND repeated in Methods, Proteomics ("The proteomics data are available at PRIDE (PXD072158)").
- Bonndata histology image deposit, reference 54, DOI 10.60507/FK2/M42B8O (reference list entry for CR54).
- Source Data files for Figs 1-4 and Extended Data Figs 2-10 (Source data section).
- No Code availability statement and no GitHub/Zenodo/figshare code link anywhere (searched). Software versions are listed in Methods instead (CellRanger 7.1.0, Seurat 5.0.1, DoubletFinder 2.0.4, clusterProfiler 4.10.0, msigdbr 7.5.1, VIPER 1.36.0, DoRothEA 1.14.1, CellChat 1.6.0, DIA-NN 1.8.1, Perseus 1.6.15, GraphPad Prism 10, FlowJo 10.8.1, ImageJ/Fiji, QuPath 0.5.1/0.6.0, LipidXplorer, AdipoQ).
- Competing interests: "The authors declare no competing interests."

## Relevance to the team
Closest high-quality MP x HFD liver study, but with important differences: polystyrene (not polyethylene), weekly gavage at 100 mg/kg, a plain 60 kJ% HFD rather than a MASH diet, and single-NUCLEUS RNA-seq of n = 2 mice per group rather than bulk liver RNA-seq of n = 4-5. It supports an MP-selective increase in hepatic lipid under both chow and HFD, with suppressed de novo lipogenesis plus increased lipid uptake and impaired VLDL export as the transcriptional pattern - a direction-of-effect template that can be compared with PE-exposed bulk liver signatures. Its causal step (KC-conditioned paracrine signalling sufficient to load hepatocytes with lipid) is validated ex vivo, whereas the HIF-1alpha/apolipoprotein axis, the TF activities and the inflammation signature are inference only; the paper reports no inflammation, no fibrosis, no oxidative-stress measurement and no circadian control, and never fits a formal particle x diet interaction.
