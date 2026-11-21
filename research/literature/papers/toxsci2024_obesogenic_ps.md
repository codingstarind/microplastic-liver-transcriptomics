# toxsci2024_obesogenic_ps

## Citation
Obesogenic polystyrene microplastic exposures disrupt the gut-liver-adipose axis. Zhao J, Adiele N, Gomes D, Malovichko M, Conklin DJ, Ekuban A, Luo J, Gripshover T, Watson WH, Banerjee M, Smith ML, Rouchka EC, Xu R, Zhang X, Gondim DD, Cave MC, O'Toole TE. 2024 (epub 2024 Jan 30). Toxicological Sciences 198(2):210-220. DOI 10.1093/toxsci/kfae013. PMID 38291899. PMCID PMC10964747. Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC10964747/ (citation fields from Europe PMC core metadata, meta_toxsci.json, and PMC header line "Toxicol Sci. 2024 Jan 30;198(2):210-220").

## Access level read
Full text (PMC HTML rendered to text: dl_B/toxsci_body.txt / toxsci.txt), including Methods, Results, figure/table legends, Discussion, Data availability. Supplementary data (kfae013_supplementary_data.zip: Supp Figs 1-9, Supp Tables 1-6) NOT read; anything that lives only in supplements is marked "NOT VERIFIED (supplement not read)".

## Biological question
Whether drinking-water polystyrene (PS) beads (0.5 and 5 um) disrupt the gut-liver-adipose axis in chow-fed mice (adipose inflammation/function, liver phenotype, hepatic transcriptome, bile acids, glutathione, ileum), and whether the anthocyanin delphinidin mitigates PS-induced weight gain/adiposity (Abstract; Introduction last paragraph).

## Experimental design (groups, controls, randomisation/blinding if stated)
- Main study: 3 groups - unsupplemented water (control), 0.5 um PS, 5 um PS in drinking water, 13 weeks (Methods "Animals and treatment").
- Allocation was by CAGE: "randomly selected cages (2 cages, n = 10 mice per treatment) were supplied with unsupplemented water ... or that containing" PS (Methods). Housed "5 per cage". Hence cage is confounded with treatment (2 cages/group).
- Delphinidin sub-study: "mice drinking PS-containing water received i.p. injections of delphinidin chloride ... (20 mg/kg) or vehicle (20% DMSO) twice weekly for 4 weeks" (Methods). Fig. 4 legend groups: unsupplemented water "U; n = 4", PS + vehicle "Veh; n = 8", PS + delphinidin "Del; n = 8"; performed with 0.5 um beads (Results). Whether this was a separate cohort from the 13-week study is not explicitly stated (Methods says body weights obtained "Prior to water supplementation and at the 4 weeks time point", implying a separate 4-week cohort) - ambiguous.
- No diet factor: all mice on chow ("the mice in the present study were fed a chow diet", Discussion). No MP x diet design.
- Blinding: not reported (searched Methods/Results). Randomisation: only at cage level as quoted above.
- Some analyses restricted to subsets: 5 um group "excluded from this analysis" of beta-catenin because it lacked significant weight gain (Results); adipokines measured only in 0.5 um vs water (Fig. 2).

## Species/strain/sex/age; tissue(s)
"Male C57BL/6 mice at 12 weeks of age were purchased from Jackson Laboratories" (Methods). Abstract says C57BL/6J. Tissues: plasma, epididymal WAT (eWAT), perivascular adipose tissue (PVAT), liver, cecum (bile acids), ileum (histology).

## Biological replicates per group (per assay - RNA-seq n may differ from phenotype n)
- Enrolled: "2 cages, n = 10 mice per treatment" (Methods).
- Liver mRNA-seq: Fig. 6 legend "Cut-off values are q-value < 0.05. n = 10."; Fig. 7 legend "n = 10" (presumably per group; not explicitly "per group"). Whether libraries were per mouse or pooled: not stated (assumed per mouse - verify).
- Adipose (eWAT) RNA-seq (delphinidin study): n not reported in main text (Supp Table 1; NOT VERIFIED). Group sizes for that study: Veh n = 8, Del n = 8 (Fig. 4 legend).
- eWAT flow cytometry: "n = 8-10" (Fig. 1 legend). Plasma adipokines: "n = 10" (Fig. 2 legend). beta-catenin western: "n = 4" (Fig. 3 legend).
- Adipose qPCR (Table 1): eWAT n = 6-9, PVAT n = 4-5 per gene (Table 1 cells, e.g. "Tnfa ... 0.92 (0.20) n = 9 | 3.02 (1.04) n = 5").
- Liver histology/liver:body weight/hepatic TG/cholesterol (Fig. 5): n not stated in legend.
- Bile acids: n per group not reported in main text; LC-MS used individual samples plus "Group-based pool samples" for MS/MS identification (Methods "Bile acid analysis").

## Polymer; particle size (nominal and measured); supplier; surface/fluorescent label; characterisation; endotoxin
- Polymer: polystyrene beads.
- Size: "'0.5 um' (size range = 0.4-0.6 um) or '5 um' (size range = 4.5-4.9 um) PS beads (Creative Diagnostics; Shirley, New York)" (Methods). These ranges appear to be supplier specifications; no in-house measured size.
- Surface chemistry / fluorescent label / catalogue number: not reported (Methods).
- Characterisation (DLS, SEM, zeta potential, FTIR): not reported in main text.
- Endotoxin testing: not reported (searched full text for "endotoxin").
- Tissue accumulation of beads not measured - an acknowledged limitation (see below).

## Dose (as reported, and units); route; vehicle; frequency; duration
- "at a concentration of 1.0 ug/ml" in drinking water, ad libitum, 13 weeks (Methods). Authors estimate intake from "average water consumption of 3-5 ml per day per mouse" (Methods), but actual water intake not reported as measured.
- Vehicle: "filtered, autoclaved water"; stocks "sonicated for 30 min prior to dilution"; "Fresh water solutions were prepared every other day"; bottles agitated on other days (Methods).
- Delphinidin: 20 mg/kg i.p. twice weekly for 4 weeks; vehicle 20% DMSO (Methods).

## Diet / metabolic-stress model (exact diet, % kcal fat, supplier, control diet)
No metabolic-stress diet. "the mice in the present study were fed a chow diet" (Discussion). Chow brand/composition/% kcal fat: not reported (Methods). Authors explicitly note this design cannot address diet-induced NAFLD (see limitations).

## Data types (RNA-seq platform/depth, proteomics, metabolomics, lipidomics, scRNA, spatial, histology, serum)
- Liver bulk mRNA-seq: "RNAseq analysis was performed by the Center for Genetics and Molecular Medicine at UofL" (Methods). Platform, read length, depth, library prep, aligner: not reported in main text.
- eWAT bulk RNA-seq (delphinidin study): performed and analysed (KEGG/GO/Reactome) by Novogene (Methods). Platform/depth not reported.
- Targeted qPCR in adipose (TaqMan assays, Hprt reference, ddCT) (Methods).
- Bile acids by LC-MS (Q Exactive HF + Ultimate 3000; ESI-) in liver and cecum (Methods).
- Hepatic GSH/GSSG/Cys/CySS by HPLC; hepatic cholesterol, TG, FFA by colorimetric kits normalised to liver weight (Methods).
- Histology: H&E of liver and ileum (Methods).
- Plasma: Luminex (adiponectin, resistin, leptin, PAI-1, ICAM-1, P-selectin, ALT, AST); E-selectin ELISA (Methods).
- Adipose flow cytometry (CD45, CD3, CD4, CD8a, CD11c, F4/80, CD301); western blot beta-catenin, Fam3b; DEXA body composition (Methods).

## Statistics & computational methods (DE tool, thresholds, multiple testing, pathway tool, models; tests for phenotypes)
- Liver DEGs: "Significance for upregulated and downregulated hepatic genes were set at a q-value <0.05 with a log2 fold change = 0" (Methods). DE software (e.g., DESeq2/edgeR) NOT reported. q-value implies FDR-type correction, method not named.
- Liver pathways: Ingenuity Pathway Analysis on DEGs per comparison; "overrepresented if -log(p-value) >1.3" (i.e., unadjusted p < 0.05) (Methods).
- Adipose RNA-seq: KEGG/GO/Reactome by Novogene; DEG thresholds not reported in main text ("DEGs based upon a fold change are listed in Supplementary Table 1", Results).
- Phenotypes: GraphPad Prism 9.5.1; mean +/- SD; normality (D'Agostino & Pearson) and Bartlett's test, then one-way ANOVA + Tukey, Kruskal-Wallis + Dunn, or Brown-Forsythe/Welch ANOVA + Dunnett; "t tests were used to compare between 2 groups. p < .05" (Methods "Statistical analysis").
- Bile acids: "One-way analysis of variance (ANOVA) ... followed by post hoc Tukey's multiple comparison" (Methods).
- Statistical unit = mouse; cage not modelled (2 cages/group). No factorial model; no interaction terms (no second factor exists in the design).

## Main findings with effect directions and statistical evidence
- Body weight/adiposity: 0.5 um PS "led to a potentiated weight gain compared with control mice after 13 weeks" (Results; data in Supp Fig. 1 - NOT VERIFIED); 5 um PS "did not show a significant, potentiated weight gain" (Results). VALIDATED (body weight, DEXA).
- eWAT F4/80+ macrophages increased with 0.5 um PS (*p < .05, Fig. 1A, n = 8-10); NOT with 5 um. M2 no difference; M1 "trending, but nonsignificant (p = .06) increase"; CD3+ no difference (Results/Fig. 1). VALIDATED (flow cytometry).
- NULL: plasma adiponectin, resistin, leptin, PAI-1 - "no differences" 0.5 um vs water (Fig. 2, n = 10). NULL: ICAM-1, P-selectin, E-selectin similar (Supp Fig. 2). VALIDATED (immunoassay).
- Adipose qPCR (0.5 um vs control, Table 1): eWAT IL-6 0.76, Mcp-1 1.04, Tnfa 0.92 ("only minor changes"); PVAT Tnfa 3.02 (SE 1.04, n = 5), Mcp-1 1.21 ("mild upregulation"). No p-values given in Table 1. VALIDATED (qPCR) but statistical significance not reported.
- NULL: beta-catenin unchanged in eWAT and PVAT (Fig. 3, n = 4). VALIDATED (western).
- Delphinidin (PS + Del vs PS + Veh): attenuated weight gain and decreased adiposity (Fig. 4; p-values not quoted in text). VALIDATED (weights/DEXA). Adipose RNA-seq: top up DEG Akr1b7, top down Fam3b "confirmed by rtPCR" and Fam3b by western (Results). KEGG: down "cAMP and PPAR signaling", up "cytokine-cytokine receptor interactions," "phagosome," "osteoclast differentiation" (Results; Supp Figs 4-5). Pathways = OMICS-ONLY; Akr1b7/Fam3b VALIDATED (qPCR/western).
- NULL liver histology: "liver histology was unchanged by either 0.5 or 5 um PS exposure ... neither lipid droplets nor necroinflammatory foci" (Results, Fig. 5A). VALIDATED (H&E).
- Liver:body weight DECREASED: 0.5 um "4.56 +/- 0.22, p = .0165"; 5 um "4.33 +/- 0.50, p = .0005" vs control "5.01 +/- 0.29" (Results, Fig. 5B). VALIDATED.
- NULL hepatic TG: "hepatic triglycerides were unchanged" (Fig. 5C). VALIDATED.
- Hepatic cholesterol INCREASED: 0.5 um "2.34 +/- 0.33 mg/g liver; p = .005"; 5 um "2.21 +/- 0.15 mg/g liver; p = .0006" vs control "1.89 +/- 0.16 mg/g liver" (Results, Fig. 5D). VALIDATED (biochemical assay).
- NULL ALT/AST: "unchanged between exposure groups" (Supp Fig. 6). VALIDATED.
- Liver DEGs (q < 0.05, n = 10): 0.5 um vs water "283 DEGs ... 203 of these upregulated and 80 downregulated"; top up Steap4, Tnfaip3, Ccl2; top down Tpm2, Trib3, Pcp4l. 5 um vs water "fewer DEGs (n = 49) ... Forty-one ... upregulated and 8 ... downregulated"; top up Firre, Fam193b, Mug-ps1; top down Slpi, Ctgf, Sgk1. Five shared DEGs: Rev1, Brpf1, Cnnt2, Cavin2, Pea15a (Results, Fig. 6). OMICS-ONLY (no qPCR validation of liver DEGs reported).
- IPA: "FXR/RXR activation" and "LXR/RXR activation" enriched by both sizes; 0.5 um only: "glutathione-mediated detoxication", "LPS/IL1 mediated inhibition of RXR function", "hepatic cholestasis", AHR/PXR/NRF2 xenobiotic metabolism, cytokine signalling, glucocorticoid receptor, linoleate metabolism; 5 um only: "sucrose degradation", "glycolysis", cytokine signalling (Results, Fig. 7). COMPUTATIONAL/OMICS-ONLY (enrichment; activation direction/z-scores not reported in text).
- Liver bile acids: 0.5 um decreased TCDCA (FC = 0.57, p = .02) and THDCA (FC = 0.33, p < .05); 5 um increased HDCA (FC = 1.93, p = .0002) and THCA (FC = 2.02, p = .03). Cecum: 0.5 um increased ACA (FC = 2.17, p = .005); 5 um increased UDCA (FC = 1.87, p = .03) and gamma-MCA (FC = 3.11, p = .03) (Results). VALIDATED (targeted LC-MS); several other changes only trends (p = .06-.13).
- NULL hepatic redox: "GSH ... GSSG ... as well as the GSH/GSSG ratio ... were unchanged by exposure to either size" (Supp Fig. 8A-C); hepatic protein content increased with 5 um (Supp Fig. 8D). VALIDATED (HPLC). Note this contrasts with the IPA NRF2/glutathione enrichment; authors say that enrichment "was driven by differential expression of glutathione S-transferases (data not shown)" (Discussion).
- NULL ileum: "no pathological changes in response to either 0.5 or 5 um PS" (Supp Fig. 9). VALIDATED (H&E).

## Pathway-theme directions (lipid, PPAR/FAO, cholesterol, inflammation, oxidative stress, mitochondria, fibrosis, circadian, interaction)
- Lipid/lipogenesis: no steatosis, hepatic TG unchanged (VALIDATED NULL). Specific lipogenic genes (Fasn, Srebf1, Acly, Pparg, Cd36): not reported in main text. Linoleate metabolism enriched (0.5 um, IPA; direction not stated).
- PPAR/FAO (liver): not reported for liver. In ADIPOSE, delphinidin (vs PS+vehicle) DOWNregulated KEGG "PPAR signaling" - a delphinidin effect, not a PS effect (OMICS-ONLY).
- Cholesterol: hepatic cholesterol UP with both sizes (VALIDATED); FXR/RXR and LXR/RXR pathways enriched (direction of activation not stated); bile acids size-dependently altered (liver mostly down with 0.5 um, some up with 5 um).
- Inflammation: liver - no necroinflammatory foci (NULL histology), but top up DEGs include Tnfaip3 and Ccl2 (0.5 um) and "LPS/IL1 mediated inhibition of RXR function"/cytokine signalling enriched (OMICS-ONLY). Adipose - F4/80+ macrophages up (0.5 um), M1 trend p = .06, eWAT cytokine genes minor; PVAT Tnfa ~3-fold (no p).
- Oxidative stress: GSH, GSSG, GSH/GSSG unchanged (VALIDATED NULL) despite NRF2/glutathione-detox pathway enrichment (OMICS-ONLY, GST-driven).
- Mitochondria: not reported.
- Fibrosis: not assessed histologically; Ctgf among top DOWN DEGs in 5 um vs water (OMICS-ONLY).
- Circadian: not reported; time of sacrifice not reported.
- Interaction: NOT APPLICABLE - no diet factor (chow only). Authors only speculate that "a second 'hit' (eg, high fat diet [HFD]) may lead to widespread adipose inflammation" and cite other HFD studies (Discussion); no MP x diet test was performed. Size (0.5 vs 5 um) compared only pairwise.

## Author-acknowledged limitations (quotes + location)
- "One limitation of these studies is that we did not determine the distribution or accumulation of ingested PS beads in murine tissues." (Discussion, penultimate paragraph)
- "It is also not known if PS consumption disrupts redox balance or promotes systemic inflammatory responses to affect these tissues." (Discussion)
- "Another limitation of these studies is the use of commercially obtained pristine beads." Natural exposure is "a mixture of plastic types and sizes" (Discussion)
- Diet: "This does not exclude the possibility that PS MPs could worsen diet-induced NAFLD, as the mice in the present study were fed a chow diet." (Discussion)
- "it is difficult to compare results between the studies given differences in MP dose, size, and mode of delivery." (Discussion)
- Mechanism of adiposity: "A clearer definition of the involvement of this mechanism requires further study." (Discussion)

## Methodological points to verify independently
- Cage confounding: 2 cages x 5 mice per group, cages randomised to treatment; was cage accounted for in any analysis (it is not mentioned)?
- Liver RNA-seq: platform, depth, library type, aligner, DE tool, normalisation - none in main text; check Dryad/supplement. Were all 30 liver libraries sequenced in one batch/run? Was library quality/depth balanced across groups?
- "log2 fold change = 0" threshold with q < 0.05: confirm the FDR method and whether any effect-size filter was applied.
- IPA pathways called at unadjusted p < 0.05 (-log p > 1.3); activation z-scores (direction) not reported in text.
- Is n = 10 for RNA-seq per group and per individual mouse (not pooled)?
- Particle characterisation (measured size, surface chemistry, aggregation in water after sonication), endotoxin, and actual water intake/dose received - not reported.
- Body-weight and liver:body-weight data: liver:body weight decrease may partly reflect increased body weight (0.5 um) - check absolute liver weights (supplement).
- Time-of-day / fasting state at euthanasia not reported (relevant to bile acids, Steap4/Tnfaip3, and clock genes).
- Delphinidin cohort: separate from 13-week cohort? Control group "U" n = 4 vs 8/8; no PS-free + delphinidin group, so delphinidin effects are not PS-specific.
- Adipose qPCR fold changes given without p-values; PVAT n = 4-5.
- Single sex (male) only.

## Data availability statement (verbatim) and all accessions/code links mentioned
- Data availability (verbatim): "Data associated with this manuscript are available at: https://doi.org/10.5061/dryad.h44j0zpsj."
- Also listed after Contributor Information: "Dryad Digital Repository DOI: https://doi.org/10.5061/dryad.h44j0zpsj".
- No GEO/SRA/BioProject accession, no code repository (GitHub/Zenodo) mentioned anywhere in the main text (grep for GSE, GEO, SRA, BioProject, github, zenodo: none). Whether raw FASTQs are in the Dryad record is not stated - verify.
- Supplementary: kfae013_supplementary_data.zip (3.8 MB) at journal/PMC.

## Relevance to the team
A polystyrene (not PE) MP liver transcriptome on CHOW only (no diet factor): 0.5/5 um PS produced liver DEGs (283/49), FXR/LXR enrichment and raised hepatic cholesterol without steatosis, inflammation histology, TG or GSH/GSSG change. Useful as a chow-only liver RNA-seq reference with cage-confounded allocation (2 cages/group) and reported NULLs; it provides no MP x diet interaction evidence.
