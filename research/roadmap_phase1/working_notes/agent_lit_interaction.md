# Methods-focused literature matrix: microplastics x metabolic liver disease

Verification pass, 2026-11-21. Six papers, priority order as requested. Everything below was re-checked against primary full text in this session (Europe PMC full-text XML, converted locally); prior session notes in `C:\Users\Soham\Desktop\microplastic\research\literature\papers\` were read first and are explicitly corrected where wrong.

Convention used throughout:
- **VALIDATED** = supported by an orthogonal, non-omics measurement (histology, qPCR, western, enzyme assay, in vivo function, genetic or pharmacological manipulation).
- **OMICS-ONLY** = exists only in the sequencing/MS data.
- **INFERENCE** = computational prediction (TF activity, enrichment, docking) with no direct measurement.
- **UNVERIFIED** = could not be checked against a primary source in this session (almost always: supplementary files not readable).

Critical cross-cutting finding: **not one of the six papers reports a tested diet x particle (or age x diet x particle) interaction term with an F/p value.** Three fit a two-way ANOVA but report only pairwise post-hoc contrasts (#2, #3, and partially #1); three use one-way ANOVA across the 2x2 cells (#4, #5, #6). Every "exacerbates"/"synergistic"/"only under HFD" claim in this set is therefore a difference-of-significance claim, not a demonstrated interaction.

---

# 1. Nature Metabolism 2026 - Size-dependent plastic exposure disrupts macrophage function and tissue-specific metabolism

Makdissi N, Viola MF, et al. (N.M. and M.F.V. equal first; E. Mass corresponding). *Nat Metab* 2026;8(10):2084-2101. DOI https://doi.org/10.1038/s42255-026-01615-8 . PMID 42811166, PMCID PMC13690586. Article type: Letter. Licence CC BY-NC-ND.

**Access level: FULL TEXT read** (Europe PMC full-text XML -> local text, 127,303 characters of body). Read in full: abstract, Main/Results, all Fig. 1-4 legends, all Extended Data Fig. 1-10 legends, all Methods subsections, Statistics, Supplementary-information inventory, Source Data inventory, Acknowledgements, Author contributions, Funding, Data availability, Competing interests. **NOT read:** Supplementary Figs. 1-8, Supplementary Tables 1-7, Supplementary Video 1, Reporting Summary, Peer Review File, Source Data files (none are in the PMC XML).

## 1.1 Exact scientific question
Whether chronically ingested micro- and nanoplastic particles (MNPs) accumulate in liver-resident macrophages, how that reprograms Kupffer cell (KC) transcription and effector function, and whether KC changes propagate to hepatocyte and systemic metabolism **in a particle-size-dependent manner**, including under a high-fat-diet (HFD) challenge. Verbatim (Abstract): "Here we identify KCs as the primary hepatic reservoir for MNPs in young male mice." Intro closing sentence frames the comparison as MPs (500 nm) vs NPs (50 nm) over a 12-week feeding model.

## 1.2 Experimental design, factors, n per group (phenotype vs omics)
Factors: **particle size** (Ctrl/PBS, NP 50 nm, MP 500 nm) x **dose** (1, 10, 100 mg/kg) x **duration** (4W, 12W, +6-month washout) x **diet** (CD vs HFD, separate cohort). Paradigms (Methods "Plastic administration model", "Chronic plastic administration with HFD feeding"):
- Acute bioavailability: single gavage 100 mg/kg fluorescent PS of 50/100/300/500 nm, organs at 24 h (ED Fig. 1a).
- Chronic: weekly oral gavage, 12 weeks, 100 mg/kg standard; titration arms 10 and 1 mg/kg; "4W" = 4-week arm (Fig. 1a).
- Washout: 12 wk exposure then 6-month washout (ED Fig. 7a).
- HFD arm (Fig. 2a): 2 wk CD acclimation -> 8 wk MNP on CD -> switch to HFD with 4 more wk MNP -> 6 further wk HFD **with MNP dosing stopped**.
- Gold-core arm: **Methods say "0.8 mg per week ... for 5 weeks"; Fig. 1d legend says "0.8 mg kg-1 twice per week, alternating NPs and MPs". CONFIRMED CONTRADICTION** (mg vs mg/kg, weekly vs twice weekly).
- Genetic models: Tnfrsf11a-Cre;Rosa26-LSL-YFP;Ms4a3-FlpO;Rosa26-FSF-tdTomato double fate-mapper; Cxcr4-CreERT2;Rosa26-LSL-tdTomato; Clec4f-DTR/+ (DT depletion, 50 ng i.p.) with Clec4f+/+ littermates; Clec4f-DTR/+;Ms4a3-Cre/+;Rosa26-LSL-tdTomato/+.
- Ex vivo: KC-hepatocyte co-culture (KCs from Ctrl or MP mice + hepatocytes from untreated 5-week-old **female** mice).

n per group, quoted from legends (all "biological replicates"):
- **snRNA-seq: n = 2 per condition** (Fig. 3d legend: "snRNA-seq of livers of mice (n = 2 per condition)"). This is the only n given; Fig. 1e legend gives none. Confirmed by GEO: 6 GSMs total.
- **Lipidomics: n = 9 (Ctrl), 8 (NP), 10 (MP)** (Fig. 1j, ED Fig. 6g,h); DBI panel n = 9/9/10 (ED Fig. 6i).
- **Proteomics: n NOT REPORTED anywhere in main text or ED legends** (ED Fig. 9a gives none). Methods: "Three samples were excluded because of 20% lower protein identifications."
- Phenotype n (examples): ORO chow 12 wk 6/10/9 (Fig. 1i); ORO 4 wk 7/8/8 (ED Fig. 5h); ORO under HFD 5/7/10 (Fig. 2i); PAS 9/9/10 (ED Fig. 5d), 10 mg/kg 4/4/4, 1 mg/kg 4/5/5; body-weight curve chow **35/32/39** (ED Fig. 4a); HFD weight gain 18/18/17 (Fig. 2b); organ weights up to 19-24 per group (ED Fig. 4b); platelet 4 wk 8/7/9 (Fig. 3j), 12 wk 7/10/12 (Fig. 3l); E. coli CFU 7/8/6 pooled from two experiments (Fig. 3k); co-culture 30 vs 40 images from n = 3 vs 4 mice (Fig. 3c); serum cytokines 5/4/4 (ED Fig. 9d); myeloid counts 8/7/8 (ED Fig. 9c); Clec4f-DTR platelet 9 vs 10 (Fig. 4d), moKC 7 vs 7 (Fig. 4e), IBA1 5 vs 6 (Fig. 4g).
- Many panels: "Data are pooled from >2 independent experiments" (Figs. 1, 2, 3 legends).

Randomisation: **NOT reported** (grep of body text: the only "randomi" hit is "250 randomizations" in the proteomics statistics). Blinding: only histopathology - "a trained pathologist blinded to experimental groups performed all assessments" (Methods, Immunohistochemistry...). No blinding statement for any other readout.

## 1.3 Datasets, accessions, deposition status - ALL THREE VERIFIED TO RESOLVE
Data availability verbatim: "The data included in the paper can be accessed as follows: for lipidomics, at Metabolomics Workbench, study ID ST003776; liver snRNA-seq, GSE291616; histology, bonndata-Image data related to publication 'Size-dependent plastic exposure disrupts macrophage function and tissue-specific metabolism'54; Proteomics, PRIDE PXD072158. Source data are provided with this paper."

- **GEO GSE291616 - EXISTS.** Verified via NCBI E-utilities esummary (db=gds, uid 200291616): title "Single nucleus (sn) RNA sequencing analysis of murine liver cells after 12 weeks of weekly treatment with either PBS, 50nm polystyrene or 500nm polystyrene"; Mus musculus; GPL30172; **6 samples GSM8838590-GSM8838595** = "Liver, PBS, Sample 1/2", "Liver, 50nm, Sample 1/2", "Liver, 500nm, Sample 1/2" -> confirms n = 2 per condition independently of the figure legend; supplementary file types MTX, TSV; public date 2026/06/09. (The GEO web page itself serves a reCAPTCHA to automated fetches; E-utilities is the primary-source route.)
- **PRIDE PXD072158 - EXISTS.** Verified via PRIDE REST v2: title "Kupffer cells exposed to 12 weeks of microplastic (500nm), nanoplastic (50nm) or PBS control."
- **Metabolomics Workbench ST003776 - EXISTS.** Verified via MW REST summary: "Investigation of Hepatic Lipid Alterations Following Micro- and Nanoplastic-Ingestion", Mus musculus, University of Bonn, analysis type MS(Dir. Inf.), **28 samples**, submitted 2025-03-06, released 2025-03-13, CC BY-NC-ND. Note 28 deposited samples vs 27 animals in the Fig. 1j lipidomics legend (9+8+10) - a one-sample discrepancy worth resolving from the deposit's sample table.
- Histology: bonndata image deposit, reference 54, DOI 10.60507/FK2/M42B8O (not independently resolved in this session - UNVERIFIED).
- **No Code availability statement and no GitHub/Zenodo link anywhere** (confirmed by full-text search). Software versions are listed in Methods instead.
- Source Data provided for Figs. 1-4 and ED Figs. 2-10.
- Re-used data: none. Niche-factor gene sets taken from ref. 49; DoRothEA/MSigDB are reference resources, not re-analysed datasets.

## 1.4 Statistical and computational models
**Phenotypes (Methods, "Statistics", verbatim-sourced):** GraphPad Prism 10; Shapiro-Wilk normality; one-way ANOVA + FDR or unpaired t-test for normal data; Kruskal-Wallis + FDR for non-normal; "When analysing data with two independent factors (for example, condition and time), a two-way ANOVA was conducted, followed by FDR correction"; Wilcoxon + FDR for ranked data; mean +/- s.d.; "Outliers were identified using Grubbs' test (alpha = 0.05) and removed."
- Two-way ANOVA appears only as **condition x TIME** (body-weight curves ED Fig. 4a, Fig. 2b; GTT/ITT curves Fig. 2e,g; adipocyte size distribution ED Fig. 4g; DBI x lipid class ED Fig. 6i; washout qPCR ED Fig. 7g), always with Geisser-Greenhouse correction. **Never as particle x diet.**
- Mixed-effects analysis used once (Fig. 3l, platelet time course).
- **ANCOVA** (regression-based, body weight as covariate for oxygen consumption) in ED Fig. 8h,i - the only covariate-adjusted model in the paper. *(This was missing from the prior notes.)*

**snRNA-seq (Methods, "snRNA-seq preparation and analysis"):** 16,000 nuclei loaded/sample, 10x 3' v3.1 Dual Index, NextSeq 2000 P4, ~30,000 reads/nucleus; CellRanger v7.1.0; Seurat v5.0.1; QC >=200 genes/nucleus, genes in >=3 nuclei, exclude <600 features or >30% mito; 20 PCs, res 0.5 then re-cluster at 0.2 after DoubletFinder v2.0.4 (assumed doublet rate 25%, pK 0.001); KC sub-clustering 20 PCs res 0.3.
- **DE: "Differential expression and ORA compared MP and NP with Ctrl using a Wilcoxon test"**, two-sided Wilcoxon rank-sum + Benjamini-Hochberg (Fig. 3h legend, Supplementary Table 1 legend). Gene filtering: per cell type drop genes outside the 0.5th-99.5th expression percentile; per sample keep genes expressed in >=1% of cells; union tested.
- **Threshold discrepancy CONFIRMED:** Methods = "adjusted P <= 0.05 with log2(FC) >= 0.75 or <= -0.75 in hepatocytes ... and adjusted P <= 0.1 with log2(FC) >= 0.25 or <= -0.25 in all other cell types"; **Fig. 1g legend = "log2(FC) above 0.75 or below -0.75 and q < 0.1"**. The hepatocyte q threshold differs between Methods (0.05) and figure legend (0.1).
- **Pseudobulk vs cell-level: DE is cell(nucleus)-level Wilcoxon; only VIPER is pseudobulked** ("Normalized expression was pseudobulked by sample and condition"). With n = 2 mice/condition and nuclei as units, hepatocyte/KC DE p-values are pseudoreplicated - this is the single biggest statistical caveat in the paper and the authors do not address it.
- ORA: clusterProfiler v4.10.0 on Hallmark/Reactome/GO-BP via msigdbr v7.5.1; GO-BP restricted to the first four hierarchy levels; gene sets 10-300 genes; FDR <= 10% with >= 3 genes; one-sided hypergeometric + BH (Supplementary Table 2 legend).
- Module scores: Seurat AddModuleScore on nine curated MSigDB lipid modules, averaged per mouse, reported as **delta vs PBS with no significance test** (Fig. 1k).
- VIPER v1.36.0 + DoRothEA mouse regulons v1.14.1 on pseudobulk; Wilcoxon + BH. Authors state plainly: **"Although no single transcription factor reached significance"**.
- CellChat v1.6.0, default parameters, genes overexpressed if detected in >=10% of cells, >=10 cells/group.
- **Proteomics:** DIA-NN v1.8.1 library-free (SwissProt mouse 2022-12-28), precursor FDR <1%, match-between-runs, MaxLFQ; Perseus v1.6.15; 3 samples excluded (QC); log2, 70% completeness in >=1 condition; imputation from a normal distribution downshifted 1.8 s.d., width 0.3 s.d.; **Welch's t-test with S0 = 0.1, FDR < 0.05, 250 randomizations**.
- **No WGCNA. No mixed model or batch/cage covariate anywhere in the omics.** Lipidomics analysed by Kruskal-Wallis + Dunn's (Fig. 1j, ED Fig. 6h).

## 1.5 Does the model suit the question?
- **Size-dependence (the paper's core claim) is well served**: all three arms (Ctrl/NP/MP) sit in one cohort and one omnibus test (one-way ANOVA or Kruskal-Wallis across three groups) with FDR-corrected post hoc, so "MP but not NP" is a within-cohort contrast, not a cross-experiment comparison. This is better than most of the set.
- **The HFD "exacerbation" claim is NOT supported by an interaction test.** The Abstract says "Under diet-induced obesity, microplastics exacerbate hepatic lipid accumulation". The HFD cohort is analysed within itself across Ctrl/NP/MP (Fig. 2i, one-way ANOVA/KW + FDR), and the chow cohort is a separate experiment with a different timeline, different ages at start (6-8 wk vs 8-12 wk) and MNP dosing that **stops 6 weeks before sacrifice** in the HFD arm. There is no model with diet and particle as crossed factors anywhere. "Exacerbates" is therefore a comparison of two separately-significant pairwise results.
- **Pseudoreplication in the snRNA-seq**: nucleus-level Wilcoxon at n = 2 mice/group means the DEG counts in Fig. 1f and the volcano in Fig. 3h are not valid biological-replicate inference. VIPER (correctly pseudobulked to n = 2) duly returns nothing significant - which is the honest answer and the authors report it.
- **Module scores with no test** (Fig. 1k) carry the "suppressed lipogenesis / increased uptake" narrative; they are descriptive deltas.
- Conversely, the **causal core of the paper does not rely on omics at all**: the co-culture, platelet-clearance, E. coli and Clec4f-DTR experiments are proper in vivo/ex vivo functional tests with reasonable n.

## 1.6 Inputs, features, outputs, targets
| Analysis | Input | Features | Output | Target |
|---|---|---|---|---|
| snRNA-seq | liver nuclei, 6 libraries (2/condition) | ~gene UMI counts, 9 cell types | DEG lists, ORA terms, module scores, UMAP | cell-type-specific transcriptional response to MNP size |
| VIPER | pseudobulk per sample per cell type | DoRothEA regulons | TF activity scores (n.s.) | upstream regulators |
| CellChat | all cell types | ligand-receptor pairs | information flow per pathway | KC->hepatocyte signalling |
| Shotgun lipidomics | 10 mg liver, internal standards, Q Exactive Plus, LipidXplorer | pmol per lipid species/class, double-bond index | class- and species-level log2FC | hepatic lipid composition |
| DIA proteomics | CD11b-enriched KCs | MaxLFQ protein intensities | differentially abundant proteins | KC proteome, validation of RNA predictions |
| Phenotype | mice | ORO area, PAS intensity, platelet %, CFU, VO2, NMR body composition, serum ELISAs/Luminex | group means | liver lipid/glycogen, phagocytosis, energy balance |

## 1.7 How results were evaluated
Thresholds as in 1.4. Evaluation is largely by **convergence across modalities** rather than by any held-out or cross-validated metric: lipidomics (direct measurement) is used to corroborate the transcriptional lipid programme; PAS glycogen corroborates the MP-specific "glycolysis" ORA term; proteomics corroborates the apolipoprotein transcript signal (APOC3, TTR up at protein level). No external dataset, no cross-validation, no permutation/robustness analysis of the snRNA-seq, no power calculation. Grubbs' outlier removal was applied to phenotype data with no per-panel accounting in the readable text (Source Data would show it).

## 1.8 Independent validation vs internal only
**VALIDATED (orthogonal/functional):**
- KCs are the hepatic MNP reservoir: flow cytometry + confocal IF + gold-core STEM/serial block-face EM + Supplementary Video 1.
- Increased hepatic neutral lipid with MP only at 12 wk (ORO, Fig. 1h,i); no difference at 4 wk (ED Fig. 5g,h).
- Decreased glycogen (PAS) in MP at 100 and 10 mg/kg, not at 1 mg/kg (ED Fig. 5c,d; 6a-f) - dose threshold.
- No injury/fibrosis: serum AST/ALT ELISA and Sirius Red negative (ED Fig. 5e,f).
- Lipidomics class/species shifts (direct MS measurement).
- **Paracrine sufficiency (strongest causal result):** KCs from MP mice + naive hepatocytes -> increased hepatocyte lipid droplet area (Fig. 3b,c; two-sided Mann-Whitney; n = 3 vs 4 mice). The mediator was **not** identified - no neutralisation or knockdown of APOC3/TTR/any ligand.
- **Phagocytosis impairment:** MP (not NP) reduced Fc-gamma-receptor-dependent clearance of opsonised platelets at 4 wk (Fig. 3j) and 12 wk (Fig. 3l) - i.e. *before* lipid accumulates.
- **Bacterial clearance:** MP mice cleared i.v. E. coli from blood less well (Fig. 3k); liver/spleen CFU unchanged (ED Fig. 9f).
- **No inflammation:** hepatic myeloid counts and 48-plex serum cytokines unchanged (ED Fig. 9c,d); blinded gWAT/BAT pathology.
- **Genetic:** KC ontogeny unchanged (double fate-mapping, Fig. 4b; Cxcr4-CreERT2, ED Fig. 10a,b); Clec4f-DTR depletion shows moKCs fail to restore phagocytosis (Fig. 4d) and are fewer at 12 wk post-DT (Fig. 4e-g) while influx/proliferation/ontogeny/niche factors are unaffected (ED Fig. 10g-j).
- Proteomic validation of APOC3 and TTR (Welch, FDR<0.05). *(Apoa1/Apoa2 protein are NOT stated as validated.)*
**OMICS-ONLY:** hepatocyte DEG counts and ORA terms; suppressed de novo lipogenesis + increased uptake programme; MP-specific VLDL-machinery and Mgll downregulation; bile-acid cascade up (Cyp7a1, Cyp8b1, Cyp27a1) with Nr0b2 down; Nr1h3/LXRa down; Lands-cycle/ceramide enzyme changes; KC ORA terms; KC cluster shift; apolipoprotein transcripts.
**INFERENCE:** all VIPER TF activities (explicitly non-significant); CellChat VCAM/OCLN/AGT/SEMA4 flow; the HIF-1alpha -> apolipoprotein -> hepatocyte-lipid axis (HIF-1alpha was never manipulated; the authors write "raising the possibility").

## 1.9 Main findings with directions
- **Lipogenesis genes (Fasn/Acly/Acaca/Srebf1): NOT NAMED ANYWHERE IN THE READABLE TEXT.** Full-text grep returns 0 hits for Fasn, Srebf, Acly, Acaca. The claim is at module/pathway level: "both particle sizes engaged a common programme of **suppressed de novo lipogenesis** together with **increased fatty acid uptake and lipoprotein-receptor expression**" (Results). Gene-level detail lives in Fig. 1l and Supplementary Table 1 - **UNVERIFIED**.
- **PPARalpha/FAO: no transcript or protein data at all.** Grep: 0 hits for Ppara, Cpt1a, Acox, Cyp4a. PPARalpha appears **only as VIPER-inferred activity** - engaged in NP hepatocytes (with HNF4a, C/EBPb) and shared in KCs of both conditions; PPARgamma NP-selective in KCs. Non-significant. No FAO flux measured.
- **Cholesterol:** "cholesterol homoeostasis" is a shared enriched hepatocyte GO term (direction of the term not stated); bile-acid synthesis genes UP with repressor Nr0b2 DOWN; **Nr1h3 (LXRa) DOWN**, which the authors use to *exclude* LXRa as the driver. Serum total cholesterol, HDL-C and TG unchanged in the HFD arm (Supplementary Fig. 8 - **UNVERIFIED numerically**).
- **Inflammation: ABSENT and this is a deliberate, well-supported negative.** No increase in hepatic KCs/monocytes/neutrophils/cDC1/cDC2, no serum cytokine elevation (IL-1b, IL-12p70, G-CSF, Ccl7, IL-18), no gWAT crown-like structures. A pro-inflammatory signature exists only as inferred TF activity in MP (STAT1, IRF2, E2F4 in hepatocytes; MYC, PU.1, IRF1 with loss of MAFB in KCs). Authors: MP "reprograms KC core functions rather than triggering a classical inflammatory response".
- **Oxidative stress: NOT MEASURED** (no ROS/MDA/SOD/GSH/Nrf2/Hmox1 anywhere). Only the inferred loss of MAFG/MAFK activity. Do not cite this paper for MNP oxidative stress.
- **Fibrosis: absent** at 12 wk and after washout (Sirius Red; fibrosis-gene RT-qPCR after washout).
- **Lipidomics:** total lipid unchanged; MP increased LPE, LPC and sphingomyelin classes; at species level "in NP, none reached significance, whereas in MP five LPC species, seven LPE species and two monogylceride species were significantly elevated". DBI increased for TG in MP and for DiHexCer/alkyl-PC/TG in NP.
- **HFD arm:** MP INCREASED and NP DECREASED hepatic ORO vs Ctrl (Fig. 2h,i); NP reduced weight gain, gWAT and BAT mass; fasting glucose higher in both NP and MP after HFD; GTT AUC increased in NP only; ITT unchanged; NP increased O2 consumption at 22 C and 4 C with reduced fat mass, no change in food/water intake or locomotion; lean mass selectively increased in NP.
- **Washout:** particles persist in KCs after 6 months but the metabolic phenotype resolves.
- **Circadian: NOT ADDRESSED.** Grep: 0 hits for circadian / Zeitgeber / time of day. Housing described only as "a 12 h light-dark cycle"; no time-of-sacrifice is stated anywhere, although PAS glycogen and bile-acid/lipogenic transcripts are strongly diurnal. GTT/ITT used 6 h fasting; no fasting statement for terminal tissue collection.

## 1.10 Limitations
**(a) ACTUAL limitations of what was done**
1. **snRNA-seq n = 2 per condition with nucleus-level Wilcoxon DE** - pseudoreplicated; the entire transcriptional narrative rests on this. VIPER, correctly pseudobulked, finds nothing significant.
2. **Threshold inconsistency** between Methods (hepatocyte q <= 0.05) and Fig. 1g legend (q < 0.1) - the reported ORA cannot be reproduced without knowing which applies.
3. **Gold-core dosing contradiction** between Methods (0.8 mg/week, 5 weeks) and Fig. 1d legend (0.8 mg/kg twice weekly, alternating).
4. **Grubbs' outlier removal applied to phenotype data** with no per-panel reporting; this is a real analytic choice affecting significance, and it is not stated to be pre-specified.
5. **No randomisation statement; blinding only for histopathology.**
6. **Proteomics n never reported**, and 3 samples were excluded post hoc on a QC criterion.
7. **MNP tissue burden could not be quantified** - Europium signal lost (~3 orders of magnitude), pyrolysis-GC-MS failed ("samples from MNP-exposed mice yielded quantitative values that were comparable to or even lower than those obtained from Ctrl tissues"; ED Fig. 3e is n = 1). So dose-at-target is unknown.
8. **The HFD arm is not a crossed design** and dosing stopped 6 weeks before sacrifice; the "exacerbation" claim is cross-cohort/pairwise.
9. **Module scores reported without significance testing** (Fig. 1k).
10. **No endotoxin, surfactant/preservative or leachate testing** (0 hits for endotoxin/LPS/surfactant); control is PBS gavage, not a formulation-matched vehicle. Commercial carboxylated PS suspensions commonly contain surfactant/azide.
11. **n varies substantially between panels drawn from the same cohort** (ORO 6/10/9 vs PAS 9/9/10 vs organ weights 19/24/23) with no accounting of which animals contributed where.
12. **Co-culture sex mismatch**: hepatocytes from 5-week-old females, all in vivo work male.
13. **KCs for proteomics and co-culture were CD11b magnetic-bead enriched**, not sorted Clec4f+/Vsig4+; purity is not reported, so contaminating myeloid cells could contribute to the co-culture effect.
14. **"Pooled from >2 independent experiments"** with no experiment/batch term in any model; cage is never mentioned as a unit.
15. **No time-of-day control** for diurnal readouts (glycogen, lipid genes, bile acids).
16. Author-acknowledged: chow-fed/short HFD models "early, mild metabolic perturbation rather than advanced steatohepatitis"; fasting insulin and NEFA not measured; mainly young male mice; mechanism not reduced to a single ligand; transcript-lipid associations "are correlative and we do not infer causation"; other nuclear-receptor inputs "we did not test functionally"; LSEC contribution "remains largely unexplored"; "behavioural experiments ... had to be removed during the revision" (Author contributions).

**(b) Simply not attempted (not flaws, but bounds on use)**
- No MASH/NASH diet (no fructose, cholesterol or CDAHFD), no fibrosis model, no longitudinal progression study.
- No female cohort; no aged mice.
- No polymer other than polystyrene; no environmentally weathered or mixed-polymer particles; no biocorona characterisation.
- No HIF-1alpha, PPARalpha, APOC3 or TTR manipulation - the proposed axis is never tested genetically.
- No tracer/flux lipid studies (the authors name this as the needed experiment).
- No oxidative-stress, mitochondrial-function or autophagy measurement in liver.
- No gut microbiome, no ABX/FMT, no intestinal barrier work.
- No circadian sampling; no human data or external dataset integration.

## 1.11 Study metadata
- **Species/strain/sex/age:** mouse, C57BL/6JRcc wild type (transgenics on C57BL/6J crossed to C57BL/6Jrcc for 2-5 generations). "Young (8-12-week-old) **male** mice were used for all experiments, with the exception of experiments with Cxcr4CreERT2/+; Rosa26LSL-tdTomato/+ mice, in which both sexes were used." HFD arm: male C57BL/6JRcc, **6-8 weeks** at diet start. Co-culture hepatocytes: 5-week-old **female**.
- **Polymer/size/supplier:** carboxylated (negatively charged) polystyrene only. NP = 50 nm, MP = 500 nm; bioavailability screen also 100 and 300 nm. Yellow-Green Fluorescent Carboxylated PS, Polysciences 15700-10, 17149-10, 16662-10, 2405-10; Red Fluorescent Carboxylated PS, Megasphere CAF050NM/CAF100NM/CAF200NM/CAF500NM; unlabelled PS Polysciences 15913-10, 09836-15; gold-core PS-coated carboxylated CD Bioparticles CPGN-500 / CPGN-50; europium MP Abvigen ABPST-1-0500-EuC.
- **Characterisation:** DLS on Zetasizer Ultra, Z-average + PDI (10 measurements), zeta potential 5x per sample; **values only in Supplementary Fig. 1 - UNVERIFIED**. No endotoxin, no leachate.
- **Dose/route/vehicle/frequency/duration:** 100 mg/kg body weight standard, titration 10 and 1 mg/kg; oral gavage; vehicle PBS (also the control gavage); once weekly; 12 weeks (4 wk, 5 wk gold-core, 24 h acute, +6-month washout).
- **Diet (catalog numbers given):** CD "13 kJ% fat; Ssniff, E15748-047"; HFD "60 kJ% fat; Ssniff, E15742-347". Units are **kJ%**, not kcal%.
- **Sampling/fasting:** GTT/ITT after 6 h fasting; glucose 2 mg/g i.p.; insulin 0.75 U/kg i.p. No fasting state or time of day stated for terminal collection.
- **Housing:** specific-pathogen-free, 12 h light-dark, food and water ad libitum; single housing only during the 3-day metabolic-cage acclimation. Mice per cage, cages per group and co-housing: **not reported**.
- **Ethics:** Directive 2010/63/EU; LANUV North Rhine-Westphalia Az 81-02.04.2019.A167 and AZ 81-02.04.2023.A255.
- **Competing interests:** none declared. Authors disclose using ChatGPT and Claude "for language correction and proofreading".

## 1.12 Prior-note reconciliation (file `natmetab2026_size_macrophage.md`)
Every load-bearing claim in the prior note was **CONFIRMED verbatim**, including: n = 2 snRNA-seq; the Methods/Fig. 1g threshold conflict; the gold-core dose conflict; Grubbs' removal; the absence of a particle x diet interaction; the "no TF reached significance" quote; the VLDL/Mgll MP-selectivity; the proteomics exclusions; all author-acknowledged limitations; and the full Data availability wording.
**Additions/corrections from this pass:**
1. **New:** ED Fig. 8h,i use a **regression-based ANCOVA with body weight as covariate** for oxygen consumption - the paper's only covariate-adjusted model; the prior note said no covariate was used anywhere, which is true for the omics but not for VO2.
2. **New:** all three accessions now **independently verified to resolve**, with GEO confirming exactly 6 libraries (2/condition) and MW listing **28** samples against 27 animals in the lipidomics legend.
3. **New:** ED Fig. 3e pyrolysis quantification is **n = 1 biological replicate** - the failed-method evidence is weaker than the prior note implied.
4. **Confirmed precisely:** prior note's "Fig. 1g q < 0.1 vs Methods 0.05" wording is exactly right.
5. Prior note listed body-weight n as 35/32/39 - **confirmed** (ED Fig. 4a).
6. Prior note's "no endotoxin" - **confirmed by grep (0 hits for endotoxin, LPS, surfactant)**.

---

# 2. J Nanobiotechnology 2026 - PS nanoplastics induce hepatic steatosis via NCoR1 / PPARalpha

Fang Y, Wei T, Zhang Z, Li Y, Yang M, Chen H, Wang Y, Hao R, Li Y, Sun Z, Wang J (corresponding). "Polystyrene nanoplastics induce hepatic steatosis by disrupting autophagic degradation of NCoR1 and suppressing PPARalpha-mediated fatty acid oxidation." *J Nanobiotechnology* 2026;24(1):437. DOI https://doi.org/10.1186/s12951-026-04304-6 . PMID 41928286, PMCID PMC13173803. Received 2025-10-28, accepted 2026-03-11.

**Access level: FULL TEXT read** (Europe PMC full-text XML -> local text, 67,294 characters of body): abstract, introduction, all Results subsections, all Fig. 1-8 legends, Table 1, all Experimental methods subsections, Discussion, Conclusion, Author contributions, Funding, Data availability, Declarations. **NOT read:** Supplementary Material 1 and 2 = Figures S1-S5 and Tables S1 (antibodies), S2 (primers). Several load-bearing results live only there: DEG volcano/counts (S2), metabolomics PCA/OPLS-DA and DEM lists (S3), the **entire lipogenesis/uptake qPCR panel (S4)**, the **entire rapamycin rescue (S5)**, and serum inflammatory factors + fibrosis staining (S1).

## 2.1 Exact scientific question
Does oral PSNP exposure exacerbate HFD-induced hepatic steatosis, and is the mechanism impaired autophagic degradation of the corepressor NCoR1 leading to suppression of PPARalpha-driven fatty acid beta-oxidation? Stated hypothesis (end of Introduction): "we hypothesized that PSNPs exacerbate lipid accumulation in HFD-fed mice."

## 2.2 Experimental design, factors, n per group (phenotype vs omics)
**True 2x2 factorial**: diet (ND vs HFD) x particle (ddH2O vs 150 mg/kg.bw PSNP), 90 days, **n = 10 per group**. Verbatim (Methods, Animals and experimental design): "Mice were randomly allocated into four experimental groups (n = 10 per group): ND: normal diet + daily oral gavage of equal volumes of ultrapure water (ddH2O); PSNPs: normal diet + 150 mg/kg.bw PSNPs; HFD: high-fat diet (D12492, Research Diets Inc., USA; 459 kcal/100 g, 60% fat, 20% protein, 20% carbohydrate) + ddH2O; HFD+PSNPs: high-fat diet + 150 mg/kg.bw PSNPs."
- Randomisation: stated ("randomly allocated"; also "two mice per group were randomly selected" for TEM, "Six randomly selected liver specimens per group" for lipid assays). **Blinding: never mentioned** (grep: 0 hits), including for NAS scoring, which is subjective.
- Separate biodistribution arm: single gavage 150 mg/kg.bw red-fluorescent PSNPs; imaging at 0, 2, 4, 6, 8, 12, 24 h; organs at 0 and 24 h. **Group sizes not reported.**
- In vitro: human L-02 hepatocytes, 100 uM oleic acid 24 h then 25 ug/mL PSNP 24 h; rapamycin arm 100 nM + 25 ug/mL PSNP 24 h. Not a 2x2.

**n by assay:**
- Body weight, food intake, fasting glucose, liver weight/ratio: **n = 10** (Figs. 2, 3a).
- Histology/NAS, ORO quantification, hepatic ALT/AST, hepatic TG/TC/FFA/FC, RT-qPCR, most westerns: **n = 6** (Fig. 3 legends; Methods "cryosections ... from six animals per group"; "Six liver samples per group were homogenized").
- **Transcriptomics: n = 4** ("hepatic transcriptomic profiling was conducted on four liver specimens"; Fig. 4a "(n = 4)").
- **Metabolomics: n = 6 per group** ("Metabolomic profiling of hepatic tissue from HFD and HFD+PSNP groups (n = 6 per group)").
- TEM: **n = 2** per group, descriptive.
- NCoR1 IHC, LC3B IF, NCoR1 nuclear IF, NCoR1 phospho-IP: **n = 3** (Figs. 5e, 6b, 7a, 7c).
- qPCR "triplicate for each biological sample (n = 3 technical replicates)" - technical, not biological.

**CRITICAL DESIGN ASYMMETRY (confirmed):** the omics cover only **two of the four cells** - HFD vs HFD+PSNP. Metabolomics says so explicitly; transcriptomics reports DEGs only "in the HFD+PSNP group compared to the HFD group". The ND and ND+PSNP arms have no transcriptome or metabolome, so **no diet x particle model is possible with these data at all.**

## 2.3 Datasets, accessions, deposition
Data availability verbatim: **"All data are available from the corresponding author upon reasonable request."**
- **NO accession of any kind.** Full-text search returns 0 hits for GSE, PRJNA, SRP, MTBLS, PXD, OMIX, CRA, Zenodo, figshare, GitHub. Neither the RNA-seq nor the metabolomics (both generated commercially by OE Biotech, Shanghai) has a public identifier. Nothing to resolve; the datasets are **effectively not deposited**.
- Supplementary Material 1 and 2 linked via the DOI; Figures S1-S5, Tables S1-S2.
- No re-used external dataset. No code availability.

## 2.4 Statistical and computational models
- **RNA-seq pipeline (Methods, Transcriptomic analysis):** OE Biotech, Illumina NovaSeq 6000, 150 bp PE; fastp **v1.0**; HISAT2; FPKM + HTSeq-count; **DESeq2, "adjusted p-value < 0.05 and |log2 Fold change| > 1"**. Result: 236 DEGs (95 up, 141 down) HFD+PSNP vs HFD. **No sequencing depth, no read counts, no mapping rate, no library-QC or batch statement anywhere** - a notable omission.
- **Metabolomics:** UHPLC-QTOF-MS/MS + GC-MS; MS-DIAL; HMDB/LipidMaps v2.3/METLIN/LuMet-Animal 3.0/LuMet-GC 5.0; PCA then OPLS-DA; **DEMs at adjusted p < 0.05, VIP > 1, |log2FC| > 1**; 102 DEMs (51 up, 51 down).
- **Enrichment: KEGG of DEGs (Fig. 4b). The enrichment TOOL and the multiple-testing correction are NOT reported** - the legend says only "color indicates the p-value". No GSEA, no Reactome, no WGCNA, no clustering beyond heatmaps, no single-cell, no ML.
- **Phenotype statistics (Methods, Statistical analysis, verbatim):** "Data are presented as mean +/- standard error of the mean (SEM). Statistical comparisons were performed using GraphPad Prism 9 software. Two-way analysis of variance (ANOVA) was applied as appropriate. Differences were considered statistically significant at p < 0.05." Every figure legend says "two-way ANOVA with a Tukey multiple comparison test" - and every legend says **mean +/- SD**. **CONFIRMED SD/SEM CONTRADICTION** between Methods and all legends.
- **INTERACTION: a two-way ANOVA is fitted, so an interaction term exists, but NO interaction F or p is reported anywhere in the paper.** All reported significance is pairwise and the legends define it explicitly: "'#' indicates statistically significant difference versus the ND control group ... '*' indicates a statistically significant difference versus the HFD control group".

## 2.5 Does the model suit the question?
- The **design** is the best in this set for an interaction question (proper 2x2, n = 10, matched vehicle gavage in controls so gavage stress is balanced). The **analysis** then throws that away: only Tukey pairwise contrasts vs ND and vs HFD are shown, so the Discussion's claim of "a synergistic amplification of hepatic lipid accumulation when HFD was combined with oral PSNP exposure, exceeding the effects of either factor alone" is asserted, not tested. The interaction F/p may exist in the authors' Prism output; it is not in the paper.
- The **omics cannot address the question at all** (two cells of four). Any statement about PSNP effects on a normal-diet background is phenotype-only.
- The strongest *honest* statement the data support: PSNPs changed the liver phenotype on the HFD background and not on ND ("this is the first study reporting that orally administered PSNPs exacerbate HFD-induced hepatic accumulation without phenotypically affecting ND-fed mice", Discussion) - consistent with, but not proof of, an interaction.
- **Mechanistic inference runs ahead of the measurements:** PPARalpha itself (Ppara transcript, PPARalpha protein) is **never measured**; the NCoR1-PPARalpha interaction central to the model is **never measured** (no co-IP, no ChIP on FAO promoters); "S6K2" in the mechanism is inferred from p-S6, which cannot distinguish S6K1 from S6K2; NCoR1 phosphorylation is detected with a pan-phospho antibody and no site mapped.
- The only causal manipulation is **rapamycin in vitro**, in Figure S5. Rapamycin both inhibits mTOR and induces autophagy, so it cannot separate the two arms of the proposed mechanism.

## 2.6 Inputs, features, outputs, targets
| Analysis | Input | Features | Output | Target |
|---|---|---|---|---|
| Bulk RNA-seq | liver, n = 4, HFD vs HFD+PSNP only | FPKM/counts | 236 DEGs; KEGG terms | PPAR pathway suppression |
| Untargeted metabolomics | liver, n = 6, same two groups | annotated metabolites | 102 DEMs; lipid bubble plot | LCFA/acylcarnitine accumulation |
| RT-qPCR | liver, n = 6 | Acsl3, Cpt1a, Acox1, Acadvl, Hmgcs1, Lc3b (+ S4 panel) | 2^-ddCt vs Gapdh | PPARa target validation |
| Western | liver and L-02, n = 6 | ACSL3, CPT1A, ACADVL, HMGCS1, NCoR1, p62, LC3B, mTOR/p-mTOR, S6/p-S6 | densitometry vs GAPDH | protein-level validation |
| IP | 30 mg liver, n = 3 | NCoR1 pulldown, pan phospho-(Ser/Thr) | band intensity | NCoR1 phosphorylation |
| IF/IHC | liver and L-02, n = 3 | LC3B area, NCoR1 nuclear signal, NCoR1/LC3B/LAMP1 colocalisation, BODIPY | quantified area/colocalisation | autophagy and NCoR1 localisation |

## 2.7 How results were evaluated
DESeq2 padj < 0.05 / |log2FC| > 1; DEM VIP > 1 + padj < 0.05 + |log2FC| > 1; phenotype p < 0.05 by two-way ANOVA + Tukey. Evaluation is by **cross-modality agreement** (transcriptome -> metabolome -> qPCR -> western -> IHC/IF), which is genuinely the paper's strength: the PPARalpha-target direction is reproduced at mRNA and protein with n = 6. There is no held-out validation, no power calculation, no replication cohort, no external dataset, and no correction across the many phenotype endpoints.

## 2.8 Independent validation vs omics-only
**VALIDATED:** hepatic PSNP accumulation (in vivo/ex vivo fluorescence imaging, labelled particles; TEM n = 2); unchanged food intake, transient weight gain, final body weight unchanged (n = 10); fasting glucose up at 90 d (vs ND); liver weight, liver/body ratio, hepatic ALT/AST unchanged (n = 6-10); steatosis worse in HFD+PSNP (H&E + ORO, n = 6); hepatic TG up, TC and FFA unchanged, FC up in HFD (enzymatic, n = 6); **Acsl3, Cpt1a, Acox1, Hmgcs1 mRNA and ACSL3, CPT1A, HMGCS1 protein all DOWN** (qPCR + western, n = 6); **NCoR1 protein UP** (western n = 6 + IHC n = 3); LC3B protein and LC3B+ area DOWN with Lc3b mRNA unchanged (IF n = 3, qPCR + western n = 6); NCoR1 phosphorylation UP (IP-western n = 3); p-mTOR and p-S6 UP (n = 6); NCoR1 nuclear accumulation UP (IF n = 3).
**OMICS-ONLY:** the 236 DEGs; KEGG "PPAR signaling pathway", "steroid biosynthesis", "terpenoid backbone biosynthesis" suppression; Fabp4, Fabp5, Cyp4a31 downregulation; the 102 DEMs and the LCFA/acylcarnitine accumulation reading.
**CORRELATIVE IMAGING ONLY:** NCoR1 colocalisation with LC3B and LAMP1 in L-02 cells (Fig. 6e) - this is colocalisation, not a degradation assay. No cycloheximide chase, no bafilomycin/chloroquine flux comparison, no Atg7 manipulation in this paper.
**CAUSAL, IN VITRO ONLY:** rapamycin rescue (Figure S5a-c) - "rapamycin treatment abolished PSNP-induced nuclear accumulation of NCoR1 and alleviated intracellular lipid deposition ... rapamycin restored the expression of key fatty acid oxidation proteins downregulated by PSNPs". Called "direct pharmacological evidence" by the authors. **Supplementary, not read - UNVERIFIED**, including whether a rapamycin-alone arm exists.
**NOT DONE ANYWHERE:** no genetic manipulation of any kind (no Ncor1 knockdown/KO, no liver-specific deletion, no Ppara KO, no Atg7 or other autophagy-gene manipulation); no PPARalpha agonist challenge; no co-IP or ChIP for NCoR1-PPARalpha; no in vivo rescue; no FAO flux or isotope tracing.

## 2.9 Main findings with directions
- **Lipogenesis (Fasn/Acly/Acaca/Srebf1): reported as UNCHANGED, but the genes are never named in the readable text** (grep: 0 hits for Fasn, Srebf, Acly, Acaca, Scd1, Cd36). Verbatim: "RT-qPCR analysis of genes involved in FFA synthesis, TG synthesis, and FFA uptake revealed no significant alterations, suggesting that PSNPs selectively impair fatty acid oxidation rather than lipid synthesis or uptake pathways (Figure S4)." Gene identities **UNVERIFIED**. Note the internal tension: Fabp4 and Fabp5 (lipid binding/transport) were among the downregulated DEGs.
- **PPARalpha/FAO: SUPPRESSED - the paper's central claim and the clearest direction in this set.** Down at mRNA: Acsl3, Cpt1a, Acox1, Hmgcs1 (n = 6); down at protein: ACSL3, CPT1A, HMGCS1 (n = 6); down in the DEG set: Cyp4a31, Fabp4, Fabp5. **Acadvl/ACADVL was assayed (Fig. 5b,c) but no change is stated - treat as not clearly altered.** KEGG "PPAR signaling pathway" suppressed. Metabolomics shows LCFA and acylcarnitine accumulation, read as a beta-oxidation block (inference, no flux measured). **Ppara itself was never measured.** The authors stress the sign explicitly: "Given the conserved nature of PPARalpha signaling and the observed suppression (rather than activation) of this pathway, our results identify a novel mechanism".
- **Cholesterol synthesis: DOWN.** "steroid biosynthesis and terpenoid backbone biosynthesis" among the suppressed KEGG enrichments; **HMGCS1 protein significantly decreased** (validated). Meanwhile hepatic TC unchanged and FC up in HFD groups vs ND. Caveat: the paper counts **HMGCS1** (cytosolic, mevalonate/cholesterol) as a PPARalpha/FAO readout; the ketogenic PPARalpha target is **HMGCS2** - primer/antibody isoform specificity is in Tables S1/S2, **UNVERIFIED**.
- **Inflammation: NO CHANGE, by design.** NAS criteria for NASH not met in any group and "no inflammation was observed among groups"; "neither serum inflammatory factors nor tissue fibrosis staining showed significant alterations (Figure S1)"; hepatic ALT/AST unchanged. No Tnf/Il1b/Il6/NLRP3/macrophage-marker data. Framed via the multiple-hit model: PSNP+HFD is a potent "first hit" (steatosis) without the second (inflammatory) hit.
- **Fibrosis: NO CHANGE** (Figure S1, **UNVERIFIED**).
- **Oxidative stress: NOT STUDIED AT ALL** - no ROS, MDA, SOD, GSH, Nrf2 or Hmox1 anywhere, in vivo or in vitro. Do not cite this paper for MNP oxidative stress.
- **Mitochondrial function: inferred only** - no membrane potential, ATP, respirometry or mitochondrial morphology scoring (TEM was used to find particles).
- **Autophagy: SUPPRESSED** (LC3B protein and LC3B+ area down; Lc3b mRNA unchanged => post-translational). **p62 was blotted (Fig. 6d) but its direction is never stated in the text** - and p62 direction is exactly what distinguishes blocked flux from reduced initiation. No flux assay (bafilomycin/chloroquine, tandem mCherry-GFP-LC3).
- **mTOR: ACTIVATED** (p-mTOR up, p-S6 up).
- **PSNPs alone (ND+PSNP) produced no hepatic phenotype** - the effect is diet-dependent.
- **Circadian: NOT STUDIED.** 0 hits for circadian/Zeitgeber. The only timing statement is "mice were fasted **overnight** and anesthetized" before collection - note that overnight fasting in mice by itself strongly induces hepatic PPARalpha and FAO genes, so all groups sit in a fasted, PPARalpha-primed state.

## 2.10 Limitations
**(a) ACTUAL limitations of what was done**
1. **Omics asymmetry** - transcriptome and metabolome cover only HFD and HFD+PSNP, so the 2x2 design is not realised in the molecular data and no interaction is computable from them.
2. **No interaction statistic reported** despite a two-way ANOVA model; "synergy" rests on the HFD+PSNP vs HFD pairwise contrast.
3. **RNA-seq QC entirely missing** - no depth, no mapping rate, no batch/library-prep statement, n = 4 selected from n = 10 with no stated rule.
4. **SEM (Methods) vs SD (every legend)** - unresolved; changes how large the effects look.
5. **No blinding** for NAS, ORO quantification or any image analysis.
6. **Mechanistic core runs on the smallest n in the study** (n = 3 for IHC/IF/phospho-IP, n = 2 for TEM).
7. **Enrichment tool and multiple-testing correction for KEGG never stated.**
8. **ALT/AST measured in liver homogenate normalised to protein, not in serum** - an unusual readout on which the "no liver injury" conclusion rests.
9. **Control diet not characterised** - "normal diet" with no brand, catalog number or composition, and no statement that it is a purified low-fat match to D12492.
10. **L-02 cell identity** - widely flagged as HeLa-contaminated in cell-line registries; no STR authentication and no mycoplasma statement (0 hits). The entire in vitro mechanism including the rapamycin rescue rests on this line.
11. **Only causal experiment is in vitro, pharmacological, and supplementary**; rapamycin cannot separate "mTOR/NCoR1-phosphorylation" from "autophagy induction".
12. **p62 direction unreported**; no autophagic-flux assay.
13. **No endotoxin, surfactant, residual monomer or leachate testing** for a 90-day, 150 mg/kg/day oral exposure (0 hits for endotoxin).
14. **Dose framing overstated** - Methods derive 150 mg/kg.bw from the **upper range** of estimated human intake (0.1-5 g/week) with a 12.3x interspecies factor, while the Conclusion calls it "environmentally relevant concentrations". These are not the same claim.
15. **Overnight fast before sacrifice** pre-activates PPARalpha/FAO in all groups - a real interpretive confound for a PPARalpha-suppression paper, never discussed.
16. **No time-of-day/ZT reporting; mice per cage not reported** (0 hits for cage).
17. Author-acknowledged (last Discussion paragraph, only two): parallel signalling pathways and the precise NCoR1 phosphorylation sites "require further mapping"; "the specific autophagic receptor responsible for NCoR1 recognition remains unidentified". Plus two conceded workarounds: no phospho-NCoR1-specific antibody (pan-phospho IP used); no S6K2-specific phospho-antibody (p-S6 used as surrogate).

**(b) Simply not attempted**
- No genetic validation anywhere (Ncor1, Ppara, Atg7).
- No co-IP/ChIP of the NCoR1-PPARalpha axis.
- No in vivo rescue (rapamycin or fenofibrate).
- No FAO flux, isotope tracing or respirometry.
- No lipidomics, proteomics, scRNA-seq or spatial data.
- No female mice; male-only, and sex is not flagged as a limitation.
- No microbiome, gut barrier or gut-liver axis work.
- No time-course or circadian sampling.
- No data deposition.

## 2.11 Study metadata
- **Species/strain/sex/age:** male C57BL/6J, **5-6 weeks** at start, Vital River Laboratories (Beijing). Male only. Housing: 12 h light/dark, 23 +/- 1 C, 60-70% RH, water ad libitum; 7-day acclimatisation. **Mice per cage not reported.**
- **Polymer/size/supplier:** PSNPs 50 nm, DK Nanotechnology Co. (Beijing), aqueous suspension 50 mg/mL in sterile deionized water. **Two formulations: red-fluorescent (ex/em 620/680 nm) for biodistribution only, and pristine (unlabelled) for the 90-day toxicity study** - this avoids the label confound that affects many MNP studies.
- **Measured characterisation (Table 1, read in full):** TEM primary diameter **49.24 +/- 5.31 nm** (500 particles, ImageJ). DLS in H2O: 0 h 53.35 +/- 0.99 nm, -30.83 +/- 4.06 mV, PDI 0.05; 24 h 53.14 +/- 0.32 nm, -25.47 +/- 1.12 mV, PDI 0.04. In DMEM: 0 h 52.68 +/- 0.35 nm, -30.70 +/- 1.61 mV, PDI 0.06; 24 h 58.26 +/- 0.12 nm, -22.70 +/- 2.14 mV, PDI 0.14. So stable in water, mild agglomeration and charge loss over 24 h in serum-free DMEM. Suspensions ultrasonicated 30 min at 160 W / 20 kHz. **No endotoxin, surfactant or leachate data.**
- **Dose/route/vehicle/frequency/duration:** 150 mg/kg.bw/day, oral gavage, vehicle ultrapure water (controls gavaged with equal volumes - gavage stress matched), daily for **90 consecutive days**.
- **Diet:** HFD **D12492, Research Diets Inc.; 459 kcal/100 g, 60% fat, 20% protein, 20% carbohydrate**. Control: "normal diet", **no brand or catalog number**.
- **Sampling/fasting:** fasting glucose biweekly after 8 h fast; terminal collection after **overnight fast** under anaesthesia; blood coagulated 60 min RT, 3000 rpm 10 min at 4 C. **No time of day / ZT stated.**
- **In vitro:** human L-02 ("normal liver", Nanjing Keygen Biotech), 100 uM oleic acid 24 h -> 25 ug/mL PSNP 24 h; rapamycin (MCE HY-1029) 100 nM. No authentication or mycoplasma statement.
- **Ethics:** GB/T 35892-2018; Animal Ethics Committee of Capital Medical University, Protocol AEEI-2024-013.
- **Funding:** Beijing Natural Science Foundation QY24342; Capital Medical University Incubation Fund PYZ24019, PYZ23017. **Competing interests: none declared.**

## 2.12 Prior-note reconciliation (file `jnanobiotech2026_ncor1_ppara.md`)
All load-bearing claims **CONFIRMED verbatim**, including: n = 10 / 6 / 4 / 6 / 3 / 2 by assay; the omics covering only HFD and HFD+PSNP; DESeq2 thresholds; missing enrichment tool and correction; missing RNA-seq depth; two-way ANOVA with **no interaction F/p anywhere** and the exact #/* legend definitions; the **SEM-vs-SD contradiction**; "All data are available from the corresponding author upon reasonable request" with **no accession**; the pristine-vs-fluorescent particle split; Table 1 DLS/zeta/PDI values (re-read and matched digit for digit); the D12492 specification and the uncharacterised "normal diet"; hepatic (not serum) ALT/AST; the overnight fast; the two Discussion limitations and the two antibody workarounds; the "first study ... without phenotypically affecting ND-fed mice" quote; and the Hmgcs1-vs-Hmgcs2 concern.
**Corrections/refinements from this pass:**
1. Prior note said "Transcriptomics and metabolomics compare ONLY HFD+PSNP vs HFD" - **confirmed for metabolomics explicitly in Methods**; for transcriptomics the Methods say only "four liver specimens" and the two-group restriction is established from the Results/Fig. 4c wording. Slightly weaker sourcing than the prior note implied, but the conclusion stands.
2. Prior note listed Fig. numbers 1-8: confirmed; the rapamycin rescue is cited in the Results as **Figure S5a and S5b** (imaging) and **S5c** (western), and again in the Discussion as a "Nuclear-cytoplasmic Western blot" in Figure S5.
3. Prior note's "no blinding, no cage information, no circadian data, no endotoxin, no STR/mycoplasma" - **all confirmed by grep (0 hits each)**.
4. No prior-note claim was found to be factually wrong.

---
# 3. Adv Sci 2026 — Xu et al., "Developmental Environmental Chronotoxicity by Polystyrene Nanoplastics Exacerbates High-Fat Diet-Induced Pediatric MAFLD Through Circadian Disruption and Mitochondrial Dysfunction"

**Access level: FULL TEXT read** (Europe PMC/PMC full-text conversion, local file `C:\Users\Soham\AppData\Local\Temp\claude\d28fb66f161d523b\5dbfd995-fe63-4ce9-8501-246b4fa39b6f\scratchpad\dl\PMC13616410.txt`, 449 lines; entire body read from the title through the Data Availability Statement, i.e. everything before and after `### References`). **Supplementary files NOT read** — Supporting File 1 `advs77840-sup-0001-SuppMat.docx` (2.6 MB) and Supporting File 2 `advs77840-sup-0002-SuppMat.docx` (58.2 KB) are .docx; Tables S1–S3 and Figures S1–S5 live there and every number sourced only from them is flagged below.

Citation: Xu P, Huang S, Zhu L, Huang J, Li K, Li G, Wu J, Zhou M, Zhou F, Ding Z, Shen Y. *Adv Sci (Weinh)* 2026; article e77840. DOI **10.1002/advs.77840** (https://doi.org/10.1002/advs.77840). PMID 42801606; PMCID PMC13616410. Received 2026-04-28; revised 2026-09-08; accepted 2026-09-10. Corresponding: Zhishan Ding (dzs@zcmu.edu.cn), Yingzhi Shen (20241160@zcmu.edu.cn). CC-BY.

---

## 1. Exact scientific question

Two questions stated verbatim at the end of the Introduction: "whether developmental stage determines susceptibility to PS-NPs-induced metabolic disturbances under dietary stress, and how PS-NPs are integrated into hepatic regulatory networks to drive MAFLD progression." Operationalised as: does chronic oral 100-nm polystyrene nanoplastic (PS-NP) exposure exacerbate HFD-induced hepatic steatosis more in juvenile than adult male mice, and is the mechanism disruption of a hepatic circadian–mitochondrial–AMPK axis (framed under DOHaD and a coined concept of "environmental chronotoxicity")?

## 2. Experimental design

**In vivo (Results 2.1, Figure 1a; Methods 4.1 Animals, 4.5 High-Fat Diet and PS-NPs Exposure)**
- Two separate age cohorts: juvenile = 3 weeks old, adult = 8 weeks old at start of exposure. Within *each* cohort a **2 × 2 factorial**: diet (standard chow vs HFD 60% kcal fat) × PS-NPs (0 vs 20 µg/mL drinking water), 10 weeks. Groups: NC, PS, HFD, HP (HFD+PS-NPs) → 8 groups total.
- Verbatim: "mice were subjected to a 2 × 2 factorial experimental design combining dietary intervention (standard chow or HFD, providing 60% of calories from fat) and PS-NPs exposure (0 or 20 µg/mL in drinking water; 100-nm GFP-labeled particles) for 10 weeks".
- **Age is NOT a modelled factor** — juvenile and adult are two parallel 2 × 2 experiments.
- n per group, phenotype/biochemistry/histology/mitochondrial assays: **n = 6** (Methods 4.23: "Animal experiments included six mice in each group"; Fig. 1 legend "n = 6/group"; Fig. 4 legend "n = 6 for c–g").
- n per group, **western blot**: n = 5 (Fig. 4 legend "n = 5 for a,b"; Fig. 9 legend "n = 5 for b").
- n per group, **RNA-seq**: n = 5 (Fig. 2 and Fig. 3 legends "(n = 5/group)"); Methods 4.15 "RNA-seq was conducted on 40 liver samples" = 8 groups × 5. **Which 5 of the 6 mice were sequenced, and how chosen, is not stated.**
- Randomisation: "Animals were randomly divided into four treatment groups" (Methods 4.5); "Group allocation was randomized" (Methods 4.23). Method of randomisation not described.
- Blinding: "histopathological assessment, image quantification, and downstream data analysis were conducted without knowledge of treatment assignment" (Methods 4.23); NAS by "A blinded pathologist" (Methods 4.13).

**In vivo pharmacological rescue arm (Results 2.7, Figure 9a; Methods 4.6)**
- Juvenile HP mice only. After 6 weeks of HFD+PS-NPs, AICAR 500 mg/kg/day or SR9009 100 mg/kg/day by oral gavage for the final 4 weeks; "controls received vehicle"; "Dosing was performed at a fixed Zeitgeber time" (the ZT is not specified). n = 5 (Fig. 9b western) / n = 6 (Fig. 9c–j). Number of arms in the panel and whether every arm was gavaged is not spelled out beyond "controls received vehicle".

**In vitro (Results 2.4–2.6; Methods 4.3, 4.7, 4.9)**
- HepG2 and Huh7 (both run in parallel for essentially every assay). 2 × 2 factorial FFA (500 µM sodium oleate + 250 µM sodium palmitate) × PS-NPs (0, 5, 20 µg/mL), 48 h.
- Circadian time course: 150 nM dexamethasone synchronisation 2 h, then sampling every 4 h from 0–48 h (13 time points: 0,4,…,48), separate cultures per time point.
- Rescue: AICAR 0.25–2 mM (Selleck S1802); SR9009 0.3125–20 µM (MedChemExpress HY-16989), applied throughout the 48 h, medium refreshed q24h; vehicle PBS or DMSO ≤ 0.1%.
- n: 3 for immunoblots/protein time course (Fig. 5c,d; 7a; 8a,b); 6 for functional assays and qPCR (Fig. 5e,f; Fig. 6; Fig. 7b–g; Fig. 8c–i). Unit of analysis = "Independent animals or independently prepared cell cultures" (Methods 4.23).

## 3. Datasets generated / re-used; accessions

- Generated: one bulk liver RNA-seq dataset, **40 libraries** (Methods 4.15). No other omics (no proteomics, metabolomics, lipidomics, scRNA-seq, spatial, microbiome). No external/public dataset is re-used anywhere.
- **Data Availability Statement, verbatim:** "The data that support the findings of this study are available from the corresponding author upon reasonable request."
- **Accessions: NONE.** Full-text search for GSE/SRA/PRJNA/BioProject/Zenodo/figshare/GitHub/"accession" returns nothing; independently confirmed by reading the PMC article page. There is therefore **no accession to resolve** — the 40-library dataset is not publicly deposited as far as the paper states. (Absence of a stated accession is not proof no record exists elsewhere.)
- Only external links in the article: the BioRender schematic credit (https://BioRender.com/o0rrsac, Agreement No. IS2A7Q9BHX, Fig. 9k) and the CC-BY licence URL.
- Supporting information: `advs77840-sup-0001-SuppMat.docx` (2.6 MB) and `advs77840-sup-0002-SuppMat.docx` (58.2 KB) — contain Tables S1–S3 and Figures S1–S5; not read.
- Ethics/licence IDs: ethics application No. 2025066; approval no. ZCMUJMJY-202503-05; animal-use licence SYXK (Zhe) 2024-0027.

## 4. Statistical / computational models

**Transcriptomics (Methods 4.15 "RNA Sequencing and Bioinformatic Analysis"; Fig. 2 legend)**
- Pipeline: fastp → HISAT2 to mouse **GRCm39** → RSEM. RSEM raw counts for DE; TPM for visualisation only. **No software version numbers given for fastp, HISAT2, RSEM, DESeq2, GOATOOLS or SciPy.**
- DE tool: **DESeq2** ("Differential expression was analyzed with DESeq2 using raw count data"). Thresholds: "**|log2 fold change| ≥ 1 and adjusted p < 0.05**" (Fig. 2 legend only). Multiple testing = DESeq2's adjusted p (BH implied, not stated).
- **Design: pairwise contrasts only.** Fig. 2b: "the 25 most significantly upregulated and 25 most significantly downregulated DEGs … in each indicated comparison". **No factorial `~ diet * PS-NPs` model, no interaction term, no interaction Wald/LRT statistic anywhere in the transcriptomic analysis. No DEG counts are reported for any contrast in the main text.**
- Enrichment: "GO and KEGG enrichment analyses were performed with **GOATOOLS** and **SciPy**, respectively, with **Bonferroni correction**". Reactome results are shown (Fig. 3c,d) but **the Reactome tool, background and significance threshold are never reported**. No GSEA, no ssGSEA.
- **Not done:** WGCNA, clustering/module analysis, cell-type deconvolution or annotation, pseudobulk vs cell-level testing (no single-cell data), mixed models, surrogate-variable/batch correction ("All libraries were constructed and sequenced within the same batch; therefore, no batch-effect correction was applied"), and **no rhythmicity statistics of any kind in vivo — no cosinor, no JTK_CYCLE, no MetaCycle, no RAIN** (impossible by design: one ZT only).
- QC reported: median 43.78 M clean reads/library (39.41–56.26 M); median overall alignment 96.93% (93.51–98.92%), unique 74.55% (**60.99–85.82%**, a wide spread; not shown to be balanced across groups).

**Phenotype / wet-lab statistics (Methods 4.23 "Statistical Analysis")**
- "Datasets containing one experimental factor were evaluated by one-way analysis of variance (ANOVA). For factorial experiments, **two-way ANOVA** was applied, with diet or FFA treatment and PS-NPs exposure or concentration entered as the two factors. Where appropriate, pairwise group differences were examined using Tukey-adjusted comparisons."
- Figures 1, 4 (in vivo) and 5a,b,f, 6 (in vitro) legends: "Two-way ANOVA (HFD × PS-NPs) with Tukey correction". **So an interaction term exists in the fitted model for all in vivo phenotypes — but no interaction F or p value is reported anywhere in the main text or figure legends; only Tukey pairwise asterisks are shown.**
- **Circadian time courses are the only place a formal interaction is explicitly reported:** "Circadian time-course experiments involved separate cultures collected at each sampling point and were therefore evaluated by two-way ANOVA incorporating treatment, sampling time, and the **treatment × time interaction**. The corresponding F and p values for the main and interaction effects are provided in **Table S3**." (Table S3 = supplementary, not read.)
- Rescue experiments (Figs 7, 8a,d–i, 9): **one-way ANOVA with Tukey**.
- Assumption checks: Shapiro–Wilk (normality), Levene (variance). Two-sided, α = 0.05. Mean ± SD throughout.
- No power calculation, no effect sizes/CIs, no pre-registration, no exclusion criteria stated.

**Legend inconsistency (verify):** Figure 6 is a pure cell-culture experiment yet its legend states "Two-way ANOVA (HFD × PS-NPs)" — the factor there is FFA, not HFD. Almost certainly a copy-paste error, but it means the factors actually modelled in Fig. 6 are unconfirmed.

## 5. Does the model suit the question? (the exacerbation / difference-of-significance problem)

- **Phenotypes: partially adequate.** The design is a true 2 × 2 and the declared model is two-way ANOVA with diet and PS-NPs as factors, so a diet × PS-NPs interaction *is* estimated. However, **every "exacerbation"/"amplification" claim in the text is worded as, and evidenced by, a pairwise contrast (HP vs HFD) with Tukey asterisks** — e.g. "serum ALT and AST levels were significantly elevated in juvenile HP mice compared with mice receiving HFD alone", "additional PS-NPs exposure further increased liver index only in juvenile mice". The interaction F/p is never shown, so the reader cannot tell whether exacerbation is supported by the interaction term or only by the pairwise step.
- **Age-dependence is NOT tested at all.** Juvenile and adult are two separate 2 × 2 models; the age claim is a comparison of significance between models — the textbook difference-of-significance error. The giveaway wording: "Comparative analysis between age groups demonstrated that alterations … were substantially more pronounced in juvenile than adult mice" (Results 2.3) and "whereas the reduction was modest and statistically insignificant in adult mice" (Fig. 4d, CPT1) — i.e. significant-in-one, not-significant-in-the-other, with no three-way age × diet × PS-NPs interaction ever fitted. **This is the paper's single biggest inferential weakness, because age-dependence is its headline claim.**
- **Transcriptomic "amplification" is not a tested interaction at all.** It rests on comparing the apparent magnitude/extent of DEGs between separate pairwise DESeq2 contrasts (and with no DEG counts even given). DEG-list size is strongly sensitive to within-group variance and to the number of libraries, both of which can differ by group; "co-exposure elicited substantially enhanced gene expression remodeling" is therefore an informal, uncalibrated comparison. A `~ diet * PSNP` DESeq2 model on the 40 libraries was feasible and was not fitted.
- **Non-monotonic results undercut a simple additive/synergy reading.** Bmal1/Npas2 mRNA and BMAL1 protein go *down* with each single stressor and *up* under co-exposure; FAO/OXPHOS transcripts go *up* with single stressors and *down* under co-exposure; TEM swelling partially reverses under co-exposure. The authors read these honestly (phase shift / compensatory feedback / adaptive-then-failing), but they mean the system is not well described by "PS-NPs makes HFD worse" on any single axis.
- **Single-time-point circadian measurement in vivo:** amplitude change, phase shift and mesor change are mathematically non-identifiable from one ZT. The authors concede this explicitly; the in vivo clock claim is therefore the weakest mechanistic link, and all genuine rhythmicity evidence is in vitro.

## 6. Inputs / features / outputs / targets per analysis

| Analysis | Input | Features | Output | Target/claim it serves |
|---|---|---|---|---|
| DESeq2 | RSEM raw counts, 40 liver libraries (8 groups × 5) | genes | DEG lists per pairwise contrast; top-25 up/down heatmaps (Fig. 2b) | "PS-NPs amplify HFD transcriptional remodeling in juveniles" |
| GO/KEGG (GOATOOLS/SciPy, Bonferroni) | DEG lists per contrast | pathway membership | enriched pathway bar/dot plots (Fig. 3a,b) | circadian + AMPK (single stressors) → OXPHOS/metabolism (co-exposure) |
| Reactome (tool unreported) | DEG lists | pathway membership | Fig. 3c,d | lipid metabolism, mitochondrial FA β-oxidation, respiratory electron transport |
| Two-way ANOVA (in vivo) | per-mouse phenotype values, n = 6 | diet, PS-NPs | group means + Tukey p | exacerbation of steatosis/dyslipidaemia/injury |
| Two-way ANOVA treatment × time (in vitro) | BMAL1 protein/mRNA at 13 time points | treatment, time, interaction | F/p in Table S3 (not read) | blunted amplitude / altered phase |
| One-way ANOVA (rescue) | per-mouse or per-culture values | drug arm | Tukey p | AICAR/SR9009 rescue |

Direction of enrichment (up- vs down-regulated pathway) is **not stated** for any KEGG/Reactome result in the text — only pathway names.

## 7. How results were evaluated

- Omics thresholds: |log2FC| ≥ 1 and adjusted p < 0.05 (DESeq2); Bonferroni for GO/KEGG. No independent validation cohort, no held-out set, no cross-validation (no ML/classifier anywhere).
- Wet-lab: α = 0.05, Tukey-adjusted; Shapiro–Wilk/Levene assumption checks; n = 5–6 biological replicates; blinded histology and image quantification; ImageJ with "a uniform threshold", ≥5 random fields/section averaged per animal (Methods 4.13); ROS = 6 random fields per independent culture averaged to one value (Methods 4.22); serum run with internal QC materials each run (Methods 4.10).
- Internal replication: every in vitro result is run in **two independent cell lines (HepG2 and Huh7)** — a real, if limited, replication layer. SR9011 is used as a second REV-ERB agonist (Figure S4a–d, not read). No independent animal cohort replication, no external dataset replication.

## 8. Independent experimental validation vs omics-only inference

**Validated by orthogonal (non-sequencing) assays:**
- Liver index, body/adiposity, serum ALT/AST/LDL-C/HDL-C/TG/T-CHO (Hitachi 7180), hepatic TG/T-CHO (Nanjing Jiancheng A110-1-1 / A111-1-1), H&E + blinded NAS (Kleiner scoring), Oil Red O area — Fig. 1b–h, Fig. 9d–j.
- BMAL1 and p-AMPKα(Thr172)/t-AMPKα **western blot** in liver (Fig. 4a,b) and cells (Fig. 5c) — Proteintech 14268-1-AP, Affinity AF3423/AF6423.
- **Functional mitochondrial assays:** ATP (luciferase, Beyotime S0026), ADP (Amplex Red, S0307S), ATP/ADP ratio; JC-1 on **isolated liver mitochondria** (Beyotime C3606/C2003S) and on cells; CPT1 activity on liver mitochondrial extracts (ml301014, DTNB/TNB, nmol/min/mg); **TEM** ultrastructure (mitochondrial area, cristae density); Seahorse XFe96 mito stress test (cells only); H2DCFDA ROS (cells only).
- **RT-qPCR** for clock, lipogenic and mitochondrial genes (liver Fig. 4h; cells Fig. 5f, 6e).
- **Pharmacological rescue in vitro and in vivo** (AICAR, SR9009; target engagement shown in Fig. 9b).

**Omics-only / inference, NOT independently validated:**
- The "substantially enhanced gene expression remodeling" under juvenile co-exposure and its attenuation in adults (Fig. 2a) — no DEG counts, no orthogonal confirmation of the global claim.
- All KEGG/Reactome pathway enrichments (Fig. 3) — used as the bridge to the circadian/mitochondrial mechanism; the *functional* assays test the conclusion, not the enrichment itself.
- Npas2 and Clock changes — mRNA only; only BMAL1 was checked at protein level.
- Hepatic particle accumulation equality across ages (Figure S2c GFP fluorescence; qualitative, "comparable", not quantified in the main text).
- The GFP-label-is-innocuous claim: asserted in the Discussion ("Comparable biological responses between GFP-labeled and unlabeled PS-NPs…") with its data only in supplementary.

**Explicitly NOT done anywhere (author-acknowledged):** any genetic loss/gain-of-function (no Bmal1, Nr1d1/Rev-erbα, Prkaa knockout, knockdown, liver-specific deletion or overexpression), ChIP, FAO flux assays, isotope tracing, acylcarnitine profiling, isoform-specific CPT1A measurement, in vivo circadian time-series sampling, temporal proteomics/PTM analysis. No ABX/FMT, no microbiome work at all.

## 9. Main findings with directions

**In vivo, juvenile (all directions vs HFD alone unless stated):**
1. Liver index ↑ further only in juveniles with PS-NPs added to HFD; PS-NPs alone no effect at either age (Fig. 1b).
2. Serum ALT ↑, AST ↑ in juvenile HP vs HFD; no further rise in adults (Fig. 1c).
3. Serum LDL-C ↑, HDL-C ↑, TG ↑, T-CHO ↑ — "PS-NPs co-exposure significantly amplified HFD-induced elevations" — juveniles only (Fig. 1d,e); hepatic TG and T-CHO ↑ (Fig. 1f).
4. Histology: more severe steatosis, hepatocellular ballooning, vacuolar degeneration, ↑NAS, ↑Oil Red O area in juvenile HP; "limited additional pathological changes" in adults (Fig. 1g,h).
5. **Not changed by PS-NPs:** body-weight gain, adiposity (Figure S3a,b), glucose tolerance and insulin sensitivity (Figure S3c,d), energy intake (Figure S1c). The liver phenotype is dissociated from obesity and systemic glycaemia.
6. p-AMPKα ↓ with either single exposure and **further ↓ under co-exposure** (monotonic; Fig. 4b). Hepatic ATP/ADP ↓ (Fig. 4c), ATP ↓ and ΔΨm ↓ under co-exposure (Fig. 4f,g).
7. **CPT1 enzyme activity ↓** in juvenile HP vs HFD alone; in adults "the reduction was modest and statistically insignificant" (Fig. 4d).
8. TEM: swelling after either single exposure; co-exposure → **cristae density ↓** with "partial alleviation of mitochondrial swelling" (Fig. 4e) — non-monotonic.

**Gene-family directions requested (liver RT-qPCR, Fig. 4h, juvenile, unless noted):**
- **Lipogenesis: UP.** *Fasn* (legend "Fas"), *Lxra*, *Srebp1c* "progressively elevated following PS-NPs and HFD co-exposure". **Acly and Acaca were NOT measured.** In cells, FAS and LXRα mRNA altered by co-exposure (Fig. 5f).
- **PPARα / FAO: Ppara itself was NEVER measured; Acox1 and Cyp4a were never measured.** The FAO axis is represented only by *Cpt1a* and *Acadl* mRNA, *Pgc-1α* mRNA and CPT1 enzyme activity. Direction: **biphasic** — "increased expression following individual stressors but decreased expression after combined exposure" (*Acadl, Pgc-1α, Cpt1a, Atp5a1, Ndufa9, Sdha*), with **CPT1 activity DOWN** in juvenile HP. Net: **FAO capacity suppressed under co-exposure.** In cells, *Acadl, Pgc-1α, Atp5a1, Ndufa9* ↓ under co-exposure (Fig. 6e).
- **Cholesterol synthesis: NOT ASSAYED.** No *Hmgcr*, *Hmgcs1*, *Sqle* or *Srebf2* anywhere. Only serum and hepatic T-CHO (both ↑ with co-exposure in juveniles).
- **Inflammation: essentially NOT STUDIED.** No *Tnf*, *Il1b*, *Il6*, NF-κB, NLRP3, F4/80/Adgre1/Cd68/Clec4f or any macrophage marker, no IHC for immune cells. The only inflammatory read-out is the lobular-inflammation component inside the composite NAS (H&E). **Treat this paper as uninformative on MP-induced hepatic inflammation.**
- **Clock genes:** *Bmal1* and *Npas2* mRNA (RNA-seq, Fig. 2b) ↓ with PS-NPs alone and ↓ with HFD alone, then "partially restored under co-exposure relative to HFD alone" — **non-monotonic**; authors refuse a restoration reading. *Bmal1, Clock, Npas2* also in the Fig. 4h qPCR panel. **BMAL1 protein** (Fig. 4a): ↓ with either single exposure, **↑ under combined exposure** at ZT4–6 — authors attribute to "altered circadian phase, compensatory feedback regulation, or temporal misalignment rather than restoration". A "persistent faster-migrating BMAL1 protein band" is noted and unexplained (Discussion).
- **Fibrosis: NOT STUDIED** (no Sirius red, no Col1a1/Acta2/Timp1, no hydroxyproline; NAS excludes fibrosis stage).
- **Oxidative stress:** ROS ↑ in hepatocytes under co-exposure (H2DCFDA, Fig. 6c), reduced by AICAR and SR9009. **No in vivo oxidative-stress measurement at all** (no MDA, SOD, GSH, Nrf2, Hmox1).

**In vitro:**
9. PS-NPs + FFA → lipid droplets ↑ and TG ↑ vs FFA alone, both cell lines (Fig. 5a,b).
10. **BMAL1 protein ↑ with PS-NPs alone and ↑ with FFA alone, most pronounced under co-exposure (Fig. 5c) — the OPPOSITE direction to the in vivo single-exposure result (↓).** p-AMPK ↓. The paper does not reconcile this discrepancy. Note also an **internal tension**: the unsynchronised 48-h blot (5c) shows BMAL1 ↑, while the synchronised time course (5d) reports "blunted oscillatory amplitude and sustained suppression of BMAL1 protein" under co-exposure.
11. Dexamethasone-synchronised 48-h time course: co-exposure blunts BMAL1 amplitude and alters phase (Fig. 5d,e; Fig. 8b,c) — **the only place a formal interaction term (treatment × time) is explicitly reported**, F/p in Table S3.
12. Mitochondria in cells: ΔΨm ↓, ATP ↓, ROS ↑, ATP-linked OCR ↓, maximal respiration ↓ (Fig. 6a–d). Caveat given: "absolute proton leak OCR was decreased rather than increased … should be interpreted cautiously", read as overall respiratory decline, not uncoupling.
13. **AICAR rescue (Fig. 7):** p-AMPK ↑ dose-dependently; ΔΨm ↑, ATP ↑, ATP-linked and maximal OCR ↑, ROS ↓, lipid droplets ↓, TG ↓.
14. **SR9009 rescue (Fig. 8):** BMAL1 ↓ dose-dependently; "partially restored both oscillatory amplitude and phase dynamics"; ΔΨm ↑, ATP ↑, ATP-linked/maximal respiration/spare capacity/coupling efficiency ↑, proton leak normalised, ROS ↓, lipid ↓, TG ↓. SR9011 reproduces mitochondrial benefit (Figure S4, not read).

**In vivo rescue (juvenile HP, Fig. 9):**
15. Target engagement confirmed (p-AMPK ↑ with AICAR; BMAL1 ↓ with SR9009, Fig. 9b). Both drugs ↓ body-weight gain. AICAR ↓ liver index, ↓ adiposity index, ↓ hepatic lipid, ↓ ALT, ↓ AST. SR9009 ↓ liver index and ↓ ALT but **not AST**, with "comparatively limited" effects on adiposity and steatosis. **Neither changed serum TG or T-CHO.** Both ↑ hepatic ΔΨm and ATP and ↓ NAS.
16. Dissociation claim: SR9009 protects mitochondria strongly while barely reducing steatosis, AICAR the reverse → authors infer "a distributed regulatory network model" (Fig. 9k schematic). **Inference from pharmacological dissociation, not a directly tested mechanism.**

## 10. Limitations

**(a) Actual limitations of what WAS done**
1. **Age-dependence — the headline claim — is never statistically tested.** Two separate 2 × 2 models; significance-in-juveniles vs non-significance-in-adults is compared informally (Results 2.3, Fig. 4d). No age × diet × PS-NPs interaction.
2. **No interaction F/p reported for any in vivo phenotype**, although two-way ANOVA was fitted. All exacerbation wording is pairwise (HP vs HFD).
3. **Transcriptomic "amplification" has no interaction model and no DEG counts** — DESeq2 pairwise contrasts only; only top-25 heatmaps per contrast shown.
4. **Single circadian time point in vivo (ZT4–ZT6)** — cannot separate amplitude, phase and mesor; authors concede this explicitly and still title the paper "chronotoxicity". All true rhythmicity data are in vitro in immortalised lines.
5. **Exposure dose is a cage-level estimate**, not an individual intake: group housing, shared drinking water. Number of mice per cage and number of cages per group are **not reported**, so treatment may be confounded with cage even though Methods 4.23 declares the individual animal the statistical unit. The estimate also covers only "the first eight weeks" of a 10-week exposure, unexplained.
6. **Control diet is unspecified "standard chow"**, not a purified low-fat control matched to D12492 (e.g. D12450H) — fibre/sucrose/micronutrient confounding for liver transcriptome and lipid phenotype.
7. **CPT1 assay is not isoform-specific** (author-acknowledged) and measures import capacity, not FAO flux.
8. **No endotoxin, surfactant, residual-monomer or leachate testing** on the particles; the GFP protein label alters the protein corona and the unlabelled-particle comparison lives only in supplementary.
9. **Single formulation, single size, single concentration** (author-acknowledged); environmental relevance of 20 µg/mL unestablished (author-acknowledged).
10. **Male mice only** (author-acknowledged).
11. **Pharmacological-only causality** — AICAR and SR9009 off-target effects "cannot be completely excluded" (author-acknowledged); no genetic validation.
12. **In vivo vs in vitro BMAL1 direction conflict** (↓ vs ↑ for single exposures) is never addressed; the in vitro 5c/5d tension likewise.
13. **Unexplained faster-migrating BMAL1 band** (author-acknowledged).
14. **qPCR normalised to β-actin alone**, no reference-gene stability check, under conditions known to alter cytoskeletal/metabolic transcripts.
15. **No software versions** for any bioinformatics tool; Reactome tool/threshold unreported; wide unique-alignment-rate spread (60.99–85.82%) not shown to be balanced by group; no batch correction (single batch asserted).
16. Figure 6 legend states the wrong factor ("HFD × PS-NPs") for a cell experiment.
17. **No data deposition** — "available from the corresponding author upon reasonable request" only; the RNA-seq cannot be reanalysed (e.g. to fit the missing interaction model).
18. Which 5 of 6 animals per group were sequenced, and by what rule, is not stated.

**(b) Things simply NOT attempted (scope, not a flaw of execution)**
- No genetic models of any kind (Bmal1, Rev-erbα/Nr1d1, Prkaa, Ppara).
- No in vivo circadian time series; no wheel-running/actigraphy; no light-phase manipulation; no measurement of day/night feeding distribution (HFD itself shifts feeding rhythm — a plausible uncontrolled proximate cause of the clock changes).
- No FAO flux, isotope tracing, acylcarnitine profiling, direct CPT1A assay (author-flagged as future work).
- No proteomics, metabolomics, lipidomics, scRNA-seq or spatial transcriptomics; no cell-type resolution.
- No inflammation or fibrosis panel; no gut microbiome, no ABX/FMT; no gut barrier work.
- No female mice; no other polymers (PE, PP, PET), sizes, aged/weathered particles or co-contaminant mixtures; no dose–response in vivo.
- No external-dataset replication, no human data.
- No sex, no dams: **this is NOT a gestational/lactational developmental exposure** — exposure starts post-weaning at 3 weeks; no dams were exposed and none are reported.

## Study-characteristics record

| Item | Value (location) |
|---|---|
| Species/strain/sex | Male C57BL/6J, Shanghai SLAC Laboratory Animal Co. (Methods 4.1) |
| Age | Juvenile 3 wk, adult 8 wk at start; 10-wk exposure → sacrificed at ~13 and ~18 wk (Results 2.1, Methods 4.1) |
| Dams / developmental window | **None — no gestational or lactational exposure.** "Developmental" here = post-weaning juvenile (3 wk) vs adult (8 wk) |
| Polymer | Polystyrene, GFP-fluorescent-labelled (PS-NPs-GFP) |
| Size | Nominal 100 nm; DLS hydrodynamic diameter **115.4 ± 1.5 nm**, PDI **0.041 ± 0.040**, zeta **−30.9 ± 0.7 mV** (Zetasizer Nano ZS90, software v7.3.0) (Methods 4.4; Results 2.1; Figure S2b) |
| Supplier/cat | **Rigor Biotechnology, FF3301G** (Methods 4.4) |
| Dose | 20 µg/mL in drinking water (in vivo); 5 and 20 µg/mL in vitro, 20 µg/mL stated as "the highest non-cytotoxic dose" (Methods 4.5) |
| Estimated internal dose | 25.53 / 21.69 / 24.17 / 24.81 µg/g bw/week (juv PS, juv HP, adult PS, adult HP); cumulative 204.24 / 173.56 / 193.33 / 198.46 µg/g bw — **first 8 weeks only, cage-level estimate, Table S1 (not read)** |
| Route / vehicle | Oral, ad libitum drinking water; vehicle = drinking water. In vitro vehicle PBS or DMSO ≤0.1% |
| Frequency / duration | Continuous; solution replaced every 3 days; bottles sonicated twice daily 30 min; stock vortexed + sonicated 40 kHz 10 min; **10 weeks** |
| Particle stability check | Daily water sampling + SEM: "No detectable aggregation or sedimentation … over the 3-day replacement cycle" (Methods 4.4; Figure S2a) |
| **Endotoxin / surfactant** | **NOT reported anywhere** (searched whole text) |
| Diet | HFD = **60% kcal fat, Research Diets D12492**; control = unspecified "standard chow" (**no catalog number**) (Methods 4.5) |
| In vitro lipotoxic stress | 500 µM sodium oleate + 250 µM sodium palmitate, 48 h |
| Sampling time / ZT | Lights on 06:00 = ZT0; 12:12 LD; 7-day acclimation; "metabolic interventions, tolerance tests, euthanasia, and terminal tissue collection were performed within the same circadian window (**ZT4–ZT6**)" — hedged by "Unless otherwise indicated" (Methods 4.8). Liver snap-frozen in LN2, −80 °C |
| Fasting state | Ad libitum except for metabolic tests: **GTT 16-h fast**, glucose 1 g/kg i.p.; **ITT 4-h fast**, insulin 0.75 U/kg i.p.; both started at ZT4 (Methods 4.12). Fasting state at terminal sampling not stated (implied fed/ad libitum) |
| Housing | 25 ± 1 °C, 50–60% RH, 12:12 LD; **group-housed**; food/water measured per cage per 3 days (100 g chow, 200 mL water per cage) (Methods 4.1, 4.2). **Mice per cage and cages per group NOT reported**; no litter-level information |
| Randomisation | Stated, method not described (Methods 4.5, 4.23) |
| Blinding | Histopathology, image quantification, downstream analysis blinded (Methods 4.13, 4.23) |
| Cell lines | HepG2, Huh7 (Cell Bank, Chinese Academy of Sciences), STR-authenticated, mycoplasma-monitored (Methods 4.3) |
| DOI | https://doi.org/10.1002/advs.77840 |

## Relevance note for the matrix

Closest published analogue to a 2 × 2 oral-nanoplastic × HFD liver-RNA-seq design (n = 5/group) with functional validation. Two directly usable contrasts: (i) sampling time **is** controlled here (ZT4–ZT6, 7-day LD acclimation) and the authors still refuse to read a single-ZT clock measurement as rhythm disruption — a useful precedent against over-reading single-time-point clock data; (ii) the exacerbation claim is modelled by two-way ANOVA for phenotypes but reported only as pairwise contrasts, and the transcriptomic amplification has no interaction model at all. Also note: PS not PE, 100 nm, 60% kcal lard HFD (not MASH/CDAHFD), male-only, group-housed with cage-level intake, no inflammation or fibrosis data, and **no public sequence accession**.

## URLs/files actually read

1. `C:\Users\Soham\AppData\Local\Temp\claude\d28fb66f161d523b\5dbfd995-fe63-4ce9-8501-246b4fa39b6f\scratchpad\dl\PMC13616410.txt` — full Europe PMC full-text conversion, read in its entirety (lines 1–449: body, all figure legends, all Methods 4.1–4.23, Statistical Analysis, Data Availability, Associated Data).
2. `C:\Users\Soham\Desktop\microplastic\research\literature\papers\advsci2026_psnp_circadian.md` — prior-session notes (re-checked, not a source).
3. https://pmc.ncbi.nlm.nih.gov/articles/PMC13616410/ — fetched to independently confirm the Data Availability Statement, the absence of any repository accession, and the supplementary file names/sizes.

Not read: https://doi.org/10.1002/advs.77840 (publisher page — not needed, PMC version is the CC-BY version of record); the two supplementary .docx files.

## UNVERIFIED items

1. **Tables S1, S2, S3 and Figures S1–S5 were not read** (.docx supplementary, 2.6 MB + 58.2 KB). Everything that lives only there is unverified: the per-group intake/cumulative-dose table (S1), qPCR primer sequences (S2), and **the F and p values for all main and interaction effects of the circadian time-course two-way ANOVAs (Table S3)** — i.e. the one formal interaction test in the paper could not be inspected.
2. Figures S1 (food/water/energy intake), S2 (SEM, DLS/zeta, hepatic GFP), S3 (body weight, adiposity, GTT, ITT), S4 (SR9011), S5 (dose-selection/toxicity) — all numbers and directions quoted from these come from the main-text narrative only.
3. **Whether an interaction F/p exists anywhere for the in vivo phenotypes** could not be settled — it is absent from the main text and legends; it may or may not be in supplementary or source data.
4. **Mice per cage and cages per group** — not reported anywhere; the cage-vs-animal unit-of-analysis question is therefore unresolvable from the publication.
5. **Which 5 of the 6 animals per group were sequenced**, and the selection rule — not reported.
6. Exact **number of arms in the in vivo rescue panel** and whether every arm (including HP control) received gavage — "controls received vehicle" is the only statement.
7. The **ZT of AICAR/SR9009 dosing** — "a fixed Zeitgeber time", value not given.
8. **Software versions** for fastp, HISAT2, RSEM, DESeq2, GOATOOLS, SciPy; the **Reactome enrichment tool, background and threshold**; the multiple-testing method behind DESeq2's "adjusted p" (BH assumed, not stated).
9. **Fasting state at terminal sacrifice** — not explicitly stated (ad libitum implied outside the metabolic tests).
10. **Control chow identity/catalog number and macronutrient composition** — never given.
11. **Endotoxin, surfactant, residual monomer, leachates** — never tested or reported; cannot be verified either way.
12. Whether the ZT4–ZT6 collection window applies specifically to the RNA-seq animals — Methods 4.8 is hedged with "Unless otherwise indicated"; no exception is stated, so it presumably does, but this is an inference.
13. The DEG counts per contrast — absent from the main text; location unknown.

## Prior-note claims checked — corrections and flags

Prior notes were accurate on essentially every load-bearing claim re-checked against the full text (design, groups, n per assay, randomisation/blinding quotes, particle characterisation, dose, diet, ZT window, DESeq2 thresholds, GOATOOLS/SciPy/Bonferroni, the pairwise-vs-interaction problem, data-availability wording, the absence of accessions, the absence of inflammation/cholesterol-synthesis/Ppara data, and all author-acknowledged limitations). Specific points:

- **Confirmed and worth strengthening:** prior notes describe the paper as "developmental exposure"-adjacent; the full text confirms there is **no gestational or lactational component and no dams** — exposure begins post-weaning at 3 weeks. Any matrix row implying maternal exposure would be wrong.
- **Confirmed:** no accession; Data Availability verbatim quote exact (re-verified against the live PMC page).
- **Confirmed:** Figure 6 legend says "Two-way ANOVA (HFD × PS-NPs)" for a cell experiment — a real legend error in the paper.
- **Minor prior-note inaccuracies (non-load-bearing):** the local file is 449 lines, not 453; Fig. 4h gene labels appear as "Fas, Lxrα" in the legend and "Fasn, Lxra" in the body; the Figure 1 legend has a typo in its significance key ("**p < 0.0001" where "****" is meant).
- **Not confirmable:** prior note stating in vitro n = 6 for "Fig. 5e,f" vs n = 3 for "5c,d" matches the legend exactly, but the legend's n assignment to panel 5e (a time-course panel) alongside 5f is internally odd — 5d and 5e are both time courses yet are assigned different n (3 vs 6). Flagged as a possible legend error in the paper, not a prior-note error.
- **New tension not in the prior notes:** Fig. 5c (unsynchronised, 48 h) reports BMAL1 protein **increased** by co-exposure, while Fig. 5d (synchronised time course) reports "sustained suppression of BMAL1 protein" under co-exposure. The paper does not reconcile these, in addition to the in vivo (↓ with single exposures) vs in vitro (↑ with single exposures) conflict the prior notes already flagged.
# 4. EHP 2024 — Zhai, Yang, Chen & Wu, "Long-Term Exposure to Polystyrene Microspheres and High-Fat Diet–Induced Obesity in Mice: Evaluating a Role for Microbiota Dysbiosis"

*Environ Health Perspect* 2024;132(9):097002. DOI [10.1289/EHP13913](https://doi.org/10.1289/EHP13913). PMID 39226184. PMCID PMC11370995. Received 2023-08-31; accepted 2024-07-23; published online 2024-09-03. Affiliations: China Agricultural University (Beijing) and Hong Kong Polytechnic University. Linked EHP Science Selection commentary: "Gut Check: Microbiota and Obesity in Mice Exposed to Polystyrene Microspheres", EHP 132(11):114001, DOI 10.1289/EHP15844 (commentary, not the research paper).

**ACCESS LEVEL: FULL TEXT.** I read the complete Europe PMC JATS body (front matter, abstract, Introduction, all Methods subsections, all Results, all seven figure legends plus the accessibility long-descriptions of every panel, Discussion, Conclusions, Acknowledgments/data statement) from the local conversion `PMC11370995.txt`, lines 1–182 (everything before `### References`). **Supplemental Material (Figures S1–S5, Tables S1–S3) was NOT read** — the PMC article page returned a reCAPTCHA interstitial and `ehp.niehs.nih.gov` did not resolve from this environment. Every numeric value that lives only in Table S3 is therefore unverified, and is marked as such below.

---

## 1. Exact scientific question

Verbatim Objective (Abstract): *"We aimed to investigate the effect of PS-MS exposure on high-fat diet (HFD)–induced obesity and underlying mechanisms."*

Verbatim framing of the gap (Introduction, final paragraph of the gap statement): *"it is still not clear whether PS-MS could disrupt gut barrier function by inducing gut dysbiosis and how PS-MS disrupts the gut microbiota and contributes to HFD-induced obesity."* The stated hypothesis is directional: *"it is highly plausible that exposure to PS-MS will induce stress on gut microbiota species, thereby compromising gut integrity and promoting obesity induced by a HFD."*

So the question has two parts: (a) does chronic (8-wk) oral 5-µm PS microsphere exposure worsen HFD-induced obesity and its liver/glucose/lipid sequelae; (b) is that worsening *mediated by* the gut microbiota (addressed by antibiotic depletion and FMT).

## 2. Experimental design

**Three independent 24-mouse experiments.** Methods, *Animals and Diets*: *"Seventy-two male C57BL/6J mice (24 mice per experiment), aged 6–8 wk"*. Figure 1A/4A: exposure started at 8 wk of age.

### Experiment 1 — main 2×2 (Figures 1, 2, 3)
- Factors: diet (ND vs HFD) × particle (water vs 10 mg/L PS-MS). Four groups: ND, ND+PS-MS, HFD, HFD+PS-MS.
- n = 6 per group for phenotype (Fig. 1 legend: *"All values are represented as mean±SD in scatter plots (n=6). Each dot represents one mouse."*). Duration 8 wk.
- Reduced n for some endpoints: MASLD activity score **n = 5** (Fig. 1G); Oil Red O — Methods/*Histology*: *"livers from 4–5 mice per group were used for staining"*; TJP1/occludin immunofluorescence **n = 3** (Fig. 3F,G, "each dot ... represents one mouse"); everything else in Fig. 3 "n=6 unless otherwise indicated".
- **16S rRNA: n = 5 per group** (Fig. 2 legend: *"Data are represented as mean±SD with scatter plots (n=5). Each dot represents one mouse."*). Feces collected from each individual mouse every 4 wk (Methods, *16S rRNA Sequencing*); the figure reports the 8-wk time point.

### Experiment 2 — ABX (antibiotic-depleted) replicate of the 2×2 (Figures 4, 5)
- Same four diet × particle groups, every mouse additionally gavaged **daily for the full 8 wk** with a broad-spectrum cocktail. Methods, *Antibiotic Treatment*: *"ampicillin (100 mg/kg BW), neomycin (100 mg/kg BW), metronidazole (100 mg/kg BW), and vancomycin (50 mg/kg BW) throughout the 8-wk ND/HFD-feeding period (n=6 per group)"* (supplier Beijing Suo Lai Bao).
- **DRUG-IDENTITY INCONSISTENCY (new, not in the prior notes):** Methods says *ampicillin*; the Figure 4 legend says *"penicillin (100 mg/kg BW), neomycin ..."* and the Discussion repeats *"penicillin, neomycin, metronidazole, and vancomycin"*. Two of three statements say penicillin, the Methods says ampicillin. Unresolved in the paper.
- Depletion efficiency: *">97% of the intestinal bacteria were eradicated (Figure S4)"* — the supporting data are in Figure S4 (not read).
- n = 6 per group; MASLD score **n = 4** (Fig. 4G); IF **n = 3** (Fig. 5F,G); goblet cells reported as **n = 9 crypts** (Fig. 5A: *"Each data represents 1 crypt, 3 slices/group, 3 crypts/slice (n=9)"*).

### Experiment 3 — FMT (Figures 6, 7)
Methods, *FMT*, in full:
- **Donors:** *"8-wk-old male donor mice were fed either the HFD or the HFD + PS-MS (10 mg/L in drinking water) (n=6 per group)"*, for 8 wk. Only HFD-background donors — no ND donors.
- **Recipients:** HFD-fed, pre-treated with antibiotics **in drinking water for 7 consecutive days before FMT**: *"vancomycin (0.5 g/L), neomycin sulfate (1 g/L), metronidazole (1 g/L), and ampicillin (1 g/L)"*. n = 6 per group; two groups, R-HFD and R-HFD+PS-MS. Fig. 6 legend: *"The mice were randomly divided into four groups (n=6)."*
- **Inoculum:** *"Fecal samples were collected from the donors ... and pooled over the course of the 8-wk experiment. Donor feces were diluted with PBS (1 g of sample in 10 mL of PBS) and homogenized for 1 min using a vortex ... centrifuged at 900×g for 3 min to remove particulate matter."* Feces pooled across donors and across time.
- **Timing/route:** *"About 200 µL of the supernatant was immediately administered to each recipient mouse by oral gavage three times a week for 8 wk"*; fresh material prepared within 10 min of gavage.
- Engraftment check: 16S PCA of recipient feces, **Figure S5** (not read).
- n = 6 per group; MASLD score **n = 4** (Fig. 6G); goblet cells **n = 9 crypts** (Fig. 7A).
- Arithmetic: 24 mice = 12 donors + 12 recipients, exactly, leaving no spares.

### Housing / randomisation / blinding
- *"Groups of mice were housed together (n=6 per group, 3/cage)"* — i.e. **2 cages per group, 3 mice per cage, analysed as 6 independent animals**. No cage-level random effect anywhere.
- Baseline standardisation: *"all groups of mice were co-housed for a period of 1 wk prior to the commencement of the experiment"*.
- Environment: 24±1 °C, 55%±5% RH, 12 h light:dark, ad libitum feed and water.
- **Randomisation** is stated only once, in the Figure 6 (FMT) legend. The word "random*" appears exactly once in the whole body text. Not stated for Experiments 1 or 2.
- **Blinding: the string "blind" does not occur anywhere in the body text.** MASLD scoring, adipocyte area, goblet-cell counts and IF intensity are all operator-scored.
- Ethics: China Agricultural University, approval No. AW72403202-1-31.

### Omics vs phenotype n, summarised
| Assay | Experiment | n |
|---|---|---|
| Body weight, organ/fat mass, OGTT, serum & liver lipids | main 2×2 | 6/group |
| MASLD activity score | main / ABX / FMT | 5 / 4 / 4 |
| Oil Red O | main | 4–5 |
| 16S rRNA (V3–V4) | main 2×2 only | **5/group** |
| Colon qPCR (Muc2, Tjp1), MUC2 ELISA, serum ELISAs | all three | 6 (unless noted) |
| TJP1/occludin IF | main, ABX | 3 |
| AB/PAS goblet cells | ABX, FMT | "n=9" crypts |
| Metagenomics / RNA-seq / metabolomics / SCFA | — | **none performed** |

## 3. Datasets

**Data availability statement, verbatim** (there is no "Data Availability" heading; it is the closing paragraph of the Acknowledgments):

> "The entirety of the data generated throughout this study has been comprehensively documented and is available within the main body of this article, as well as its supplementary information. The 16S rRNA raw reads are available at the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA) repository (BioProject PRJNA990701)."

The same accession is repeated in the Figure 2 legend.

**Verification of PRJNA990701 — it resolves.** Checked at NCBI BioProject and at ENA:
- Title: "High-throughput sequencing of mouse fecal microorganisms"; description "The study aims to investigate the relationship between intestinal flora and obesity"; center/submitter **China Agricultural University**; registered 3-Jul-2023; first public **2023-07-26**; public status.
- **20 read runs**, all `AMPLICON`, Illumina MiSeq, sample title "Metagenome or environmental sample from feces metagenome": SRR25114495–SRR25114514 (19 runs) plus **SRR29975743** (a later-added run, accession block consistent with mid-2024, i.e. around revision/acceptance).
- **FLAG: 20 runs = exactly 4 groups × n=5 of the main 2×2 experiment.** The 16S data underlying **Figure S4** (ABX depletion check) and **Figure S5** (FMT engraftment PCA) are therefore almost certainly **not deposited**, even though those two supplementary figures carry the paper's entire mediation argument. (Strictly: 20 runs cannot cover main + ABX + FMT cohorts; I could not inspect per-run sample metadata to confirm which cohort each belongs to.)

**No other accessions exist.** No GEO/GSE, no SRA PRJNA beyond this one, no PRIDE, no MetaboLights, no Zenodo/figshare, no GitHub, no code repository, no analysis scripts. Nothing was re-used from an external dataset.

Software cited in place of code: GraphPad Prism 9.0; R 3.5.3 (correlations) and R 3.6.0 (PCA); QIIME 1.8.0; Pear 0.9.6; Vsearch 2.7.1 (chimera removal and Uparse-algorithm OTU clustering); BLAST; Illumina Analysis Pipeline 2.6; ImageJ; Image-Pro Plus.

Supplementary files referenced: **Table S1** (qPCR primers), **Table S2** (diet compositions), **Table S3** (all numeric source data — cited in *every* figure legend as "Corresponding numeric data can be found in Table S3"), **Figure S1** (TEM of PS-MS), **Figure S2** (water intake), **Figure S3A** (observed species), **Figure S3B** (Shannon), **Figure S4** (ABX depletion), **Figure S5** (FMT engraftment PCA). Supplemental Material is hosted at https://doi.org/10.1289/EHP13913.

## 4. Statistical / computational models

**Statistical Analysis subsection, verbatim and complete:**

> "Analysis of the data was performed using GraphPad software (version 9.0). Unless otherwise indicated, the parametrically distributed data underwent Student's t-test for comparisons between two groups, whereas nonparametrically distributed data were analyzed using the Mann–Whitney test, and differences in mean values among more than two groups were determined using one-way analysis of variance (ANOVA) with a Tukey multiple comparison test. The Brown–Forsythe and Welch ANOVA test was performed when the variances were not equal. Multiple unpaired t-tests and two-way repeated measures ANOVA with a Tukey multiple comparison test were performed as indicated in figure legends. Correlation statistical analyses are performed by Pearson's correlation by using R (version 3.5.3; R Development Core Team). Statistical significance was considered when p<0.05."

### Was a diet × particle interaction term fitted? **NO.**
- Figures 1, 2, 3, 4, 5 legends all state: *"Data were analyzed using one-way ANOVA with a Tukey multiple comparison test"* (Fig. 2 adds that panels H, J, K, M, N used *"the Brown–Forsythe and Welch ANOVA test"*). The 2×2 is collapsed into a **one-factor, four-level** ANOVA.
- The only "two-way" models are **two-way repeated-measures ANOVA on the body-weight-over-time curves** (Fig. 1B, Fig. 4B). Their two factors are **group and time**, not diet and particle; the legends report per-time-point Tukey contrasts (Fig. 1B: `###p<0.001, ####p<0.001, HFD + PS-MS vs. HFD`; `**/***/****`, HFD vs. ND).
- **No interaction F statistic or interaction p-value is reported anywhere in the paper.** The word "interaction" occurs four times in the body, all in a biological sense ("biological interactions", "mucus–bacterial interactions"), never statistical.
- The ABX arm is a *separate cohort*, not an ABX factor crossed into the main design, so there is likewise no diet × PS × ABX three-way model and no formal mediation analysis (no causal-mediation, no bootstrapped indirect effect).

### 16S pipeline and diversity statistics
- Reads merged with Pear 0.9.6; chimeras removed with Vsearch 2.7.1; OTUs by the Uparse algorithm of Vsearch at **97% similarity**; taxonomy by BLAST against NCBI; rarefaction curves and α-diversity (**Shannon, Simpson, Chao1**) in QIIME 1.8.0; **PCA (not PCoA) on OTU tables in R 3.6.0**; relative abundances at phylum→species; differential abundance by **Metastats and LEfSe**.
- **No β-diversity hypothesis test of any kind.** The strings PERMANOVA, ANOSIM, adonis, PCoA, ANCOM do not appear. Group separation is asserted from a PCA scatter (Fig. 2B, PC1 33.89%, PC2 18.93%) with no accompanying statistic — the Results claim *"microbiota composition in both ND and HFD mice were considerably different in mice exposed to PS-MS compared with those not (Figure 2B,C)"* rests on eyeballing an ordination.
- **LEfSe thresholds are not stated** (no LDA score cut-off, no α for the Kruskal–Wallis/Wilcoxon steps). Metastats thresholds not stated. **No multiple-testing correction is reported anywhere**: "FDR", "false discovery", "Benjamini" return zero hits in the body text.
- α-diversity (Chao1, observed species, Shannon) is compared by the same one-way ANOVA + Tukey, n=5.

### Differential-expression tool
Not applicable — there is **no transcriptomics**. Gene expression is SYBR-based RT-qPCR (ABI 7500), 2^−ΔΔCT, normalised to **GAPDH only**, on **two genes in colon only (Muc2, Tjp1)**; primers "designed and saved by our laboratory" / "commonly used primers in our laboratory" (Table S1); 4,000 ng RNA reverse-transcribed.

### Correlation analysis
- **Internal contradiction (confirmed):** Methods says *"Correlation statistical analyses are performed by Pearson's correlation"*; Figure 3M legend says *"Pearson correlation values were used for the matrix"*; but the Discussion says *"we conducted a correlation analysis between the relative abundance of PS-associated gut microbiota and obesity-related indicators using **Spearman's** correlation."* Two of three say Pearson. Unresolved.
- Scale: top 15 bacterial families × 17 obesity/inflammation parameters (Fig. 3M x-axis: weight gain, liver weight, iWAT, pWAT, eWAT, LDL-C, HDL-C, TC, LDL-C/HDL-C, colon MUC2 protein, ZO-1 relative expression, fasting glucose, liver TG, liver TC, IL-1β, IL-6, TNF-α) ≈ **255 coefficients, with no multiplicity correction**, on n=5 mice per group (20 animals total for the microbiota side).

### ML / random forest
**None.** No machine learning, no random forest, no classifier, no predictive model, no cross-validation anywhere in the paper.

### Other absences
No power calculation, no pre-registration, no exclusion criteria (despite n varying 6→5→4→3 across endpoints), no cage random effect, no covariate adjustment.

## 5. Does the model suit the biological question?

**No — the central claim is a textbook difference-of-significance error.** The design is explicitly factorial (diet × particle), the headline claims are explicitly interaction claims — *"PS-MS exacerbated HFD-induced obesity"*, *"the lipid accumulation in the liver and white adipose tissue after PS-MS exposure was more pronounced in HFD-fed mice but not in ND-fed mice"*, *"a superimposed effect between PS-MS and HFD on intestinal integrity and flora"* — yet the analysis is a one-way ANOVA over four groups. The evidence offered for "exacerbation" is *significant in the HFD stratum, non-significant in the ND stratum*. The difference between those two simple effects was never tested, and no interaction F/p exists in the paper. With n=6 per cell, the study is in any case badly underpowered to detect an interaction even if one were fitted.

Three further mismatches between model and question:

1. **The mediation argument is structurally the same error, one level up.** ABX mediation is inferred from PS effects being *absent* in a separate antibiotic-treated cohort. An absence of significance in cohort B is not evidence that cohort A's effect was microbiota-dependent; that requires an ABX × PS interaction within one randomised experiment, or a formal mediation model. Neither was done, and batch/time differences between the three 24-mouse experiments cannot be excluded.
2. **The picture is not uniformly diet-conditional, which the one-way framing hides.** Chao1/observed species rose with PS-MS **in both ND and HFD**; goblet cells, Muc2/MUC2, Tjp1/TJP1 and occludin fell with PS-MS **in both ND and HFD**. So the barrier mechanism is *not* diet-dependent while the obesity phenotype is claimed to be — a dissociation that an interaction model would have had to confront explicitly.
3. **Cage is not modelled.** 3 mice/cage, 2 cages/group, coprophagy homogenises cagemates; for a microbiota-mediated phenotype the effective n is plausibly closer to 2 per group than 6. Treating mice as independent inflates every p-value in the 16S analysis, which is the evidentiary core of the paper.

Minor but in the same direction: the AB/PAS goblet-cell "n=9" is **9 crypts (3 slices × 3 crypts)**, not 9 mice — explicit pseudoreplication in Figs. 5A and 7A.

## 6. Inputs, features, outputs, targets per analysis

| Analysis | Input | Features | Output / target |
|---|---|---|---|
| Phenotype (Figs. 1, 4, 6) | per-mouse weights, calipered/weighed tissues, glucometer, colorimetric kits | BW, BW gain, food & energy intake, liver/iWAT/eWAT/pWAT mass, fasting glucose, OGTT AUC, liver TC/TG, serum TC/TG/LDL-C/HDL-C, LDL-C:HDL-C | group means compared by one-way ANOVA + Tukey (FMT: t-test) |
| Histology | paraffin H&E (5 µm), frozen Oil Red O (8 µm), AB/PAS | MASLD activity score (Kleiner: steatosis 0–3 by % hepatocytes with droplets, ballooning 0–2, lobular inflammation 0–3 by foci per ×200 field), Oil Red O % area, iWAT adipocyte area, goblet cells/crypt | ordinal/continuous scores, Image-Pro Plus quantified |
| Colon qPCR | colon total RNA | Muc2, Tjp1 Ct vs Gapdh | 2^−ΔΔCT relative expression |
| ELISA | colon homogenate supernatant; serum | MUC2 (ng/g tissue); serum LPS, TNF-α, IL-6, IL-1β, IL-17A (pg/mL) | concentrations from standard curves (all Jiangsu Meimian kits) |
| Immunofluorescence | frozen distal colon, 5 µm | TJP1/ZO-1 and occludin, Cy3 secondary, Hoechst nuclei, Nikon Ti2 confocal | relative fluorescence intensity (ImageJ), n=3 |
| 16S (Fig. 2) | per-mouse fecal DNA (Omega Bio-tek Stool DNA Kit), V3–V4, 338F/806R, barcoded, MiSeq PE300 (Allwegene, Beijing) | 97% OTU table; taxa at phylum→genus | Chao1/observed species/Shannon/Simpson; PCA; F:B ratio; family and genus relative abundances; LEfSe cladogram (Fig. 2E) |
| Correlation (Fig. 3M) | family relative abundances (n=5/group) + 17 host parameters | 15 families × 17 phenotypes | correlation matrix heat map, r from −0.4 to +0.8 |
| ABX arm | same panel in microbiota-depleted mice | as above | null results used as mediation evidence |
| FMT arm | pooled donor fecal slurry → ABX-pretreated HFD recipients | same phenotype + barrier + cytokine panel | transfer of phenotype; engraftment by PCA (Fig. S5) |

## 7. How results were evaluated

Entirely by **frequentist null-hypothesis significance testing at p<0.05 on pairwise group contrasts**, with stars (`*p<0.05, **p<0.01, ***p<0.001, ****p<0.0001`) and no effect-size estimates, no confidence intervals, and no adjusted p-values beyond Tukey within each one-way ANOVA. There is no hold-out, no cross-validation, no replication cohort, no external validation dataset, no sensitivity analysis.

Effect sizes appear only as ad hoc percentages in the Results prose: *">17% more BW (BW relative to week 0: 158.17%±8.25 vs. 175.21%±5.35)"*; *"76% greater inguinal white adipose tissue (iWAT) mass, 64% greater perirenal white adipose tissue (pWAT) mass, and 39% greater epididymal white adipose tissue (eWAT) mass"*; *"The Firmicutes-to-Bacteroidetes (F:B) ratio was less than one-half"*; *">97% of the intestinal bacteria were eradicated"*. The FMT illustration (Fig. 6A long description) gives recipient BW gains of **10.87 g (R-HFD) vs 13.95 g (R-HFD+PS-MS)**.

All underlying numbers are deferred to Table S3, which I could not read.

## 8. Independent validation vs omics-only

**Validated by independent, non-sequencing measurement (strong):**
- Obesity phenotype — direct weighing, repeated over 8 wk, plus OGTT, plus enzymatic lipid kits, plus H&E/Oil Red O histology, plus adipocyte morphometry. Multiple orthogonal modalities.
- Colon barrier — **four independent modalities agreeing**: AB/PAS goblet-cell histology, Muc2 mRNA (qPCR), MUC2 protein (ELISA), TJP1/occludin protein (immunofluorescence), plus Tjp1 mRNA. This is the best-corroborated mechanism in the paper.
- Systemic inflammation/endotoxemia — ELISA only (LPS, TNF-α, IL-6, IL-1β, IL-17A); single modality, single supplier, no western blot or multiplex cross-check.

**Causally supported by the ABX arm:** that the PS-MS-associated differences in obesity outcomes, barrier proteins and serum cytokines **require an intact microbiota** — stated as *"the PS-induced obesity outcomes were not detected in mice fed with either the ND or the HFD after receiving ABX treatment"* and, for barrier, *"this lower expression depended on the presence of gut flora"*. Caveats: depletion was incomplete (authors' own limitation), the cohort is separate, and the inference is from a null.

**Causally supported by the FMT arm (the strongest causal claim):** that the **PS-MS-shaped microbiota is itself sufficient** to transfer, into HFD-fed ABX-pretreated recipients, greater BW/BW gain, liver weight, MASLD score, iWAT/eWAT/pWAT mass, adipocyte area, serum TC, LDL-C and LDL-C:HDL-C, lower goblet cells / MUC2 / Muc2 / Tjp1, and higher serum TNF-α and IL-17A — *without a food-intake difference* (Fig. 6D), which rules out hyperphagia as the transferred mechanism. Caveats: recipients are ABX-conditioned, not germ-free (authors concede *"Transplanting the gut microbiota to germ-free mice would provide more compelling evidence"*); donor feces were **pooled across donors and across the whole 8 wk**, so recipients within a group share one inoculum and are not independent; and there is **no ND-recipient arm**, so microbiota-mediated effects on a normal-diet background were never tested.

**Correlative / omics-only (16S) and NOT causally supported:**
- Every specific taxon attribution — Lachnospiraceae, Oscillospiraceae, Marinifilaceae, Akkermansiaceae, Deferribacteres, Desulfovibrio, Bifidobacterium, Parabacteroides, Atopobiaceae. No gnotobiotic mono-colonisation, no taxon-specific intervention, no supplementation experiment. The FMT transfers a *whole community*, so it cannot attribute the effect to any named taxon.
- The entire LPS-translocation mechanism. **Serum LPS was measured only in the main experiment (Fig. 3H); it is not among the ABX panel (Fig. 5H–K: TNF-α, IL-6, IL-1β, IL-17A) nor the FMT panel (Fig. 7E–G: TNF-α, IL-17A, IL-6).** So the mechanistic link "barrier loss → LPS translocation → inflammation" is never tested in either causal arm.
- The Desulfovibrio-as-LPS-source hypothesis is pure literature inference (*"exhibiting endotoxin activity that is 1,000-fold higher than that of other gut bacteria"*, citing ref. 93).
- **SCFAs were never measured** (no occurrence of SCFA, butyrate, acetate, propionate or "short-chain" in the body text), despite the paper's Firmicutes/Lachnospiraceae framing.
- No external validation dataset of any kind.

## 9. Main findings with directions

**Body weight / adiposity.** ND vs ND+PS-MS: no difference (*"mice treated with PS-MS gained a similar amount of weight under ND"*). HFD+PS-MS vs HFD: BW **↑** >17% (158.17%±8.25 → 175.21%±5.35 of week-0 BW; Fig. 1B,C), *"despite consuming less high-fat feed and calories"* (Fig. 1D) — i.e. not hyperphagia-driven. iWAT **↑76%**, pWAT **↑64%**, eWAT **↑39%**; iWAT adipocyte area **↑** (Fig. 1M,N). No adipose differences under ND. Water intake unchanged (Fig. S2).
*(Minor panel-labelling slip: the Results text cites "(Figure 1J–L)" for iWAT/pWAT/eWAT in that order, but the Fig. 1 legend assigns J = perirenal, K = epididymal, L = inguinal.)*

**Hepatic steatosis.** Liver weight **↑** in HFD+PS-MS vs HFD, unchanged ND vs ND+PS-MS (Fig. 1E). MASLD histological activity score **↑** and Oil Red O lipid area **↑** (Fig. 1F–I). Liver TC **↑** and liver TG **↑** (Fig. 1P,Q). Steatosis-dominant: the Kleiner-derived score as defined in Methods is **steatosis + ballooning + lobular inflammation only — no fibrosis stage**, and no Sirius red/trichrome was done.

**Lipogenesis genes (Fasn, Acly, Acaca, Srebf1/SREBP-1c): NOT MEASURED.** Zero hits for any of these in the body text. There is **no liver molecular work of any kind in this paper** — no hepatic RNA-seq, no hepatic qPCR, no hepatic western blot. The only qPCR in the study is colonic Muc2 and Tjp1. Lipogenic direction is inferred solely from lipid *content*.

**PPARα / fatty-acid oxidation (Ppara, Cpt1a, Acox1, Cyp4a): NOT MEASURED** (zero hits; the single "PPAR" string match is inside the word "apparent").

**Cholesterol synthesis:** no mechanistic measurement (no Hmgcr, Srebf2, Ldlr). Phenotypically, liver TC **↑**, serum TC **↑**, serum LDL-C **↑**, LDL-C/HDL-C **↑** (Fig. 1S,T,V). **NULL for serum TG and serum HDL-C** — *"There was no difference observed in TG and HDL-C with PS-MS supplementation when compared with the HFD group (Figure 1R,U)"*.

**Glucose.** Fasting glucose **↑** and OGTT AUC **↑** in HFD+PS-MS vs HFD; *"no significant difference in fasting blood glucose levels and blood glucose AUC between the ND group and the ND + PS-MS group"* (Fig. 1O). OGTT: 6-h fast, 2 g/kg BW oral glucose, tail vein at 0/15/30/60/90/120 min.

**Inflammation (Tnf, Il1b, Il6, NLRP3, LPS/endotoxemia).** Serum **LPS, IL-1β, IL-6, TNF-α and IL-17A all ↑**, but note the *exact* contrast in the Results: *"serum levels of LPS, IL-1β, IL-6, TNF-α, and IL-17A were significantly higher in the **HFD + PS-MS group than in the ND group** (Figure 3H–L)"* — **the text compares HFD+PS-MS against ND, not against HFD**, so the prose does not establish an added PS effect on top of HFD for these analytes. The Discussion nevertheless speaks of *"the disparities in serum inflammatory cytokine levels between the HFD and the HFD + PS-MS groups"*, implying a HFD contrast that the Results text never states. Whether HFD+PS-MS exceeds HFD can only be settled from Table S3 / the plotted significance bars. **NLRP3 was not measured** (zero hits); no NF-κB, no hepatic inflammatory markers, no macrophage markers, no TLR4 measurement (TLR4/NF-κB appears only as cited literature).

**Gut barrier genes.** Goblet cells per crypt **↓**, Muc2 mRNA **↓**, MUC2 protein **↓**, Tjp1 mRNA **↓**, TJP1 protein **↓**, occludin protein **↓** — and critically *"regardless of whether they were fed a ND or a HFD"* / *"in both the ND and HFD groups"* (Fig. 3A–G). **This core effect is diet-independent.** No Cldn1/Cldn2/Muc1/Tff3 measured; no intestinal permeability assay (no FITC-dextran).

**SCFAs: NOT MEASURED AT ALL.**

**Specific taxa (16S, n=5/group).**
- α-diversity: **Chao1 ↑ and observed species ↑ with PS-MS in BOTH ND and HFD** (Fig. 2A; Fig. S3A). Shannon: *"no significant difference among the ND, ND + PS-MS, and HFD groups"* (Fig. S3B) — a null the authors foreground, concluding the effect is on "abundance rather than ... diversity".
- **F:B ratio ↓** in HFD+PS-MS (less than half that of HFD; Fig. 2D) — **opposite to the classic obesity-associated high F:B**, and the authors concede the obesity was *"independent of the F:B ratio"*.
- Family **↑**: Lachnospiraceae, Oscillospiraceae, Marinifilaceae, Akkermansiaceae, Deferribacterales (Fig. 2E–G,I,J); **↓**: Bifidobacteriaceae, Atopobiaceae (Fig. 2H,K).
- Genus **↑**: Desulfovibrio (*"in the HFD, but not the ND groups"*); **↓**: Bifidobacterium, Parabacteroides (Fig. 2L–N).
- **Abstract/Results mismatch (new, not in the prior notes):** the Abstract lists **Bacteroidaceae** among the taxa increased by HFD+PS-MS, but no Figure 2 panel for Bacteroidaceae is described in the Results — it appears in the paper only in the Fig. 3M correlation heat map and the Fig. 2E cladogram key. The Abstract also writes "Deferribacteres" where the Results write "family Deferribacterales" (an order name used as a family); taxonomic ranks are handled loosely throughout.
- Correlations (Fig. 3M): Streptococcaceae, Lachnospiraceae, Oscillospiraceae, Desulfovibrionaceae, Marinifilaceae **positively** associated with BW gain, fat depots, serum/liver lipids, fasting glucose and cytokines and **negatively** with MUC2; Prevotellaceae, Muribaculaceae, Erysipelotrichaceae, Bifidobacteriaceae, Atopobiaceae in the opposite direction; Bacteroidaceae positive with IL-1β and IL-6.
- **Internal contradiction (confirmed):** Akkermansiaceae was reported **increased** by PS-MS under HFD (Fig. 2J) yet *"f__Akkermansiaceae showed significant negative associations with BW gain, liver weight, epididymal adipose, fasting glucose, and serum TC"* — the taxon moves opposite to the phenotype it supposedly predicts. The Discussion devotes a long paragraph to Akkermansia but never reconciles the sign.

**ABX arm (Figs. 4–5) — the deliberate nulls.** >97% bacterial eradication (Fig. S4). No PS effect on BW gain, liver weight, liver fat, MASLD score, iWAT/eWAT/pWAT, adipocyte area, or serum lipids (Fig. 4C–R); no PS effect on goblet cells, Muc2 mRNA/protein, Tjp1 mRNA, TJP1 or occludin protein (Fig. 5A–G); serum cytokine differences abolished (Fig. 5H–K). **Two awkward in-arm results the Discussion only partly addresses:** (i) food and energy intake were significantly **HIGHER** in ABX HFD+PS-MS than ABX HFD while BW did not differ (Fig. 4D) — the authors attribute this to absent microbial fat absorption; (ii) *"the blood glucose AUCs in the ND + PS-MS and HFD + PS-MS groups were significantly **lower** than in the ND and HFD groups, respectively"* (Fig. 4M) — i.e. under ABX, PS-MS looks *protective* for glucose tolerance, a direction reversal the Discussion never mentions. Also, ABX did not abolish the HFD effect itself: *"ABX treatment was not completely resistant to HFD-induced weight gain"*, which the authors read as microbes being contributory but *"not the sole determining factor"*.

**FMT arm (Figs. 6–7).** R-HFD+PS-MS vs R-HFD: BW/BW gain **↑** (10.87 g → 13.95 g), liver weight **↑**, MASLD score **↑**, iWAT/eWAT/pWAT **↑**, adipocyte area **↑**, serum TC **↑**, LDL-C **↑**, LDL-C/HDL-C **↑**; food and energy intake **unchanged**; fasting glucose **↑** but **glucose AUC NULL** (*"no significant difference in blood glucose AUC"*, Fig. 6M — explicitly flagged by the authors). Colon: goblet cells **↓**, MUC2 protein **↓**, Muc2 mRNA **↓**, Tjp1 mRNA **↓** (Fig. 7A–D). Serum TNF-α **↑**, IL-17A **↑**, *"However, no difference was found for IL-6 (p=0.086)"* (Fig. 7E–G) — a rare explicit non-significant p. Statistics here: Student's t-test by default, Fig. 6B by **multiple unpaired t-tests**, Fig. 6N serum TC by **Mann–Whitney**.

## 10. Limitations

### (a) Actual limitations of what WAS done

1. **No interaction test for an explicitly interaction-shaped claim.** One-way ANOVA + Tukey across four cells; "exacerbation"/"superimposed effect" rests entirely on significant-in-HFD vs non-significant-in-ND. No F, no p, no effect size for diet × particle.
2. **Mediation inferred from nulls in separate cohorts.** ABX is a different 24-mouse experiment, not a crossed factor; no three-way model, no formal mediation. Batch/time confounding between the three experiments cannot be excluded and is never discussed.
3. **Cage ignored.** 3 mice/cage, 2 cages/group, analysed as 6 independent animals, in a study whose entire thesis is microbiota transfer between co-housed, coprophagic animals. The 1-wk pre-study co-housing helps the baseline but does nothing for the 8-wk exposure period.
4. **Pseudoreplication in the goblet-cell data**: "n=9" is 9 crypts (3 slices × 3 crypts), not 9 mice (Figs. 5A, 7A).
5. **Unexplained, undeclared n attrition**: 6 → 5 (16S, MASLD main) → 4 (MASLD in ABX and FMT; Oil Red O 4–5) → 3 (immunofluorescence). No exclusion criteria, no statement of which animals were dropped or why.
6. **No β-diversity statistic.** Compositional separation is asserted from a PCA scatter alone; no PERMANOVA/ANOSIM/adonis. PCA (not PCoA) on OTU abundances is also a questionable ordination choice for compositional microbiome data.
7. **No multiple-testing correction anywhere**: ~255 correlations in Fig. 3M, plus LEfSe/Metastats across hundreds of taxa, with no FDR, no stated LDA cut-off and no stated α for LEfSe.
8. **Pearson vs Spearman contradiction** between Methods/Fig. 3M and the Discussion.
9. **Antibiotic identity contradiction**: ampicillin (Methods) vs penicillin (Fig. 4 legend and Discussion).
10. **Incomplete microbiota depletion** — authors' own: *"this investigation had certain limitations, most notably our inability to fully eradicated the intestinal microbial population."*
11. **FMT into ABX-conditioned, not germ-free, recipients** — authors' own: *"Transplanting the gut microbiota to germ-free mice would provide more compelling evidence."*
12. **FMT inoculum pooled** across donors and across all 8 wk, destroying donor-level variation and making recipients within a group non-independent.
13. **FMT has no ND-recipient arm**, so the microbiota-mediated effect was never tested on a normal-diet background — exactly the stratum where the main experiment found no phenotype.
14. **Serum LPS, the linchpin of the proposed mechanism, was not re-measured in either the ABX or the FMT arm.**
15. **No blinding** for MASLD scoring, adipocyte area, goblet-cell counts or fluorescence intensity — all operator-dependent.
16. **Randomisation stated only for the FMT experiment.**
17. **Particle characterisation is minimal and internally inconsistent**: nominal "5 µm" but stated range "4.0 to 4.9 µm"; TEM morphology only (Fig. S1); **no measured size distribution/DLS, no zeta potential, no FTIR/Raman polymer confirmation, no surface chemistry (plain vs COOH vs NH2 unstated), no additive/leachate analysis, and no endotoxin testing** — the last being a direct threat to the study's primary mechanistic readout (serum LPS) given that commercial bead suspensions can carry endotoxin.
18. **Delivered dose unverified.** 10 mg/L nominal, bottles changed every 2 d, "shaken to ensure even distribution", but no measurement of actual suspended concentration over the interval — 5-µm beads sediment and adhere to bottle walls. The headline 36 µg/d/mouse figure is *derived* from an assumed 3.6 mL/d water intake, not from measured consumption of particles.
19. **Slightly mismatched vehicle handling**: controls got *"filtered and autoclaved water"*, PS groups got *"ultrapure sterile water"* plus 10 min sonication. No particle-free sonicated-vehicle control.
20. **Single dose, single particle size, single polymer, single sex (male), single strain** — no dose–response, which the authors concede matters because cited work reports *"relatively low doses of PS-MS were found to trigger obesity, whereas higher doses did not."*
21. **Gapdh as the sole reference gene** under HFD — authors' own caveat, mitigated only by protein-level agreement.
22. **Serum cytokine contrast as written is HFD+PS-MS vs ND**, not vs HFD, which if taken literally makes the inflammation result largely a diet effect.
23. **Direction reversal in the ABX arm** (better glucose AUC with PS-MS) left entirely unaddressed.
24. **Deposition is partial**: PRJNA990701 holds 20 amplicon runs, consistent with the main 2×2 only; the ABX-depletion and FMT-engraftment 16S data behind Figures S4 and S5 appear not to be deposited, and no analysis code is shared.
25. **Direct conflict with the companion 2023 EHP study** (0.45–0.53 µm PS-COOH), which the authors themselves cite as reporting **weight reduction** in HFD-fed mice, versus weight gain here. The authors rationalise this as a size/property difference without testing it.

### (b) Things simply not attempted (scope, not flaws)

- Any liver molecular work: no hepatic RNA-seq, microarray, proteomics, qPCR or western blot; hence no Fasn/Acly/Acaca/Srebf1, no Ppara/Cpt1a/Acox1/Cyp4a, no cholesterol-synthesis genes, no NLRP3/NF-κB.
- Shotgun metagenomics, metatranscriptomics, fecal or serum metabolomics, bile acids.
- SCFA quantification of any kind.
- Fibrosis staging, liver ALT/AST, insulin/HOMA-IR, insulin tolerance test, indirect calorimetry, energy expenditure, fecal energy loss.
- Intestinal permeability assay (FITC-dextran), mucus layer thickness measurement.
- Particle biodistribution/tissue burden (no fluorescent label, unlike the 2023 companion study), so no evidence on whether 5-µm beads reach the liver at all.
- Germ-free mice, gnotobiotic mono-colonisation, taxon-specific supplementation or depletion (e.g. Akkermansia, Bifidobacterium, Parabacteroides, Desulfovibrio) — i.e. no causal test of any named taxon.
- Female mice; other strains; other polymers (PE, PP, PVC); other sizes; nanoplastics; weathered/environmentally realistic particles; mixtures.
- Recovery/reversibility, transgenerational, or post-exposure follow-up.
- Dose–response design.
- External/independent-cohort validation.

---

## Study characteristics (requested checklist)

| Item | Value (verbatim where quoted) |
|---|---|
| Species/strain/sex/age | Mouse, **C57BL/6J, male only**, *"aged 6–8 wk"* at purchase (Vital River Laboratory Animal Center); exposure *"started at 8 wk of age"* (Fig. 1A/4A); FMT donors *"8-wk-old male"* |
| Polymer | Polystyrene microspheres (PS-MS); **surface chemistry not reported**; **no fluorescent label** |
| Particle size | Nominal **5 µm**; supplier-stated range *"4.0 to 4.9 µm"*; **no independently measured size distribution** (TEM morphology only, Fig. S1) |
| Supplier / lot | **Shanghai Macklin Biochemical Technology Co., Ltd.**, stock *"25 mg/mL"* aqueous suspension, *"5.0% wt/vol"* in deionized water, lot/catalogue *"identified as L815942"* |
| Dose & expression | **10 mg/L in drinking water** ("10-mg/L PS-MS ... at a dose of 10 mg/L per day"); authors' derived per-animal dose *"each control mouse consumed 3.6 mL of water/d ... resulting in a PS-MS intake of 36 µg/d per mouse"* (≈1.8 mg/kg BW/d for a 20 g mouse), equated to *"the minimum potential exposure dose for humans (0.1 g of MPs/wk)"*. **Not expressed as particles/L or particles/d anywhere.** Single dose level |
| Route / vehicle | Oral, ad libitum drinking water; vehicle = *"ultrapure sterile water"* for PS groups, *"filtered and autoclaved water"* for controls |
| Frequency / duration | Continuous in water; *"PS-MS water was replaced every 2 d"*; **8 wk** (all three experiments) |
| Dispersion | *"sonication for 10 min"*; *"Bottles were shaken to ensure an even distribution ... to prevent sedimentation"*; stock *"underwent thorough sonication, followed by dilution with water prior to their use"* |
| Diets (catalog nos.) | ND = **D12450B**, 3.85 kcal/g, 10% kcal fat, Research Diets Inc.; HFD = **D12492**, 5.24 kcal/g, 60% kcal fat, Research Diets Inc. Compositions in **Table S2**. (D12450B is the proper purified low-fat match for D12492 — a genuine methodological strength.) Ad libitum, **not pair-fed** |
| Sampling / fasting state | *"the mice underwent an overnight fasting period before cervical dislocation euthanasia"*; blood by **eyeball sampling**, centrifuged 3,000 rpm 15 min 4 °C for plasma, −80 °C. OGTT at week 8 after a **6-h** fast. Feces collected **every 4 wk** from each individual mouse. **Time of day of sacrifice not stated** |
| Tissues | Liver, distal colon, inguinal (iWAT), epididymal (eWAT) and perirenal (pWAT) adipose — excised and weighed; one piece in paraformaldehyde, one in O.C.T., remainder snap-frozen in liquid N₂ at −80 °C |
| Housing / cage | *"Groups of mice were housed together (n=6 per group, 3/cage)"*; 24±1 °C, 55%±5% RH, 12 h light:dark; 1-wk co-housing of all groups before study start; **cage not used as a unit of analysis, no cage random effect** |
| Randomisation | Stated **only** for FMT (Fig. 6 legend: *"The mice were randomly divided into four groups (n=6)"*) |
| Blinding | **Not reported anywhere** ("blind" has zero occurrences) |
| Endotoxin / surfactant checks | **None.** No endotoxin assay, no surfactant/dispersant disclosure, no leachate analysis. ("surfactant", "zeta", "FTIR", "Raman" all zero occurrences; the single "endotoxin" occurrence is a literature claim about *Desulfovibrio* LPS) |
| Ethics | China Agricultural University Laboratory Animal Welfare and Animal Experimental Ethical Committee, approval No. **AW72403202-1-31** |
| Funding | National Key R&D Program of China (2022YFF1100102, 2022YFF1100104); NSFC No. 31625025; 2115 Talent Development Program, China Agricultural University |
| Competing interests | *"The authors affirm that they have no known financial interests or personal relationships that could have potentially influenced or biased the findings"* |
| DOI | https://doi.org/10.1289/EHP13913 |

---

## Prior-note claims re-checked

**Confirmed as written** (spot-checked against the body text): the three-experiment / 72-mouse structure; all group sizes and the n=5/4/3 reductions; housing 3/cage and the 1-wk co-housing; the water-handling asymmetry; randomisation only in Fig. 6; no blinding; particle supplier/size/lot and the 5 µm vs 4.0–4.9 µm mismatch; the absence of endotoxin, zeta, DLS, FTIR/Raman and surface-chemistry data; the 10 mg/L dose and 36 µg/d derivation; both diet catalogue numbers; the verbatim Statistical Analysis paragraph; **no interaction term and no interaction p anywhere**; the two-way RM ANOVA being group × time only; the 16S pipeline list and the missing LEfSe thresholds/FDR; the Pearson-vs-Spearman contradiction; the >17% BW figure and the 158.17%/175.21% values; the 76/64/39% fat-depot percentages; the serum TG and HDL-C nulls; Chao1 up in both diets and the Shannon null; the F:B reversal; all named taxa directions; the Akkermansiaceae sign contradiction; barrier loss in **both** diets; the serum-cytokine contrast being stated as "HFD+PS-MS vs ND"; every ABX null plus the higher-food-intake and lower-glucose-AUC reversals; every FMT result including the AUC null and IL-6 p=0.086; the absence of any liver molecular work, of SCFAs, of RNA-seq/metagenomics/metabolomics, of NLRP3/PPARα/lipogenesis genes; all five author-acknowledged limitations; and the verbatim data-availability statement.

**Prior-note items that were incomplete, and are corrected here:**
1. **ABX drug identity** — the prior notes recorded only "ampicillin" (Methods). The Figure 4 legend and the Discussion both say **penicillin**. This internal contradiction was missed.
2. **PRJNA990701 was marked "existence NOT verified"** — it is now verified to resolve, with the metadata and **20 amplicon runs** listed in §3; and the run count indicates the **ABX and FMT 16S data are not deposited**, which the prior notes did not raise.
3. **Abstract-only taxon** — the prior notes' taxon list omitted **Bacteroidaceae**, which the Abstract claims was increased but which has no corresponding Figure 2 panel in the Results.
4. **FMT statistics** — the prior notes said "Student's t-test"; more precisely, Fig. 6B used **multiple unpaired t-tests** and Fig. 6N (serum TC) used **Mann–Whitney**, with Student's t-test as the default elsewhere.
5. **Figure 2 panel-specific test** — panels H, J, K, M, N used the **Brown–Forsythe and Welch ANOVA**, not Tukey one-way; the prior notes mentioned Brown–Forsythe only as a generic Methods sentence.
6. **Fig. 1 panel-letter slip** in the Results text for the WAT depots (text "Figure 1J–L" in iWAT/pWAT/eWAT order vs legend J=perirenal, K=epididymal, L=inguinal) — not previously noted.
7. **Serum LPS absent from the ABX and FMT panels** — a materially important gap for the proposed mechanism that the prior notes did not flag.
8. **No β-diversity statistic at all** (no PERMANOVA/ANOSIM/adonis; PCA not PCoA) — the prior notes listed the pipeline but did not state this explicitly as an evidentiary gap.

**No prior-note claim was found to be factually wrong.** The eight items above are omissions or under-specifications, not errors.

---

## URLs/files actually read

1. `C:\Users\Soham\AppData\Local\Temp\claude\d28fb66f161d523b\5dbfd995-fe63-4ce9-8501-246b4fa39b6f\scratchpad\dl\PMC11370995.txt` — Europe PMC JATS full text converted to text, lines 1–182 (entire body up to `### References`). **Primary source, read in full.**
2. `C:\Users\Soham\Desktop\microplastic\research\literature\papers\ehp2024_ps_hfd_obesity.md` — prior session notes (re-checked, not a source).
3. https://www.ncbi.nlm.nih.gov/bioproject/PRJNA990701 — BioProject record (title, submitter, registration date, 20 linked SRA experiments).
4. https://www.ebi.ac.uk/ena/browser/api/summary/PRJNA990701 — ENA project summary (public status, first-public date 2023-07-26, centre).
5. https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA990701&result=read_run&fields=run_accession,sample_title,library_strategy,instrument_model,read_count — ENA run listing (20 runs, AMPLICON, Illumina MiSeq, accessions SRR25114495–SRR25114514 and SRR29975743).

**Attempted and failed (so NOT read):**
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11370995/ — returned a reCAPTCHA browser-verification interstitial; no article or supplementary content retrieved.
- https://ehp.niehs.nih.gov/doi/10.1289/EHP13913 — DNS resolution failed (`ENOTFOUND ehp.niehs.nih.gov`) from this environment.
- Consequently **Supplemental Material (Figures S1–S5, Tables S1–S3) was never accessed.**

## UNVERIFIED items

1. **All of Table S3** (per-animal numeric source data for every figure). In particular: whether serum LPS/IL-1β/IL-6/TNF-α/IL-17A in HFD+PS-MS are significantly above **HFD** (not merely above ND); the exact p-values and plotted significance bars for every panel; and whether the data would support a diet × particle interaction if re-analysed.
2. **Table S1** — qPCR primer sequences for Muc2, Tjp1, Gapdh (primer specificity/efficiency unassessable).
3. **Table S2** — full diet compositions (sucrose, cholesterol, fibre, protein source matching between D12450B and D12492).
4. **Figure S1** — TEM images of the PS-MS (actual morphology, aggregation state).
5. **Figure S2** — water intake data supporting "no significant difference in daily water consumption".
6. **Figure S3A/S3B** — observed-species and Shannon panels.
7. **Figure S4** — the ">97% of the intestinal bacteria were eradicated" evidence, including what was quantified (16S qPCR? culture? sequencing depth?) and in how many mice.
8. **Figure S5** — the FMT engraftment PCA; the number of recipient samples sequenced and whether donor and recipient communities genuinely overlap.
9. **Which cohort each of the 20 deposited SRA runs belongs to** — I verified the count and library type but could not inspect per-run sample attributes, so "ABX/FMT 16S not deposited" is a strong inference from arithmetic (20 = 4 × 5), not a direct confirmation.
10. **Whether the three 24-mouse experiments were run concurrently or sequentially** — never stated; batch confounding cannot be assessed.
11. **Time of day of euthanasia** — not stated; only "overnight fasting".
12. **How many mice the "n=9" crypts in Figs. 5A and 7A came from** ("3 slices/group" implies possibly 3 mice, but this is not stated).
13. **Which animals were excluded** to produce n=5, n=4 and n=3 subsets, and on what criteria.
14. **The ampicillin-vs-penicillin discrepancy** cannot be resolved from the paper.
15. **The Pearson-vs-Spearman discrepancy** cannot be resolved from the paper.
16. **Whether the Abstract's "Bacteroidaceae" increase corresponds to any tested panel** — no Figure 2 panel for it is described.
17. **Delivered (as opposed to nominal) PS-MS concentration in the drinking bottles** — never measured by the authors, so unverifiable in principle from this publication.
18. **Endotoxin content of the Macklin bead suspension** — never measured; unverifiable from this publication.
# 5. IJBS 2025 — Ahn et al., "Polystyrene Nanoplastics Exacerbate HFD-induced MASLD by Reducing Cathepsin Activity and Triggering Large Vacuole Formation via Impaired Lysosomal Acidification"

Int J Biol Sci 2025; 21(9): 3867–3885. DOI 10.7150/ijbs.108268. PMID 40607257. PMCID PMC12210384.
Authors: Ahn J, Ryu K, Kim H, Seo HW, Jang M, Kim SH, Park Y, Son MJ, Lee HJ, Kwon OS, Chung KS (KRIBB, Daejeon). Received 2024-12-05; Accepted 2025-05-11. "The authors have declared that no competing interest exists."

**ACCESS LEVEL: FULL TEXT read** — complete Europe PMC full-text body (Abstract, Introduction, Materials and Methods, Results incl. all main figure legends Fig. 1–7, Discussion, Conclusions, Acknowledgments, Availability of data and materials, Ethics, Author contributions, Abbreviations), i.e. everything up to the "### References" marker (lines 1–188 of PMC12210384.txt). **Supplementary figures/tables (Fig. S1–S5, Tables S1–S4) were NOT read** — they are a separate file not present in the download directory. Every S-figure/S-table statement below is therefore only as *described in the main text*, and is marked as such.

---

## 1. Exact scientific question

Stated in Introduction, final paragraph: "we aimed to explore the combined effects of obesity and NPs on human health by the concurrent administration of a HFD and PS-NPs to mice and human liver cell lines." Operationally: (a) does daily oral 50 nm carboxyl-PS nanoplastic exposure aggravate HFD-induced MASLD in mice, over 8 vs 20 weeks; and (b) what is the cellular mechanism — the stated hypothesis being PS-NP endocytosis → lysosomal accumulation → impaired cathepsin maturation/activity → blocked autophagic flux → large vacuoles → hepatic lipid accumulation.

Prior note on the question: **CONFIRMED verbatim.**

## 2. Experimental design

### In vivo (Methods, "Animal studies and experimental design"; Fig. 1A legend)
2 × 2 factorial *layout* (diet: CD vs HFD; particle: vehicle vs PS-NPs), run at two durations (8 or 20 weeks).
Verbatim: "Following a 2-week acclimation period, 4-week-old mice were grouped as follows: i) control diet (CD, n=6, Teklad Global 18% protein rodent diet; Teklad 2018S; Envigo, Madison, WI, USA); ii) control diet with PS-NPs (0.5 mg) per day (CD + PS-NPs, n=6); iii) HFD (n=6) with 45 kcal% fat (D12451; Research Diet, New Brunswick, NJ, USA), or iv) HFD with PS-NPs (HFD + PS-NPs, n=6) for 8 or 20 weeks."

- n = 6 per group (stated in Methods and repeated in Fig. 1A legend for all four groups). Whether the 8-week and 20-week cohorts are separate animals is **not stated**; "for 8 or 20 weeks" implies separate cohorts (so 24 + 24 mice) but is never made explicit, and no total animal number is given anywhere.
- Randomisation: "Mice were randomly assigned to four groups" (Results §1); Fig. 1A legend: "Mice were randomly divided into four groups". No method of randomisation (sequence generation, allocation concealment) given.
- Blinding: **not reported anywhere** (checked Methods incl. Histochemical evaluation and Statistical analysis, and all figure legends).
- Housing: "To ensure that the PS-NPs were dosed accurately, mice were individually housed and weighed before and after the study" — individual housing, so cage = animal (no cage-level pseudoreplication). No room temperature, light cycle, or bedding information.
- Vehicle control: "control jelly cube" (Fig. 1A legend), i.e. the same 2.5% agarose cube containing 0 mg PS-NPs — a well-matched vehicle.

### Omics sample size (Fig. 2; Methods "Total RNA extraction, real-time PCR, and transcriptome analysis")
**n per group for RNA-seq is NEVER STATED.** Methods say only "mRNA was sequenced by ebiogen Inc. (Seoul, Republic of Korea). Detailed information on the methods is provided in Table S4." Fig. 2B legend says only "Each column represents a sample", with no count. Nothing in Results gives a number. It is also not stated whether all four groups were sequenced, from which duration (the Fig. 2 title says "over 20 weeks", so 20-week livers), nor whether all n=6 per group or a subset were used, nor batch/time-of-sacrifice information.
Prior note item on this: **CONFIRMED (RNA-seq n is genuinely absent).**

### In vitro (Methods "Cell culture", "Crystal violet staining", "Assessment of neutral lipid…", "Determination of CYP3A4 and cathepsin B activity", "Triglyceride quantification assay", "Labeling with organelle-specific tracers and uptake", "Immunoblot analysis and immunostaining"; Results Figs. 3–7)
- Cell-line screen for vacuolisation/lipid accumulation: crystal violet in "indicated hepatocyte cell lines" (Fig. 3A) and BODIPY in **four** lines — "HepaRG, HepG2, SK-HEP-1, and 293T" at 100 µg/mL (Fig. 3B). Main text Results names "HepaRG, HepG2, and SK-Hep1, as well as the 293T kidney cell line". Methods additionally list Hep-G2/2.2.15, Hep3B, Huh7 and PLC/PRF/5 as cultured lines; **which of these actually appear in Fig. 3A cannot be verified from text** (figure images not readable here).
- Main model: HepaRG only. "The formation of large vacuoles and a substantial increase in lipid accumulation were only observed in HepaRG cells (Fig. 3); therefore, we performed further analyses using HepaRG cell line."
- Dose–response: 0, 10, 20, 50, 100 µg/mL PS-NPs for 24 or 48 h (Methods, BODIPY/Nile Red section); crystal violet at 0/50/100 µg/mL 24 h; organelle tracers at 0 or 50 µg/mL 24 h; flow cytometry at 0 or 100 µg/mL 24 h; TG assay 24 h. Viability: "the IC50 value exceeded 1,000 μg/mL (Fig. S3)" (main-text statement only).
- Oleate co-exposure: "HepaRG cells were pre-treated with 0.2 mM oleate for 3 h and then treated with PS-NPs for 24 h" (Methods; Fig. 5H–J, Fig. 4F).
- Endocytosis-inhibitor panel (Fig. 6A legend, verbatim): "5 mM 3-methyladenine (3-MA), a PI3K inhibitor; 1 μM cytochalasin B (CB), an inhibitor of actin polymerization; 5 μM 5-(N-ethyl-N-isopropyl) amiloride (5-AL), a micropinocytosis inhibitor; 5 μM chlorpromazine hydrochloride (CPZ), a clathrin-mediated endocytosis inhibitor; or 20 μM genistein, a caveolae-mediated endocytosis inhibitor, followed by treatment with 0 or 50 μg/mL PS-NPs for 24 h." (Wortmannin is listed among purchased chemicals in Methods but does **not** appear in the Fig. 6A legend — prior note listed wortmannin among inhibitors used; it is purchased, not demonstrably used in Fig. 6A. **Prior-note item partially wrong/unverifiable.**)
- Autophagy modulators: rapamycin 2 µg/mL pre-treatment (Fig. 7G, 7H) and 1 µg/mL (Fig. 7I); bafilomycin A1 10 nM for 3 h (Fig. 7H, 7J). Prior note's "rapamycin (1–2 ug/mL) and bafilomycin A1 (10 nM)": **CONFIRMED.**
- **In vitro replicate numbers (n) are not reported in any legend or in Methods.** No statement of independent biological repeats.

## 3. Datasets generated or re-used; accessions; deposition

Generated: one bulk liver mRNA-seq dataset (mouse, 20 weeks, four groups). Re-used: none — no public dataset, no external validation cohort, no re-analysis of GEO/TCGA data anywhere in the paper.

**Data Availability statement, verbatim (section "Availability of data and materials"):**
> "All data generated or analyzed during this study are included in this published article and its additional files."

**Accessions: NONE.** No GEO, SRA, BioProject, ENA, ArrayExpress, Zenodo, figshare or GitHub identifier appears anywhere in the full text (searched the whole body). Therefore **no repository page could be fetched or verified** — there is nothing to resolve. The transcriptome is effectively not deposited; the raw/processed RNA-seq data are not obtainable, and the DEG lists cannot be re-derived or re-thresholded by a third party.
Sequencing vendor: ebiogen Inc. (Seoul). Acknowledgments mention "the assistance of the KOBIC Research Support Program" (Korea Bioinformation Center) — a support programme, not a deposit.
Funding: KRIBB Research Initiative KGM5362521; NRF 2020R1I1A2072074 and RS-2024-00352135; NST GTL24021-000. Ethics: KRIBB-AEC-21279.
Prior note on data availability and accessions: **CONFIRMED verbatim.** (The prior note's extra claim that "The PMC record also carries an empty Data Availability Statement heading" could not be checked in the converted text — **UNVERIFIED**, and immaterial.)

## 4. Statistical / computational models

### Classical statistics (Methods, "Statistical analysis", verbatim)
> "Statistical significance was assessed using Tukey's multiple comparison tests or Student's t-tests following one-way analysis of variance using GraphPad Prism 9 software (San Diego, CA, USA). The error bars represent the standard error of the mean. Asterisks denote p values less than 0.05, 0.01, 0.001, or 0.0001."

- Model = **one-way ANOVA + Tukey** (or t-tests). Figure legends 1, 4, 5, 7 all say "one-way ANOVA and Tukey's multiple comparison tests". Fig. 5D legend is internally odd: "* p < 0.05 by ANOVA two-tailed unpaired t-test" — ambiguous as to which test was actually applied. (Prior note flagged this; **CONFIRMED.**)
- Reference group: Fig. 1E,F legend — "Comparisons are with the CD group unless otherwise noted." Fig. 5I/5J legend — "Comparisons are made with the 0 μg/mL PS-NPs in with or without 0.2 mM oleate group."
- **TWO-WAY ANOVA: not used. INTERACTION TERM (diet × particle): NOT FITTED ANYWHERE. No F statistic, no p value, no effect size for an interaction is reported anywhere in the main text, in any figure legend, or in Methods.** The four-group in vivo design is analysed only as a one-factor, four-level design with pairwise post-hoc contrasts. Error bars are SEM, not model-based CIs.
- The words "synergistically"/"synergistic effect" appear in the Abstract ("synergistically exacerbated the symptoms of steatosis") and Results ("HFD exerts a synergistic effect with PS-NPs-induced inflammation"), and the Methods BODIPY section says "To explore the synergistic effects between PS-NPs and oleate" — with **no formal synergy, interaction, or combination-index statistic behind any of them.**

### Differential expression (Results §2; Fig. 2C legend)
- Thresholds, verbatim: "DEGs were selected under the following conditions: P value < 0.05, |log2 fold change| > 0.5"; Fig. 2C legend repeats "(p value < 0.05 & |log2 FC| > 0.5)".
- **Multiple-testing correction at the DE step: NONE.** The threshold is a raw/unadjusted p-value. No FDR, no q-value, no adjusted-p is mentioned for DEG calling.
- **DE tool and version: NOT NAMED.** No DESeq2/edgeR/limma/NOISeq is mentioned; the method is deferred to Table S4 ("Detailed information on the methods is provided in Table S4"), which was not read. Likewise aligner, genome build, quantifier, normalisation method, platform, read length and depth are all absent from the main text.
- Design: **three separate pairwise contrasts** — CD+PS-NPs vs CD; HFD+PS-NPs vs HFD; HFD+PS-NPs vs CD+PS-NPs — plus a Venn overlap (Fig. 2D). **The HFD vs CD contrast is never reported**, so the magnitude of the diet main effect is unknown and no interaction contrast can be reconstructed from the published numbers.

### Enrichment
- "Differentially expressed gene sets were then used for over-representation analysis with the ClusterProfiler 4.0 package in R, using parameters set as adjusted p < 0.05 (Benjamini-Hochberg method). The biological functionality of the gene sets was investigated through Gene Ontology and Kyoto Encyclopedia of Genes and Genomes (KEGG) analyses." (Methods)
- So: ORA (not GSEA), clusterProfiler 4.0, BH-adjusted p < 0.05, GO + KEGG. Note Fig. 2F legend says "The dotted line represents the p value threshold of 0.05", which is a *nominal* p threshold on the bar plot — inconsistent with the BH-adjusted criterion in Methods. Universe/background gene set not specified.
- **No GSEA, no WGCNA, no co-expression modules, no deconvolution, no machine learning/classifier, no single-cell or snRNA-seq (hence no pseudobulk-vs-cell-level question), no spatial data, no molecular docking or MD simulation anywhere.**
- Unsupervised analysis present: PCA (Fig. 2A) and a heatmap of the 2,180-DEG union with z-scores −4 to +4 (Fig. 2B). Clustering method/distance/linkage for the heatmap is not stated. Note the heatmap is drawn on genes *selected* by the same data (selection-biased display of group separation).
- Colocalisation quantification (Fig. 5K, 6C–G, 7L) is reported descriptively; no Pearson/Mander's coefficients are given in text or legends.

Prior note on statistics: **CONFIRMED in full**, including the absence of an interaction term, the uncorrected DEG p-threshold, the unnamed DE tool, and the missing HFD-vs-CD contrast.

## 5. Does the model suit the biological question?

**No — the central "exacerbate/synergistic" claim rests entirely on pairwise tests and is exposed to the difference-of-significance fallacy.**

- The design is a genuine 2 × 2 (plus duration), which *is* the right design to ask whether PS-NPs modify the diet effect. But it is analysed with a one-factor ANOVA over four levels. The quantity of interest — (HFD+PS-NPs − HFD) − (CD+PS-NPs − CD) — is never estimated or tested.
- The headline Results sentence compares across *different reference groups*: "the HFD group demonstrated a 25% increase in body weight and a 20% increase in liver weight compared with the CD group, whereas the HFD+PS-NPs group displayed a 37% increase in body weight and a substantial 61% increase in liver weight compared with the CD+PS-NPs group" (Results §1, Fig. 1B). 20% vs 61% are HFD-vs-CD and HFD+PS-NPs-vs-CD+PS-NPs ratios respectively; the difference between those two ratios is exactly an interaction, asserted but not tested.
- The Abstract's "no effects were observed with PS-NPs on a normal chow diet" is **contradicted by the paper's own Results**: PS-NPs increased hepatic macrophages and neutrophils in CD-fed mice (Fig. S1C–F as described), increased collagen deposition in CD-fed mice (Fig. S1H–I as described), and produced 361 DEGs in CD+PS-NPs vs CD. So a "significant in HFD, not significant in CD" reading is not even internally consistent, let alone a tested interaction.
- Transcriptomics: three independent pairwise contrasts at uncorrected p < 0.05 cannot support a claim about modification of the diet effect. Comparing 361 DEGs (CD arm) with 628 DEGs (HFD arm) is a difference of counts, not a test; and the 16-gene Venn overlap has no null model attached.
- The exposure confound compounds this: the dose is a fixed 0.5 mg/animal/day while HFD animals end 25–37% heavier, so the mg/kg dose systematically diverges between diet arms over 20 weeks — a diet-correlated exposure gradient that an untested interaction cannot separate from a true biological interaction.
- The in vitro analogue (oleate + PS-NPs, Fig. 5H–I) is likewise pairwise: "Concurrent treatment with PS-NPs and oleate resulted in greater lipid accumulation compared to treatment with PS-NPs alone" — again no interaction test, despite Methods calling it an exploration of "synergistic effects".
- Secondary concerns: SEM error bars for inference; p<0.05 significance stars across many panels with no per-figure multiplicity accounting beyond Tukey within a panel; no power/sample-size justification; n=6 for a four-group factorial is low for detecting an interaction even if one were fitted.

## 6. Inputs, features, outputs, targets of each analysis

| Analysis | Input | Features | Output | Target/comparison |
|---|---|---|---|---|
| Body/liver weight (Fig. 1B) | Individually housed mice, 20 wk | body wt, perfused liver wt, liver:body ratio | group means ± SEM | 4 groups, one-way ANOVA + Tukey vs CD |
| Serum chemistry (Fig. 1E,F) | Serum collected "1 week before euthanization" (Methods); Fig. 1E,F legend says "after 19 weeks of treatment" | ALT, AST (340 nm), CHO (540/600 nm), TG (660/800 nm) | concentrations | same 4-group contrast |
| Histology (Fig. 1C,D; Fig. S1B,H,I) | 4 µm paraffin (H&E, Picro-Sirius Red); 10 µm OCT cryosections (Oil Red O, IF) | % Oil Red O-positive area by ImageJ; steatosis pattern; collagen; F4/80, MPO/neutrophil IF | area ratios, representative images | 4 groups |
| Bulk liver RNA-seq (Fig. 2) | 20-wk liver RNA (RNAlater + RNeasy), sequenced by ebiogen | gene-level expression (quantifier unstated) | PCA, heatmap of 2,180 DEGs, DEG counts, volcano, KEGG/GO bars | 3 pairwise contrasts; **no HFD vs CD** |
| ORA (Fig. 2E,F; S2) | the three DEG lists | GO BP / KEGG term membership | enriched terms, BH-adj p < 0.05 | per contrast |
| qRT-PCR | mouse liver and HepaRG RNA | 2^−ΔΔCT, reference gene **RPL13A** | relative expression | dose or group |
| Western blot | HepaRG lysates | band intensity / β-actin, ImageJ | fold change | dose or time |
| Cathepsin B activity (Fig. 4C) | HepaRG homogenate (Bioruptor) | fluorometric substrate, Ex/Em 400/505 nm | activity normalised to protein | ± PS-NPs |
| CYP3A4 activity (Fig. 5G) | HepaRG, 24 h | P450-Glo luminescence / ATP (CellTiter-Glo 3D) | normalised activity | dose |
| Intracellular TG (Fig. 5D) | HepaRG lysate, 5% NP-40 | OD 570 nm / total protein | TG level | dose |
| BODIPY/Nile Red imaging + flow (Fig. 5A–C) | HepaRG 24/48 h | mean fluorescence intensity (Cytation 5 Gen5; FACSCalibur + Flow Software 2) | lipid load | dose |
| Uptake inhibitor panel (Fig. 6A,B) | HepaRG + inhibitor, 50 µg/mL PS-NPs | dextran and LysoTracker intensity normalised to cell count | % inhibition | vs no-inhibitor |
| Organelle identity (Fig. 6C–G) | HepaRG, YG/GFP-PS-NPs, FITC-dextran, LysoSensor Y/B DND-160, LysoTracker, EEA1/RAB7/LAMP1/LAMP2/LGP85 | confocal colocalisation (qualitative) | vacuole identity and acidity | ± PS-NPs |
| Autophagy flux (Fig. 7I–K) | HepaRG ± rapamycin ± BafA1 | p62/β-actin ratio; LC3-II | flux inference from the BafA1-induced *difference* | ± PS-NPs |

## 7. How results were evaluated

Null-hypothesis significance testing with one-way ANOVA + Tukey (GraphPad Prism 9), SEM error bars, asterisk thresholds 0.05/0.01/0.001/0.0001. Omics evaluated by DEG counts at uncorrected p<0.05 with |log2FC|>0.5, plus BH-corrected ORA. Internal consistency across orthogonal assays is the main evaluation device (e.g. mRNA + protein + enzymatic activity for CYP3A4; mRNA + protein for MYC and C/EBPα; LysoTracker intensity + LAMP1 western). **No cross-validation, no held-out or external dataset, no independent replication cohort, no pre-registration, no formal effect-size or power reporting.** Histology quantification sampling: "Three spots from three lobes in each tissue section were randomly chosen to assess the relative area of fatty changes, and the ratio of the area of fatty changes to the whole tissue was calculated using ImageJ" (Methods, Histochemical evaluation) — selector and blinding unstated.

## 8. Independent experimental validation vs omics-only — and which mechanistic steps are causally tested

**Causally tested (perturbation experiments, all in vitro, all HepaRG):**
1. *Route of entry* — pharmacological endocytosis inhibition (3-MA, cytochalasin B, 5-AL, CPZ, genistein): "Each inhibitor reduced the uptake of 50 nm PS-NPs by more than 50%. Inhibition of actin-, caveolin-, and clathrin-mediated endocytosis significantly decreased the internalization of 50 nm PS-NPs by more than 80%" (Results; Fig. 6A,B). This is a genuine loss-of-function test of uptake, though with the usual off-target caveats of these inhibitors and no genetic confirmation.
2. *mTORC1/autophagy as driver of vacuoles and lipid* — rapamycin rescue: "Rapamycin markedly reduced the occurrence of large vacuoles induced by PS-NPs (Fig. 7G), curbing the increase in phosphorylated-S6K1 and perilipin2 (PLIN2) expression (Fig. 7I and S5)" and "rapamycin treatment impeded the accumulation of p62 and LC3 proteins". A real rescue, but pharmacological and pleiotropic.
3. *Autophagic flux* — BafA1 flux assay: "the decrease in the difference in the p62 ratio based on Baf A1 treatment in the presence of PS-NPs indicate that PS-NP treatment impairs autophagic flux (Fig. 7J and 7K)". This is the standard flux design and is the strongest flux evidence in the paper.
4. *Lipid × particle combination in vitro* — oleate co-treatment (Fig. 5H,I), descriptive/pairwise, not an interaction test.

**Measured but not causally manipulated (correlative, though biochemically direct):**
- Cathepsin D maturation: "The protein levels of CtsD were proportionally increased, but predominantly existed in immature forms, pro-CtsD and single chain CtsD, along with a concurrent decrease in the levels of mature form of CtsD (Fig. 4B)" — western blot only.
- Cathepsin B activity: "the activity of the lysosomal enzyme CtsB was also reduced by more than 80% (Fig. 4C)" — a direct enzymatic assay, i.e. functional evidence, but still an *observed* consequence, not a tested cause.
- **Critically: the titular causal step — that reduced cathepsin activity causes the vacuoles / lipid accumulation / MASLD — is NEVER tested.** There is no cathepsin knockdown/knockout, no cathepsin overexpression or re-expression, no cathepsin inhibitor arm (e.g. CA-074, pepstatin) shown to phenocopy PS-NPs, and no rescue by restoring cathepsin activity. The Discussion itself only says such therapies "are under investigation".
- Likewise **"impaired lysosomal acidification" (in the title) is never demonstrated**: no ratiometric or quantitative lysosomal pH measurement. LysoSensor Yellow/Blue DND-160 is used qualitatively and shows the opposite-signed observation — "a live imaging assay using LysoSensor DND-160, which produces yellow fluorescence in an acidic environment, demonstrated the acidic pH of these large vacuoles" (Fig. 6E) — and LysoTracker signal *increases*. The acidification claim rests on analogy to BafA1 and on blocked CtsD maturation, not on measured de-acidification. V-ATPase is invoked throughout the Discussion/Conclusions but **never measured** (no V-ATPase subunit western, qPCR, or activity assay). **This is the single largest title-vs-evidence gap in the paper.** Prior note flagged it; **CONFIRMED.**
- No genetic perturbation of any kind anywhere (no siRNA, shRNA, CRISPR, dominant-negative, knockout mouse).
- **No in vivo mechanistic intervention**: no rapamycin, no autophagy modulator, no cathepsin manipulation in mice. The in vivo and in vitro arms are linked purely by inference.
- **No measurement of hepatic PS-NP burden in mice** — no biodistribution, no particle quantification in liver, no Py-GC/MS, no imaging of particles in tissue. So it is never shown that mouse liver ever saw a particle concentration comparable to 50–100 µg/mL.
- Note also one internal citation slip: the Discussion refers to "Figure 6J-6K" for the autophagic-flux/BafA1 result, which is in Fig. 7J–7K; and "(Fig. 5F and 5G)" and "(Fig. 6C and 6D)" in the Discussion appear to point to Fig. 6F/6G and Fig. 7B–7D respectively.

**Omics-only (no orthogonal confirmation):** the DEG counts and all enrichment themes (PPAR signalling, lipid metabolism, AMPK, fatty acid biosynthesis, cell–cell adhesion, actin organisation); the inflammation gene enrichment (Fig. S2C,D as described); the fibrosis-gene null result (Fig. S2E as described); the 16-gene overlap gene list (Ctse, Il7r, Cd46, Mt1, Ophn1, Myom3 and the further names Coro6, Tnfrsf25, Racgap1, Mecom, Antxr1, Plekhh2 given in Discussion). **Only MYC was carried from the transcriptome to orthogonal in vivo validation** — "MYC … was significantly increased in mouse liver tissue obtained from HFD+PS-NPs group (Fig. S2 and 4D)" (liver qPCR, Fig. 4D) and then in HepaRG by qPCR and western (Fig. 4E,F). Ctse was singled out in text but **never validated by qPCR or protein** in either species.

## 9. Main findings with directions

**In vivo, 20 weeks (Fig. 1, Fig. S1 as described):**
- Body weight ↑, liver weight ↑, liver:body ratio ↑ with HFD, and more so with HFD+PS-NPs (Fig. 1B; percentages quoted above, different reference groups).
- Duration: at 8 weeks "only body weight significantly increased (Fig. S1A)" and "8-week of oral PS-NPs administration was not enough to affect the levels of ALT, AST, and CHO (Fig. S1G)". A clear time dependence.
- **Serum ALT ↑, AST ↑, cholesterol ↑** in HFD+PS-NPs "than those in CD, CD+PS-NPs, and HFD alone (Fig. 1E)"; **serum TG unchanged** — "no notable differences were observed in serum triglyceride (TG) levels between each group (Fig. 1F)". This ALT/AST/CHO panel is the clearest HFD+PS-NPs > HFD contrast in the paper.
- **Steatosis ↑**: CD and CD+PS-NPs "displayed no dramatic changes in fat accumulation"; "20 weeks of HFD exhibited zonal aggregates of mainly microvesicular hepatocytes, whereas HFD+PS-NPs group had diffuse lesions containing a mixture of macrovesicular and microvesicular hepatocytes (Fig. S1B)". Oil Red O quantified in Fig. 1D.
- **Hepatic TG: never measured.** There is no liver triglyceride assay in the mice — steatosis in vivo is histological (Oil Red O area) only, and the only TG numbers are serum (unchanged) and in vitro HepaRG.
- **Inflammation ↑ on BOTH diets**: "oral administration of PS-NPs induced hepatic inflammation in CD-fed mice as well as HFD-fed mice, as evidenced by increased macrophages (Fig. S1C and S1D) and neutrophils (Fig. S1E and S1F)" — F4/80 and MPO immunofluorescence.
- **Fibrosis — internally inconsistent**: histology ↑ ("The exposure to PS-NPs showed increased collagen deposition in CD or HFD fed mice compared to those in the CD-alone and HFD-alone groups (Fig. S1H and S1I)", Picro-Sirius Red) vs transcriptome null ("fibrosis-related genes were not significantly altered (Fig. S2E), suggesting that the concentration and duration of HFD exposure applied in this study were insufficient to fully induce liver fibrosis"). The paper never reconciles these. **CONFIRMED as in prior note.**
- Transcriptome: 361 DEGs (CD+PS-NPs vs CD; 169 up / 192 down), 628 (HFD+PS-NPs vs HFD; 226 up / 402 down), 1501 (HFD+PS-NPs vs CD+PS-NPs; 775 up / 726 down); union 2,180; 16 DEGs shared between the two PS-NP-effect contrasts. Enriched: "PPAR signaling pathway, lipid metabolism (biosynthetic and catabolic processes), cell-cell adhesion, and actin organization" (Fig. 2E,F); Discussion adds "PPAR signaling, AMPK signaling, and fatty acid biosynthesis pathways".
- Myc ↑ in HFD+PS-NPs liver (qPCR, Fig. 4D).

**Against the specific gene list requested:**
- **Lipogenesis — Fasn, Acly, Acaca, Srebf1/SREBP-1c: NOT REPORTED anywhere** (neither in mouse liver nor in HepaRG). Lipogenesis is represented only by PPARG, CD36, PLIN2, PGC-1A, C/EBPA in HepaRG. This is a notable gap for a steatosis paper.
- **PPARα/FAO — Ppara, Cpt1a, Acox1, Cyp4a: NOT REPORTED.** A PPARα antibody is listed in Methods ("Antibodies against … PPARα were obtained from Santa Cruz Biotechnology") but **no PPARα result appears anywhere in the main text or any figure legend.** The "PPAR signaling" enrichment is never resolved into α vs γ vs δ in text; the only PPAR measured is PPARγ, which is ↑. Direction of the PPAR pathway bar in Fig. 2F is conveyed by colour only (pink = up, blue = down per legend) and is not stated in words — **UNVERIFIED from text.**
- **Cholesterol synthesis genes (Hmgcr, Sqle, Srebf2): NOT REPORTED.** Only serum cholesterol ↑.
- **Inflammation — Tnf, Il1b, NLRP3, NF-κB: NOT REPORTED individually.** Only IL-6 ↑ (HepaRG mRNA, Fig. 5L) and the cell-count immunofluorescence (F4/80 ↑, MPO ↑). Inflammatory gene enrichment is omics-only (Fig. S2C,D as described).
- **Autophagy/lysosome: p62 ↑ (mRNA and protein, dose- and time-dependent, Fig. 7B–D), LC3-II/LC3B ↑ (Fig. 7B–D), LAMP1 ↑ (protein, dose-dependent, Fig. 5N and 6G), LAMP2 ↑ and LGP85 ↑ (Fig. 6G), TFEB → nuclear translocation ↑ (Fig. 7F), mTORC1 activity ↑ (p-S6K1 ↑, Fig. 7A), autophagic flux ↓ (Fig. 7J,K), mature cathepsin D ↓ with pro-/single-chain CtsD ↑ (Fig. 4B), cathepsin B activity ↓ >80% (Fig. 4C).** Cathepsin L is **not measured** (mentioned only when citing others: "Patients with MASLD display abnormal autophagy and reduced liver expression of CtsB, CtsD, and CtsL"). Cathepsin E appears only as a mouse DEG, never validated. All of the above are **HepaRG in vitro except the CtsE DEG and Myc qPCR** — there is **no in vivo western or IF for LC3, p62, LAMP1, TFEB or any cathepsin.**
- Other in vitro directions: PPARG ↑ (mRNA, dose-dependent, Fig. 4A, 5L); CD36 ↑, PLIN2 ↑, PGC-1A ↑, C/EBPA ↑ (mRNA and protein); **FABP1 ↓**; C/EBPB unchanged; HO-1 ↑ and IL-6 ↑; HNF4A ↓, ALB ↓, CYP3A4 ↑ (mRNA, protein and activity); intracellular TG "modestly increased" (Fig. 5D); lipid droplets colocalise mostly with ER-Tracker, some with acidic organelles and mitochondria (Fig. 5K).
- Vacuole identity: EEA1-negative, RAB7/LAMP1/LAMP2/LGP85-positive, dextran-filled, LysoSensor-acidic (Fig. 6D–G).
- Cell-type specificity: the phenotype is **HepaRG-only** among HepaRG/HepG2/SK-HEP-1/293T (Fig. 3) — a major caveat for generalising the mechanism.

## 10. Limitations

### (a) Actual limitations of what was done
1. **The core "exacerbate/synergistic" claim is not supported by the statistics used** — one-way ANOVA + Tukey across four groups, no interaction term, no two-way model; the headline 20%-vs-61% comparison uses two different reference groups.
2. **DEGs called at an uncorrected p < 0.05** with |log2FC| > 0.5, no FDR; DEG counts (361/628/1501) are therefore not FDR-controlled and are not comparable with FDR-based counts from other studies. The 16-gene Venn overlap has no null expectation computed.
3. **No DE tool, version, aligner, genome build or normalisation named** in the main text (deferred to Table S4); **RNA-seq n per group is never stated anywhere.** Both make the transcriptome non-evaluable as reported.
4. **No data deposition and no accession** — the RNA-seq cannot be re-analysed or re-thresholded by anyone else; the Data Availability statement is boilerplate that does not actually provide the data.
5. **The HFD vs CD DEG contrast is omitted**, so the diet main effect cannot be sized and no interaction contrast can be reconstructed.
6. **Dose is per animal (0.5 mg/day), not per kg body weight**, while HFD animals end 25–37% heavier — the mg/kg exposure systematically diverges between diet arms over 20 weeks, confounding the very diet × particle comparison being claimed. Consumption of each jelly cube is asserted to be controlled by individual housing but **actual consumption is never verified or reported**.
7. **No particle characterisation beyond the supplier's nominal 50 nm** — no DLS, no measured hydrodynamic diameter, no zeta potential, no PDI, no TEM, no aggregation check in the agarose cube, in gut-relevant media, or in serum-containing culture medium. For a carboxyl-modified anionic particle at 100 µg/mL in 10% FBS this is a material omission (protein corona/agglomeration unknown).
8. **No endotoxin test and no surfactant/preservative control** for the Polysciences suspension; given that inflammation is a reported endpoint on both diets, LPS or surfactant contamination is an unexcluded alternative explanation. The commercial suspension's stabiliser/azide content is never addressed. No particle-free "spent suspension" control.
9. **Control and HFD diets are not ingredient-matched** — grain-based Teklad 2018S chow vs purified D12451 (45 kcal% fat). Diet-composition differences (fibre, phytoestrogens, micronutrients) are confounded with fat content, a well-known problem in HFD studies.
10. **Blinding not reported** for any outcome, including the ImageJ quantification of Oil Red O and Picro-Sirius Red and the histopathological grading; randomisation stated but method unspecified. No exclusion criteria, no animals-lost reporting, no ARRIVE-style reporting.
11. **Males only, single strain (C57BL/6J), exposure starting very young.** "Following a 2-week acclimation period, 4-week-old mice were grouped" is internally ambiguous — whether the mice were 4 weeks old at arrival or at grouping is not reconciled, so the actual age at dosing onset (4 vs 6 weeks) is unclear. No female arm; sex-specific metabolic responses cannot be assessed.
12. **Title overreach on two counts**: (i) "Impaired Lysosomal Acidification" is never measured — no quantitative/ratiometric lysosomal pH, and the LysoSensor data shown indicate the vacuoles *are* acidic with LysoTracker signal increased; V-ATPase is invoked repeatedly but never assayed. (ii) "Reducing Cathepsin Activity" is measured but never causally tested — no cathepsin knockdown, inhibitor-phenocopy, or activity-restoration rescue.
13. **Mechanism rests on a single cell line that the paper itself shows is exceptional** — "The formation of large vacuoles and a substantial increase in lipid accumulation were only observed in HepaRG cells (Fig. 3)". The vacuolar phenotype is absent in HepG2, SK-HEP-1 and 293T, which argues the mechanism may be HepaRG-specific rather than a general hepatocyte response.
14. **In vitro concentrations (50–100 µg/mL) are far above any demonstrated in vivo hepatic concentration**, and hepatic particle burden in the mice was never measured, so the in vivo/in vitro bridge is inference only.
15. **No in vivo mechanistic confirmation**: not one of LC3, p62, LAMP1, TFEB, p-S6K1, cathepsin B activity or cathepsin D maturation was measured in mouse liver. The entire lysosome/autophagy mechanism is in vitro.
16. **Internal contradictions left unresolved**: fibrosis histology (↑) vs fibrosis transcriptome (null); the Abstract's "no effects … on a normal chow diet" vs the Results' CD+PS-NPs increases in macrophages, neutrophils, collagen and 361 DEGs; macropinocytosis described as minor in Results but "clathrin- and caveolin-mediated" endocytosis given as the primary pathway in Discussion; wortmannin purchased but absent from Fig. 6A; several mis-pointed figure citations in the Discussion (Fig. 6J–6K for a Fig. 7 result, etc.).
17. **No replicate counts for any in vitro experiment**, and Fig. 5D's "ANOVA two-tailed unpaired t-test" leaves the test actually applied unclear. SEM used throughout rather than SD or CIs.
18. **Serum sampling timing is stated inconsistently** — Methods: "collected 1 week before euthanization"; Fig. 1E,F legend: "after 19 weeks of treatment". **Fasting state at sampling and at sacrifice is never stated**, nor is time of day — both materially affect ALT/AST/CHO/TG and hepatic lipid/transcriptome readouts.

### (b) Things simply not attempted (absences, not flaws in what was done)
- Female mice; other strains; other ages; other species.
- Other polymers (PE, PP, PVC, PET), other sizes, pristine vs environmentally weathered particles, and unmodified vs carboxyl- vs amine-modified surface chemistry.
- Lower, environmentally realistic doses, or a dose–response in vivo (only one dose tested).
- Hepatic triglyceride/cholesterol biochemical quantification in liver tissue.
- Biodistribution and hepatic particle quantification; gut barrier, microbiome, or portal LPS measurement (despite the Introduction citing LPS findings).
- Proteomics, lipidomics, metabolomics, single-cell or single-nucleus RNA-seq, spatial transcriptomics, ATAC/epigenomics.
- GSEA, WGCNA, or any network/module analysis; any ML model; molecular docking or MD simulation.
- Genetic perturbation of any node (cathepsin B/D/L, TFEB, mTOR, Atg genes), in vitro or in vivo.
- In vivo therapeutic rescue (rapamycin, TFEB activation, cathepsin restoration).
- Primary human hepatocytes, differentiated-HepaRG vs progenitor comparison, organoids, or human liver samples/biopsies; no human exposure or epidemiological anchor.
- Direct lysosomal pH quantification; V-ATPase expression or activity; lysosomal membrane permeabilisation assays (galectin puncta, cathepsin release to cytosol).
- Oxidative-stress measurement in vivo (ROS, MDA, SOD, GSH, Nrf2) — only HO-1 mRNA in HepaRG.
- Mitochondrial function (respirometry, mtDNA copy number, OXPHOS complexes) — only MitoTracker colocalisation imaging and PGC-1A mRNA.
- Insulin resistance/glucose metabolism (no GTT/ITT, no insulin, no HOMA-IR) despite an obesogenic model.
- Circadian/clock genes; time-of-day control.
- Recovery/washout arm; transgenerational or developmental exposure.
- Reporting-guideline compliance (ARRIVE), pre-registration, or a protocol.

---

## Additional recorded details

- **Species/strain/sex/age**: C57BL/6J male mice (DBL, Eumseoung-gun, Chungcheongbuk-do, Republic of Korea); males only; "Following a 2-week acclimation period, 4-week-old mice were grouped" (age at dosing onset ambiguous, see limitation 11).
- **Polymer**: carboxyl-modified polystyrene (PS-COOH), anionic.
- **Nominal size**: 50 nm. Verbatim: "Carboxyl-modified polystyrenes (PS-COOH NPs) with a size of 50 nm (2.5% w/v aqueous suspensions including 3.64 × 10^14 particles/mL) were purchased from Polysciences, Inc. (Warrington, PA, USA)." **Measured size / DLS / zeta / PDI: NOT REPORTED.** Catalogue number: not given.
- **Fluorescent particles**: "Yellow Green-conjugated PS-NPs" (Results) / "GFP-conjugated PS-NPs (PS-NP-GFP)" (Fig. 6C legend) — supplier, size and characterisation of the labelled particles are **not separately specified**, and the two names are used interchangeably for what is presumably one reagent.
- **Dose/route/vehicle/frequency/duration**: 0.5 mg PS-NPs per mouse per day (absolute, not mg/kg; no mg/kg figure given anywhere), oral via voluntary jelly cube — "Oral administration involved the use of PS-NP-jelly cubes prepared as previously described 28, containing either 0 or 0.5 mg PS-NPs with 2.5% agarose"; vehicle = identical 2.5% agarose cube with 0 mg; daily; 8 or 20 weeks concurrent with diet. Not gavage.
- **Diets with catalogue numbers**: control "Teklad Global 18% protein rodent diet; Teklad 2018S; Envigo, Madison, WI, USA"; HFD "45 kcal% fat (D12451; Research Diet, New Brunswick, NJ, USA)".
- **Sampling time / fasting state**: serum "collected 1 week before euthanization" (Methods) / "after 19 weeks of treatment" (Fig. 1E,F legend); liver at 8 or 20 weeks; perfused liver weights taken "under anesthesia with isoflurane". **Fasting state and time of day: NOT REPORTED.**
- **Housing/cage**: individually housed ("mice were individually housed and weighed before and after the study"); no environmental details.
- **Randomisation**: stated ("randomly assigned/divided into four groups"), method not described. **Blinding**: not reported.
- **Cell lines with identifiers**: HepaRG (HPRGC10, Gibco); HepG2 (HB-8065, ATCC); Hep-G2/2.2.15 (SCC249, Sigma-Aldrich); Hep3B (88064, KCLB); Huh7 (60104, KCLB); SK-HEP-1 (HTB-52, ATCC); 293T (CRL-3216, ATCC); PLC/PRF/5 (CRL-8024, ATCC). **No STR authentication and no mycoplasma testing reported.** HepaRG culture: William's E + 1% Glutamax, 10% FBS, 1% Pen-Strep, 5 ng/mL insulin, 50 µM hydrocortisone hemisuccinate; others in DMEM + 10% FBS + 1% Pen-Strep; 37 °C, 5% CO2. **HepaRG differentiation state (progenitor vs DMSO-differentiated) is never stated** — an important omission, since HepaRG phenotype depends heavily on it and may explain the line's uniqueness in Fig. 3.
- **Endotoxin / surfactant checks**: **NONE reported.**
- **Reference gene for qPCR**: RPL13A; 2^−ΔΔCT. **Loading control for western**: β-actin.
- **Software**: GraphPad Prism 9; ImageJ; clusterProfiler 4.0 (R); Gen5 (BioTek); Flow Software 2; instruments Cytation 5, ZEISS LSM 800, FACSCalibur, SpectraMax i3x, VERSAmax, GloMax Navigator, UVITEC MINI H9, Bioruptor, Leica microscope, Olympus BX1, Leica cryostat.
- **DOI link**: https://doi.org/10.7150/ijbs.108268
- **Ethics**: KRIBB IACUC, approval number KRIBB-AEC-21279.

---

## URLs/files actually read

1. `C:\Users\Soham\AppData\Local\Temp\claude\d28fb66f161d523b\5dbfd995-fe63-4ce9-8501-246b4fa39b6f\scratchpad\dl\PMC12210384.txt` — Europe PMC full-text XML converted to text; **lines 1–188 read in full** (the entire body up to the "### References" marker at line 189). This is the primary source for every quotation above.
2. `C:\Users\Soham\Desktop\microplastic\research\literature\papers\ijbs2025_psnp_lysosome.md` — the prior-session notes, read first and then checked claim by claim against (1).

**No web URL was fetched in this session.** The planned fetches of https://doi.org/10.7150/ijbs.108268, https://www.ijbs.com/v21p4017.htm and https://pmc.ncbi.nlm.nih.gov/articles/PMC12210384/ for the supplementary material were **not performed**, so the supplement remains unread (see UNVERIFIED items). Note also that the journal page for this article is at v21p3867 (pages 3867–3885), not v21p4017 as given in the task brief — the v21p4017 URL is likely a different article and should not be used.

---

## UNVERIFIED items

Items I could NOT confirm against the primary source, plus prior-note claims that are wrong or unsupported:

1. **All supplementary material (Fig. S1–S5, Tables S1–S4) — NOT READ.** Everything attributed to S-figures above is only what the main text asserts about them. In particular these remain unverified: the 8-week results (S1A, S1G), steatosis pattern (S1B), macrophage/neutrophil quantification (S1C–F), Picro-Sirius Red quantification (S1H–I), GO terms and the inflammation/fibrosis gene panels (S2A, S2C–E), the IC50 > 1,000 µg/mL viability curve (S3), dextran trafficking (S4), the rapamycin westerns (S5), and **Table S4, which is the only place the RNA-seq pipeline (platform, read length, depth, aligner, quantifier, DE tool and version) is documented.** Tables S1–S3 (reagents, primary and secondary antibodies) likewise unread — so antibody catalogue numbers, dilutions and primer sequences are unverified.
2. **RNA-seq n per group** — genuinely absent from the main text; may or may not be in Table S4. Currently unknown.
3. **DE tool, version, and whether any FDR was computed upstream** — absent from main text; deferred to Table S4.
4. **Direction (up vs down) of the PPAR signalling and other KEGG bars in Fig. 2F** — encoded only in bar colour per the legend; cannot be read from text. Same for the GO terms in Fig. S2A.
5. **Which figure panels used Tukey vs Student's t-test** — Methods give both as alternatives without mapping; Fig. 5D's "ANOVA two-tailed unpaired t-test" is self-contradictory.
6. **Whether the 8-week and 20-week experiments used separate animals**, and the total number of mice used — never stated.
7. **Exact composition of the Fig. 3A cell-line panel** — Results name four lines; Methods culture eight; the legend says "indicated hepatocyte cell lines" and the image cannot be read.
8. **Whether the PS-NP jelly cubes were actually fully consumed**, and the resulting realised dose — asserted to be controlled by individual housing, never measured or reported.
9. **Any repository record** — no accession exists to resolve, so nothing could be fetched or verified; absence of deposition is itself confirmed from the verbatim statement.
10. **Prior-note claim that is WRONG/overstated**: the prior notes list **wortmannin** among the endocytosis inhibitors used. Wortmannin appears only in the Methods list of purchased chemicals; the Fig. 6A legend names 3-MA, cytochalasin B, 5-AL, CPZ and genistein only. Wortmannin is not shown to have been used in any reported experiment.
11. **Prior-note claim NOT CONFIRMABLE**: "The PMC record also carries an empty 'Data Availability Statement' heading with no further text." The converted text contains only the "Availability of data and materials" section; no separate empty heading is present in the file read. Immaterial to conclusions.
12. **Prior-note claim not confirmed as written**: the prior notes describe Fig. 5 lipid accumulation as measured by "BODIPY, Nile Red, flow cytometry; Fig. 5A-5C" — correct — but also state in vitro replicate numbers are "not reported in the legends read"; I confirm they are reported **nowhere**, including Methods.
13. **Prior note's in vitro dose list for crystal violet** is 0/50/100 µg/mL (Methods), narrower than the 0/10/20/50/100 range used for BODIPY/Nile Red; the prior note gave a single combined range. Minor correction.
14. Everything else load-bearing in the prior notes — the question, the 2×2 design, n=6, randomisation, individual housing, no blinding, particle identity and the absence of DLS/zeta/endotoxin, dose/route/vehicle/duration, diets and catalogue numbers, the verbatim statistics and DEG-threshold text, the absence of any interaction term, the DEG counts and splits, the serum/histology/fibrosis findings, the cathepsin and autophagy results, the HepaRG-only phenotype, the data availability statement and absence of accessions, and the title-vs-acidification gap — **was CONFIRMED verbatim against the full text.**
# 6. EHP 2023 — Okamura et al., "Oral Exposure to Polystyrene Microplastics of Mice on a Normal or High-Fat Diet and Intestinal and Metabolic Outcomes"

**Citation.** Okamura T, Hamaguchi M, Hasegawa Y, Hashimoto Y, Majima S, Senmaru T, Ushigome E, Nakanishi N, Asano M, Yamazaki M, Sasano R, Nakanishi Y, Seno H, Takano H, Fukui M. *Environ Health Perspect* 2023;131(2):027006. DOI https://doi.org/10.1289/EHP11072. PMID 36821708, PMCID PMC9945580. Received 2022-02-08; revised 2023-01-05; accepted 2023-01-13.

**Access level: FULL TEXT read** (complete Europe PMC JATS body, i.e. everything from Abstract through the end of the Discussion and funding notes, lines 1–278 of the converted text, before `### References`). **Supplemental Material NOT read** — only two PDFs exist (`ehp11072.s001.acco.pdf`, `ehp11072.smcontents.508.pdf`) and they were not downloaded; Table S1 (primer list), Table S2 (all summary data), Figures S1–S4 are therefore unverified. Note the paper's exact title uses "Intestinal and Metabolic Outcomes", not "Intestinal/hepatic effects" as in the task brief.

---

## 1. Exact scientific question

Stated hypothesis (Introduction, final paragraph of the MP/LGS section): *"we hypothesized that the toxicity of MPs might be more marked in LGS [leaky gut syndrome] induced by a high-fat diet than under a normal diet."*

Stated objective (Abstract, Objectives; repeated at the end of the Introduction): *"We aimed to evaluate the pathophysiology of intestinal outcomes associated with consuming a high-fat diet and simultaneous intake of MPs, focusing on endocrine and metabolic systems."*

So the question is explicitly an **effect-modification question** (does HFD-induced barrier failure amplify MP toxicity?), with the small intestine (jejunum) as the primary target organ and liver/glucose metabolism as downstream readouts. The liver is secondary, justified in the Introduction via portal LPS influx and NAFLD/NASH.

## 2. Experimental design

- **2 × 2 factorial**, diet (ND vs HFD) × MP (absent vs present). One particle size only, one dose only. Methods "Mice": *"The mice were divided into four groups (n=10 per group) according to diet and whether or not they were given MPs."* Groups: ND, ND+MP, HFD, HFD+MP.
- **Timeline** (Figure 1A): purchased at 7 wk; diet + MP from 8 wk to 12 wk of age = **4 wk exposure**; euthanasia at 12 wk after an overnight fast.
- **Phenotype n = 10/group** for: body weight, food/water intake (Fig. 1B–D), FITC-dextran permeability (Fig. 1I; Methods), serum ALT/TG/NEFA/LDL-C (Methods, Blood Biochemistry), serum/fecal/hepatic palmitic acid (Methods, GC-MS), jejunal histology + Muc2 IF + GFP+ area (Fig. 2B–G), lamina propria flow cytometry (Fig. 2H–K), fecal SCFAs (Fig. 3), liver weight/NAS/Oil Red O/hepatic palmitate (Fig. 4).
- **n = 5/group** for iPGTT and for ITT, and Methods state *"iPGTT and ITT were performed in different mice"* — so the two functional tests used disjoint subsets of the 10 mice.
- **"Omics" n = 6/group**: 16S rRNA sequencing (Figure 5B–G legend, "n=6" for OTUs, Chao1, Shannon, Gini-Simpson and both LEfSe comparisons). How the 6 of 10 were selected is not stated.
- **qPCR n is internally inconsistent**: Methods (Gene Expression Analysis) say *"Expression levels from six mice from each group were determined, and RT-PCR was performed in triplicate for each sample (n=10)"*; the Table 1 footnote says "(n=10)" for both jejunum and liver panels. 6 vs 10 unresolved.
- **In vitro arm**: MODE-K murine intestinal epithelial cells (sex: female; gift from R. Blumberg), 1×10⁴ cells/well in 96-well plates, treated on day 5 for 24 h with MP 1,000 µg/L, palmitic acid 200 µM, recombinant mouse IL-22 50 ng/mL (NBP2-35122, Funakoshi). Methods say n = 6; **Figure 6 legend says n = 10** — second internal inconsistency.
- **Housing**: *"One mouse was kept per cage."* Singly housed, so cage = animal (no pseudo-replication from cage effects, but also no cage-level microbiome coprophagy control across groups).
- **Pair feeding**: *"equal amounts of feed and water were supplied for pair feeding"*; intake measured every 3 d and *"normalized based on the lowest intake (group fed a normal diet with MPs added)"*. All groups were therefore mildly restricted to the ND+MP group's intake.
- **Randomisation: NOT REPORTED** (no occurrence of "random" anywhere in the full text). **Blinding: NOT REPORTED** (no occurrence of "blind"), which matters for NAS scoring, villus morphometry and flow gating.
- Ethics: Committee for Animal Research, Kyoto Prefectural University of Medicine, approval M2021-61.

## 3. Datasets generated or re-used; accessions; deposition

- Data generated: 16S rRNA V3–V4 amplicon sequencing (Illumina MiSeq, Macrogen) on n = 6/group; targeted GC-MS (palmitic acid in serum/feces/liver; SCFAs in feces); TaqMan RT-qPCR panels; histology/IF images; flow cytometry; serum chemistry. **No RNA-seq, microarray, proteomics or untargeted metabolomics anywhere in this paper.**
- Data re-used: none; no public dataset is re-analysed.
- **Data Availability statement: THERE IS NONE.** The article has no Data Availability / Availability of Data section. Verified three ways: (i) grep of the Europe PMC JATS XML for "data availab", "accession", "PRJNA", "PRJDB", "DRA", "GSE", "SRA", "GitHub", "Zenodo", "figshare" → zero hits; (ii) the XML `<sec sec-type="associated-data">` block contains only *"Associated Data — Supplementary Materials — Click here for additional data file. Click here for additional data file."*; (iii) independent fetch of https://pmc.ncbi.nlm.nih.gov/articles/PMC9945580/ returned no data-availability sentence and no accession.
- The **only** verbatim data-sharing-type sentence in the whole article is the front-matter line: *"Supplemental Material is available online (https://doi.org/10.1289/EHP11072)."* Figure legends repeatedly add: *"Summary data can be found in Table S2."*
- **Accessions: NONE.** The 16S raw reads are not deposited, or at least no deposit is declared anywhere in the article or on the PMC record. There is therefore **no accession to resolve**. (Caveat: one string resembling an accession, "MA516870", appears in the flow-cytometry antibody list — it is a Thermo/eBioscience catalogue number for FITC-CD206, not a data accession.)
- Supplemental files (names from the XML `<media xlink:href>` attributes): `ehp11072.s001.acco.pdf` (Tables S1–S2, Figures S1–S4) and `ehp11072.smcontents.508.pdf` (contents listing). Both are PDFs — summary data only, **no machine-readable per-animal data file**.

## 4. Statistical / computational models

**Core statistics (Methods, "Statistical Analysis", verbatim):** *"The data were analyzed using JMP software (version 13.0; SAS Institute, Inc.). One-way ANOVA with Holm–Šídák's multiple comparisons test was used to compare the results of different groups. Statistical significance was set at p<0.05."* Figures drawn in GraphPad Prism 9.0.

**The critical point for the team: the 2 × 2 design was analysed as a ONE-WAY ANOVA across four independent groups followed by pairwise Holm–Šídák contrasts. No two-way ANOVA was fitted. No diet × particle interaction term was fitted. No interaction F statistic or interaction p-value is reported anywhere** — not in Methods, Results, Table 1, Table 2, or any figure legend. Grep of the full text for "two-way" returns zero hits; the only two hits for "interaction" are (a) "biological interactions" in the Introduction and (b) nothing else of statistical relevance. Every figure legend without exception reads *"Data were analyzed using one-way ANOVA with Holm-Šídák's multiple comparisons test."*

- Pairwise contrast coding (Table 1 footnote, verbatim): *"The statistical analyses among four groups were performed by one-way ANOVA with Holm-Šídák's multiple comparisons test: *, p<0.05 vs. ND; †, p<0.05 vs. ND+MP; ‡, p<0.05 vs. HFD group."* Table 2 and the figures use \*p<0.05, \*\*p<0.01, \*\*\*p<0.001, \*\*\*\*p<0.0001.
- **Multiple testing**: only the within-ANOVA Holm–Šídák step-down across the 4-group comparisons. No family-wise or FDR correction across the >50 outcomes spanning histology, flow cytometry, qPCR (15 genes), 8 phyla, metabolites and functional tests.
- **16S pipeline** (Methods, "16s rRNA Sequencing"): QIIME v1.9.1; barcodes/primers with scores <75% excluded; OTU clustering by **UCLUST at 97% similarity**; taxonomy by **BLAST against UNITE 2017** with "the UNITE sequence set of the Greengenes core set aligned with UCLUST and ITS" (as written — UNITE is a fungal ITS database, so this sentence is internally incoherent for a bacterial 16S study; flagged below). Alpha diversity: Chao1, Shannon, Gini–Simpson. No beta-diversity ordination (no PCoA/UniFrac/PERMANOVA) is reported. No rarefaction depth, no read counts, no DADA2/ASVs.
- **Differential abundance**: phylum level by the same one-way ANOVA + Holm–Šídák (Table 2). Genus/clade level by **LEfSe** (http://huttenhower.sph.harvard.edu/lefse/), *"A p-value threshold of 0.05 (Wilcoxon rank sum test) and an effect size threshold of 2 were used for all biomarkers"* (LDA log10 score > 2). **LEfSe was run as two separate two-class comparisons (ND vs ND+MP; HFD vs HFD+MP), never as a diet × MP design.** No FDR within LEfSe is stated.
- **No DE tool in the transcriptomic sense** (no DESeq2/edgeR/limma) because there is no transcriptomics. qPCR quantified by 2^−ΔΔCt normalised to *Gapdh*, ND set to 1.0.
- **Enrichment analysis: none** (no GO/KEGG/GSEA). **Correlation analyses: none** (no microbiome–phenotype, microbiome–metabolite or gene–phenotype correlations are computed; all cross-domain links are verbal inference in the Discussion). **No ML of any kind** beyond LEfSe's LDA, which is used as an effect-size filter, not as a classifier with held-out evaluation.
- No power calculation, no pre-registration, no exclusion criteria, no mixed/hierarchical model, no covariate adjustment, no batch variable.

## 5. Does the model suit the question? (difference-of-significance problem)

**No — this is the paper's central methodological weakness, and it is exactly the failure mode the team is auditing.**

The scientific question is an interaction question ("is MP toxicity greater under HFD than under ND?"). The statistical model is a one-way ANOVA over four nominal groups. An interaction is never estimated or tested. Every diet-dependence claim in the paper is built from the pattern "significant in HFD+MP vs HFD, not significant in ND+MP vs ND", i.e. the **difference-of-significance fallacy** — a difference in significance is not evidence of a significant difference.

The claims that rest on this fallacy include the headline sentences:
- Abstract/Discussion: *"Our findings suggest that in this study, MP induced metabolic disturbances, such as diabetes and NAFLD, only in mice fed a high-fat diet."*
- Conclusion: *"the present study is the first to suggest that MPs cause metabolic disturbances under HFD intake conditions"*.
- Results, liver qPCR: *"On the other hand, their expressions were not different between ND mice and mice fed the ND with MPs."*

Three further points sharpen this:
1. **Two readouts actually do show MP effects under ND**, contradicting a strict "only under HFD" reading: iPGTT glucose and AUC (*"blood glucose levels were higher in mice fed a ND with MPs added and mice fed HFD with MPs added than in both ND and HFD mice"*, Fig. 1E,F) and goblet cell number / Muc2 mucus thickness (*"The number of goblet cells in mice fed the ND containing MPs and the those fed HFD containing MPs was significantly lower than in ND and HFD mice, respectively"*, Fig. 2E,F). LEfSe also found taxa shifted under ND (Fig. 5F,H).
2. The authors themselves repeatedly describe **non-significant ND trends** — serum lipids *"although the difference was not statistically significant in the present sample size"* and *"Red O staining area of the liver was higher in the group fed the ND with MPs than in the ND-only group, although not significantly, suggesting that the ND group was also affected by MP-induced metabolic disturbances to a small extent"* (Discussion). This is an explicit power-limitation admission that undercuts "only in HFD".
3. Because the ND control is **chow** and the HFD is a **purified lard diet**, even a properly fitted interaction would confound "fat" with "diet matrix" (see §10).

A reader could in principle recompute a two-way ANOVA from group means and SDs with n per group if Table S2 reports them — but Table S2 is a PDF and was not read here, and the qPCR n (6 vs 10) is ambiguous, so even a reconstructed interaction test would be approximate.

## 6. Inputs, features, outputs, targets

- **Inputs (experimental)**: diet (ND chow vs D12492 60% kcal lard HFD) and PS-MP in drinking water (1,000 µg/L vs none), 4 wk.
- **Features / measured variables**: body weight; food and water intake; iPGTT and ITT glucose curves and AUCs; plasma FITC-dextran minus baseline GFP signal (permeability surrogate); serum ALT, TG, NEFA, LDL-cholesterol, palmitic acid; jejunal villus height and width, crypt depth, PAS+ goblet cells per area, Muc2 mucus-layer thickness, GFP+ area fraction (MP deposition); lamina propria ILC1/CD45+, T-bet+ILC3/CD45+, ILC3/CD45+, M1:M2 ratio; fecal palmitic acid, acetic, propanoic, butanoic acid; liver wet weight, liver:body ratio, NAS, Oil Red O area, hepatic palmitic acid; 16S OTU counts, phylum relative abundances, Chao1/Shannon/Gini–Simpson, LEfSe taxa; qPCR ΔΔCt for 9 jejunal and 6 hepatic genes; MODE-K *Muc2* expression and intracellular GFP fluorescence.
- **Outputs/targets of inference**: the metabolic phenotype (glucose intolerance, NAFLD severity) and the intestinal barrier/immune phenotype, as a function of MP exposure *conditional on diet*. There is no predictive model, no train/test split, no held-out set, no classifier — "targets" exist only in the statistical-hypothesis sense.

## 7. How results were evaluated

Purely by null-hypothesis significance testing: one-way ANOVA + Holm–Šídák pairwise contrasts, α = 0.05, reported as star annotations in figures and dagger symbols in Table 1. Effect sizes are reported only implicitly as fold-change (qPCR, ND = 1.0) and percentage abundance (Table 2); LDA score > 2 is the only explicit effect-size threshold. **No confidence intervals, no estimation plots, no cross-validation, no replication cohort, no sensitivity analysis, no reporting of exact p-values** for most endpoints (Table 2 is the exception — it gives per-phylum p-values). Data presented as mean ± SD throughout.

## 8. Independent experimental validation vs omics-only

This paper is overwhelmingly **assay-based, not omics-based**; the only "omics" is 16S amplicon sequencing.

**Independently measured / multi-modal (not omics-dependent):**
- Glucose phenotype: iPGTT + ITT, physiological tests (n = 5).
- Barrier: FITC-dextran permeability (functional) **plus** PAS goblet-cell counts **plus** Muc2 immunofluorescence **plus** jejunal *Muc2* qPCR — four independent modalities agreeing on mucus-barrier loss.
- MP deposition: GFP+ tissue area (imaging) plus plasma GFP signal (Fig. S3D,E).
- Liver: NAS (H&E histology) **plus** Oil Red O area **plus** hepatic palmitic acid by GC-MS **plus** hepatic lipogenic qPCR — four modalities agreeing on increased hepatic lipid.
- Inflammation: lamina propria flow cytometry **plus** jejunal and hepatic cytokine qPCR.
- SCFA/palmitate: targeted GC-MS, an orthogonal chemistry readout to the 16S data.

**Omics-only / unvalidated:** every microbiome claim. The 16S results (Bacteroidetes down, Proteobacteria up, *Desulfovibrio* up, diversity down) have **no gnotobiotic, fecal-transplant, antibiotic, qPCR-of-taxa or culture validation**. The mechanistic chain "dysbiosis → less SCFA → less IL-22 → thinner mucus → MP uptake → inflammation" is argued, not demonstrated: there is **no in vivo IL-22 or SCFA intervention**.

**Mechanistic in vitro arm (MODE-K)** partially substitutes: palmitic acid (± MP) lowered *Muc2*, MP alone did not (*"The expression of Muc2 was not different between control group and MP alone group"*, Fig. 6A); IL-22 restored *Muc2* (Fig. 6B); MP accumulation in cells was higher with PA and reduced by IL-22 (Fig. 6C,D). This is a cell line only, female-derived, and does not test the mouse mechanism.

## 9. Main findings with directions

All contrasts below are **pairwise** (one-way ANOVA + Holm–Šídák), never interaction tests.

**Body weight**: **no MP effect in either diet.** *"We observed no difference in the body weight of mice fed a ND and a ND with MPs added and of mice fed a HFD and a HFD with MPs added"* (Results; Fig. 1B). Food and water intake also not different among the four groups (Fig. 1C,D) — but intakes were pair-fed/normalised by design, so this is partly constructed. **Liver**: absolute wet liver weight was **LOWER** in HFD+MP than HFD (Fig. 4A); liver:body ratio **not significantly different** across all four groups (Fig. 4B). The lower absolute liver weight alongside higher steatosis is unexplained by the authors.

**Hepatic lipid**: NAS **UP** in HFD+MP vs HFD; *"The NAS was zero in ND-fed mice, with and without exposure to MPs, whereas NAS in mice fed the HFD containing MPs was higher than that in mice fed the HFD without MPs"* (Fig. 4D). Oil Red O area **UP** in HFD+MP vs HFD, unchanged ND vs ND+MP (Fig. 4E; with a non-significant ND trend admitted in the Discussion). Hepatic palmitic acid **UP** in HFD+MP vs HFD, unchanged under ND (Fig. 4F, GC-MS).

**Lipogenesis genes (liver, Table 1; mean (SD), ND = 1)** — note the paper measured **Scd1, Elovl6 and Fasn only**:
| Gene | ND | ND+MP | HFD | HFD+MP |
|---|---|---|---|---|
| Scd1 | 1 (0.12) | 1.16 (0.09) | 1.22 (0.20)\*† | 3.03 (0.27)\*†‡ |
| Elovl6 | 1 (0.12) | 1.16 (0.09) | 1.44 (0.24)\*† | 4.28 (0.95)\*†‡ |
| Fasn | 1 (0.12) | 1.16 (0.09) | 2.05 (0.09)\*† | 9.82 (0.90)\*†‡ |
Direction: **UP with MP on HFD (‡ = p<0.05 vs HFD), unchanged with MP on ND.** **Acly, Acaca and Srebf1 were NOT measured.** *Pparg* not measured.
**Data-integrity oddity (confirmed against Table 1):** *Scd1*, *Elovl6* and *Fasn* have **byte-identical ND and ND+MP means and SDs** — 1 (0.12) and 1.16 (0.09) — across three different genes; likewise liver *Il6* and *Il1b* both read 1 (0.19) and 1.01 (0.24). This could be a transcription/copy error in the table and should be checked against Table S2.

**PPARα / fatty-acid oxidation (Ppara, Cpt1a, Acox1, Cyp4a): NOT MEASURED.** There is no FAO readout of any kind in this paper. Equally **cholesterol synthesis was not measured** in liver (no *Hmgcr*, *Sqle*, *Srebf2*); the only cholesterol-related datum is serum LDL-cholesterol, reported as higher in HFD+MP vs HFD in Figure S4D (figure not read — value unverified).

**Inflammation (liver, Table 1)**: *Tnfa* 1 (0.21) → 1.05 (0.28) → 2.44 (0.83)\*† → **5.93 (1.18)\*†‡**; *Il6* 1 (0.19) → 1.01 (0.24) → 2.22 (0.36)\*† → **5.51 (0.50)\*†‡**; *Il1b* 1 (0.19) → 1.01 (0.24) → 1.66 (0.23)\*† → **6.17 (1.72)\*†‡**. All **UP** with MP on HFD, unchanged on ND. No NF-κB, NLRP3, F4/80 or CD68 liver readouts; no liver protein/western confirmation.

**Gut barrier**: permeability **UP** — *"plasma GFP expression was significantly higher in the mice fed an HFD that included MPs than in HFD-fed mice"* (Results; Fig. 1I, Fig. S3D,E); no ND vs ND+MP difference. Villus height and width **DOWN**, crypt depth **UP** in HFD+MP vs HFD only (Fig. 2B–D). Goblet cells and Muc2 mucus thickness **DOWN with MP in BOTH diets** (Fig. 2E,F). GFP+ deposition area **UP only** in HFD+MP vs HFD (Fig. 2G). Jejunal qPCR (Table 1): *Tnfa* 1.66→2.53, *Il6* 2.58→3.93, *Il1b* 2.65→4.04, *Ffar2* 2.43→6.26, *Ffar3* 2.45→3.28, *Cd36* 1.47→7.63, *Sglt1* 2.31→6.11 all **UP**; *Il22* 0.47→0.27 and *Muc2* 0.48→0.18 **DOWN** (HFD→HFD+MP); all ND vs ND+MP contrasts null.

**Innate immunity (Fig. 2H–K)**: ILC1/CD45+, T-bet+ILC3/CD45+ and M1:M2 ratio **UP**, ILC3/CD45+ **DOWN** in HFD+MP vs HFD; no ND vs ND+MP differences.

**Metabolites (Fig. 3)**: fecal palmitic acid **DOWN** and acetic, propanoic and butanoic acid all **DOWN** in HFD+MP vs HFD; no ND changes. Serum palmitate/ALT/TG/NEFA/LDL-C **UP** in HFD+MP vs HFD only (Fig. S4A–E, not read).

**Microbiome (Table 2, phylum means ± SD with per-phylum ANOVA p)**: HFD+MP vs HFD — Bacteroidetes **DOWN** 49.77→29.12% (p=0.026), Proteobacteria **UP** 14.66→44.98% (p=0.001), Deferribacteres **UP** 4.90→8.65% (p=0.035); Firmicutes NS (30.30→19.90%, p=0.220). ND vs ND+MP: nothing significant at phylum level (Bacteroidetes p=0.859, Proteobacteria p=0.995, Firmicutes p=0.996). Alpha diversity (Chao1, Shannon, Gini–Simpson) **DOWN** in HFD+MP vs HFD, unchanged under ND; OTU number unchanged in both comparisons (Fig. 5B–E). LEfSe: ND+MP enriched for Proteobacteria, *Parasutterella*, Sutterellaceae, Betaproteobacteria and depleted for Bacteroidales/Bacteroidia/Bacteroidetes/*Tannerella* (6 over-, 5 under-represented taxa); HFD+MP enriched for *Desulfovibrio* and depleted for Prevotellaceae, Bacteroidales, Bacteroidia, Bacteroidetes, *Alloprevotella*, Ruminococcaceae (7 under-, 1 over-represented).

**Oxidative stress, mitochondrial function, fibrosis staining, circadian genes: NOT MEASURED.** Fibrosis appears only as a 0–4 sub-score inside NAS; no Sirius red/Masson stain, no *Col1a1*/*Acta2*.

## 10. Limitations

### (a) Actual limitations of what was done

1. **Design–analysis mismatch (most serious).** A 2 × 2 factorial analysed by one-way ANOVA with pairwise contrasts; the central "MP only harms on HFD" claim is never tested as an interaction. See §5.
2. **No multiplicity control across outcome families.** >50 endpoints, only within-ANOVA Holm–Šídák.
3. **Unmatched control diet.** ND is a grain-based CLEA chow (345 kcal/100 g, 4.6% kcal fat); HFD is purified D12492 (459 kcal/100 g, 60% kcal lard). These differ in fibre, sucrose, phytoestrogens and micronutrients, not just fat — the correct control would be D12450B/J. Any "diet-dependent MP effect" is therefore confounded with diet matrix, and the microbiome contrast especially so.
4. **Internal inconsistencies in reported n**: qPCR 6 mice (Methods) vs n=10 (Table 1 footnote); MODE-K n=6 (Methods) vs n=10 (Fig. 6 legend). If technical triplicates were counted as n, the qPCR p-values are inflated.
5. **ITT dose discrepancy**: 0.5 U/kg (Methods, Analytical Procedures) vs **0.75 U/kg** (Figure 1 legend).
6. **Permeability assay is acknowledged to be imperfect**: *"It was difficult to separate the FITC signal from FITC-dextran and GFP signal from GFP-MP using absorbance spectrometry"* — the index is a crude after-minus-before subtraction, and the same GFP signal doubles as the MP-deposition measure, so dye leaching from beads is an unexcluded artefact (the authors note only ~50% fluorescence retention at 3 d, ref. 58). Timing is also inconsistently reported: blood 3 h after gavage in Methods and Results vs *"plasma dextran 4 h after oral administration"* in the Discussion.
7. **No randomisation or blinding reported**, despite subjective endpoints (NAS, villus morphometry, flow gating).
8. **Small and unexplained omics subsample**: 16S on 6 of 10 mice per group with no stated selection rule.
9. **16S pipeline is dated and partly incoherent as written**: QIIME 1.9.1, 97% UCLUST OTUs (not ASVs), taxonomy via *"BLAST (UNITE 2017) … with the UNITE sequence set of the Greengenes core set aligned with UCLUST and ITS"* — UNITE/ITS is a fungal reference, so the bacterial taxonomy assignment is not reproducibly described. No read depths, no rarefaction, no beta diversity/PERMANOVA.
10. **Inconsistent sampling site for "feces"**: SCFAs are introduced as small-intestinal contents (*"Samples obtained during euthanasia from the small intestines"*) but then measured in *"rectal fecal samples (20 mg)"*; 16S DNA came from *"frozen appendicular fecal samples"*. Three different gut compartments are pooled under the word "feces", which weakens the SCFA↔microbiome link.
11. **No raw data and no sequence deposit** (§3) — the 16S reads cannot be reanalysed, and no per-animal data file exists (Table S2 is a PDF).
12. **Table 1 numeric anomaly**: identical ND/ND+MP means and SDs across three lipogenic genes and across two cytokines (§9).
13. **Single dose, single particle size, single sex, single time point**; pair feeding to the lowest-intake group means all animals were mildly restricted, and the actual ingested particle dose (µg/kg/d or particles/d) is never derived despite water intake being measured.
14. **No particle characterisation beyond dispersion stability**: no DLS/SEM size distribution, no zeta potential, no agglomeration assessment, **no endotoxin assay and no check for surfactant/azide preservative** (grep confirms "endotoxin" appears only as background in the Introduction; "surfactant", "azide", "zeta" appear nowhere). COOH-modified sub-micron PS suspensions are commonly supplied with such additives, and the beads were not reported as washed.
15. **Author-acknowledged** (Discussion, penultimate paragraph, verbatim): *"There was a limitation of this study. In iPGTT, we did not measure the plasma insulin levels. If we had had the data, more accurate assessment of insulin resistance would have been possible. In addition, the evaluation of dose response and verification with other size and concentration combinations of microplastics has not been conducted, and this is a future research topic."* Plus the implicit power admission: ND+MP serum lipids were higher *"although the difference was not statistically significant in the present sample size"*.
16. **Time of sacrifice / circadian phase not stated**; all mice fasted 16 h before liver sampling, which itself suppresses *Fasn*/*Scd1*/*Elovl6* — fold-changes are therefore measured in a fasted, lipogenesis-suppressed state.
17. **Particle size is sub-micron (0.45–0.53 µm)** yet the paper is framed throughout as "microplastics"; this is at the nanoplastic/MP boundary and is not the size range of most dietary MP exposure.

### (b) Things simply not attempted (scope, not flaws)

- No liver (or intestinal) transcriptomics, proteomics or untargeted metabolomics — the molecular liver readout is a 6-gene qPCR panel.
- No PPARα/FAO axis, no cholesterol-synthesis genes, no SREBP/ACLY/ACACA, no *Pparg*.
- No oxidative-stress, mitochondrial-function or ER-stress measurements.
- No dedicated fibrosis histology or fibrogenic gene expression.
- No protein-level confirmation (western/ELISA) for any qPCR finding.
- No germ-free/antibiotic/FMT causal test of the microbiome arm; no in vivo IL-22 or SCFA supplementation.
- No female mice, no second strain, no recovery/washout arm, no dose–response, no alternative polymer (PE/PVC/PET) or size comparison.
- No MP quantification in liver tissue (only intestinal GFP+ area), so hepatic effects are not tied to hepatic particle burden.
- No insulin, GLP-1, bile acid or LPS/endotoxaemia measurements in serum.
- No correlation or mediation analysis linking microbiome → SCFA → barrier → liver.

## Reported methods details (checklist)

| Item | Value (source) |
|---|---|
| Species/strain | Mouse, C57BL/6J wild type, littermates (Methods, Mice) |
| Sex | **Male only** |
| Age | Purchased at 7 wk; exposure 8→12 wk; euthanised at 12 wk |
| Supplier | Shimizu Laboratory Supplies |
| Polymer | Polystyrene, carboxyl (COOH)-modified, green-fluorescent |
| Product code | F-K1 050, **Merck** |
| Size | **0.45–0.53 µm** (nominal; no measured distribution) |
| Dose | **1,000 µg/L in drinking water**; ad libitum; no mg/kg/d derived |
| Route | Oral, via drinking water |
| Vehicle | Drinking water, sonicated 20 kHz 15 min; jugs foil-wrapped; water changed every 3 d |
| Duration | 4 wk (8→12 wk of age) |
| ND | 345 kcal/100 g, 4.6% kcal fat, **CLEA** (no catalog number given) |
| HFD | 459 kcal/100 g, 20% protein, 20% carbohydrate, 60% fat (lard), **D12492, Research Diets Inc.** |
| Feeding | Pair-fed; food and water normalised to lowest-intake group (ND+MP), measured every 3 d |
| Sampling / fasting | Euthanasia at 12 wk **after an overnight fast**; jejunum and liver from mice **fasted 16 h**; iPGTT after 16-h fast, 2 d before sacrifice; ITT after 5-h fast, 3 d before sacrifice; FITC-dextran after 4-h fast, 5 d before sacrifice; SI pathology evaluated 2 d after sacrifice |
| Anaesthesia/euthanasia | 4.0 mg/kg midazolam + 0.3 mg/kg medetomidine + 5.0 mg/kg butorphanol |
| Housing | **One mouse per cage**; pathogen-free controlled environment |
| Randomisation | **Not reported** |
| Blinding | **Not reported** |
| Endotoxin / surfactant checks | **Not reported / not performed** |
| Particle QC performed | Only dispersion stability in water (count + fluorescence intensity on days 0/1/2, ×3; Figure S3) |
| Ethics | Kyoto Prefectural University of Medicine, approval M2021-61 |
| Funding | JST-CREST JPMJCR19H3 (extensive unrelated pharma disclosures) |
| DOI | https://doi.org/10.1289/EHP11072 |

## Prior-note reconciliation (every load-bearing claim re-checked)

**Confirmed verbatim against the full text:** 2×2 design and n=10/group; single-housing; pair feeding and normalisation to ND+MP; littermates; 7-wk purchase and 8→12 wk exposure; the F-K1 050 / 0.45–0.53 µm / COOH / Merck particle description; 1,000 µg/L in drinking water; CLEA ND and D12492 HFD specs; iPGTT/ITT n=5 in different mice; FITC-dextran n=10; 16S n=6; the "six mice … (n=10)" qPCR ambiguity; the MODE-K n=6 vs n=10 conflict; the ITT 0.5 vs 0.75 U/kg conflict; the exact Statistical Analysis paragraph; **no two-way ANOVA / no interaction term / no interaction p anywhere**; QIIME 1.9.1, UCLUST 97%, UNITE 2017/Greengenes wording, LEfSe Wilcoxon p<0.05 + LDA>2 run within each diet; absence of randomisation and blinding; absence of endotoxin testing; all Table 1 values and symbols; all Table 2 values and p-values; the body-weight, liver-weight, NAS, Oil Red O, hepatic palmitate, SCFA, flow-cytometry and LEfSe directions; the identical ND/ND+MP means for *Scd1*/*Elovl6*/*Fasn*; the iPGTT and goblet-cell/Muc2 exceptions under ND; both author-acknowledged limitations verbatim; the permeability-assay caveat verbatim; the "Red O … although not significantly" quote verbatim; absence of any Data Availability statement or accession.

**Prior-note errors / imprecisions found:**
1. **Title.** The prior note's citation title is correct ("…Intestinal and Metabolic Outcomes"); the *task brief*'s "Intestinal/hepatic effects" is not the published title. (Not a prior-note error, but worth flagging.)
2. **Jejunal gene list.** The prior note writes "Sglt1/'Sglt2' in one sentence". Precisely: the Methods gene list says ***Sglt2***, while the Results text and Table 1 both say ***Sglt1***. This is a genuine inconsistency in the paper, and Table S1 (primer list, not read) would settle which assay was run.
3. **"Serum LDL-cholesterol higher in HFD+MP vs HFD (Fig. S4D)".** The main text says only that serum ALT, TG, NEFA, LDL-C and palmitate *"in mice fed the HFD containing MPs were significantly higher in HFD-fed mice"* — the per-analyte panel letters and values are in Figure S4, which was **not read**. Treat the per-analyte attribution as unverified.
4. **"Fecal palmitic acid lower in HFD+MP vs HFD".** Correct, but the prior note omits the stronger Discussion claim: *"significantly lower in mice fed the HFD with MPs in comparison with that in the other three groups"*, and omits the Results observation that HFD+MP fecal palmitate returned to ND-like levels.
5. **SCFA/16S sampling site.** The prior note says "small-intestinal feces (palmitic acid, SCFA)". The SCFA Methods actually say the samples came from the small intestine *and then* describe homogenising *"rectal fecal samples"*; the 16S used *"frozen appendicular fecal samples"*. Three compartments, inconsistently described — a point the prior note missed.
6. **Permeability timing.** The prior note does not flag that the Methods/Results say blood was taken **3 h** after gavage whereas the Discussion says **4 h**.
7. **Microbiome LEfSe direction sentence in the Discussion.** The Discussion contains an apparent error of its own: *"The phylum Proteobacteria, the genus Parasutterella, the family Sutterellaceae, and the class Betaproteobacteria were more frequently present in mice fed the HFD with MPs than in **ND-fed mice**"* — this repeats the ND+MP sentence with a mismatched comparator and does not correspond to the HFD LEfSe result in Figure 5G/I. Not noted previously.
8. **Figure 5I caption text.** The accessibility description for Figure 5I labels the HFD/HFD+MP cladogram taxa as "enriched in normal diet mice / normal diet plus microplastics mice" — a caption transcription error in the published figure description. Minor, but noted.
9. **"No power calculation / no exclusion criteria / no covariates".** Confirmed (nothing of the kind appears), as is the prior note's statement that there are no correlation analyses — I add explicitly that there is **no enrichment analysis and no beta-diversity test** either.

Everything else in the prior note that I could check was accurate; I found no claim in it that contradicts the primary source other than items 2–6 above, which are refinements rather than outright errors.

## URLs/files actually read

1. `C:\Users\Soham\Desktop\microplastic\research\literature\papers\ehp2023_ps_hfd.md` — prior session notes (read in full).
2. `C:\Users\Soham\AppData\Local\Temp\claude\d28fb66f161d523b\5dbfd995-fe63-4ce9-8501-246b4fa39b6f\scratchpad\dl\PMC9945580.txt` — Europe PMC full-text body, lines 1–278 (everything before `### References`), read in full.
3. `C:\Users\Soham\AppData\Local\Temp\claude\d28fb66f161d523b\5dbfd995-fe63-4ce9-8501-246b4fa39b6f\scratchpad\dl\PMC9945580.xml` — the matching JATS XML, queried by grep/python for: data-availability, accession identifiers (PRJNA/PRJDB/DRA/GSE/SRA), `<sec sec-type>` listing, `<supplementary-material>` and `<media xlink:href>` entries, and the full `associated-data` section text.
4. `https://pmc.ncbi.nlm.nih.gov/articles/PMC9945580/` — fetched via a subordinate fetch agent solely to confirm the absence of a Data Availability statement/accession and to list the supplemental file names. (`https://ehp.niehs.nih.gov/doi/10.1289/EHP11072` was attempted and failed with a DNS error; the publisher page was therefore NOT read.)

**Not read:** the Supplemental Material PDFs `ehp11072.s001.acco.pdf` and `ehp11072.smcontents.508.pdf` — hence Table S1 (TaqMan primer IDs), Table S2 (all per-group summary data), Figure S1 (ILC gating), Figure S2 (macrophage gating), Figure S3 (MP dispersion stability, plasma GFP), Figure S4 (serum ALT/TG/NEFA/LDL-C/palmitate) are unverified.

## UNVERIFIED items

1. **Table S2 contents** — all figure-panel numeric values (body weight, AUCs, villus metrics, flow ratios, SCFA and palmitate concentrations, NAS, Oil Red O area, diversity indices). Everything I report from figures comes from the Results narrative and legends, not from the plotted numbers.
2. **Figure S4A–E** — individual serum ALT, TG, NEFA, LDL-cholesterol and palmitic acid values and which analytes were individually significant; only the combined Results sentence was verified.
3. **Figure S3A–E** — MP dispersion-stability data and the plasma GFP panels underpinning the deposition claim.
4. **Figures S1–S2** — flow-cytometry gating strategies; ILC and macrophage definitions could not be audited.
5. **Table S1** — TaqMan assay IDs; cannot resolve the *Sglt1* vs *Sglt2* discrepancy.
6. **Whether the qPCR biological n is 6 or 10**, and whether triplicates were treated as replicates in the ANOVA. Unresolvable from the article.
7. **Whether MODE-K n is 6 or 10.** Unresolvable from the article.
8. **Which ITT dose (0.5 or 0.75 U/kg) was actually used.** Unresolvable from the article.
9. **Whether a 16S deposit exists anywhere** — no accession is declared, so none could be tested. I searched the XML and the PMC record and found none; I did not search SRA/ENA by author name, and the publisher's own page could not be loaded (DNS failure), so a statement on the EHP site cannot be fully excluded (the PMC record normally mirrors it).
10. **Any interaction statistic** — none exists to verify; I could not recompute one because Table S2 was not read and the qPCR n is ambiguous.
11. **Whether the ND/ND+MP identical means for *Scd1*/*Elovl6*/*Fasn* (and *Il6*/*Il1b*) are a typesetting error or real** — requires Table S2.
12. **Exact p-values** for everything except the Table 2 phylum abundances; the paper reports only star thresholds.
13. **Particle characterisation data** (size distribution, zeta potential, endotoxin, surfactant) — not merely unread, but not generated by the authors.

---

# Cross-paper synthesis

## A. Interaction testing - the matrix the team asked for
| # | Paper | Design | Model fitted | Interaction term fitted? | Interaction F/p reported? | Basis of the "exacerbation" claim |
|---|---|---|---|---|---|---|
| 1 | Nat Metab 2026 | 3-arm (Ctrl/NP/MP) x separate chow and HFD cohorts | one-way ANOVA / Kruskal-Wallis + FDR; two-way ANOVA only for condition x TIME; ANCOVA (body weight) for VO2 | **No** (particle x diet never crossed in a model) | No | pairwise within the HFD cohort, plus a separate chow cohort |
| 2 | J Nanobiotech 2026 | true 2x2, n = 10 | **two-way ANOVA + Tukey** | **Yes (implicitly)** | **No - no F/p anywhere** | pairwise HFD+PSNP vs HFD; omics only in 2 of 4 cells |
| 3 | Adv Sci 2026 | two separate 2x2s (juvenile, adult) | **two-way ANOVA (HFD x PS-NP) + Tukey**; DESeq2 pairwise contrasts | **Yes for phenotypes** | **No for any in vivo phenotype** (only in vitro treatment x time, F/p in unread Table S3) | pairwise; age-dependence claim compares "significant in juveniles" vs "not significant in adults" across two separate models |
| 4 | EHP 2024 | 2x2 plus separate ABX and FMT cohorts | **one-way ANOVA + Tukey** across four cells; RM-ANOVA (group x time) for weight curves | **No** | No | significant-in-HFD vs non-significant-in-ND; mediation inferred from nulls in a separate ABX cohort |
| 5 | IJBS 2025 | 2x2 layout, n = 6 | **one-way ANOVA + Tukey only** | **No** | No | pairwise; headline "20% vs 61%" uses two different reference groups |
| 6 | EHP 2023 | true 2x2, n = 10 | **one-way ANOVA over four groups + Holm-Sidak** | **No** (0 hits for "two-way") | No | pairwise; two readouts (iPGTT, goblet cells/Muc2) actually changed under ND too |

**Nothing in this set tests the interaction the team cares about.** Papers 2 and 3 are the only ones whose fitted model even contains the term; neither prints it. The defensible sentence for the team's own report: across six studies, diet-conditional microplastic effects are claimed exclusively on the difference-of-significance pattern, and no published diet x particle interaction statistic exists in this literature.

## B. Data deposition
| # | Paper | Omics generated | Accession | Resolves? |
|---|---|---|---|---|
| 1 | Nat Metab 2026 | snRNA-seq, lipidomics, proteomics | GSE291616 / ST003776 / PXD072158 (plus bonndata histology 10.60507/FK2/M42B8O) | **All three VERIFIED this session** (GEO 6 GSMs; PRIDE title matches; MW 28 samples). bonndata not checked. |
| 2 | J Nanobiotech 2026 | bulk RNA-seq, untargeted metabolomics | **none** ("available from the corresponding author upon reasonable request") | n/a - not deposited |
| 3 | Adv Sci 2026 | bulk RNA-seq (40 libraries) | **none** ("available from the corresponding author upon reasonable request") | n/a - not deposited |
| 4 | EHP 2024 | 16S V3-V4 | **PRJNA990701** | **VERIFIED** (NCBI + ENA; 20 MiSeq amplicon runs = 4 groups x n=5, so ABX and FMT 16S appear NOT deposited) |
| 5 | IJBS 2025 | bulk liver RNA-seq | **none** ("included in this published article and its additional files") | n/a - not deposited |
| 6 | EHP 2023 | 16S | **none; no Data Availability statement at all** | n/a - not deposited |

Four of six deposited nothing. No paper in the set provides a code-availability statement or analysis repository.

## C. Direction summary for the pathway themes
| Theme | #1 Nat Metab | #2 J Nanobiotech | #3 Adv Sci | #4 EHP 2024 | #5 IJBS | #6 EHP 2023 |
|---|---|---|---|---|---|---|
| Hepatic lipid / steatosis | **UP (MP only)** chow and HFD; NP DOWN under HFD | **UP** (HFD only) | **UP** (juveniles only) | **UP** (HFD only) | **UP** (ORO area; hepatic TG never measured) | **UP** (HFD only) |
| Lipogenesis (Fasn/Acly/Acaca/Srebf1) | **DOWN** (module level; genes never named) | **unchanged** (qPCR, Fig. S4; genes never named) | **UP** (Fasn/Fas, Srebp1c, Lxra; Acly/Acaca not measured) | not measured | not measured | **UP** (Fasn 2.05->9.82, Scd1 1.22->3.03, Elovl6 1.44->4.28; Acly/Acaca/Srebf1 not measured) |
| PPARalpha / FAO | VIPER-inferred only and non-significant; no transcript or protein | **DOWN** (Cpt1a, Acox1, Acsl3, Cyp4a31 mRNA; ACSL3/CPT1A protein). Ppara itself not measured | FAO biphasic, **DOWN under co-exposure** (Cpt1a, Acadl, Pgc-1a) with CPT1 activity DOWN; Ppara/Acox1/Cyp4a not measured | not measured | not measured (PPARalpha antibody listed, no result) | **not measured at all** |
| Cholesterol synthesis | GO "cholesterol homoeostasis" (direction unstated); bile acids UP, Nr1h3 DOWN | **DOWN** (steroid/terpenoid KEGG; HMGCS1 protein DOWN) | not assayed | not measured (serum TC/LDL-C UP) | not measured (serum CHO UP) | not measured |
| Inflammation | **absent** (cell counts + 48-plex cytokines) | **absent** (NAS, serum factors) | essentially not studied | serum LPS/cytokines UP vs ND only | macrophages and neutrophils UP on BOTH diets | **UP** (Tnfa 2.44->5.93, Il6 2.22->5.51, Il1b 1.66->6.17) |
| Fibrosis | absent | absent | not studied | not studied | collagen UP by histology, fibrosis genes null | not studied |
| Oxidative stress | **not measured** | **not measured** | ROS measured (in vitro / mitochondria) | not measured | HO-1 UP in vitro | not measured |

The PPARalpha direction conflict the team is tracking is real and sits between papers 2 and 3 (suppressed FAO) and the PE literature. Papers 1, 4, 5 and 6 contribute nothing to it, because none of them measured Ppara, Cpt1a, Acox1 or Cyp4a.

## D. Shared blind spots across all six
- **No endotoxin testing** in any of the six; no surfactant, leachate or residual-monomer characterisation; no formulation-matched vehicle control (all use water, PBS or a jelly cube).
- **Male-only** in all six (paper 1 uses both sexes in one fate-mapping line and female donor hepatocytes for the co-culture).
- **Cage/litter is never the unit of analysis**, and cage composition is reported in only one paper (EHP 2024: 3 mice/cage, 2 cages/group, n = 6 treated as independent - in a microbiota-transfer study).
- **Time of sacrifice / Zeitgeber time** is reported only by paper 3 (ZT4-ZT6). Papers 2 and 5 fast mice overnight before collection, which by itself pre-activates hepatic PPARalpha/FAO.
- **Blinding** is stated only by papers 1 (histopathology only) and 3 (histology/quantification).
- **No power calculation** anywhere, and no multiplicity control across the dozens of phenotype endpoints in any paper.

# URLs and files actually read (sections 1-2; per-paper lists are inside sections 3-6)
- https://www.ebi.ac.uk/europepmc/webservices/rest/search (six DOI lookups, to obtain PMIDs/PMCIDs)
- https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13690586/fullTextXML -> local dl/PMC13690586.xml and natmetab_body.txt (FULL TEXT read)
- https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13173803/fullTextXML -> local dl/PMC13173803.xml and jnano_body.txt (FULL TEXT read)
- https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds (GSE291616) and esummary.fcgi?db=gds&id=200291616 - GEO verification; the GEO web page itself returns a reCAPTCHA to automated fetches, so E-utilities is the primary-source route
- https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD072158 - PRIDE verification
- https://www.metabolomicsworkbench.org/rest/study/study_id/ST003776/summary - Metabolomics Workbench verification
- Prior-session notes read before verifying: C:\Users\Soham\Desktop\microplastic\research\literature\papers\natmetab2026_size_macrophage.md and jnanobiotech2026_ncor1_ppara.md (papers 3-6 notes were read by the per-paper verification agents, listed in their sections).

# Global UNVERIFIED items
1. **All supplementary material for all six papers** (figures and tables are not in the Europe PMC XML). Load-bearing instances: Nat Metab Supplementary Figs. 1-8 (particle DLS/zeta values, serum lipids) and Tables 1-7 (DEG counts, ORA results, TF activity, proteomics, antibodies, primers); J Nanobiotech Figs. S1-S5 (serum inflammation/fibrosis, DEG volcano, metabolomics PCA, the lipogenesis/uptake qPCR panel, the entire rapamycin rescue) and Tables S1-S2; Adv Sci Tables S1-S3 (Table S3 is the only place any interaction F/p appears) and Figs. S1-S5; EHP 2024 Figs. S1-S5 and Tables S1-S3; IJBS Figs. S1-S5 and Tables S1-S4 (Table S4 is the only documentation of its RNA-seq pipeline); EHP 2023 Tables S1-S2 and Figs. S1-S4.
2. **Nat Metab:** proteomics n per group; which animals contributed to which panel; how many points Grubbs' test removed; whether the 6 snRNA-seq libraries were a single batch; measured particle size/PDI/zeta; the 28-vs-27 lipidomics sample-count discrepancy; the bonndata histology deposit (10.60507/FK2/M42B8O) was not resolved.
3. **J Nanobiotech:** which 4 of 10 animals were sequenced and by what rule; RNA-seq depth/mapping/batch (never reported, unverifiable in principle); SEM vs SD; p62 direction; whether a rapamycin-alone arm exists; Hmgcs1-vs-Hmgcs2 primer and antibody specificity; control-diet composition.
4. **Adv Sci:** whether any in vivo interaction F/p exists (Table S3); mice per cage; which 5 of 6 animals were sequenced; software versions; control chow identity; fasting state at sacrifice.
5. **EHP 2024:** all of Table S3 (every per-animal value and p); whether serum LPS and cytokines exceed HFD rather than only ND; which cohort each of the 20 deposited SRA runs belongs to (the "ABX/FMT 16S not deposited" conclusion is an arithmetic inference); the ampicillin-vs-penicillin and Pearson-vs-Spearman contradictions.
6. **IJBS:** RNA-seq n per group (never stated anywhere); the DE tool and version; whether the 8-wk and 20-wk cohorts are separate animals; up/down direction of the Fig. 2F KEGG bars; HepaRG differentiation state.
7. **EHP 2023:** true biological n for qPCR (Methods "six mice" vs Table 1 "n = 10") and MODE-K (6 vs 10); the ITT insulin dose (0.5 vs 0.75 U/kg); whether the byte-identical Table 1 values for Scd1/Elovl6/Fasn and for Il6/Il1b are a typesetting error; Sglt1 vs Sglt2; whether any 16S deposit exists anywhere (none declared; SRA was not searched by author name).
8. **Endotoxin, surfactant and leachate content of the particles in all six studies** - never measured by any of the authors, so unverifiable in principle rather than merely unread.
