# References (Phase 1 roadmap)

**Verification date:** 2026-11-21. Each entry gives its **access level**, meaning what was actually read: FT = full text, SI = supplementary material read, ABS = abstract only, GEO = the accession record was read. "Key locations" point to the figures, tables or sections that support the claims we rely on. **UNVERIFIED** marks anything not checked against the primary source.

Detailed per-paper notes, with quotes and the URLs read, are in the session scratchpad files listed at the end. The earlier per-paper records are in `research/literature/papers/`.

---

## A. Microplastic–liver original research

1. **Jung W, Abushukair H, … Joshi AD, Oh TG (2026).** Spatial transcriptome mapping identifies Ppara-Anxa2 cross-talk in microplastic-induced hepatotoxicity. *Sci Adv* 12(25):eaec8681. https://doi.org/10.1126/sciadv.aec8681 · PMC13274612 · **FT + SI (Tables S1–S18)**
   - **Design and n:** Methods "Animals and in vivo treatments" and "Bulk RNA-seq and analysis"; Table S1 (phenotype n = 5/group).
   - **Pairwise DE:** Tables S5/S6. **Interaction table:** Table S7 (76 genes padj < 0.05; mostly skeletal-muscle genes; not mentioned in the text).
   - **Lipogenic genes down:** Fig. 2D. **Ppara:** Fig. 2F legend (vs CD+Veh only). **Spatial findings:** Figs 3D, 4E, 5E–F, 6C. **ChIP and in vitro:** Figs 6H, 6M–N, Table S18.
   - **Re-used data:** GSE245069 (fig. S2F–G), GSE8295, GSE282998, SRX5142528/SRX5028180.
2. **Djouina M et al. (2023).** Oral exposure to polyethylene microplastics induces inflammatory and metabolic changes and promotes fibrosis in mouse liver. *Ecotoxicol Environ Saf* 264:115417. https://doi.org/10.1016/j.ecoenv.2023.115417 · **FT** (figure images not re-viewed this session)
   - **Key locations:** Methods 2.1–2.9; Fig. 1B–G (qPCR and biochemistry); Fig. 7 (CCl4); Results 3.1–3.6.
   - **Fold changes read from bar charts in the earlier session: UNVERIFIED.**
3. **Li Q, Zhu K, et al. (2025).** Polystyrene microplastics induce liver fibrosis and lipid deposition in mice through three hub genes revealed by the RNA-seq. *Sci Rep* 15:2583. https://doi.org/10.1038/s41598-025-86810-5 · PMC11747203 · **FT + SI**
   - **Key locations:** Fig. 6 (lipid-gene qPCR); Fig. 8E–H (hub qPCR; Fmo2 not confirmed); Suppl. Fig. 1 (oxidative stress null).
   - **Data:** GSE245069.
4. **Chiu HW, Chu CW, Huang CC, Chia ZC, Wang YL, Lee YH (2025).** Polystyrene microplastics induce hepatic lipid metabolism and energy disorder by upregulating the NR4A1-AMPK signaling pathway. *Environ Pollut* 369:125850. https://doi.org/10.1016/j.envpol.2025.125850 · PMID 39956511 · **ABS + GEO (GSE277772)**
   - **Missing source: the full text is needed.** Particle size, statistics and validation are all UNVERIFIED.
5. **Xu R, Cao JW, Lv HL, Geng Y, Guo MY (2024).** Polyethylene microplastics induced gut microbiota dysbiosis leading to liver injury via the TLR2/NF-κB/NLRP3 pathway in mice. *Sci Total Environ* 917:170518. https://doi.org/10.1016/j.scitotenv.2024.170518 · **ABS only. Missing source.**
6. **Makdissi N et al. (2026).** Size-dependent plastic exposure disrupts macrophage function and tissue-specific metabolism. *Nat Metab* 8(10):2084–2101. https://doi.org/10.1038/s42255-026-01615-8 · **FT** (supplement not readable)
   - **Key locations:** Fig. 1g (threshold conflict q ≤ 0.05 vs q < 0.1); ED Fig. 8h,i (ANCOVA); Methods statistics.
   - **Data:** GSE291616, PXD072158, ST003776.
7. **Fang Y et al. (2026).** Polystyrene nanoplastics induce hepatic steatosis by disrupting autophagic degradation of NCoR1 and suppressing PPARα signalling. *J Nanobiotechnol* 24:437. https://doi.org/10.1186/s12951-026-04304-6 · **FT** (Figs S4/S5 unread)
   - **Key locations:** Table 1 (particle characterisation); omics only for HFD vs HFD+PSNP; two-way ANOVA with no interaction reported.
8. **Xu P et al. (2026).** Developmental environmental chronotoxicity by polystyrene nanoplastics exacerbates high-fat diet-induced pediatric MASLD. *Adv Sci* e77840. https://doi.org/10.1002/advs.77840 · **FT** (Table S3 unread)
   - **Key locations:** ZT4–ZT6 sampling (Methods).
   - **Design:** two separate 2 × 2 cohorts, juvenile and adult. Exposure starts at weaning, so this is not a gestational design.
9. **Zhai Z et al. (2024).** Long-term exposure to polystyrene microspheres and high-fat diet-induced obesity in mice: evaluating a role for microbiota dysbiosis. *Environ Health Perspect* 132(9):97002. https://doi.org/10.1289/ehp13913 · **FT**
   - **Causal arms:** antibiotic depletion (ABX) and faecal transplant (FMT).
   - **Data:** PRJNA990701 (16S, 20 runs).
10. **Ahn J et al. (2025).** Polystyrene nanoplastics exacerbate HFD-induced MASLD by reducing cathepsin activity and triggering large vacuole formation. *Int J Biol Sci* 21(9):3867–3885. https://doi.org/10.7150/ijbs.108268 · **FT** (Table S4 unread)
    - DEGs were called at uncorrected p < 0.05. The abstract's "no effects on chow" contradicts the paper's own results.
11. **Okamura T et al. (2023).** Oral exposure to polystyrene microplastics of mice on a normal or high-fat diet and intestinal and metabolic outcomes. *Environ Health Perspect* 131(2):27006. https://doi.org/10.1289/ehp11072 · **FT**
    - **Key location:** Table 1. Several genes have byte-identical ND and ND+MP values (data-integrity question).
12. **Zhao J et al. (2024).** Obesogenic polystyrene microplastic exposures disrupt the gut-liver-adipose axis. *Toxicol Sci* 198(2):210–220. https://doi.org/10.1093/toxsci/kfae013 · PMC10964747 · **FT** (supplement not read)
13. **Park S et al. (2026).** Chronic PET-microplastic exposure: disruption of gut-liver homeostasis and risk of hepatic steatosis. *Adv Sci* 13(1):e12030. https://doi.org/10.1002/advs.202512030 · PMC12766997 · **FT + SI**
    - **Key location:** Table S3 (PET not detected in liver).
14. **Han J et al. (2025).** Chronic nanoplastic exposure promotes the development and progression of MASH. *Liver Int* 45(8):e70224. https://doi.org/10.1111/liv.70224 · PMC12257899 · **FT**
    - **No transcriptomics.** Unbalanced 2 × 2 (6/12/6/12) analysed by one-way ANOVA.

**Context only (from the earlier review, not re-verified this session):** Roh Y et al. 2024 *IJMS* 25:4964 (https://doi.org/10.3390/ijms25094964; PS × leptin knockout, sign reversal); Luo G et al. 2026 *Liver Int* 46:e70635 (https://doi.org/10.1111/liv.70635; ABS, "low dose protective"); Horvatits T et al. 2022 *EBioMedicine* 82:104147 (https://doi.org/10.1016/j.ebiom.2022.104147; human cirrhotic liver; **which polymers were detected is UNVERIFIED**).

**New candidates found but not reviewed:**
- Maurya et al. 2024, Research Square preprint rs-5196445: HDPE in Wistar rats; **not peer reviewed**; data on request.
- Lyu et al. 2026, *J Hazard Mater*, https://doi.org/10.1016/j.jhazmat.2026.143575: pristine vs chlorinated PE with liver transcriptomics; **methods UNVERIFIED**.
- PET single-gavage rat study, *Drug Chem Toxicol* 2026, PMID 41853959: GSE303044, 141 samples; ABS + GEO only.

## B. Statistical and modelling methods

15. **Wu D, Smyth GK (2012).** Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res* 40:e133. https://doi.org/10.1093/nar/gks461 · **FT** · Table 3: error control at residual df 6–27.
16. **Wu D, Lim E, Vaillant F, Asselin-Labat ML, Visvader JE, Smyth GK (2010).** ROAST: rotation gene set tests for complex microarray experiments. *Bioinformatics* 26:2176–2182. https://doi.org/10.1093/bioinformatics/btq401 · **FT** · Section 5.1: 2 × 2 factorials with 3 or 5 per cell. (fry is limma's fast approximation to ROAST and has no separate paper.)
17. **Schurch NJ et al. (2016).** How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? *RNA* 22:839–851. https://doi.org/10.1261/rna.053959.115 · **FT** · Figs 1 and 4: with 3 replicates, 20–40% of DE genes are found.
18. **Ramasamy A, Mondry A, Holmes CC, Altman DG (2008).** Key issues in conducting a meta-analysis of gene expression microarray datasets. *PLoS Med* 5:e184. https://doi.org/10.1371/journal.pmed.0050184 · **FT** (a guide, not original research) · Issue 6, Table 6.
19. **Sweeney TE, Haynes WA, Vallania F, Ioannidis JP, Khatri P (2017).** Methods to increase reproducibility in differential gene expression via meta-analysis. *Nucleic Acids Res* 45:e1. https://doi.org/10.1093/nar/gkw797 · **FT** · Figs 3–4, Table 3.
20. **Haynes WA, Vallania F, Liu C, et al. (2017).** Empowering multi-cohort gene expression analysis to increase reproducibility. *Pac Symp Biocomput* 22:144–153. https://doi.org/10.1142/9789813207813_0015 · **FT**
21. **Hong F, Breitling R, McEntee CW, Wittner BS, Nemhauser JL, Chory J (2006).** RankProd: a bioconductor package for detecting differentially expressed genes in meta-analysis. *Bioinformatics* 22:2825–2827. https://doi.org/10.1093/bioinformatics/btl476 · **ABS**
22. **Leek JT, Scharpf RB, Bravo HC, et al. (2010).** Tackling the widespread and critical impact of batch effects in high-throughput data. *Nat Rev Genet* 11:733–739. https://doi.org/10.1038/nrg2825 · **FT** (a Perspective) · Fig. 1, Table 1. When batch is confounded with biology, the two cannot be separated.
23. **Zhang Y, Parmigiani G, Johnson WE (2020).** ComBat-seq: batch effect adjustment for RNA-seq count data. *NAR Genom Bioinform* 2:lqaa078. https://doi.org/10.1093/nargab/lqaa078 · **FT** · Discussion: "may not work well" with confounded designs.
24. **Teufel A et al. (2016).** Comparison of gene expression patterns between mouse models of nonalcoholic fatty liver disease and liver tissues from patients. *Gastroenterology* 151:513–525. https://doi.org/10.1053/j.gastro.2016.05.051 · **ABS** (section and figure references UNVERIFIED)
25. **Vacca M et al. (2024).** An unbiased ranking of murine dietary models based on their proximity to human MASLD. *Nat Metab* 6:1178–1196. https://doi.org/10.1038/s42255-024-01043-6 · **FT** · Figs 3–5: direction agreement on human-significant genes and pathways.
26. **Bernau C, Riester M, Boulesteix AL, et al. (2014).** Cross-study validation for the assessment of prediction algorithms. *Bioinformatics* 30:i105–i112. https://doi.org/10.1093/bioinformatics/btu279 · **FT** · Sections 3.1–3.4, Figs 2–4.
27. **Patil P, Parmigiani G (2018).** Training replicable predictors in multiple studies. *PNAS* 115:2578–2583. https://doi.org/10.1073/pnas.1708283115 · **ABS**
28. **Ambroise C, McLachlan GJ (2002).** Selection bias in gene extraction on the basis of microarray gene-expression data. *PNAS* 99:6562–6566. https://doi.org/10.1073/pnas.102102699 · **ABS**
29. **Varma S, Simon R (2006).** Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics* 7:91. https://doi.org/10.1186/1471-2105-7-91 · **FT** · Table 1, Figs 1–4. The abstract (38%) and Table 1 (22.2%) disagree; the discrepancy is in the paper itself.
30. **Fielden MR et al. (2007).** A gene expression biomarker provides early prediction and mechanistic assessment of hepatic tumor induction by nongenotoxic chemicals. *Toxicol Sci* 99:90–100. https://doi.org/10.1093/toxsci/kfm156 · **ABS**
31. **Natsoulis G et al. (2005).** Classification of a large microarray data set: algorithm comparison and analysis of drug signatures. *Genome Res* 15:724–736. https://doi.org/10.1101/gr.2807605 · **ABS**

**Tools already used in our analysis (from the final report; citations still to be checked by the team):** Love MI, Huber W, Anders S 2014 *Genome Biol* 15:550 (DESeq2); Soneson C, Love MI, Robinson MD 2015 *F1000Research* 4:1521 (tximport); Patro R et al. 2017 *Nat Methods* 14:417 (Salmon). The edgeR-QL, limma-voom, fgsea, ashr and MSigDB citations are still marked [TEAM] in the final report.

## C. Datasets

See `dataset_feasibility_screen.csv`. All accession records were read directly on GEO on 2026-11-21 unless marked otherwise.

## D. Missing sources (team action)

| Item | Why it matters | Route |
|---|---|---|
| Chiu 2025 full text | Particle size, statistics, Nr4a1 validation | Library login, save the PDF to Downloads |
| Xu 2024 STOTEN full text | Only other PE oral liver study with a mechanism | Library login |
| Djouina Fig. 1B fold-change values | Ppara direction per bead arm | Re-read the figure by eye from the PDF (it is already in Downloads) |
| Supplements: Fang S4/S5, Xu P Table S3, Ahn Table S4, Makdissi Suppl. Fig. 1 and Tables 1–5 | Pipeline details, the only reported interaction statistics | Browser download |
| Horvatits 2022 polymer list | Jung's claim that PE is the dominant polymer in human liver | Read Results/Table |
| Paper behind GSE261893/GSE261725 (if any) | Design details (particle size, dose) | Search again later, or contact the submitters |

## E. Working notes produced this session

The detailed per-paper and per-dataset verification notes, with quotes, locations and the URLs read, are copied into `research/roadmap_phase1/working_notes/`:

- `agent_lit_core.md`: Jung, Djouina, Li, Chiu, Xu 2024
- `agent_lit_interaction.md`: Makdissi, Fang, Xu P, Zhai, Ahn, Okamura
- `agent_lit_independent.md`: Toxicol Sci 2024 section only. The file was cut off mid-write; the remaining findings appear in the matrix and here.
- `agent_methods.md`: methods papers
- `agent_datasets.md`: GEO screen

AI agents produced these notes. Each records what was read, but the team has not checked them.
