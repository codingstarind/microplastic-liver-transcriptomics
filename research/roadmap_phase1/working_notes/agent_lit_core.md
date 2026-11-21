# Methods-focused literature matrix: microplastics and liver transcriptomics (5 core papers)

Verification session: 2026-11-21. Earlier local notes (`research/literature/papers/*.md`, `published_gene_evidence.csv`) were treated as leads and re-checked against the primary sources below.

**Access levels**

| # | Paper | What was read this session |
|---|---|---|
| 1 | Jung 2026 Sci Adv | FULL TEXT (Europe PMC XML, PMC13274612) + supplementary PDF `sciadv.aec8681_sm.pdf` (figs S1–S6 legends) + `sciadv.aec8681_tables_s1_to_s18.xlsx` (all 18 tables parsed) + GEO GSE294180 SOFT record |
| 2 | Djouina 2023 EES | FULL TEXT (publisher PDF `1-s2.0-S0147651323009211-main.pdf`, which the user downloaded on 2026-11-20; text extracted). Figure *images* could not be rendered this session (no PDF rasteriser), so bar heights or fold changes read from figures earlier are NOT re-verified. Suppl. Fig. 1 not read. |
| 3 | Li 2025 Sci Rep | FULL TEXT (Europe PMC XML, PMC11747203) + supplement `41598_2025_86810_MOESM1_ESM.pdf` (primer table, Suppl. Figs 1–4 legends, full blots) + GEO GSE245069 SOFT + GSM7835567 records |
| 4 | Chiu 2025 Env Pollut | ABSTRACT ONLY (Europe PMC core record, Crossref) + GEO GSE277772 record, which is linked to PMID 39956511 and gives the design. Publisher full text was unavailable: closed access per Unpaywall; ScienceDirect returned 403 and then a CAPTCHA, which I did not attempt; no PMC copy; Claude in Chrome was not connected. |
| 5 | Xu 2024 STOTEN (TLR2) | ABSTRACT ONLY (Europe PMC core record incl. MeSH). Closed access per Unpaywall; ScienceDirect 403; no PMC copy; no GEO/SRA link found. |

---

## 1. Jung W et al. 2026. "Spatial transcriptome mapping identifies *Ppara-Anxa2* cross-talk in microplastic-induced hepatotoxicity." *Sci Adv* 12(25): eaec8681

DOI: https://doi.org/10.1126/sciadv.aec8681 · PMID 42308314 · PMC13274612 · Access: FULL TEXT + full supplement.

### Exposure record
| Field | Value (location) |
|---|---|
| Species/strain/sex/age | Mouse, C57BL/6 (Jackson), **male**, 8 weeks old at start (Methods "Animals and in vivo treatments") |
| Polymer / size | Polyethylene, "10 to 150 μm in diameter, Cospheric" (in vitro Methods give the catalogue number: CPMS-0.96). Size, zeta potential and endotoxin were **not measured** |
| Dose / route / vehicle | 2 mg PE per mouse per day, oral gavage, in "0.1% Tween 20 + 99.9% molecular-grade water". Vehicle groups were gavaged daily with the same solution without PE (Methods) |
| Duration | 8 weeks (diet and gavage together) |
| Diet | MD = D09100310 (40 kcal% fat, 20 kcal% fructose, 2% cholesterol); CD = D09100304 (10 kcal% fat, matching sucrose), Research Diets, ad libitum |
| Housing | 22 ± 1 °C, 25 ± 5% humidity, 12 h light/dark; "Groups of mice (n = 4 or 5 per group) were housed together", i.e. **one cage per group**, so cage is confounded with treatment |
| Sampling time / fasting | **Not reported** (no time of sacrifice, no fasting statement) |
| Randomisation / blinding | Randomisation not reported. Only histology scoring was blinded (Methods "Histology") |

### Q1. Scientific question
How does oral PE alter liver gene expression, cell composition and zonation, alone and on a MASH diet? Which nuclear-receptor (NR) programme mediates the effect? (Abstract; last two paragraphs of the Introduction.)

### Q2. Design and n
- 2 × 2 factorial: diet (CD/MD) × gavage (Veh/PE). Groups are named G1 CD+Veh, G2 CD+PE, G3 MD+Veh, G4 MD+PE (table S-list sheet).
- **Phenotype n**: Table S1 holds raw data for 5 mice in each of the 4 groups (body mass, liver weight, glucose, ALT, TG, NAS components). This matches the Fig. 1B legend ("n = 5 within each of the four groups"); the Methods say "4 or 5".
- **Bulk RNA-seq n = 4 / 5 / 5 / 5 (19 libraries)**. Sources: GEO GSE294180 overall design ("Chow Diet (n = 4) … (n = 5)…"); column headers of table S8 (CD_Veh_1–4, CD_PE_1–5, MD_Veh_1–5, MD_PE_1–5). The main text never states it. One CD+Veh mouse therefore has phenotype data but no RNA-seq library.
- **Xenium**: "A total of 19 samples from two Xenium slides" (Methods "Spatial transcriptomic data processing"), with 4–5 biological replicates per group. How samples were allocated to the two slides is not reported.
- **RT-qPCR (liver)**: "four to five samples" (Methods). Table S18 shows n = 5 per group for Anxa2 (fig. S6B).
- **ChIP-qPCR**: n = 5 per group (table S18, Fig. 6H).
- **In vitro**: 3 replicates per condition for agonist, antagonist and siRNA (table S18). MTT viability: 6 per polymer (table S1). Whether these are wells or independent preparations is **not stated**.

### Q3. Datasets
- Generated: GSE294180, bulk RNA-seq raw gene-level counts (`GSE294180_MP_raw_counts_gene_level.csv.gz`, GPL30172 NextSeq 2000); Xenium data on Zenodo record 19603085 (Data availability).
- Re-used: GSE245069 (Li 2025 PS liver data, reprocessed through the authors' own pipeline, fig. S2F–G); GSE282998 (GW7647 dose series, B6129SF2/J); GSE8295 (WY14643 in WT vs Ppara-null 129S1/SvImJ mice, microarray, limma, probe 1419091_a_at); ChIP-seq SRX5142528 and SRX5028180 (from ref. 58, Liang et al. 2019).

### Q4. Models implemented
- **Bulk DE**: fastp 0.23.4 → Salmon 1.10.3 (GENCODE GRCm39) → tximport → DESeq2. Thresholds: |log2FC| ≥ 0.38 and FDR < 0.05 (Methods "Bulk RNA-seq and analysis"). Contrasts were pairwise across all 6 group pairs (tables S5, S6).
- **Interaction term: YES, but only in the supplement.** Table S7 "BULK_RNA-seq_Interaction_Analysis" gives DESeq2 columns PE_in_CD, Diet_Effect and Interaction (log2FC, p, padj). **76 genes have padj_Interaction < 0.05.** Table S12 gives a sample-based spatial interaction analysis per cell type (505 rows). Neither analysis is described in the Methods or Results, and every main-text claim ("exacerbates", "synergistic", "MD-specific PE genes") rests instead on pairwise contrasts or set subtraction. *Correction to the earlier local note, which said no interaction model existed.*
- **Content of the S7 interaction list**: dominated by skeletal-muscle transcripts (Acta1, Myh1, Myh4, Ttn, Ckm, Tnnt3, Tnnc2, Mb, Ryr1, Mylpf, Xirp2 …), several with interaction log2FC between −9 and −28. This pattern suggests muscle-tissue carry-over in one or a few libraries rather than biology (my inference; the paper does not discuss it). Ppara, Fasn, Acly and Anxa2 are **not** in the interaction list.
- **GSEA**: fgsea on ranked statistics, minimum set size 1, maximum 500, 10,000 permutations, custom set LIPID_METABOLISM_AND_TOXICITY. NES and FDR values are not given in the text.
- **Spatial**: Seurat 5.2.0, SCTransform, Harmony ("across experimental groups and samples"), clustering at resolution 0.8 → 32 clusters → 15 cell types. Two kinds of DE are reported (tables S10, S11):
  - cell-based FindMarkers: |log2FC| ≥ 0.38, minimum expressing fraction 0.33, FDR < 0.05;
  - sample-based pseudobulk DESeq2.

  Cell proportions: betareg with emmeans and Holm correction (sample level), and glmmTMB beta GLMM `Proportion ~ condition + (1|sample)`. Shannon diversity: lme4 with `(1|sample)` and Tukey (table S13).
- **WGCNA**: on central-hepatocyte pseudobulk (19 samples), signed network, soft power 8, minModuleSize 15, mergeCutHeight 0.25. Module–trait correlation is Pearson with an asymptotic p-value; no correction across modules or traits is stated. Blue module vs MD+PE: r = 0.68, p = 0.001 (fig. S6D legend).
- **Anxa2 co-expression**: cell-level Spearman correlation on SCT data, cut at R > 0.273, followed by Enrichr with WikiPathways (Fig. 6E–F).
- **Phenotype and qPCR statistics**: ordinary **one-way ANOVA + Tukey** on the 4 groups (Fig. 1, Fig. 6 and fig. S1 legends). This is not a two-way ANOVA. Raman band test: Mann-Whitney (Fig. 5L).
- **PPI/hub analysis**: none (the only "hub" is the WGCNA intramodular hub).

### Q5. Fit of model to question
- DESeq2 on Salmon/tximport counts suits n = 4–5 bulk libraries.
- The 2×2 question ("does PE worsen MASH?") requires an interaction contrast. The authors computed one (S7/S12) but built their narrative on set subtraction ("DEGs uniquely induced by PE under the MD condition, excluding those also altered by PE in CD", Results Fig. 4C). Set subtraction confuses differences in power with biology: there are more MD DEGs because MD effects are large.
- One-way ANOVA + Tukey on 4 groups answers "which groups differ", not "is there synergy".
- The mixed models (glmmTMB/lme4 with sample as a random effect) are the right unit-of-analysis choice for spatial proportions and diversity. Cell-level FindMarkers and Spearman correlations treat thousands of cells per mouse as independent, which is pseudoreplication. The supplement does also report sample-based DESeq2 alongside.
- Because there is one cage per group, cage is fully confounded with treatment and no model can separate them.
- The GSEA minimum set size of 1 is permissive.

### Q6. Inputs / features / outputs
- Bulk: raw counts → DEG lists (S5–S7), GSEA NES, NR heatmap (S8).
- Xenium: ~479 targeted genes (379-gene Mouse Tissue Atlassing panel + 100 add-on genes, table S2) → cell types, zonation, proportions, Shannon index, inflammation-score tertiles (0.33/0.67 quantiles), pseudobulk.
- WGCNA: pseudobulk expression → modules and eigengenes; targets = binary group traits.
- In vitro targets: Anxa2 and Acox1 mRNA.

### Q7. How results were evaluated
- FDR thresholds as above.
- Concordance of pseudobulk with bulk: Pearson correlation (Fig. 4A; coefficient not given in text). This is a platform-concordance check, not a test of the model.
- No held-out data, no permutation of group labels, no sensitivity analysis.
- I recounted DEGs from S5/S6 at |log2FC| ≥ 0.38 and padj < 0.05: CD+PE vs CD+Veh = 615; MD+Veh vs CD+Veh = 1458; MD+PE vs MD+Veh = 287. The text reports 600, 1437 and 279 ("uniquely altered"). The counts are close but not identical; the definition of "uniquely" may explain the difference.

### Q8. Independent validation
- **Liver RT-qPCR of Anxa2** (fig. S6B; 5 per group). Whether the RNA is from the same mice as RNA-seq is not stated; it is probably the same cohort, so this is technical, not biological, replication.
- **ChIP-qPCR** for PPARα at 3 Anxa2 peaks in the authors' own livers (Fig. 6H; antibody PA5-85125; 5 per group). No IgG or negative-region control and no %-input normalisation is described (Methods; table S3 lists only the antibody).
  - **Reporting problem:** in table S3 the forward and reverse primers for Peak3 (CTCACACTGCTCCAGCTGGT / ACCAGCTGGAGCAGTGTGAG) are exact reverse complements of each other. That pair cannot amplify a genomic interval as listed, so the table probably contains an error. This is the peak with the largest reported enrichment.
- **Protein**: none (no western blot or IHC).
- **Knockout**: only in external data (GSE8295, Ppara-null, WY14643, different strain). No Ppara-null mouse was exposed to PE.
- **In vitro loss of function** (Fig. 6M–N; S18): GW6471 antagonist and Ppara siRNA. Knockdown efficiency can be read from S18: Ppara mRNA about 0.24–0.61 of NT control.
  - **Design gap:** the S18 columns are "Control+Veh, GW6471, GW6471+PE" and "NT-siRNA, Ppara-siRNA, Ppara-siRNA+PE". There is **no PE-alone or NT-siRNA+PE arm**. The claim that inhibition "abolished PE-induced Anxa2 up-regulation" therefore lacks the PE-induced reference condition within those experiments. Fig. 6I likewise has no PE-only arm.
- **External cross-polymer check**: GSE245069 reprocessed; reported as Ppara up and Esr down in PS-exposed mice (fig. S2F–G).
- **Particle presence in tissue**: micro-Raman on post-Xenium slides (MD+PE vs MD+Veh only, fig. S5K–L) and O-PTIR on serial sections ("two independent tissue sections, four to seven regions per section").

### Q9. Main findings (with directions; bulk values from S5/S6, log2FC = second group vs first)
- **Phenotype (S1, 5 per group)**:
  - ALT up with PE in both diets (Fig. 1E);
  - liver TG higher in MD+PE than MD+Veh (Fig. 1F);
  - steatosis score up;
  - no liver-weight change (Fig. 1C, fig. S1B);
  - body weight up in CD+PE from week 7 (Fig. 1B).
- **De novo lipogenesis: DOWN with PE.**
  - CD (PE vs Veh): Fasn −1.51, Acly −1.92, Acaca −1.44, Me1 −1.31, Pnpla3 −3.06, Scd1 −1.02, Srebf1 −0.71 (all padj < 0.05).
  - MD (PE vs Veh): Fasn −1.42, Acly −1.10, Acaca −0.64, Me1 −0.83, Pnpla3 −2.35, Scd1 −1.25 (padj < 0.05); Srebf1 not significant.
  - Text: shared PE DEGs "down-regulation (e.g., Pnpla3, Acly, Me1, and Fasn)" (Fig. 2D).
- **PPARα: modest UP, and not significant within MD in bulk.**
  - Ppara: CD+PE vs CD+Veh +0.53 (padj 0.011); MD+PE vs CD+Veh +0.75 (padj 2e-5); **MD+PE vs MD+Veh +0.34 (padj 0.23, ns)**; MD+Veh vs CD+Veh +0.41 (padj 0.053).
  - The Fig. 2F legend compares only with CD+Veh, which mixes the diet and PE effects.
  - Nr1i3: MD+PE vs CD+Veh +0.91 (padj 0.0035); MD+PE vs MD+Veh ns.
  - Targets: Cyp4a10 up with PE in CD (+2.44, padj 4e-4); Acox1 not significant for PE in CD (+0.22); Cpt1a no change.
  - Spatial: Ppara and Nr1i3 up in Kupffer cells in MD+PE (Fig. 4E). WGCNA: Ppara is the hub of the MD+PE Blue module.
- **Cholesterol synthesis: DOWN with PE in CD.** Hmgcr −1.38 (padj 0.007), Sqle −1.45 (padj 0.04); Hmgcs1 −0.55 (ns); Ldlr ns. Under MD, PE further lowered Tm7sf2 (−0.80) and Ldlr (−0.58) and raised Mylip (+0.68), all padj < 0.05 (Fig. 2E genes; directions taken from S6 because the text gives none).
- **Lipid droplet / uptake: UP.** Plin4 +1.01 (CD), +0.98 (MD); Cd36 up only with diet.
- **Anxa2**: bulk PE effect in CD −0.38 (ns), in MD +0.19 (ns). It is up with diet: MD+Veh vs CD+Veh +1.6. The PE-specific Anxa2 rise is spatial (central hepatocytes, Fig. 6C) and by qPCR (fig. S6B), **not in bulk RNA-seq**.
- **Spatial**: midlobular hepatocyte proportion up with PE (beta regression, Fig. 3D); Kupffer cells, stellate cells and lymphocytes down in MD+PE vs CD+Veh (Fig. 5F); inflammation signature vs Shannon diversity R = −0.52 (Fig. 5E).
- **Cross-polymer**: Ppara up in PS (GSE245069) as well (fig. S2G).

### Q10. Limitations
**(a) Actual methodological/biological limitations of what they did**
1. Cage is fully confounded with treatment (one cage of 4–5 per group).
2. The synergy claims rest on one-way ANOVA, pairwise contrasts and set subtraction. Their own interaction table (S7) is neither reported in the text nor used, and it is dominated by apparent muscle-contamination genes.
3. The Ppara induction headline uses CD+Veh as the comparator. Within MD, the PE effect on Ppara is not significant in bulk.
4. Cell-level spatial DE and correlations are pseudoreplicated (partly mitigated by the sample-based tables).
5. Harmony was run "across experimental groups and samples", which risks absorbing true group effects.
6. The in vitro antagonist and siRNA experiments lack a PE-alone arm; n = 3; independent preparations are not stated.
7. ChIP-qPCR has no IgG or input control described, and the Peak3 primer pair is self-complementary as listed.
8. The polymer screen is not size-matched (55–5000 µm at 100 µg/ml).
9. PE is 10–150 µm, with no size, endotoxin or zeta characterisation.
10. Dose is per mouse, not per kg.
11. Raman was run only on MD+PE vs MD+Veh slides, with no procedural blanks and no quantification.
12. Unit inconsistency: "PE (100 μM)" in the siRNA paragraph vs 100 μg/ml elsewhere.
13. The RNA-seq n is not stated in the paper, and one CD+Veh mouse lacks a library.

**(b) Not attempted**
- Females; a dose–response; a time course.
- Non-plastic particle control; sham (non-gavaged) arm.
- In vivo Ppara knockout or hepatocyte-specific deletion under PE.
- Protein-level validation; lipidomics; gut or microbiome analysis.
- Particle biodistribution quantification.
- Reporting of time of sacrifice or fasting.

---

## 2. Djouina M et al. 2023. "Oral exposure to polyethylene microplastics induces inflammatory and metabolic changes and promotes fibrosis in mouse liver." *Ecotoxicol Environ Saf* 264:115417

DOI: https://doi.org/10.1016/j.ecoenv.2023.115417 · PMID 37651791 · Access: FULL TEXT (publisher PDF, CC BY-NC-ND); figure images and Suppl. Fig. 1 not viewed this session.

### Exposure record
| Field | Value (location) |
|---|---|
| Species/strain/sex/age | C57BL/6 (Janvier), **female**, "7 to 12 weeks old" (Methods 2.2). CCl4 experiment: 8-week-old females |
| Polymer / size | PE microbeads (Cospheric), FTIR-confirmed in a previous study. RB (red fluorescent, UVPMS-BR-1.090) 10–45 µm, median 36 µm; GB (green fluorescent, UVPMS-BG-1.00) 106–125 µm, median 116 µm (Methods 2.1). Supplied by ANSES (Acknowledgments) |
| Dose / route / vehicle | Diet spiked with 100 µg/g food (RB, GB) or 100 µg each (RB+GB, i.e. 200 µg/g total). Ad libitum feeding, no vehicle. Food intake and ingested dose were not measured |
| Duration | Exp. 1: 6 or 9 weeks. Exp. 2: 8 weeks + CCl4 0.2 ml/kg i.p. (olive oil 1:4). Frequency is "twice a week" in Methods 2.2 but "once a week" in the Fig. 7 legend, an unresolved internal inconsistency |
| Diet | Standard food (composition not reported); no metabolic diet |
| Housing / sampling | "housed under standard conditions". Cage numbers, time of sacrifice and fasting are not reported. CCl4 mice were sacrificed 24 h after the last injection |
| Randomisation / blinding | "randomly assigned" (method not given); blinding not reported (image quantification was done in QuPath) |

### Q1. Question
"explore the hepatic effects of sub-chronic PE MP exposure in mice" (end of Introduction). Secondary question: does dietary PE worsen CCl4-induced fibrosis? (Section 3.6.)

### Q2. Design and n
- Exp. 1: 4 groups (CT, RB, GB, RB+GB), n = 10 per group, at each of 2 time points (6 and 9 weeks), i.e. separate cohorts (Methods 2.2; every figure legend says n = 10).
- Exp. 2: 5 groups (CT/olive oil, CCl4, CCl4+RB, CCl4+GB, CCl4+RB+GB), n = 8 per group (Fig. 7 legend).
- No omics, so phenotype n equals assay n. Per-assay n is not reported separately.

### Q3. Datasets
None generated or re-used. Data availability: "Data will be made available on request." No accession of any kind.

### Q4. Models
- Mann-Whitney U test vs control (or vs the CCl4 group, marked #), GraphPad Prism, p < 0.05 (Methods 2.9).
- **No multiple-testing correction** across roughly 60 qPCR targets × 3 groups × 2 times.
- No ANOVA, no interaction term (PE × CCl4 or bead size × time), no DE tool, no enrichment, PPI or WGCNA.
- qPCR: SYBR Green, 2^−ΔΔCt, normalised to Polr2a. Primers "available upon request" (Methods 2.7).

### Q5. Fit
- Nonparametric pairwise tests with n = 10 are acceptable for single endpoints. The study is a multi-arm, multi-time design, however, and correction would be needed across arms, times and genes.
- "Exacerbation" of CCl4 fibrosis requires a PE × CCl4 interaction (a PE effect larger in CCl4 than in healthy liver). It was tested only as CCl4+PE vs CCl4.
- RB+GB has double the particle mass, so bead combination and dose are confounded.
- A targeted panel answers "do these chosen genes change", not a transcriptome-wide question.

### Q6. Inputs / outputs
Liver mRNA (targeted panels: FA uptake, β-oxidation, FA synthesis, TG synthesis, xenobiotic metabolism, oxidative stress, cytokines, cell-type markers, fibrosis markers); liver cholesterol, TG and FFA; MDA; serum ALT/AST; H&E foci (QuPath); Sirius red; Ki-67, α-SMA IHC; TUNEL.

### Q7. Evaluation
Raw p < 0.05 per comparison. No effect-size thresholds, no model validation.

### Q8. Independent validation
- Mostly orthogonal endpoints within the same mice:
  - Ki-67 mRNA confirmed by Ki-67 IHC (Fig. 4A–B);
  - HSC markers corroborated by Sirius red and α-SMA IHC in Exp. 2 (Figs 6–7);
  - lipogenic mRNA paired with liver TG and cholesterol biochemistry (Fig. 1E–G).
- No external dataset, no knockout, no protein for lipid genes.

### Q9. Main findings (text-verified, Results 3.1–3.6)
- **6 weeks**
  - Ppara: "A two-fold increase … after exposure to both types of PE microbeads". *Earlier figure reading (not re-verified this session): in Fig. 1B the rise was significant only in RB+GB, and GB alone was lower than control.*
  - Acox1, Cpt1a, Acadm unchanged.
  - "**Fasn, Acly and Acaca were significantly upregulated in all groups**".
  - Pparg and Scd1 up in RB+GB; Cd36 and Fabp1 up in RB+GB; Fabp1 down in RB; Gpam up in GB.
- **9 weeks**
  - Ppara and Pparg up in GB and RB+GB.
  - Cd36, Acaca, Mttp up in all PE groups.
  - RB: Fabp1 and Acox1 down, Dgat1 up. GB: Fatp2, Cpt1a, Gpam up.
  - *Earlier figure reading: Fasn not significant at 9 weeks.*
- **Biochemistry**: liver cholesterol up (GB); TG up (RB+GB); FFA unchanged (Fig. 1E–G).
- **Xenobiotic metabolism**: Cyp2a4 and Cyp2b10 up in all groups at 6 weeks; Cyp2a4 up at 9 weeks; Gsta1 up (6 weeks); Ugt1a1 up (RB+GB, 9 weeks).
- **Oxidative stress**: Nrf2 battery up; MDA down (RB, 6 weeks) but up (RB+GB at 6 weeks, RB at 9 weeks), i.e. the direction depends on time.
- **Inflammation**: foci larger (RB+GB, 6 weeks) and more numerous (GB, RB+GB, 9 weeks). Cytokines are mixed: Il1b down in single-bead groups, Il6 strongly up in RB+GB at 9 weeks.
- **Proliferation**: up (Mki67 and Ki-67 IHC); TUNEL unchanged.
- **Fibrosis**: Sirius red up (GB at 6 weeks; GB and RB+GB at 9 weeks). On CCl4: Sirius red and α-SMA areas higher in GB and RB+GB than CCl4 alone; Acta2, Pdgfa, Col1a1, Vim enhanced in RB+GB.
- **Cholesterol-synthesis genes**: not assayed.

### Q10. Limitations
**(a) Actual**
1. No multiple-testing correction.
2. No interaction test for the "exacerbation" claim.
3. RB+GB carries double the dose.
4. Fluorescent dyes, with no unlabelled-bead control.
5. Ingested dose unmeasured.
6. Wide starting age range (7–12 weeks).
7. Female only, with sex-biased targets: the authors note Cyp2a4 is "strongly female-biased".
8. CCl4 frequency contradiction between Methods and the Fig. 7 legend.
9. Primer sequences unpublished.
10. Bulk-tissue cell-type markers cannot separate per-cell expression changes from changes in cell composition.
11. The authors themselves state "it seems unlikely that PE MPs reached the liver" (Discussion), so the mechanism is presumed to be indirect (via the gut) but is not tested.

**(b) Not attempted**
- Transcriptomics; males; a metabolic diet.
- Particle tracking in liver.
- Gut or microbiome measurements in this study (cited from Djouina 2022).
- Dose–response; protein for lipid genes; endotoxin testing.

---

## 3. Li Q, Zhu K, et al. 2025. "Polystyrene microplastics induce liver fibrosis and lipid deposition in mice through three hub genes revealed by the RNA-seq." *Sci Rep* 15:2583

DOI: https://doi.org/10.1038/s41598-025-86810-5 · PMID 39833454 · PMC11747203 · Access: FULL TEXT + supplement.

### Exposure record
| Field | Value (location) |
|---|---|
| Species/strain/sex/age | ICR (outbred), **male**, "Five-week-old" (Methods "Animal model"). **Discrepancy:** GEO GSE245069 says "Six-week-old"; MeSH indexes C57BL, which is erroneous because both the text and GEO say ICR |
| Polymer / size | Polystyrene, 80 nm, 0.5 µm, 5 µm (Baseline ChromTech, Tianjin). Nominal sizes only; no size, surface, label or endotoxin characterisation |
| Dose / route / vehicle | 1 mg/L in sterilised drinking water, ad libitum; ingested dose not measured |
| Duration | 12 weeks |
| Diet | Unspecified standard diet |
| Housing / sampling | SPF, 12 h light/dark, 21–23 °C. Cage, time of sacrifice and fasting not reported |
| Randomisation / blinding | "equally and randomly divided"; blinding not reported ("ten random fields" per mouse for IHC) |

### Q1. Question
What mechanism links chronic oral PS nano/microplastics of three sizes to liver fibrosis, and which hub genes appear in liver RNA-seq? (Introduction; Abstract.)

### Q2. Design and n
- 4 groups (H2O, 80 nm, 0.5 µm, 5 µm), n = 10 allocated per group.
- **Phenotype n is 5–8**:
  - liver/body weight, ALT, AST n = 5–8 (Fig. 1 legend);
  - H&E n = 8; α-SMA WB n = 5; α-SMA IHC n = 8; α-SMA qPCR n = 6 (Fig. 2);
  - CD68 IHC n = 8; cytokine qPCR n = 6; WB n = 5 (Fig. 3);
  - lipid-gene qPCR n = 5–6 (Fig. 6);
  - oxidative-stress kits n = 7 (Suppl. Fig. 1 legend).
- **Sequencing n = 4 per group** (16 libraries; "Four independent biological replicate samples were sequenced for each group"; GEO has 16 GSMs). Platform: GPL28330 **DNBSEQ-T7**, at BGI.
- Hub-gene qPCR used "expanded" samples (n = 5–6, Fig. 8E–H). It is not stated whether these include the sequenced mice.
- Attrition from 10 to 5–8 is unexplained.

### Q3. Datasets
Generated: GSE245069 (`GSE245069_gene_expression.txt.gz`, `GSE245069_rna_expression_.txt.gz`; GEO says these hold "raw expression", aligned to GRCm38.p6). No re-used datasets.

### Q4. Models
- DE with limma ("linear models for the microarray data") on RNA-seq. It is not stated whether voom or limma-trend was used, or on which input (raw counts or FPKM). Thresholds: |log2FC| > 1 and adjusted p < 0.05; adjustment method not named. Three pairwise contrasts vs H2O.
- Venn and UpSetR overlap gives 87 shared DEGs (46 up, 41 down).
- Enrichment: DAVID (GO BP/CC/MF, KEGG, Reactome); background and FDR cut-off not stated.
- PPI: STRING 11.5 (score ≥ 0.4) → 85 nodes, 51 edges → Cytoscape 3.4.0 **CytoHubba "12 algorithms"** → hubs Acot3, Abcc3, Nr1i3, Fmo2. How the 12 rankings were combined is not stated. GeneMANIA was used for neighbours.
- Phenotypes: one-way ANOVA + **LSD** post hoc (SPSS 23), mean ± SEM. No WGCNA, no ML, no size-trend test.

### Q5. Fit
- limma (with voom) is valid for n = 4; without voom on raw counts it would be inappropriate, and the paper does not say which was done.
- |log2FC| > 1 is a conservative cut-off.
- **CytoHubba hub ranking is a weak inferential tool here.** It ranks topological centrality in a sparse network (51 edges among 85 nodes) built from prior-knowledge interactions. It does not measure the genes' relevance to fibrosis, and the "hub" label carries no statistical uncertainty.
- LSD post hoc does not control family-wise error across 3 sizes.
- Each size vs control is a reasonable design, but no dose (size) trend model was used.

### Q6. Inputs / outputs
Liver RNA-seq counts → DEG lists → 87 shared DEGs → PPI → 4 hub genes. Targets: fibrosis (Masson, α-SMA), inflammation (CD68, cytokines), lipid deposition (Oil Red O, TEM), lipid-gene qPCR.

### Q7. Evaluation
Adjusted p and FC thresholds; enrichment p-values; hub genes "validated" by qPCR in 5–6 mice. No external replication, permutation or network robustness check.

### Q8. Independent validation
- qPCR of hubs: Acot3, Abcc3 and Nr1i3 significantly up; **Fmo2 not significant** (Results; Fig. 8E–H).
- This is mRNA only, probably in an overlapping set of mice (not stated). There is no protein, knockdown or overexpression for any hub gene.
- Fibrosis and inflammation phenotypes are supported by histology, IHC and WB (α-SMA, IL-1β, IL-6), but these are not causally linked to the hub genes.
- No external dataset.

### Q9. Main findings
- Liver weight and liver/body ratio decreased (more with 0.5 and 5 µm); ALT, AST and body weight unchanged (Fig. 1).
- Masson area up with 0.5 and 5 µm; α-SMA protein up but mRNA unchanged (Fig. 2).
- Inflammation up, depending on size (Fig. 3).
- Oxidative stress not changed: "levels of GPX, MDA, and SOD were not increased" (Suppl. Fig. 1, n = 7).
- Lipid deposition up, largest at 80 nm (Fig. 4).
- DEGs: 333, 602 and 758 for 80 nm, 0.5 µm and 5 µm.
- **Lipid genes by qPCR (Fig. 6, text)**:
  - PPARα/FAO up: Acox1, Ppara, Cpt1a, Mcad up (80 nm, 0.5 µm); Ppara up (5 µm); Cpt1b up (80 nm).
  - Synthesis up: Fas and Apoc4 (0.5 µm), Pparg (0.5, 5 µm), Scd1 (80 nm). Acc1 is in the primer table but no change is described in the text.
  - Transport up: Fatp1 (80 nm, 0.5 µm), Fatp2 (0.5 µm).
- Hub genes Acot3, Abcc3, Nr1i3 up.
- Cholesterol genes: not reported.

### Q10. Limitations
**(a) Actual**
1. limma usage on RNA-seq is unspecified, and the alignment and quantification pipeline is not described in the paper (GEO says GRCm38.p6).
2. Hub genes come from network topology on a sparse PPI, so the hub status is a prediction. Fmo2 failed qPCR.
3. LSD post hoc, no multiplicity control.
4. Attrition from 10 to 5–8 is unexplained.
5. Copy errors: oxidative stress measured "in ovaries" in male mice; GPX vs GSH; the Discussion says 0.5 and 5 µm "resulted in a reduction in … liver fibrosis"; "four RNA-seq datasets" for three contrasts.
6. The age discrepancy (paper 5 weeks vs GEO 6 weeks).
7. No particle characterisation, and the dose actually ingested is unknown.
8. Liver accumulation of particles was not measured.

**(b) Not attempted**
- Females; metabolic diet; dose–response.
- Functional tests of hub genes.
- Protein for lipid genes; liver TG/cholesterol biochemistry.
- Time course; particle tracking.

---

## 4. Chiu HW, Chu CW, Huang CC, Chia ZC, Wang YL, Lee YH. 2025. "Polystyrene microplastics induce hepatic lipid metabolism and energy disorder by upregulating the NR4A1-AMPK signaling pathway." *Environ Pollut* 369:125850

DOI: https://doi.org/10.1016/j.envpol.2025.125850 (PII S0269749125002234; Crossref, published March 2025) · PMID 39956511 · no PMC copy · **ACCESS: ABSTRACT ONLY** (plus a GEO record linked to the paper). Not open access (Unpaywall: closed). ScienceDirect returned 403 to automated fetch and then a CAPTCHA, which I did not attempt. Author order per Crossref and PubMed: Chiu HW first; a library catalogue lists Chia ZC first, which is unverified.

### Exposure record
| Field | Value | Source |
|---|---|---|
| Species/strain/sex/age | C57BL/6, male, 6 weeks | **GEO GSE277772** sample characteristics (submitted by H-W Chiu, Y-H Lee et al., National Cheng Kung University; PubMed ID 39956511). MeSH also lists "Mice, Inbred C57BL" and "Male" |
| Polymer / size | Polystyrene. The abstract mentions "different sizes (0.08-100 μm)" only in background citing prior research. **The particle size used is UNVERIFIED** |
| Dose / route / frequency | "orally dosed with 0.4 mg/day PS-MPs twice weekly for eight weeks"; controls "orally with water twice weekly" | GEO. Gavage vs drinking water: UNVERIFIED ("orally dosed") |
| Duration | 8 weeks | GEO |
| Diet, housing, sampling time, fasting, vehicle | UNVERIFIED | — |

### Q1. Question (abstract)
Does NR4A1 mediate PS-MP disruption of hepatic lipid and energy metabolism?

### Q2. Design and n
- Sequencing n = **3 PS-MP vs 3 water control** (GEO GSE277772, 6 GSMs, NovaSeq 6000).
- Phenotype n, number of groups (one or several sizes or doses) and in vitro replicates: **UNVERIFIED**.
- In vitro: "NR4A1 knockdown in hepatocytes" (abstract); cell line and method UNVERIFIED.

### Q3. Datasets
Generated: GSE277772 (raw counts matrix and FPKM in `GSE277772_RAW.tar`; GPL24247). GEO processing: Trimmomatic → HISAT2 → featureCounts v1.6.0 → DESeq2 v1.22.1 (RLE). The GEO record states "Assembly: GRCh38" and "e.g., H. sapiens", which is probably template text, but the mouse reference actually used is unverified. Re-used datasets: UNVERIFIED.

### Q4. Models
- Abstract: "Ingenuity Pathway Analysis (IPA) software and Gene Set Enrichment Analysis (GSEA)".
- GEO: DESeq2 (and DEGseq for any comparisons without biological replicates).
- DE thresholds, multiple-testing method and statistical tests for phenotypes: **UNVERIFIED**.

### Q5. Fit (from what is verifiable)
- n = 3 vs 3 gives low power for transcriptome-wide DE.
- IPA upstream-regulator and GSEA analyses are hypothesis-generating.
- Selecting NR4A1 as the driver from such data needs functional follow-up, which the abstract says was done in vitro by knockdown.

### Q6. Inputs / outputs
Liver RNA-seq → IPA/GSEA → NR4A1. Phenotypes: liver TG (down), ROS (up). Proteins: AMPK, autophagy and lipid-metabolism proteins (abstract).

### Q7. Evaluation
UNVERIFIED.

### Q8. Validation (abstract)
- In vitro NR4A1 knockdown under PS-MP "reduced the expression of AMPK and lipid metabolism-related proteins", i.e. protein-level and loss-of-function evidence in vitro.
- In vivo knockout: none mentioned.
- External dataset: none mentioned.
- Whether NR4A1 was confirmed by qPCR or western blot in liver: UNVERIFIED.

### Q9. Main findings (abstract)
- "PS-MPs reduced triglyceride levels and significantly increased reactive oxygen species (ROS) in liver tissue".
- "significant upregulation in NR4A1 gene expression".
- "PS-MP accumulation in the liver activated NR4A1 and the AMPK-autophagy pathway, reducing lipid biosynthesis". So the direction for lipogenesis is DOWN, and AMPK and autophagy are UP.
- PPARα and cholesterol genes: UNVERIFIED.
- The magnitude of the Nr4a1 change is UNVERIFIED. It could be checked in the GSE277772 counts, which was not done here.

### Q10. Limitations
**(a) Actual (verifiable)**
1. n = 3 per group for RNA-seq.
2. Male only.
3. Ambiguous reference assembly in GEO.
4. "Accumulation in the liver" is asserted, but the detection method is unverified.

**(b) Not attempted (per abstract)**: no in vivo NR4A1 loss of function; no metabolic diet mentioned.

Everything else is UNVERIFIED pending the full text.

---

## 5. Xu R, Cao JW, Lv HL, Geng Y, Guo MY. 2024. "Polyethylene microplastics induced gut microbiota dysbiosis leading to liver injury via the TLR2/NF-κB/NLRP3 pathway in mice." *Sci Total Environ* 917:170518

DOI: https://doi.org/10.1016/j.scitotenv.2024.170518 (PII S0048969724006557) · PMID 38286276 · no PMC copy · **ACCESS: ABSTRACT ONLY** (Unpaywall: closed; ScienceDirect 403; no GEO/SRA link found).

### Exposure record
- Mice (title). MeSH lists no sex or strain tag. Strain, sex, age: **UNVERIFIED**.
- Polymer: polyethylene (title, MeSH).
- Size, dose values, route, vehicle, duration, diet, housing, sampling: **UNVERIFIED**. The abstract implies **two concentrations** plus a control ("both concentrations").

### Q1. Question
"investigate the potential of MPs to cause liver injury by disrupting the balance of gut microbiota" (abstract).

### Q2. Design and n
Two PE concentrations plus a control (inferred from the abstract). n: UNVERIFIED.

### Q3. Datasets
None stated in the abstract. The microbiota profiling platform (presumably 16S) is UNVERIFIED. No transcriptome is mentioned.

### Q4. Models
UNVERIFIED. The abstract mentions a correlation of "fibril-related proteins" with concentration.

### Q5. Fit
The causal claim (dysbiosis leads to liver injury via TLR2) needs an intervention such as antibiotic depletion, faecal microbiota transfer, germ-free mice or TLR2 inhibition/knockout. **Whether any was performed is UNVERIFIED.** From the abstract alone, the link is associative.

### Q6. Inputs / outputs
Liver organ coefficient, liver function, fibrosis proteins, microbiota composition, oxidative stress, inflammation, apoptosis/necrosis, TLR2/NF-κB/NLRP3 (abstract).

### Q7. Evaluation
UNVERIFIED.

### Q8. Validation
UNVERIFIED (assays for TLR2/NF-κB/NLRP3 are presumably protein- or mRNA-level but not stated).

### Q9. Findings (abstract)
- Liver organ coefficient DOWN.
- Fibril-related proteins UP, concentration-dependent.
- Desulfovibrio, Clostridia, Enterorhabdus, Bacteroides and Gemella UP; Dubosiella DOWN (spelled "Dubosoella").
- Oxidative stress, inflammation, apoptosis and necrosis UP.
- TLR2/NF-κB/NLRP3 activated.
- No lipid, PPARα or cholesterol findings mentioned.

### Q10. Limitations
- (a) Not assessable without the full text.
- (b) Per the abstract: no transcriptomics, no diet factor.

---

## Cross-paper notes relevant to the team (verified)
- **Lipogenic genes**: DOWN with PE in Jung (male, gavage, 10–150 µm, purified diets: Fasn, Acly, Acaca, Me1, Pnpla3 in both diets) vs UP in Djouina (female, dietary 36/116 µm, chow: Fasn, Acly, Acaca at 6 weeks). PS studies: Li reports Fas and Scd1 up by qPCR; Chiu's abstract reports lipid biosynthesis reduced. Directions are inconsistent across studies that differ in sex, route, size, diet and polymer.
- **Ppara**: up in all three studies that report it (Jung, which is significant vs CD+Veh but ns within MD; Djouina; Li qPCR). Jung's own reprocessing of GSE245069 also shows Ppara up.
- **Only Jung provides a reanalysable 2×2 design.** Its supplement already contains a DESeq2 interaction table (S7), and the table looks contaminated by muscle genes. Worth checking against our own reanalysis of GSE294180.

## URLs / sources actually read
- https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13274612/fullTextXML (Jung full text)
- https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13274612/supplementaryFiles (→ sciadv.aec8681_sm.pdf, sciadv.aec8681_tables_s1_to_s18.zip/.xlsx)
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294180 (SOFT brief) and GPL30172
- https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11747203/fullTextXML (Li full text)
- https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11747203/supplementaryFiles (→ 41598_2025_86810_MOESM1_ESM.pdf)
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245069, GSM7835567, GPL28330 (SOFT)
- Local file C:\Users\Soham\Downloads\1-s2.0-S0147651323009211-main.pdf (Djouina publisher PDF, text extracted)
- https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:10.1016/j.envpol.2025.125850 (Chiu abstract, MeSH)
- https://api.crossref.org/works/10.1016/j.envpol.2025.125850 ; https://api.unpaywall.org/v2/10.1016/j.envpol.2025.125850
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277772 (series + GSM SOFT; linked to PMID 39956511)
- https://www.sciencedirect.com/science/article/pii/S0269749125002234 (403 / CAPTCHA, not attempted)
- https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:10.1016/j.scitotenv.2024.170518 (Xu abstract, MeSH)
- https://api.unpaywall.org/v2/10.1016/j.scitotenv.2024.170518 ; https://doi.org/10.1016/j.scitotenv.2024.170518 → linkinghub → ScienceDirect (403)
- https://aminer.org/pub/660cf4d113fb2c6cf61f5fe6/... (header only, no content)
- NCBI E-utilities GEO search (polystyrene AND liver)
