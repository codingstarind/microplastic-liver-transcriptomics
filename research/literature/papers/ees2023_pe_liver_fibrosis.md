# ees2023_pe_liver_fibrosis (extracted by coordinator from the user-supplied publisher PDF)

## Citation
Djouina M, Waxin C, Dubuquoy L, Launay D, Vignal C, Body-Malapel M. "Oral exposure to polyethylene microplastics induces inflammatory and metabolic changes and promotes fibrosis in mouse liver." Ecotoxicology and Environmental Safety 2023; 264:115417. DOI 10.1016/j.ecoenv.2023.115417. PMID 37651791. No PMCID (not in Europe PMC full text). Open access, CC BY-NC-ND. Received 20 June 2023; accepted 26 August 2023. Univ. Lille, Inserm, CHU Lille, U1286-INFINITE, France.

## Access level read
Full text of the publisher PDF (`1-s2.0-S0147651323009211-main.pdf`, downloaded by the user 2026-11-20), read as extracted text: Abstract, Methods, Results, all figure legends, Discussion, Conclusion, Data availability. **Supplementary Fig. 1 (body weight, serum ALT/AST at 6 wk) NOT read** — not included in the PDF text.

## Biological question
Hepatic effects of sub-chronic dietary exposure to polyethylene (PE) microplastics in mice: "the purpose of this study was to explore the hepatic effects of sub-chronic PE MP exposure in mice" (Introduction, end). Second question: does dietary PE worsen fibrosis in a chemically induced fibrosis model?

## Experimental design (groups, controls, randomisation/blinding if stated)
**Experiment 1 (healthy mice).** "mice were randomly assigned to 4 experimental groups (n = 10/group): 1) control group, receiving control food; 2) RB group... 3) GB group... 4) RB+GB group" (Methods 2.2). Treatment periods "either 6 or 9 weeks" (so 6-wk and 9-wk cohorts). Control = uncontaminated food. **No vehicle/inert-particle control** (no non-plastic bead of matched size).
**Experiment 2 (fibrosis susceptibility).** CCl4 model, "Mice were randomly divided into 5 groups (n = 8/group)": control (control food + olive oil), CCl4, CCl4+RB, CCl4+GB, CCl4+RB+GB. Sacrificed 24 h after last CCl4 injection.
- Randomisation: stated ("randomly assigned"/"randomly divided"); method not described. Blinding: **not reported** (quantification was software-assisted in QuPath, which mitigates but does not replace blinding).
- **Note: this is not a 2 × 2 metabolic-diet design.** The second factor is a hepatotoxin (CCl4), not diet.
- Cage/housing density: "housed under standard conditions" only; cage not reported as a unit.

## Species/strain/sex/age; tissue(s)
- "C57BL/6 mice were purchased from Janvier Labs... **Female** mice from 7 to 12 weeks old were used in the study." (Methods 2.2). **Female only** — opposite sex to Jung 2026 (male). Age range 7–12 wk is wide and not reported per group.
- Exp. 2: "8-week old C57BL/6 female mice".
- Tissues: liver (histology, IHC, TUNEL, qPCR, MDA, lipids), serum (ALT/AST).

## Biological replicates per group (per assay)
- Exp. 1: n = 10/group per time point (Methods and every figure legend: "n = 10/group... or control food (CT; n = 10)").
- Exp. 2: n = 8/group (Methods 2.2; Fig. 7 legend).
- Replicate = one mouse (liver sampled at necropsy per mouse). No pooling stated. Attrition/exclusions: none reported.

## Polymer; particle size; supplier; label; characterisation; endotoxin
- Polymer: **polyethylene**, two sizes used singly and in combination.
- "PE microbeads were acquired from Cospheric (Santa Barbara, USA). Two categories... red fluorescent beads (RB; item# UVPMS-BR-1.090; 10–45 µm) and green fluorescent beads (GB; item# UVPMS-BG-1.00; 106–125 µm)." "The median size of microbeads was 36 and 116 µm for RB and GB, respectively. Both were composed of a pure PE polymer as confirmed by Fourier transform infrared spectra." (Methods 2.1)
- Characterisation: size range and polymer composition characterised **in a previous study** (Djouina et al. 2022a) — FTIR-confirmed; no DLS/zeta (not applicable at this size).
- **Fluorescently labelled beads** (red/green dye) — note dye leaching is not addressed.
- Endotoxin: **not reported**.
- Beads provided by ANSES Laboratoire de Sécurité des Aliments (Acknowledgments).

## Dose; route; vehicle; frequency; duration
- "food spiked with 100 µg RB microbeads/g of food" (and the same for GB; RB+GB group received 100 µg **each**/g of food — i.e. 200 µg/g total). Rationale: "this amount is within the realistic range of human ingestion levels" (Methods 2.2).
- Route: **dietary** (spiked chow), ad libitum. No gavage, no vehicle/surfactant.
- Duration: 6 or 9 weeks (Exp. 1); 8 weeks (Exp. 2).
- **Food intake / per-mouse ingested dose (mg/kg/day) not reported.**

## Diet / metabolic-stress model
**No metabolic-stress diet.** Standard food spiked with beads; chow composition/supplier not reported. The second-hit model is **CCl4-induced fibrosis** (0.2 ml/kg, in olive oil 1:4, i.p.), not a high-fat/MASH diet.
- **Internal inconsistency to flag:** Methods 2.2 says CCl4 "i.p. injection **twice a week** for 8 weeks"; the Fig. 7 legend says "by i.p. injection **once a week** for 8 weeks". Unresolved in the text.

## Data types
- **RT-qPCR only — there is NO RNA-seq / transcriptome in this study.** Targeted panels: fatty-acid uptake, β-oxidation, FA synthesis, TG synthesis, xenobiotic enzymes, oxidative-stress genes, cytokines/chemokines, cell-type markers, fibrotic markers. Normalised to *Polr2a*, 2^-ΔΔCt. "The primer sequences were designed using Primer Express 3... and are available upon request." (Methods 2.7) — **primer sequences not published**.
- Histology: H&E, sirius red (QuPath quantification, whole-slide, Zeiss AxioScan.Z1, 20×).
- IHC: Ki-67, α-SMA. TUNEL (In Situ Cell Death Kit, TMR red).
- Biochemistry: liver cholesterol, triglyceride, free fatty acids; liver MDA (normalised to protein); serum ALT/AST.

## Statistics & computational methods
- "Results are expressed as mean ± standard error of the mean. The statistical significance of differences between experimental groups was calculated using the **Mann-Whitney nonparametric U test** (GraphPad Prism). Statistical significance was defined as p < 0.05." (Methods 2.9)
- **No multiple-testing correction** anywhere, despite dozens of qPCR targets × 3 exposure groups × 2 time points. All comparisons are **pairwise vs control** (Exp. 1) or vs control and vs CCl4 (Exp. 2, "# p < 0.05... compared to the CCl4 group").
- **No interaction term is fitted** for PE × CCl4; "exacerbation" rests on PE-exposed CCl4 groups differing from the CCl4-alone group by Mann-Whitney.
- No DE tool / pathway analysis / multivariate model (no omics).

## Main findings with effect directions and statistical evidence
- **No mortality, no clinical signs, no body-weight change** at 6 or 9 wk (Results opening; Suppl. Fig. 1A — not read). VALIDATED-null.
- **Serum AST and ALT increased at 9 wk in the RB (36 µm) group** (Suppl. Fig. 1B — not read; stated in Results text). VALIDATED (serum enzymes).
- **Lipid metabolism (6 wk, Fig. 1):** *Cd36* and *Fabp1* up in RB+GB; *Fabp1* **down** in RB; *Slc27a2/Fatp2* unchanged. "**A two-fold increase in expression of the peroxisome proliferator-activated receptor alpha (Ppara) gene was detected after exposure to both types of PE microbeads**", while *Acox1*, *Cpt1a*, *Acadm* unchanged. *Pparg* and *Scd1* up in RB+GB; "**Fasn, Acly and Acaca were significantly upregulated in all groups of PE-exposed mice**". TG-synthesis genes unchanged except *Gpam* up in GB. VALIDATED (qPCR).
- **Lipid metabolism (9 wk):** RB down-regulated *Fabp1*, *Acox1*; up *Dgat1*. GB increased *Fatp2*, *Cpt1a*, *Gpam*. *Ppara* and *Pparg* higher in GB and RB+GB. *Cd36*, *Acaca*, *Mttp* up in **all** PE groups. Summary (authors): "disturbed genes associated with fatty acid metabolism **in favor of increased mRNA levels** of genes associated with fatty acid uptake, β-oxidation, and synthesis." VALIDATED (qPCR).
- **Liver lipids:** cholesterol **up** in GB; triglycerides **up** in RB+GB; free fatty acids unchanged (Fig. 1E–G). VALIDATED (biochemistry).
- **Xenobiotic metabolism (Fig. 2A):** *Cyp2a4*, *Cyp2b10* up in all three PE groups at 6 wk; *Cyp2a4* also up at 9 wk; *Gsta1* up at 6 wk; *Ugt1a1* up in RB+GB at 9 wk. VALIDATED (qPCR).
- **Oxidative stress (Fig. 2B–C):** *Nos2*, *Cat*, *Nfe2l2*, *Nqo1*, *Hmox1* up in RB+GB at 6 wk; *Cat* **down** in RB. *Sod1*, *Sod2*, *Gpx1* unchanged. At 9 wk: *Nfe2l2* up (GB), *Nqo1* up (RB, RB+GB), *Hmox1* up in all three. **MDA decreased** in RB at 6 wk but **increased** in RB+GB at 6 wk and in RB at 9 wk — explicitly time-dependent and direction-inconsistent. VALIDATED (qPCR + MDA assay).
- **Inflammation (Fig. 3):** no necrosis/ballooning/edema on H&E. Inflammatory foci **size** up in RB+GB at 6 wk; foci **number** up in GB and RB+GB at 9 wk. Cytokines at 6 wk: *Tnf*, *Ifng*, *Il1b*, *Il6*, *Ccl2*, *Ccl5* up in RB+GB; but *Il1b* **down** in GB and *Il6* **down** in RB. At 9 wk: *Il1b* **down** in RB and GB; *Il6* strongly up in RB+GB; *Ccl2* up in GB. Authors' reading: "co-exposure to RB and GB PE microbeads induced **low grade** inflammation". VALIDATED (histology + qPCR); **direction is group-dependent, not uniform**.
- **Proliferation (Fig. 4):** *Mki67* mRNA up (RB+GB at 6 wk; GB and RB+GB at 9 wk); Ki-67 IHC confirmed in RB+GB at both times. **TUNEL unchanged** (no apoptosis). VALIDATED (qPCR + IHC + TUNEL).
- **Cell-type markers (Fig. 5):** *Alb*, *Hnf4a* up at 9 wk in RB+GB; cholangiocyte markers *Hnf1b*/*Sox9* up in RB and down in GB at 6 wk, both down in RB+GB at 9 wk; HSC markers *Acta2*, *Pdgfa* up at 6 wk in RB+GB, *Pdgfa* up in all groups at 9 wk; portal-fibroblast markers unchanged; Kupffer markers *Cd68* down / *Clec4f* up in GB at 9 wk. VALIDATED (qPCR). **Note: bulk-tissue markers cannot separate per-cell expression from cell-composition change.**
- **Fibrosis in healthy mice (Fig. 6):** sirius red area up in GB at 6 wk; in GB and RB+GB at 9 wk. VALIDATED (histology).
- **Fibrosis susceptibility (Fig. 7):** all CCl4 groups raised ALT/AST vs control; "elevated AST and ALT activity levels were significantly **worsened** in mice fed RB+GB". Sirius red and α-SMA positive areas "significantly increased in mice exposed to GB and RB+GB compared to fibrotic mice exposed to uncontaminated food". *Acta2*, *Pdgfa*, *Col1a1*, *Vim* up in all CCl4 groups; RB+GB "significantly enhanced" all four; GB enhanced *Vim*. Conclusion: "dietary PE exposure exacerbated CCl4-induced liver fibrosis." VALIDATED (histology, IHC, qPCR) — but **by pairwise contrast vs the CCl4 group, not by an interaction test**.

## Pathway-theme directions
- **Lipid / lipogenesis: UP.** *Fasn*, *Acly*, *Acaca* up in all PE groups at 6 wk; *Acaca* up at both time points; *Scd1*, *Pparg* up; liver TG and cholesterol up. **This is the opposite direction to Jung 2026's PE-associated down-regulation of Fasn/Acly/Me1 in male mice on CD/MASH diets** — note the many differing conditions before treating it as a contradiction (sex, route, dose, bead size, diet).
- **PPARα / FAO: UP at the receptor, mixed downstream.** *Ppara* ~2-fold up at 6 wk (both bead types) and up at 9 wk (GB, RB+GB) — **same direction as Jung 2026's reported Ppara induction**. But downstream *Acox1*, *Cpt1a*, *Acadm* unchanged at 6 wk, and *Acox1* **down** in RB at 9 wk, so receptor induction is not accompanied by consistent FAO-target induction.
- **Cholesterol:** liver cholesterol up (GB). No cholesterol-synthesis gene panel.
- **Inflammation: UP overall but group-dependent** (several *Il1b*/*Il6* down-regulations in single-bead groups).
- **Oxidative stress: altered, direction time-dependent** (Nrf2 battery up; MDA down at 6 wk in RB, up at 9 wk).
- **Mitochondria:** not assayed.
- **Fibrosis: UP**, weak in healthy liver, clear exacerbation on a CCl4 background.
- **Circadian:** not reported; sacrifice time of day not stated.
- **Interaction:** **not formally tested.** PE × CCl4 "exacerbation" is pairwise.

## Author-acknowledged limitations (quotes + location)
- Particle fate: "**it seems unlikely that PE MPs reached the liver in our experimental conditions. Therefore, the effects are more likely to be mediated indirectly, in particular through the effects of MPs on dysbiosis and intestinal immunity**" (Discussion). Also: "our data do not allow to determine whether the administered MPs were degraded by digestion".
- Cross-study comparability: "It is not possible at this time to determine whether these dissimilarities on fibrosis outcomes were due to differences in experimental design or chemical differences in the polymers." (Discussion)
- Cell-specific immune effects: "The MP induced cell-specific immune modulation remains to be elucidated." (Discussion)
- Human relevance: "More research is necessary to determine the human relevance of these results." (Discussion)
- Sex-specificity is raised for *another* study's results (Yang et al. 2022), not for their own female-only design.

## Methodological points to verify independently
- **Sex:** female only, where Jung 2026 and most other MP-liver studies use males. The authors themselves note *Cyp2a4* is "strongly female-biased (female/male ratio of ~1000)" — several xenobiotic findings may be sex-specific.
- **No multiple-testing correction** across a large qPCR panel; many single-group significances should be read cautiously.
- **CCl4 frequency contradiction** (Methods "twice a week" vs Fig. 7 legend "once a week").
- RB+GB group received **double the total particle mass** (100 µg/g of each); several effects are strongest in RB+GB, so dose and bead-type combination are confounded.
- Age range 7–12 weeks at start — not reported per group; could add variance.
- Fluorescent dye in both bead types; no unlabelled-bead control for dye effects.
- Food intake and therefore delivered dose per mouse not measured.
- Cage/housing not reported as a unit of analysis; sacrifice time of day and fasting state not reported.
- Primer sequences "available upon request" — not published, so the qPCR panel cannot be independently reproduced as published.
- Suppl. Fig. 1 holds the body-weight and 6-wk ALT/AST data: not read here.

## Data availability statement (verbatim) and accessions
- Verbatim: "**Data will be made available on request.**"
- **No accession of any kind** (no GEO/SRA/BioProject/Zenodo/figshare/GitHub). Consistent with there being no omics data in the study.
- Supporting information: "Supplementary data associated with this article can be found in the online version at doi:10.1016/j.ecoenv.2023.115417."

## Relevance to the team
The only other **polyethylene** oral-exposure mouse liver study in this set, and the closest independent comparison to Jung 2026 for RQ1. It agrees with Jung on *Ppara* induction but reports **lipogenic genes (Fasn, Acly, Acaca) up** where Jung reports them down, under materially different conditions (female mice, dietary beads at 100 µg/g food, 36/116 µm, no MASH diet). Its second-hit model is CCl4, not diet, so it speaks to RQ2 only by analogy: a pre-injured liver showed worse fibrosis with PE, tested pairwise rather than as an interaction. Its authors argue PE effects are likely **indirect** (gut-mediated) because particles of this size probably never reached the liver — a mechanistic alternative the team should keep in view when interpreting liver-only bulk data.
