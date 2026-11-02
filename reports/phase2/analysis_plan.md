# Phase 2 analysis plan (fixed 2026-11-02, before pathway, concordance or robustness results were seen)

This plan records the analytic choices made **before** looking at the corresponding results, so
later readers can tell planned analyses from data-driven ones. Any departure will be listed in
the relevant report as a deviation.

The plan was written after the provisional DESeq2 contrast summary was seen (numbers of
significant genes per contrast) but before any gene-level, pathway-level or concordance output.

## Count input

- **Primary input:** decided by the re-quantification pilot (`reports/quantification_validation/pilot_protocol.md`).
- **Until then:** every result uses the deposited matrix without an offset and is labelled **PROVISIONAL**.

## Model (Direction B core)

`~ diet + microplastic + diet:microplastic` with reference levels CD and Veh (checked in code). Contrasts:

| Name | Definition |
|---|---|
| PE_in_CD | `microplastic_PE_vs_Veh` |
| PE_in_MD | `microplastic_PE_vs_Veh + dietMD.microplasticPE` |
| Diet_in_Veh | `diet_MD_vs_CD` |
| Interaction | `dietMD.microplasticPE` |

- Each contrast is its own BH-FDR family (DESeq2 `results()`, α = 0.05, independent filtering on).
- Shrunken LFCs (ashr) are reported alongside MLEs, and are used for ranking and plotting only, not for testing.

## Pathway analysis

**Gene sets.** MSigDB 2026.1.Mm via msigdbr 26.1.1, matched by mouse gene symbol.

- **Primary family:** all 50 MH Hallmark sets. BH-FDR across the 50 sets, separately for each contrast.
- **Focus themes (subset of the primary family, reported first):**

  | Theme | Sets |
  |---|---|
  | Fatty-acid metabolism | HALLMARK_FATTY_ACID_METABOLISM |
  | Cholesterol | HALLMARK_CHOLESTEROL_HOMEOSTASIS |
  | Peroxisome / bile acids | HALLMARK_PEROXISOME, HALLMARK_BILE_ACID_METABOLISM |
  | Lipid storage | HALLMARK_ADIPOGENESIS |
  | Inflammation | HALLMARK_INFLAMMATORY_RESPONSE, HALLMARK_TNFA_SIGNALING_VIA_NFKB, HALLMARK_IL6_JAK_STAT3_SIGNALING, HALLMARK_INTERFERON_GAMMA_RESPONSE |
  | Xenobiotic | HALLMARK_XENOBIOTIC_METABOLISM |
  | Oxidative / cellular stress | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY, HALLMARK_UNFOLDED_PROTEIN_RESPONSE, HALLMARK_P53_PATHWAY, HALLMARK_HYPOXIA |
  | Mitochondrial | HALLMARK_OXIDATIVE_PHOSPHORYLATION |
  | Fibrosis / remodelling | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, HALLMARK_TGF_BETA_SIGNALING |
  | Contamination flag | HALLMARK_MYOGENESIS (expected to track muscle carry-over, not treatment) |

- **PPAR (secondary family, not in Hallmark):** REACTOME_PPARA_ACTIVATES_GENE_EXPRESSION and WP_PPAR_SIGNALING_PATHWAY, with their own BH family.
- **Exploratory only:** Reactome pathways (M2 CP:REACTOME). Reported as supporting detail, never as headline findings.

**Method 1: pathway scores (primary, unit of analysis = mouse).**

1. For each set, take the genes present after filtering and z-score their VST expression across the 19 mice.
2. The score is the mean z per mouse.
3. Fit an ordinary linear model `score ~ diet * microplastic` (classical two-way factorial regression, 15 residual df).
4. Report each contrast's estimate, SE, 95% CI and *t*-test p-value, with BH within family per contrast.

Score caveats:

- The mean z assumes genes in a set move in the same direction. Mixed-direction sets dilute.
- Hallmark sets are curated mainly from up-regulated responses.

**Method 2: preranked GSEA (fgsea, supporting).**

- Ranking statistic: DESeq2 Wald statistic for each contrast.
- Settings: `fgseaMultilevel`, `minSize = 15`, `maxSize = 500`, `eps = 0`.
- Reported as NES and padj. The **gene** is the sampling unit here, and gene–gene correlation is ignored, so its p-values are treated as less reliable than Method 1's.

**Categorising the PE response across diets** (Hallmark, using Method 1 estimates and 95% CIs). The rules use the CIs of the PE effect in CD (β_CD), the PE effect in MD (β_MD) and the interaction (β_int):

| Category | Rule |
|---|---|
| Similar in both diets | β_CD and β_MD both exclude 0, same sign; β_int CI includes 0 |
| Amplified under MD | Same sign; β_int CI excludes 0 in the direction of the effect |
| Attenuated under MD | Same sign; β_int CI excludes 0 against the direction of the effect |
| Opposing | β_CD and β_MD exclude 0 with opposite signs |
| CD-only / MD-only, interaction supported | One PE CI excludes 0; β_int CI excludes 0 |
| CD-only / MD-only, interaction not supported | One PE CI excludes 0; β_int CI includes 0. **Not** evidence of diet specificity. |
| No clear PE effect | Neither PE CI excludes 0 (not evidence of no effect) |

CIs are unadjusted 95% intervals; BH-adjusted interaction p-values are reported next to them.

## Direction A: PE vs MASH-diet similarity

**Shared-reference problem.** PE_in_CD and Diet_in_Veh both subtract the CD+Veh mean, so their estimation errors are positively correlated. With equal per-group variance and n = 4/5/5, the null correlation of the two estimates for a gene is 0.25 / 0.45 ≈ 0.56 of the noise component. Naive concordance is therefore biased upward.

**Primary approach: split reference (sample-disjoint).**

1. Split the 4 CD+Veh mice into two pairs. Use one pair as the reference for the PE contrast (vs CD+PE) and the other as the reference for the diet contrast (vs MD+Veh).
2. All 6 ordered assignments of the 3 possible pair partitions are run. The two contrasts then share **no** mice, so their errors are independent.
3. Concordance statistics are computed for each split and summarised by their range.
4. Each contrast is estimated with DESeq2 on its 7 or 7 mice. The cost is precision (reference n = 2).

**Statistics:**

- Spearman correlation of MLE log2FCs over genes with baseMean > 50;
- sign agreement among the top 500 genes ranked by |diet Wald statistic| in the diet half;
- Hallmark NES correlation (fgsea on each half).

**Calibration:** the same statistics computed when the "PE" group is replaced by an equally sized random subset of the *other* CD+Veh … not possible (n = 4). Instead the null is approximated analytically, and the naive (shared-reference) value is shown next to the split values to quantify the bias.

**Secondary:** the same comparison with PE_in_MD against Diet_in_Veh. These share no mice, but MD+Veh appears in PE_in_MD, and MD+Veh vs CD+Veh is the diet contrast, so MD+Veh is shared. Split MD+Veh (n = 5) into 2 + 3 in the same way.

**External signature.** Candidate independent MASH datasets are listed with their comparability. They are implemented only if time allows, and labelled accordingly.

## Direction B: interaction

- Gene level: the interaction coefficient with SE, CI and padj (from the model above).
- Pathway level: Method 1 interaction estimates, plus fgsea on the interaction Wald statistic.

## Direction C: robustness (all against the primary model)

| Check | Implementation |
|---|---|
| Drop MD_1 | Refit without MD_1 |
| Leave-one-mouse-out | 19 refits; top-gene and pathway-category stability |
| Muscle covariate | Standardised log2 muscle-marker CPM added to the design |
| Depth covariate | Standardised log2 read pairs added (collinear with PE; interpreted as a stress test, not a correction) |
| Quality-proxy covariate | Standardised gene-counts-per-read-pair proxy added |
| Alternative method | edgeR quasi-likelihood and limma-voom with the same design |
| Gene-class checks | Fraction of headline genes that are muscle-tracking, mitochondrial, circadian/clock-controlled or among the 50 most abundant |

**Headline findings** (fixed before the checks) are the pathway categories of the focus themes and the sign/magnitude of the PE–diet concordance. The robustness matrix reports which survive which check.

## Deviations and additions (recorded as they happened)

| # | When | What changed | Why | Effect on interpretation |
|---|---|---|---|---|
| D1 | After the diet pathway results were seen, before the PE pathway results were examined | Added limma **fry** (self-contained, mouse-level) and **camera** (competitive) as pathway methods. PE-response categories now use fry; mean-z categories are still reported. | Mean-z scores failed a positive control: no strong diet effect passed FDR because Hallmark sets mix up- and down-regulated genes. | fry is more appropriate for mixed sets and more conservative than fgsea |
| D2 | While building the gene-set lists | REACTOME_PPARA_ACTIVATES_GENE_EXPRESSION is called REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA in MSigDB 2026.1.Mm; SANDERSON_PPARA_TARGETS (mouse liver) added to the PPAR family | Release naming; a more liver-specific PPARα target set | The PPAR family has 3 sets instead of 2 |
| D3 | After the primary interaction genes were seen | Added a **clock-phase covariate** variant and a per-mouse clock score | Clock genes shifted in PE_in_CD and Diet_in_Veh in the same direction | New confounder identified (CD+Veh circadian/feeding state). Reported as a sensitivity analysis, not as a correction. |
| D4 | Leave-one-out concordance | The relabelling null was not recomputed for each of the 19 removals (cost ~25 min each); the observed split-reference ρ is reported per removal | Compute time | The LOO check covers stability of the estimate, not of its p-value |
| D5 | Bug fix | The pathway script originally ignored the covariate of covariate variants (design `~ diet * PE` was hard-coded). Fixed to use the DESeq2 fit's design; all covariate variants re-run. | Correctness | Earlier covariate-variant pathway numbers were wrong and are superseded |

The **primary model was not changed** (no covariate). Whether to make the muscle-covariate model primary is left to the team (see `interpretation.md`, decision 2).

## Team decisions

| Date | Decision | Consequence |
|---|---|---|
| 2026-11-02 | **Primary model stays as pre-specified** (`~ diet + microplastic + diet:microplastic`, no covariate). The muscle-covariate model remains a pre-specified sensitivity analysis. | Headline results are reported from the primary model. Wherever the muscle-covariate fit changes a conclusion (PE_in_MD gene counts, the gene-level interaction list), both are shown and the dependence on contamination is stated. |
