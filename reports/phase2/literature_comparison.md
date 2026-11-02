# Comparison with the original publication: reproduction vs extension

**Paper:** Jung W. *et al.* (2026) *Sci. Adv.* 12(25): eaec8681 (PMC13274612). Full text read via Europe PMC.

**Supplementary tables S1–S18:** not yet obtained. PMC serves a proof-of-work bot challenge to scripted requests, and a manual browser download has been requested. Until they arrive, the comparison uses the genes and numbers stated in the main text and figure legends, and **no supplementary content is described or used**.

> Agreement with the paper is **reproducibility evidence** (same mice, same reads), not
> independent biological validation.

## Methods differences that can change results

| Aspect | Paper | Ours |
|---|---|---|
| Count input | tximport → DESeq2; setting not reported | Deposited integer matrix, no length offset (**provisional**); re-quantification pilot under way |
| Model | Pairwise condition comparisons (4-level factor implied) | 2 × 2 factorial `~ diet * PE`, explicit references, interaction coefficient |
| Thresholds | \|log2FC\| ≥ 0.38 and FDR < 0.05 | FDR < 0.05 per contrast; the 0.38 cut reported for comparison only |
| Pre-filter | not reported | ≥ 10 counts in ≥ 4 samples |
| Contamination handling | not reported | Diagnosed; muscle-covariate and MD_1 sensitivity |
| Pathways | fgsea with a custom "LIPID_METABOLISM_AND_TOXICITY" NR gene set | MSigDB 2026.1.Mm Hallmark/PPAR; fry, camera, fgsea, mouse-level scores |

## Reproduced (provisional)

| Paper statement | Our estimate (primary model; muscle-covariate model in brackets where it differs) |
|---|---|
| PE-induced DEGs: 503 CD-only, 182 MD-only, 97 shared (\|log2FC\| ≥ 0.38, FDR < 0.05) | PE_in_CD: 513 genes pass both thresholds. PE_in_MD: 295 (119 with the muscle covariate). The order of magnitude matches for CD; for MD it depends strongly on contamination handling. |
| Shared PE-down lipid genes: Pnpla3, Acly, Me1, Fasn | All four lower with PE in **both** diets (log2FC −0.8 to −3.1; FDR ≤ 0.002 where tested). **Reproduced.** |
| Shared PE-up: Plin4, mt-Nd6, Cyb5a | Plin4 up in both (≈ +1.0; FDR ≈ 0.02). mt-Nd6 up (FDR 0.02). Cyb5a up under CD (FDR < 0.001), weaker under MD. **Reproduced**; Fubp3 not reproduced (n.s.). |
| CD-specific PE-down: Angptl8, Acp3, Lepr, Hmgcr | Angptl8, Acp3, Hmgcr down under CD (FDR ≤ 0.009); Lepr n.s. **Mostly reproduced.** |
| CD-specific PE-up: Cyp4a12a, Cyp26b1, Gadd45g, Prkcb | Cyp4a12a and Cyp26b1 up (FDR 0.02–0.03); Gadd45g borderline (FDR 0.05); Prkcb not testable (NA). **Partly reproduced.** |
| MD+PE-specific: Tm7sf2, Ldlr, Mylip (cholesterol), Pask, Oit3, Nkap | Tm7sf2, Ldlr, Mylip, Oit3 significant under MD (FDR ≤ 0.02); Pask NA; Nkap n.s. **The formal interaction is not significant for any of them.** |
| Ppara induced in CD+PE or MD+PE vs CD+Veh | Within-diet PE effect on Ppara is small and n.s. (+0.16, FDR 0.62 under CD). The paper's comparison against CD+Veh mixes diet and PE. **Not reproduced as a within-diet PE effect.** |
| Nr1i3 (CAR) induced in MD+PE vs CD+Veh | Within-diet PE effect +0.47 (CD) and +0.56 (MD), n.s. **Not reproduced as a within-diet effect.** |
| GSEA: PE enriches a lipid metabolism/toxicity NR set in both diets | Not directly comparable (the custom set may be in the supplement). With public sets: PE lowers cholesterol homeostasis in both diets and raises PPARα targets under CD. |
| PCA shows distinct clustering by group | Diet separation reproduced. Phase 1 showed PC1 of the primary PCA is skeletal-muscle carry-over. |

## Genuine extensions (checked against the paper's full text)

1. **Formal diet × PE interaction** with effect sizes and CIs, at gene and pathway level. The paper reports pairwise comparisons and Venn overlaps only.
2. **Quantitative PE-vs-diet concordance with a shared-reference-corrected design** (split reference plus exhaustive relabelling null). Not in the paper.
3. **QC and confounding analysis:**
   - skeletal-muscle carry-over;
   - the library depth/quality confound with PE (same flowcell, unequal pooling);
   - the circadian-phase difference of the CD+Veh group;
   - Cook's-distance influence;
   - leave-one-out and covariate sensitivity;
   - edgeR/limma comparison.

   None is reported for the bulk data in the paper.
4. **Quantification audit and re-quantification pilot** (tximport setting, length offset).

Novelty has been checked against the paper only. A literature search for other microplastic–liver bulk transcriptomic studies with factorial modelling has **not** been done yet. That search is needed before calling (1) and (2) novel in the final report.

## Bulk vs spatial findings

- The paper's mechanistic claim (PPARα → Anxa2 in pericentral hepatocytes, inflammatory hotspots) rests on **Xenium spatial data, ChIP and in-vitro work**, not on bulk RNA-seq.
- In bulk, Anxa2 shows no PE effect (CD −0.46, MD +0.14; n.s.) and inflammatory Hallmark sets show no PE change.
- This is **not a contradiction**: a change confined to a subset of central-zone hepatocytes or to rare immune niches is diluted in whole-tissue RNA.
- Bulk data neither support nor refute the spatial claims.
