# Direction B: does the PE response depend on diet? (diet × PE interaction)

> **Data status:** deposited GEO matrix used without a tximport length offset, **pilot-verified under assumptions A1–A3** (`reports/quantification_validation/quantification_validation_report.md` §5.1; team decision 2026-11-02). Mitochondrially encoded genes and processed pseudogenes are less reliable. All analyses remain **exploratory** until the team adopts hypotheses.

**Code:**

- `analysis/deseq2/run_deseq2_factorial.R` (coefficient `dietMD.microplasticPE`);
- `analysis/interaction/interaction_patterns.py`;
- `analysis/pathways/pathway_analysis.R`.

**Outputs:** `results/interaction/<label>/`, `results/pathways/<label>/`
**Figure:** `figures/phase2/provisional_deposited/fig18_interaction_vs_muscle.png`

## Question and estimand

Is the PE effect under the MASH diet different from the PE effect under the control diet? The formal quantity is the **interaction coefficient**:

> β_int = (log2 MD+PE − log2 MD+Veh) − (log2 CD+PE − log2 CD+Veh)

This is a difference of log fold changes, i.e. a **ratio of fold changes**. It is defined on the log scale: an effect that is additive in log2 has β_int = 0 even if absolute count changes differ. A large change in MD+PE relative to CD+Veh is **not** evidence of synergy, because it includes the diet effect.

**Pattern categories** (rules fixed in `reports/phase2/analysis_plan.md`, applied to each PE effect's and the interaction's 95% CI):

- similar;
- amplified or attenuated under MD;
- opposing;
- one-diet-only (with or without interaction support).

## Gene level

| Model | Interaction genes (FDR < 0.05) | of which muscle-tracking | Patterns |
|---|---|---|---|
| Primary (`~ diet * PE`) | 87 | **23** (26% vs 1% of tested genes) | CD-only 31, MD-only 31, opposing 22, attenuated 3 |
| + muscle covariate | 16 | 0 | CD-only 6, opposing 6, MD-only 4 |
| without MD_1 | 49 | 3 | CD-only 24, opposing 17, MD-only 6, attenuated 2 |
| 19 leave-one-out refits | 26–131 (median 66) | – | – |
| edgeR quasi-likelihood (same design) | **1** | – | – |
| limma-voom (same design) | **1** | – | – |

**Interpretation:**

1. **A large share of the primary interaction signal is skeletal-muscle carry-over**, not biology.
   - The muscle genes (Myh1, Acta1, Ttn, Atp2a1, Tnnt3 and others) show an apparent PE *increase* under CD, driven by contaminated CD_MP_1, and a strong apparent *decrease* under MD, because the MD+Veh reference contains contaminated MD_1, MD_2 and MD_3.
   - Several estimates are implausible (|log2FC| > 10) because of near-zero counts in some groups.
2. **The remaining interaction genes are method-fragile.** DESeq2 reports 16–87 depending on contamination handling; edgeR-QL and limma-voom, which are more conservative with outlying samples, report 1. Gene-level interaction calls therefore **cannot be treated as established**.
3. **Genes with a clear PE effect in both diets mostly show no detectable difference.**
   - 640 genes have same-sign PE effects with both CIs excluding 0. Only 3 have interaction FDR < 0.05, and 37 have an interaction CI excluding 0.
   - With a median interaction CI width of 0.80 log2, differences up to about ±0.4 log2 (1.3-fold) cannot be excluded. **"No significant interaction" does not mean "identical response".**
4. **Power:** detecting an interaction needs a difference in fold change of roughly 0.5–0.7 log2 for typical genes (feasibility report). That is consistent with few robust hits.

## Pathway level

Using fry (self-contained, mouse = unit) on the Hallmark and PPAR families:

- **No set has a directional interaction at FDR < 0.05** in the primary model.
- **Cholesterol homeostasis:** PE effect −0.19 (CD) and −0.20 (MD), mean shrunken LFC; interaction fry p = 0.84; mean-z interaction −0.003 (95% CI −0.41 to 0.40). The most coherent PE pathway signal is **similar under both diets**: not detectably amplified or attenuated, with the CI allowing moderate differences.
- **PPAR / peroxisome / bile-acid sets, UPR, IFN-γ:** competitive tests (camera, fgsea) suggest diet-dependence; fry does not.
  - PPARα targets: PE up under CD (fry FDR 0.041), weaker under MD (fry p = 0.35); camera interaction FDR 0.003.
  - Read as **"possibly attenuated under MD, exploratory"**.

## Summary of Direction B (provisional)

| Pattern | Evidence |
|---|---|
| Similar PE response in both diets | **Cholesterol-homeostasis program** (down); most genes with clear PE effects in both diets |
| Attenuated under MD | PPARα targets / PPAR signalling: exploratory (competitive tests only) |
| Amplified under MD | No robust example |
| Opposing | Gene-level calls exist but are dominated by contamination and are method-fragile |
| One-diet-only | Many genes reach FDR only in one diet. Without interaction support this reflects power, not established specificity. |

## Differences from the paper

- The paper's MD-specific PE genes (Tm7sf2, Ldlr, Mylip, Pask, Oit3, Nkap) and its description of a "compounded" MD+PE effect come from **pairwise** contrasts and Venn overlaps.
- In our factorial model, Tm7sf2 and Ldlr are significant under MD (FDR ≤ 0.001), and not under CD in either the primary or the muscle-covariate fit (FDR ≥ 0.07). Their **interaction** is not significant.
- So "specific to MD" is not supported as a formal interaction (see literature comparison).

## Course connections

- Interaction = effect modification in a two-way factorial model, as in two-way ANOVA.
- An interaction needs about twice the standard error of a simple effect, so it needs more data.
- Scale dependence: a log-scale (multiplicative) model defines "no interaction" as equal *fold* changes.
- An influential observation (MD_1) can create an apparent interaction. That is why the diagnostics and the leave-one-out analysis matter.
