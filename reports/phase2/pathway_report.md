# Pathway-level analysis: which hepatic programs respond to PE, and does diet change the response?

> **Data status:** deposited GEO matrix used without a tximport length offset, **pilot-verified under assumptions A1–A3** (`reports/quantification_validation/quantification_validation_report.md` §5.1; team decision 2026-11-02). Mitochondrially encoded genes and processed pseudogenes are less reliable. All analyses remain **exploratory** until the team adopts hypotheses.

**Code:** `analysis/pathways/pathway_analysis.R`
**Outputs:** `results/pathways/provisional_deposited/`
**Gene sets:** MSigDB **2026.1.Mm** (mouse-native) via msigdbr 26.1.1, matched by mouse gene symbol.

| Family | Sets | Role |
|---|---|---|
| Hallmark (MH) | 50 | Primary |
| PPAR | WP_PPAR_SIGNALING_PATHWAY, REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA, SANDERSON_PPARA_TARGETS | Separate BH family |
| Reactome | 911 sets with ≤ 500 genes | Exploratory only |

## Biological question

Single genes are noisy and numerous. Do coordinated programs, such as lipid synthesis, cholesterol homeostasis, PPARα signalling, inflammation, xenobiotic handling, mitochondria and fibrosis, shift with PE, and is the shift the same under both diets?

## Why three methods, and what each assumes

| Method | Unit of analysis | Question asked | Main weakness here |
|---|---|---|---|
| **Mean-z score + linear model** (two-way factorial `lm`) | mouse | Does the *average* standardised expression of the set change? | Dilutes sets whose genes move in opposite directions. Hallmark sets do: e.g. 30 cholesterol-homeostasis genes go down and 10 go up with diet. |
| **limma fry** (rotation test) | mouse | Are the set's genes, *taken together*, changed (directional or mixed)? | The mixed-direction test flags almost every set for every contrast, so it is uninformative; we use the directional test. |
| **limma camera** | gene, correcting for inter-gene correlation | Are the set's genes *more* changed than other genes (competitive)? | Competitive tests answer a different question from self-contained tests. |
| **fgsea** (preranked on the Wald statistic) | gene | Are set genes concentrated at the top or bottom of the ranking? | Ignores correlation between genes, so p-values are **anti-conservative**. Supporting evidence only. |

**Deviation from the pre-specified plan** (`reports/phase2/analysis_plan.md`):

- The plan named the mean-z score as the primary method. A positive-control check showed it detected none of the strong diet effects at FDR < 0.05, because of dilution.
- We therefore added fry and camera (mouse-level, correlation-aware) and **base the PE-response categories on fry**. The mean-z categories are still reported (`pe_response_categories.csv`).
- The decision was made after seeing the diet results, but before inspecting PE pathway results for this purpose.
- Two set changes: the Reactome PPARα set has a new name in the 2026 release, and SANDERSON_PPARA_TARGETS was added from the set list.

## Results (focus themes; fry directional test unless stated)

| Program | PE under control diet | PE under MASH diet | Diet (vehicle) | Interaction |
|---|---|---|---|---|
| **Cholesterol homeostasis** | Down (fry p = 0.004, FDR 0.11; camera FDR < 0.001; fgsea NES −2.3) | Down (fry p = 0.005, FDR 0.09; camera FDR < 0.001; NES −2.6) | Down (fry FDR < 0.001) | Nothing: fry p 0.84, camera FDR 0.93, mean-z estimate −0.003 (95% CI −0.41 to 0.40) |
| **PPARα targets** (Sanderson) | **Up (fry FDR 0.041**; camera FDR < 0.001) | Up, weaker (fry p = 0.35; camera FDR 0.075) | Up (fry FDR < 0.001) | camera FDR 0.003 (down); fry p = 0.21 |
| WP PPAR signalling | Up, n.s. | Down (fry p = 0.055) | Up (fry FDR 0.018) | camera FDR < 0.001; fry p = 0.08 |
| Adipogenesis | n.s. | Down (fry p = 0.017, FDR 0.14; camera FDR 0.002) | n.s. | n.s. by fry |
| Oxidative phosphorylation | Down by camera/fgsea only (camera FDR 0.014) | Down (fry p = 0.039; camera FDR < 0.001) | n.s. by fry | n.s. |
| Fatty-acid metabolism | n.s. by fry/camera (fgsea NES +1.6) | n.s. | Up (camera FDR 0.001) | n.s. by fry |
| Xenobiotic metabolism | n.s. (fgsea NES +1.6 only) | n.s. | n.s. | n.s. |
| Inflammation (4 sets) | n.s. (all fry p > 0.5) | n.s. (fry p > 0.39) | n.s. by fry | n.s. by fry; camera flags IFN-γ response (FDR 0.033) |
| Oxidative / cellular stress | n.s. | n.s. | ROS up (fry FDR 0.04) | UPR: camera FDR 0.02 only |
| Fibrosis / EMT, TGF-β | n.s. | TGF-β up (fry p = 0.025, FDR 0.18) | EMT up (camera FDR 0.001) | n.s. by fry |

(n.s. = not significant at FDR < 0.05 in the stated test.)

**PE-response categories (fry, Hallmark + PPAR, 53 sets):**

- 52: "no clear PE effect" at FDR < 0.05;
- 1: "CD-only, interaction not supported" (SANDERSON_PPARA_TARGETS).

Read this as **limited power at the pathway level after multiplicity correction**, not as an absence of effects. The figure `figures/phase2/provisional_deposited/fig16_pathway_pe_response_forest.png` shows the estimates with 95% CIs.

## What we can say (provisional)

1. **Lower cholesterol-homeostasis expression with PE, under both diets, with similar size and no detectable interaction.**
   - Mean shrunken gene LFC: −0.19 (CD) and −0.20 (MD).
   - Same direction as the diet effect (−0.61).
   - This is the most consistent PE pathway signal across methods: fry p < 0.01 in both diets, camera and fgsea strongly significant. It does not reach fry FDR < 0.05.
   - It fits the gene-level reproduction of the paper (Hmgcr, Ldlr, Tm7sf2 lower).
   - "Similar under both diets" is a *lack of evidence for difference*. The interaction CI still allows moderate differences.
2. **Higher PPARα target expression with PE under the control diet**, also in the diet's direction. Diet-dependence is suggested by camera but not by fry; **unresolved**.
3. **No coherent inflammatory, fibrotic or xenobiotic program change with PE** at this sample size. This contrasts with the paper's emphasis on inflammation, which came from histology and *spatial* data, not bulk pathways.
4. **Interaction at the pathway level:** competitive tests (camera, fgsea) flag PPAR/peroxisome/bile-acid sets (less PE-induced under MD) and UPR. The self-contained fry test does not. Because camera and fgsea answer a competitive question and fgsea ignores gene correlation, we classify these as **exploratory**.

## Alternative explanations

- **Feeding and circadian state.** Cholesterol synthesis (SREBP2 targets) and PPARα targets are among the most feeding- and time-of-day-responsive programs in liver; fasting lowers SREBP targets and raises PPARα targets. The CD+Veh mice have a distinct clock-gene phase (see interpretation), which could produce exactly the PE_in_CD pattern. Under MD, where the clock does not differ, the cholesterol signal is still present, which partly argues against a pure timing artefact.
- **Library quality.** Adding depth or quality-proxy covariates removes most PE_in_CD gene-level signal (see the robustness report), but these covariates are collinear with PE, so this cannot discriminate.
- **Reduced food intake under PE gavage** would also produce a fasting-like profile. Food intake was measured by the authors but is not available per mouse.

## Course connections

- **Unit of analysis:** gene-level tests (fgsea) treat ~13,000 correlated genes as independent observations and so overstate certainty. Mouse-level tests (fry, mean-z regression) respect that the experiment has only 19 independent units.
- **Self-contained vs competitive hypotheses** are different null hypotheses, not different power levels of one test.
- An **enrichment score is not a pathway activity measurement**: it summarises mRNA changes of member genes, not flux or protein activity.
