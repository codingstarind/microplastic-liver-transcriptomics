# Direction C: robustness and sensitivity analysis

> **Status: PROVISIONAL** (all variants use the deposited counts). A second layer of
> robustness, re-quantified counts vs deposited counts, will be added after the pilot.

**Code:**

- `analysis/robustness/run_sensitivity.sh` (variants);
- `analysis/robustness/summarise_robustness.py`;
- `analysis/deseq2/run_alternative_methods.R`;
- `analysis/concordance/pe_vs_diet_concordance.py --drop=`.

**Outputs:** `results/robustness/provisional_deposited/`, `results/deseq2/*__*`, `results/deseq2/loo/`

## Design of the checks (fixed in `reports/phase2/analysis_plan.md`)

| Check | What it probes | Notes |
|---|---|---|
| Leave-one-mouse-out (19 refits) | Dependence on any single animal | No group falls below 3 mice |
| Drop MD_1 | Heaviest skeletal-muscle contamination | Pre-specified sensitivity, not a QC exclusion |
| Muscle covariate (z of log2 muscle-marker CPM) | Carry-over of diaphragm tissue | |
| Depth covariate (z of log2 read pairs) | Sequencing-depth difference | **Collinear with PE** (R² = 0.78 on the design); a stress test, not a correction. Size factors already normalise depth, so this is not double correction of library size but a test of depth-linked composition effects. |
| Quality-proxy covariate (gene counts per read pair) | Library quality/composition | **Collinear with PE** (R² = 0.80) |
| Clock-phase covariate | Circadian/feeding state at sacrifice | Nearly an indicator of the CD+Veh group |
| edgeR quasi-likelihood, limma-voom | Software/model dependence | Same design and contrasts. Agreement ≠ validation. |
| Gene-class checks | Muscle-tracking, mitochondrial, clock, top-abundant genes among FDR hits | |

## Robustness matrix (headline findings × checks)

Legend:

- ✔ survives: same direction, and at the primary's significance level or better;
- ~ weakened: same direction, weaker evidence;
- ✖ lost: no longer distinguishable from 0;
- – not applicable.

| Headline finding (provisional) | 19 LOO | −MD_1 | +muscle | +clock | +depth | +quality | edgeR/limma |
|---|---|---|---|---|---|---|---|
| **G1** PE under CD changes many genes (858 FDR) | ✔ 545–909 genes; stat ρ ≥ 0.92 | ✔ 836 | ✔ 928 | ~ 575 | ✖ 116 (ρ = 0.70) | ✖ 66 (ρ = 0.50) | ~ 544 / 517 (LFC ρ ≥ 0.96) |
| **G2** PE under MD changes many genes (472 FDR) | ✔ 273–622 | ~ 286 | ✖ 135 | ✔ 522 | ✖ 75 | ✖ 133 | ~ 211 / 178 |
| **G3** Diet effect (2,262 FDR) | ✔ 1,699–2,537 | ✔ 1,924 | ~ 1,483 | ~ 747* | ✔ 2,228 | ✔ 2,265 | ✔ 2,080 / 2,101 |
| **G4** Gene-level interaction (87 FDR) | ~ 26–131 | ~ 49 | ✖ 16 | ~ 22 | ✔ 89 | ~ 67 | ✖ **1 / 1** |
| **P1** PPARα targets up with PE under CD (fry FDR 0.04) | ✔ 19/19 same direction, p < 0.05 (11/19 FDR < 0.05) | ✔ | ✔ (FDR 0.03) | ✔ (FDR 0.03) | ✖ (p = 0.18) | ✖ (p = 0.22) | – |
| **P2** Cholesterol homeostasis down with PE under CD (fry p = 0.004) | ✔ 19/19, p < 0.05 | ✔ | ✔ | ✔ | ~ (p = 0.035) | ✖ (p = 0.55) | – |
| **P3** Cholesterol homeostasis down with PE under MD (fry p = 0.005) | ✔ 19/19, p < 0.05 | ✔ | ~ (p = 0.025) | ✔ | ~ (p = 0.025) | ~ (p = 0.14) | – |
| **P4** No pathway interaction detected (fry) | ✔ 18/19 refits detect none; 1 refit flags 1 set | ✔ | ✔ | ✔ | ✔ | ✔ | – |
| **A1** PE(CD)–diet gene-level resemblance (split ρ = 0.33) | ✔ median ρ 0.25–0.36; every split ≥ 0.21 | ✔ ρ = 0.35 | – | not interpretable** | – | – | – |
| **A2** No PE(MD)–diet resemblance (ρ = 0.05) | ✔ median ρ −0.04 to 0.13 | ✔ ρ = 0.10 | – | – | – | – | – |

\*The clock score is nearly an indicator of CD+Veh, so it absorbs much of the diet contrast.
\*\*The concordance from the clock-adjusted model is not interpretable: both contrasts share the adjusted reference.

**Leave-one-out concordance (A1, A2):** `results/concordance/provisional_deposited/leave_one_out_split_rho.csv`.

- This gives the observed split-reference ρ for each of the 19 single-mouse removals.
- The relabelling null was not recomputed per removal: about 25 min each, against an expected null centre near 0, as in the full analysis.
- Removing CD_MP_2 gives the lowest PE(CD)–diet median ρ (0.25); removing CD_MP_5 the highest (0.36).

## Gene-class checks (primary model, FDR < 0.05 genes)

| Contrast | Muscle-tracking | Mitochondrial (mt-) | Clock genes | Top-30 abundant |
|---|---|---|---|---|
| PE_in_CD | 1 (0.1% vs 0.95% background) | **11 (1.3% vs 0.11%)** | 3 | **9 (1.0% vs 0.22%)** |
| PE_in_MD | **31 (6.6% vs 0.95%)** | 2 | 0 | 1 |
| Diet_in_Veh | 55 (2.4%) | 1 | 11 | 7 |
| Interaction | **23 (26%)** | 0 | 3 | 0 |

**Readings:**

- PE_in_MD and the interaction are contaminated by muscle genes.
- PE_in_CD is enriched for mitochondrially encoded and very abundant transcripts. Those are the transcript classes most sensitive to library quality and depth, and PE libraries are shallower with a lower counts-per-read-pair proxy.
- Clock genes: Bmal1, Clock and Npas2 are higher and Bhlhe40, Dbp and Per3 lower in PE_in_CD, the same direction as Diet_in_Veh, because the CD+Veh reference has a distinct clock phase.

## Confounders that cannot be separated from treatment

| Confounder | Aligned with | Evidence | Consequence |
|---|---|---|---|
| Library pooling depth + quality proxy | PE (both diets) | Perfect separation of read pairs; lower gene counts per read pair; same flowcell | Every PE effect is partly or wholly explainable by library properties. Adjusting removes most PE gene-level signal and the P1/P2 pathway signals. |
| Library-prep order | Group (S39–S57 in group order) | Index order | Any prep-batch effect aligns with group |
| Clock phase / feeding state | CD+Veh group | All CD+Veh mice below every other mouse on the clock score | Inflates PE_in_CD and Diet_in_Veh in the same direction for clock-controlled genes (relevant to A1) |
| Cage | Group (group-housed) | Methods | Cage effects are inseparable from treatment |

## What is robust, what is fragile

- **Robust to single mice, contamination, clock adjustment and software:**
  - the diet effect;
  - the direction of the PE effect on lipogenic/cholesterol genes and on PPARα targets under CD;
  - the absence of detectable pathway-level interaction.
- **Fragile:**
  - the gene-level interaction list (method- and contamination-dependent);
  - the PE_in_MD gene count (contamination-dependent);
  - every PE finding under depth/quality adjustment, which is collinear with PE and so **unresolvable with this design**.
