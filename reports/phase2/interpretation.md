# Phase 2 scientific interpretation (provisional, for team review)

> Everything below rests on the **deposited count matrix used without a tximport length
> offset**, which is still being validated (`reports/quantification_validation/`). These are
> computational results for the team to examine. They have not been independently validated
> by students, and they are not final claims. Mechanisms are not inferred: bulk mRNA changes
> do not measure protein, flux or pathway activity.

## Central question

How does oral polyethylene (PE) microplastic exposure relate to hepatic transcriptional programs, and does a MASH-inducing diet change that relationship?

## What the data currently show

### 1. A reproducible PE-associated transcriptional signature exists, but it is confounded with library properties

- Under the control diet, PE is associated with lower expression of **lipogenic and cholesterol-synthesis genes** and higher **PPARα target** expression.
  - Lipogenic: Fasn, Acly, Me1, Pnpla3, about 2–8-fold lower.
  - Cholesterol synthesis: Hmgcr; Hallmark cholesterol homeostasis.
  - PPARα targets: Sanderson set; Cyp4a12a, Cyp4a14.
- The lipogenic and cholesterol pattern is also seen under the MASH diet.
- These directions reproduce the paper's gene-level statements.
- The patterns survive leave-one-mouse-out, the muscle-contamination covariate, dropping MD_1 and the clock-phase covariate.
- **However:**
  - every PE library is shallower and has a lower gene-counts-per-read-pair proxy than every vehicle library, all on the same flowcell;
  - adjusting for these properties removes most of the PE signal;
  - because the properties are almost perfectly aligned with PE, **the data cannot tell whether the PE signature is biological, technical, or both**;
  - PE_in_CD hits are over-represented among mitochondrial and very abundant transcripts, the classes most sensitive to library quality.

**Status:** reproducible within this dataset; causal attribution to PE **unresolved**.

### 2. Under the control diet, the PE response partly resembles the diet response; under the MASH diet it does not

- Gene-level fold changes for PE under CD correlate moderately with the diet effect: Spearman ρ ≈ 0.33 after removing the shared-control artefact, which accounted for about 0.2 of the naive 0.54.
- The resemblance is concentrated in **cholesterol (down) and PPARα (up) programs**, not transcriptome-wide pathway structure.
- Under the MASH diet, PE changes do not align with the diet axis (ρ ≈ 0.05).
- **Major alternative explanation:** the four CD+Veh mice, the shared reference for both comparisons, have a distinct **circadian/feeding state**. All four score below every other mouse on a pre-defined clock-gene phase score. Sacrifice time or feeding differences in this group would produce exactly this resemblance for clock- and feeding-controlled liver genes, and would explain its absence under MD.

**Status:** exploratory. "PE mimics metabolic dysfunction" is **not supported** as a claim.

### 3. Diet-dependence of the PE response is not established

- The pathway-level interaction is not detected by self-contained (fry) tests in any model or leave-one-out refit, except one refit that flags a single set.
- The most coherent PE program (cholesterol homeostasis down) has nearly identical PE effects under both diets. Its interaction estimate is ≈ 0, but the CI is wide (±0.4 in score units).
- Gene-level interaction hits are dominated by skeletal-muscle carry-over: 23 of 87 are muscle-tracking. They shrink to 16 with a muscle covariate and to **1** with edgeR or limma.
- The competitive tests hint that the PE induction of PPAR/peroxisome programs may be attenuated under MD (exploratory).

**Status:** "similar response under both diets" for cholesterol homeostasis is the best-supported description. Amplification or synergy under MASH is **not supported**. Absence of an interaction is **not demonstrated**; the power is low.

### 4. No bulk evidence for PE-associated inflammatory, fibrotic or xenobiotic program changes

- Hallmark inflammation, IL6/JAK/STAT3, TNF-α, IFN-γ, EMT/TGF-β and xenobiotic sets show no PE effect in self-contained tests.
- The paper's inflammatory and Anxa2 findings come from spatial data, which bulk RNA dilutes.

**Status:** consistent with spatially restricted effects; neither supports nor refutes them.

## Statistical significance, magnitude, relevance

- Lipogenic genes change by large factors (2–8-fold).
- Pathway-level shifts are modest: about −0.2 mean log2FC for cholesterol homeostasis.
- Whether such changes matter physiologically cannot be judged from mRNA alone. The authors' phenotypes (ALT, triglycerides, histology) were not deposited per mouse.

## Decisions for the research team

| # | Decision | Options and consequences |
|---|---|---|
| 1 | **Count input for confirmatory analysis** | Pending the pilot. If the deposited matrix needs an offset: full re-quantification (about 1 day of downloads at observed speeds, ~35 GB) or accept the provisional status. |
| 2 | **Primary model: include the muscle covariate?** | Pre-specified primary = no covariate. Evidence shows contamination drives PE_in_MD and interaction hits. Making the muscle covariate primary is defensible but is a post-hoc change; it must be declared as such. |
| 3 | **How to present the library-depth/quality confound** | It cannot be removed. It should be stated as a principal limitation of any PE claim. |
| 4 | **Direction A framing** | Keep as exploratory, or invest in an independent MASH signature (candidate GEO series listed in the similarity report) to remove the same-data circularity. |
| 5 | **Circadian confound** | Ask the authors for sacrifice times and food intake; otherwise treat clock-controlled programs as confounded. |
| 6 | **Which pathway method is primary** | fry (mouse-level, conservative) vs camera/fgsea (gene-level/competitive, more sensitive). Recommendation: fry primary, others supporting. |
| 7 | **Hypotheses to adopt** | Candidate testable hypotheses: (a) PE lowers hepatic cholesterol-synthesis program expression under both diets; (b) the PE effect on PPARα targets differs by diet. Both must be stated as conditional on the unresolved library confound. |

## Exploratory vs confirmatory

All analyses so far are **exploratory**: provisional counts, and hypotheses not yet adopted by the team. A confirmatory analysis would:

- fix the hypotheses and model (including the covariate decision) **before** re-running on validated counts;
- report only those pre-specified tests as confirmatory.
