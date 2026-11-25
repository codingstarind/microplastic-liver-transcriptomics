# Transcript: Published gene replication in microplastic-liver transcriptomics

- Claude Code session ID: `eb0c5c48-cf91-4643-bd2e-12cb718a1785`
- Models in this export: claude-opus-5-5 (129 replies)
- Prompts and written replies only; tool calls and outputs omitted.

### User (2026-11-21T09:48:09.732Z)

BT3041 — Cross-Paper Gene Replication and Course-Methods Analysis
We have identified a promising research lead in our microplastic-liver transcriptomics project.
Published mouse experiments report alterations in lipid-metabolism genes following microplastic exposure. In our GSE294180 dataset, our preliminary analysis suggests that some of these genes show changes in the opposite direction.
Our immediate objective is to determine whether these published findings reproduce in our dataset and whether the differences are statistically convincing.
I want this analysis to explicitly incorporate appropriate methods from our BT3041 course, rather than relying only on standard differential-expression analysis.
Our research questions and biological interpretations will remain team-owned. Do not independently invent or validate hypotheses.
Work in our existing `microplastic-liver-transcriptomics` repository, starting from `phase2/factorial-analysis`. Preserve existing analyses and create a separate research directory.
Part 1 — Published gene replication
Start with Fasn, Acly, Acaca and Ppara from Djouina et al. (2023).
Verify their exact claims and experimental contrasts directly from the original paper.
As a separate exploratory panel, examine Acot3, Abcc3, Nr1i3 and Nr4a1 from published polystyrene studies.
Use the existing DESeq2 results for:

* PE under control diet
* PE under MASH diet
* Diet effect
* Diet × PE interaction

For every gene, report log2 fold change, 95% confidence interval, p-value, adjusted p-value, expression level and agreement with the published direction.
Do not interpret matching signs alone as statistically demonstrated replication.
Part 2 — BT3041 Topic 2: Hypothesis Testing
Use the existing factorial analysis to examine the prespecified published gene panel.
Present the null hypotheses, effect estimates, uncertainty and statistical evidence.
Distinguish:

* Same direction, statistically supported
* Opposite direction, statistically supported
* Same direction, statistically inconclusive
* Opposite direction, statistically inconclusive
* Not assessable

Explain the role of multiple-testing correction and why nonsignificant results do not prove absence of an effect.
Produce a gene-level forest plot showing effect sizes and confidence intervals.
Part 3 — BT3041 Topic 4: Regression and Interactions
Use the existing factorial design:
Gene expression ~ Diet + PE + Diet × PE
Explain what each term represents.
For the published gene panel, estimate whether the PE-associated expression change differs between the control and MASH diets.
Report the interaction estimates with uncertainty rather than only significant-gene counts.
Compare the existing additive and factorial formulations where statistically appropriate, making clear that model complexity is not itself evidence of biological interaction.
Part 4 — BT3041 Topic 5: Model Validation and Robustness
Perform leave-one-mouse-out sensitivity analysis for the selected genes.
For each omitted mouse, refit the appropriate model and record changes in estimated PE effects, confidence intervals and direction.
Summarize how frequently each gene retains the original direction and whether any single mouse drives the conclusion.
Where feasible, perform a clearly labelled exploratory stratified bootstrap within experimental groups. Preserve the factorial structure and acknowledge the limitations caused by the small sample size.
Compare the primary analysis against existing depth-adjusted and muscle-adjusted results.
Be especially careful because some technical covariates are strongly associated with treatment. Adjustment may not identify a causal PE effect.
Part 5 — BT3041 Topic 9: Hierarchical Clustering
Prepare a heatmap of the published candidate genes using appropriately normalized, transformed and gene-standardized expression values.
Cluster genes by similarity of expression patterns and show mouse annotations for diet and PE treatment.
Compare reasonable similarity measures and explain which one is appropriate.
Investigate whether lipid-metabolism genes form coherent expression groups.
Do not claim that a small cluster establishes a biological pathway or causal mechanism.
Part 6 — BT3041 Topic 10: PCA
Perform PCA restricted to the predefined published gene panel as an exploratory visualization.
Show whether gene-expression variation aligns more with diet, PE exposure or individual-mouse differences.
Report explained variance and relevant gene loadings.
Because the panel contains very few genes, explicitly state that this is a descriptive visualization, not an independent validation or classification test.
Part 7 — Biological interpretation and literature comparison
Compare our gene-level results with the original experimental literature.
For every discrepancy, document relevant differences in:

* Mouse strain and sex
* Polymer type and particle size
* Exposure dose, route and duration
* Dietary and metabolic conditions
* RNA-seq versus qPCR or protein-based measurements
* Statistical modelling and validation methods

Do not assume that differences between studies are caused by diet or exposure conditions without supporting evidence.
Do not conduct pathway enrichment yet. We'll use the gene-level findings to decide which enrichment analysis is biologically justified.
Deliverables
Save everything under `research/cross_study_replication/`:

1. `published_gene_evidence.csv` — verified literature claims and experimental contexts
2. `gene_replication_results.csv` — effect sizes, intervals, p-values and directional agreement
3. `interaction_results.csv` — diet × PE interaction estimates
4. `leave_one_out_results.csv` — sample-level robustness analysis
5. `gene_replication_forest.png` — effect-size comparison
6. `gene_expression_heatmap.png` — hierarchical clustering
7. `gene_panel_pca.png` — exploratory PCA
8. `replication_report.md` — methods, findings, uncertainty and literature comparison
9. Reproducible scripts and instructions

For every analysis, explain:

* Which BT3041 concept is being applied
* Why it is appropriate
* What biological question it addresses
* What the result establishes
* What the result cannot establish

Prioritize scientific correctness over generating more plots.
The final report should distinguish verified results, exploratory findings, unsupported claims and limitations.
Our team will evaluate the evidence and independently decide the next research questions.

*[7 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T09:50:04.360Z)

Reading the existing DESeq2 outputs and the source papers now; the Djouina PDF is in your Downloads folder, so I can check its claims against the original text.

*[7 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T09:51:23.252Z)

Checking whether the R/DESeq2 environment in WSL is reachable so I can refit models; no additive model exists in the repo yet, so I'll fit one.

*[6 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T09:52:47.751Z)

Literature claims for Djouina are verified from the PDF; a background agent is checking the polystyrene papers (Acot3/Abcc3/Nr1i3 and the Nr4a1 source). Now writing the model-fitting script.

*[10 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T09:55:08.251Z)

Writing the Python script that builds the tables, forest plot, heatmap and PCA while the bootstrap runs.

*[17 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T10:00:30.773Z)

Fixed the label overlaps in the forest plot and heatmap. The full 2000-replicate bootstrap is still running; meanwhile I'm reading the summary tables.

*[6 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T10:02:02.838Z)

Drafting the replication report now; the bootstrap is still running in WSL.

*[8 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T10:05:31.179Z)

The bootstrap has been running about 10 minutes. Checking it's still computing and not stuck.

*[4 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T10:06:20.256Z)

Everything except the bootstrap section is built: evidence table, replication, interaction, leave-one-out and covariate tables, all three figures, the report and the README. The 2000-draw bootstrap is still running in WSL, about 11 minutes in. When it finishes I'll fill in the report's bootstrap section (§6.2) and regenerate the outputs that use it.

*[6 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T10:17:46.569Z)

The Djouina 2023 lipogenic findings do not replicate in our mice: with PE, Fasn, Acly and Acaca are clearly lower, not higher, under both diets. Ppara shows no detectable PE effect. Everything is in `research/cross_study_replication/`, and existing analyses are untouched. Read [replication_report.md](research/cross_study_replication/replication_report.md) first.

**Published claims, checked against the originals**
- **Djouina 2023:** read from the original PDF, including Fig. 1 at high resolution.
  - The text says Ppara rose "two-fold after exposure to both types" of bead. The figure shows the rise only in the group given both bead sizes together, which got double the particle mass. The 116 µm beads alone gave a significant *decrease* (about 0.74×).
  - The Fasn rise was gone by 9 weeks. Only Acaca was up in every group at both time points.
  - Fold changes come from reading the bar charts by eye (about ±0.1), because the paper prints no values.
- **Li 2025 (Acot3, Abcc3, Nr1i3):** full text read. These three were all reported up.
- **Nr4a1:** the main source (Chiu 2025) is paywalled, so I could only read its abstract. One full-text qPCR study supports the "up" direction. Someone on the team should get Chiu 2025 in full before relying on it.

**Results**

| Gene | PE effect, control diet (log2FC [95% CI]) | Category |
|---|---|---|
| Fasn | −1.58 [−2.17, −0.99] | Opposite direction, statistically supported |
| Acly | −1.82 [−2.28, −1.36] | Opposite direction, statistically supported |
| Acaca | −1.61 [−1.96, −1.26] | Opposite direction, statistically supported |
| Ppara | +0.16 [−0.14, +0.45] | Same direction, inconclusive |
| Abcc3, Nr1i3 | about +0.5 | Same direction, inconclusive |
| Acot3, Nr4a1 | very wide intervals | Data uninformative |

- **Lipogenic genes:** the decrease is the same under the MASH diet (FDR ≤ 0.004). It holds in all 19 leave-one-mouse-out refits, all 2000 bootstrap draws, and under muscle and clock adjustment.
- **Ppara:** the interval rules out a two-fold rise, so the data argue against Djouina's headline Ppara effect, not just fail to find it.
- **Abcc3 and Nr1i3:** pooled over both diets they reach FDR < 0.05. They lose significance after adjusting for sequencing depth or yield.
- **Diet × PE interaction:** seven of eight genes show none. Acaca passes the genome-wide threshold (FDR 0.038), but significance is lost in 13 of 19 refits and under the muscle covariate. The diet alone already lowers Acaca strongly, so the "interaction" could be a floor effect.
- **Additive vs factorial model:** for six genes the penalised criteria prefer the simpler additive model.
- **Clustering:** the three lipogenic genes always group together, and Ppara/Abcc3/Nr1i3/Acot3 form a second group, under every distance measure tried. The negative correlation between the two groups disappears once group means are removed: it came from opposite responses to treatment, not from the genes tracking each other mouse to mouse.
- **PCA:** the first component explains 58% of variance; diet accounts for 43% of that and PE for 38%. This is descriptive only.

**What the data cannot show**
- **PE is not shown to be the cause.** Every PE library is shallower and lower-yield than every vehicle library, so adjusting for depth cannot separate a PE effect from a library artefact.
- **Clock state is a partial caveat.** The CD+Veh group has a distinct circadian state, and these lipogenic genes are strongly diurnal. The MASH-diet comparison doesn't use CD+Veh at all and shows the same decrease, so the clock state doesn't explain the result on its own.
- **Selection effect.** The panel was examined after a preliminary look at these same data. That's why I used genome-wide FDR as the main test rather than a 4-gene correction.

The literature comparison lists the differences for every discrepancy (sex, polymer and particle size, dose and route, purified high-sucrose diet vs chow, qPCR vs RNA-seq, statistics, timing) without claiming any of them as the cause. Hypotheses and pathway enrichment are left to the team; the report ends with open items for you to decide.

The R step that refits the models takes about 22 minutes (mostly the bootstrap); the README gives the commands for all three scripts. The new scripts pass `ruff`.
