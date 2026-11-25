# Transcript: BT3041 microplastic-liver research strategy

- Claude Code session ID: `5dbfd995-fe63-4ce9-8501-246b4fa39b6f`
- Models in this export: claude-opus-5-5 (88 replies)
- Prompts and written replies only; tool calls and outputs omitted.

### User (2026-11-21T15:26:32.879Z)

Sairam. I want to restart our BT3041 microplastic-liver project in a more systematic way, because I think we've been doing a lot of analysis without a clear enough scientific contribution.
I don't want us to simply keep analysing GSE294180, adding plots, or trying models because we have learned them in class. The point is to first understand the existing research, identify a meaningful gap, examine how other researchers approached similar problems, and then decide how we can improve on their methods and compare our results fairly.
This chat is for Phase 1 only: Research Problem, Literature Gaps, and Model Audit.
Background and existing work
Our repository is:
https://github.com/codingstarind/microplastic-liver-transcriptomics
Our original dataset is GSE294180, containing 19 mouse liver RNA-seq samples in a 2×2 experimental design:

* Control diet + vehicle: 4 mice
* Control diet + PE: 5 mice
* MASH-inducing diet + vehicle: 5 mice
* MASH-inducing diet + PE: 5 mice

Our original two research questions were:

1. Do polyethylene-associated liver gene-expression changes resemble the changes induced by a MASH-inducing diet?
2. Does an already metabolically stressed liver respond differently to PE exposure?

We can retain these questions and add new ones, provided they remain in the same broad area: microplastics, hepatic metabolism, metabolic dysfunction, transcriptomics and molecular responses.
We've already done extensive analysis, particularly in the `phase2/factorial-analysis` branch. This includes preprocessing, QC, PCA, clustering, DESeq2, factorial regression, pathway analysis, PE–diet similarity testing, model sensitivity analysis and a re-quantification pilot.
Please read the existing project reports, scripts and results rather than starting from scratch.
There is also newer work from a previous Claude session in `research/cross_study_replication/`, which reportedly includes cross-paper comparisons of Fasn, Acly, Acaca, Ppara and other genes against Djouina, Li and Chiu, with additive/factorial models, leave-one-out testing, bootstrapping, PCA, clustering and plots.
Check whether these files are actually available locally. Don't assume they are committed to GitHub.
Some important limitations we've already encountered are technical differences strongly associated with PE exposure, muscle-associated markers in certain samples, uncertainty over the deposited count matrix's transcript-length offsets, and fragile gene-level interaction findings.
These issues must remain visible in our future analysis.
Task 1 — Consolidate our existing work
Before doing anything new, inspect what the project actually contains.
Read the Phase 1 and Phase 2 reports, analysis plan, quantification validation, pathway and robustness reports, final report draft, and the newer cross-study replication outputs if available.
I want a concise inventory separating:

* What we've completed computationally.
* Which results are reasonably stable.
* Which results are fragile or exploratory.
* Which published findings we've already reproduced using the same data.
* What we haven't independently validated.
* Which analyses are missing or unfinished.

Don't rerun the original analysis just to make this inventory.
Task 2 — Revisit the literature with a modelling perspective
We previously reviewed roughly 28 distinct papers and resources, but I feel that we focused too much on summarizing their biological findings.
This time, I want a focused review of around 12–18 important original research papers, with particular attention to their methods.
Start from Jung et al. 2026, Djouina et al. 2023, Li et al. 2025, Chiu et al. 2025, relevant independent PE/PS liver transcriptomic experiments, metabolic-disease interaction studies and appropriate transcriptomic modelling or meta-analysis papers.
Verify the primary sources. Use full texts and supplementary materials whenever accessible.
For each paper, determine:

1. What exact scientific question was it trying to answer?
2. What was the experimental design and sample size?
3. Which datasets were used?
4. Which statistical or machine-learning models were implemented?
5. Why was that model suitable for the paper's biological question?
6. Which inputs, features, outputs and targets were used?
7. How were the model's results evaluated?
8. Was there independent experimental validation, external validation or only an internal analysis?
9. What were the paper's main findings?
10. What methodological or biological limitations remain?

Please distinguish actual limitations from things a paper simply did not attempt to investigate.
Provide DOI links and precise figure, table, or methods-section references for important claims.
Task 3 — Identify research gaps properly
I don't want a vague list of possible future directions.
For each candidate gap, explain:

* What has already been demonstrated.
* What remains uncertain.
* Why that uncertainty matters biologically.
* Which published methods are relevant.
* Whether our available data could realistically resolve it.
* What would count as meaningful new evidence.
* What would make the proposed investigation infeasible or inconclusive.

Also check whether other studies have already attempted the proposed improvement. Don't assume novelty from one paper's limitations.
Task 4 — Identify possible modelling approaches
Please examine three broad directions, without committing to any of them yet.
A. Factorial and pathway-level statistical modelling
Compare conventional pairwise differential-expression approaches with factorial models and gene-set-level inference. Assess whether our improvements are genuinely useful given the limitations of GSE294180.
B. Cross-study molecular-response modelling
Investigate whether cross-study effect-size synthesis, meta-analysis, or pathway-level reproducibility methods could improve on isolated DEG and hub-gene analyses.
Look carefully at study heterogeneity, biological independence and whether the effect sizes are comparable.
C. Predictive machine learning
Investigate whether we can build a scientifically meaningful predictive model, potentially using regularized regression, published gene panels or pathway features.
But first establish a meaningful prediction target and check whether enough independent biological samples exist.
If classification would merely learn which laboratory produced a sample, or if the sample size makes validation unreliable, say so and reject the option.
For every candidate model, show the link:
Published paper → existing method → demonstrated limitation → proposed improvement → appropriate data → baseline → evaluation metric → external benchmark.
We need the reason for selecting a model before implementing it.
Task 5 — Preliminary external dataset screen
Examine the public records for GSE245069, GSE277772, GSE261893, GSE261725 and any other highly relevant dataset identified from the literature.
Check whether processed count matrices, group labels and suitable controls exist.
Pay attention to shared controls, overlapping samples, differences in polymers, particle sizes, exposure duration, animal species and diet.
Do not pool studies simply because all involve microplastics.
At this stage, a feasibility screen is enough. We'll do detailed dataset acquisition and harmonization in Phase 2.
Task 6 — Incorporate our BT3041 course
We have covered hypothesis testing, regression, model selection, bootstrapping, cross-validation, classification, clustering and PCA.
Identify where these techniques could contribute to a meaningful model or evaluation procedure.
Don't force every technique into the project.
We need the scientific question to guide the method, rather than the other way around.
Deliverables
Create a new folder `research/roadmap_phase1/` and prepare:

1. `existing_work_inventory.md` — what we've done, with important limitations.
2. `paper_methods_matrix.csv` — detailed comparison of published studies and modelling methods.
3. `dataset_feasibility_screen.csv` — candidate datasets and their suitability.
4. `candidate_model_tracks.md` — at most three coherent research/model directions, with baselines, strengths, weaknesses and evaluation procedures.
5. `phase1_decision_memo.md` — a concise recommendation, no more than two pages, explaining what the team should decide next.
6. `references.md` — verified primary sources and specific figure/methods references.

The final decision memo should answer in simple language:
What problem should we investigate? Why is it worth studying? What have researchers already done? What could we improve? What would we compare against? Which data would let us establish whether our approach is actually better?
Please clearly mark uncertain claims and missing sources.
Our team will own the research questions, biological hypotheses, interpretations and final scientific decisions. Treat your suggested directions as candidates for us to evaluate, not conclusions.
Don't change or overwrite existing analyses. Don't start a new model or proceed to Phase 2 automatically.
Once you've completed this phase, give me a concise findings summary, the file locations, your strongest two or three candidate research directions, and the decisions we need to make as a team.
I want us to finish this foundation properly before we start the next phase.

*[7 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T15:27:12.106Z)

I've read the final report and the cross-study replication outputs, which exist only locally and aren't on GitHub. Next I'm reading the earlier literature review and dataset audit.

*[9 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T15:29:15.512Z)

Meanwhile, reading the remaining Phase 2 reports for the inventory.

*[3 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T15:29:41.380Z)

I have enough for the inventory, so I'll write it while the agents run.

*[3 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T15:30:38.054Z)

The existing-work inventory is written. Five research agents are still running in the background, verifying papers, methods and the GEO datasets. I'll write the remaining five deliverables once their results arrive.

### Claude (claude-opus-5-5, 2026-11-21T15:44:04.898Z)

The methods agent has finished. Still waiting on the four literature and dataset agents.

### Claude (claude-opus-5-5, 2026-11-21T15:50:21.937Z)

The core-paper agent finished, and two findings change our earlier claims. Jung's Table S7 already contains a DESeq2 interaction analysis that the paper never mentions, so "nobody tested the interaction" is no longer true. Jung also reports a significant Ppara rise with PE under control diet, while our analysis found none. I'm updating the inventory now.

*[1 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T15:50:29.686Z)

The inventory now reflects both corrections. Three agents are still running: the dataset screen and two literature reviews.

### Claude (claude-opus-5-5, 2026-11-21T15:51:02.740Z)

The independent-transcriptomes agent finished with partial file output. Its message has what I need:

- **GSE261893:** a polyethylene-diet study with only 3 mice per group and no linked paper.
- **GSE261725:** the matching polystyrene series from the same lab.
- **GSE277772:** the data behind Chiu 2025, also 3 vs 3.

Still waiting on the dataset screen and the diet-interaction papers.

### Claude (claude-opus-5-5, 2026-11-21T15:52:36.358Z)

The diet-interaction agent is done; only the dataset screen is still running. Before writing the methods matrix I'm pulling exact figure and table locations from the core-paper notes.

*[3 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T15:53:05.902Z)

I have enough to write the methods matrix now, generating it with a script so the CSV quoting is safe.

*[2 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T15:55:37.658Z)

The methods matrix is written, with 18 rows: 13 microplastic studies and 5 methods papers. Next, `references.md`.

*[1 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T15:56:50.610Z)

`references.md` is written. I'll check whether the dataset agent has written anything yet before drafting the model tracks, which depend on which MASH reference datasets exist.

*[3 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T15:58:01.056Z)

The dataset screen is written. Now the main analytical document: the research gaps (Task 3), then the three candidate tracks with their full evidence chains (Task 4), and the course mapping (Task 6).

*[1 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T15:59:20.338Z)

Now the decision memo, kept under two pages.

*[4 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T16:00:27.602Z)

Phase 1 is done. The six files are in `research/roadmap_phase1/`, with the verification notes from the five research agents in `working_notes/`. I didn't change, re-run or commit anything existing, and I haven't started Phase 2.

## Findings summary

- **Your earlier work is on disk, but not all on GitHub.** The Phase 2 branch is pushed. `research/` (the literature review and the cross-study replication) and `presentation/` exist only on this laptop. The replication folder is complete: scripts, outputs and a record of the run.
- **Two corrections from reading Jung 2026's supplement:**
  - **An interaction was computed after all.** Their Table S7 is a DESeq2 interaction analysis with 76 significant genes, never mentioned in the paper and mostly skeletal-muscle genes. That's the same artefact we found. So "nobody tested the interaction" has to become "nobody reported or interpreted it".
  - **Our Ppara result doesn't match theirs on the same counts.** Their tables give +0.53 log2 (padj 0.011) for PE under control diet; ours gives +0.16 (FDR 0.62). The difference is unexplained, so we should resolve it before making any Ppara claim.
- **The literature disagrees on direction.** Fat-synthesis genes go down with plastic in Jung, Makdissi and Chiu, but up in Djouina, Okamura and Xu P. PPARα is suppressed in the one study that fasted its mice overnight.
- **No other public dataset combines microplastics with a MASH diet in liver.** The only other polyethylene liver dataset, GSE261893, has 3 vs 3 mice over 71 weeks and no paper. Its controls are byte-identical to those in three sibling series from the same lab.
- **Useful find: GSE316023 uses the same MASH diet as Jung (D09100310).** That gives us an independent MASH reference. GSE246328 is a second one.
- **Uncertain points are marked.** Chiu 2025 and Xu 2024 were read as abstracts only. Two checks disagreed on whether GSE277772 includes raw counts. No fasting or time-of-day reference dataset has been searched for yet.

## Files

- [existing_work_inventory.md](existing_work_inventory.md)
- [paper_methods_matrix.csv](paper_methods_matrix.csv): 13 microplastic studies and 5 methods papers
- [dataset_feasibility_screen.csv](dataset_feasibility_screen.csv): 14 entries
- [candidate_model_tracks.md](candidate_model_tracks.md): the research gaps, the three tracks with their full evidence chains, and the course mapping
- [phase1_decision_memo.md](phase1_decision_memo.md)
- [references.md](references.md)

## Strongest candidate directions

1. **Track A (primary).** Use regression and model selection to ask whether the PE response is closer to independent MASH, PPARα or fasting/clock reference signatures. It has to pass a positive control (our own diet effect must line up with the external MASH axis), plus a depth-downsampling control and a muscle-contamination control. It doesn't solve the problem that PE libraries were sequenced more shallowly; no analysis can.
2. **Track B (complementary).** Re-analyse the roughly four independent labs' data with one pipeline and compare whole programmes (fat synthesis, cholesterol, PPARα, inflammation) rather than DEG lists or hub genes, including a check of each dataset's sequencing depth against treatment. This can only be descriptive, because the studies are too few, too small and too different to pool.
3. **Track C (machine learning).** Exposure classifiers are rejected: they would learn sequencing depth or which lab produced the sample, and the sample sizes are far too small. The one sensible version, a MASH score trained on external data, is really part of Track A.

## Decisions for the team

- **Deadline:** the draft report says it's due 2026-11-25. Is this restart for that submission or for later work?
- **Track A:** adopt it as primary, and choose which reference signatures to use.
- **Track B:** include it in full, or only a small direction-of-effect check?
- **Track C:** formally drop exposure classification as a model aim?
- **Novelty claim:** reword it to "first reported and interpreted diet × plastic interaction, with contamination and depth checks".
- **Ppara:** resolve the gap between Jung's +0.53 and our +0.16 first.
- **Housekeeping:** commit `research/` and `presentation/`, get the Chiu and Xu full texts through the library, and decide whether to email Jung's authors about sacrifice times, food intake and the spatial-data slide layout.

### User (2026-11-21T16:05:28.421Z)

Yes, new branch

*[3 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T16:05:45.933Z)

The only hits were the word "secretion", so there are no secrets. `literature.zip` is just an archive of `research/literature/`, so I'll leave it out of the commit as a duplicate.

*[2 tool calls]*

### Claude (claude-opus-5-5, 2026-11-21T16:06:00.765Z)

Everything is committed on a new branch, `research/phase1-roadmap` (commit `c285f6d`), which starts from `phase2/factorial-analysis`. That covers `research/literature/`, `research/cross_study_replication/`, `research/roadmap_phase1/` (including the working notes) and `presentation/`.

- **One file left out:** `research/literature.zip` is still uncommitted. It's just a zip of `research/literature/`, which is now committed, so committing it too would add nothing.
- **Not pushed:** the branch exists only on this laptop, so it isn't backed up on GitHub yet.
- **Nothing to worry about:** Git printed line-ending warnings (LF → CRLF) when committing; they don't affect the files.

### User (2026-11-21T16:07:18.772Z)

Keep local

### Claude (claude-opus-5-5, 2026-11-21T16:07:23.057Z)

I'll leave it there: `research/phase1-roadmap` stays on this laptop and isn't pushed. Until it is, the commit exists only here, so it isn't backed up if the laptop fails. To publish it later:

```powershell
git push -u origin research/phase1-roadmap
```
