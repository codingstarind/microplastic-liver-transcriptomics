# AI Usage Declaration

**Project:** Polyethylene microplastic exposure and hepatic transcriptional programs (re-analysis of GEO GSE294180)
**Course:** BT3041 – Analysis and Interpretation of Biological Data, IIT Madras
**Repository:** https://github.com/codingstarind/microplastic-liver-transcriptomics
**Prepared:** 2026-11-02. Drafted with Claude Code from the repository and the recorded Claude Code session, for the team to check and complete.

> **How to read this document.** Statements are based on the repository contents, its commit
> history and the recorded Claude Code session. Anything those records cannot confirm, in
> particular independent student review of results, is marked **[TEAM TO CONFIRM]**. The team
> must complete or correct these items before submission. The BT3041 AI-use guidelines were not
> available when this was drafted, so **no claim of compliance with them is made here**.

---

## 1. Summary

AI was used **extensively** in this project, mainly **Claude Code** (Anthropic; model Claude Opus 5.5, run through the Claude Code desktop application). It:

- wrote, ran and debugged nearly all of the project's code;
- produced the figures;
- ran the statistical analyses;
- drafted most of the written reports, including the final report draft.

The student team set the research problem and the overall scope of each phase, made the recorded methodological decisions listed in Section 4, and retains responsibility for the scientific claims submitted.

AI-generated outputs, including p-values, figures, interpretations and draft conclusions, are **computational results and drafts**. They do not become accepted scientific claims until the team has reviewed and approved them.

## 2. Tools used

| Tool | How it was used | Evidence |
|---|---|---|
| **Claude Code** (Claude Opus 5.5) | Code writing and execution; data acquisition and audit; statistical analysis; figures; repository management (commits, branches, pull requests); drafting reports | Recorded session. All 15 commits in the repository were created through Claude Code under the team's Git identity: 14 carry a `Co-Authored-By: Claude Opus 5.5` line, and the 15th is the PR #1 merge, which Claude Code performed at the team's request. |
| **ChatGPT** (OpenAI) | The team stated that it reviewed the PR #2 methodology proposal with ChatGPT before requesting corrections | Team message in the session. **[TEAM TO CONFIRM extent]** |
| Other tools (e.g. for slides, literature search, or the FASTQ download tool in `D:\mpliver_work\fastq_pilot`) | Not recorded in the Claude Code session | **[TEAM TO CONFIRM]** whether other AI tools were used, and for what |
| GitHub Copilot reviewer | Requested on PR #1 and PR #2, but produced no review (quota) | PR comments |

## 3. AI assistance by research activity

### 3.1 Research problem, dataset and scope

| Activity | Who | Detail |
|---|---|---|
| Choosing the research problem and dataset (GSE294180, Jung *et al.* 2026) | **Team** | Provided in the team's instructions; the project had already been approved |
| Phase scopes and task lists (Phase 1 audit/EDA, feasibility review, Phase 2 framework) | **Team** | Detailed written instructions from the team for each phase |
| Research directions A–C (PE–MASH similarity, diet × PE interaction, robustness) | **AI proposed, team selected** | Candidate questions Q1–Q4 were first drafted by Claude Code in `reports/preliminary_observations.md` (Phase 1). The team then set Directions A–C as the Phase 2 framework. |
| Formal hypotheses | **Not yet adopted** | The candidate hypotheses H1 and H2 in `reports/phase2/interpretation.md` and the final report draft were **suggested by Claude Code**. **[TEAM TO CONFIRM]** which hypotheses the team adopts. |

### 3.2 Literature and background

| Activity | Who | Detail |
|---|---|---|
| Reading the GEO records and the source paper; summarising design, methods and published findings | AI | Claude Code retrieved the paper's full text (Europe PMC) and GEO records, and wrote `reports/dataset_understanding.md` and `reports/phase2/literature_comparison.md` |
| Checking which analyses the paper already reported (reproduction vs extension) | AI drafted | Based on the paper's main text only; the supplementary tables were not obtained |
| Searching for external MASH datasets and gene sets | AI | GEO/MSigDB queries; candidate list in `reports/phase2/pe_mash_similarity_report.md` |
| Explaining statistical and biological concepts (e.g. tximport offsets, shared-reference bias, Cook's distance) | AI | In reports and session replies |
| Broader literature search on microplastics and liver | **Not done** | Flagged as outstanding in the final report draft. **[TEAM TO CONFIRM]** any literature work done outside this session. |
| Assessing source relevance and forming conclusions from the literature | **Team responsibility** | **[TEAM TO CONFIRM]** that the team has checked the cited sources |

### 3.3 Data acquisition, audit and preprocessing

| Activity | Who | Detail |
|---|---|---|
| Repository setup (structure, environment, tests, CI) | AI | Created by Claude Code at the team's request |
| Downloading GEO files | **Team (manual) + AI** | NCBI blocked scripted downloads, so the team downloaded the count matrix and SOFT file in a browser; Claude Code verified and registered them with checksums |
| Downloading reference annotation (GENCODE) | AI | Claude Code |
| Downloading pilot FASTQs (6 runs) | **Team** | The team downloaded the files and placed them on the machine; Claude Code verified them against ENA MD5 checksums. The download tool used is **[TEAM TO CONFIRM]**. |
| Data audit, sample metadata, data dictionary | AI | Code and reports by Claude Code |
| Preprocessing choices (filtering, normalisation, transformation) | AI proposed and implemented | Reasoning in `reports/preprocessing_decisions.md`. **[TEAM TO CONFIRM]** that the team reviewed and approved these choices. |

### 3.4 Programming and computational implementation

| Activity | Who | Detail |
|---|---|---|
| Python and R code (`src/`, `scripts/`, `analysis/`, `pipelines/`) | **AI wrote essentially all code** | Claude Code generated, ran and revised the code |
| Choice of statistical methods (DESeq2 factorial model, limma fry/camera, fgsea, split-reference concordance, sensitivity analyses) | Team set requirements; AI designed specifics | The team's Phase 2 instructions asked for the DESeq2 factorial model, reference levels, contrast types, pathway themes, a shared-reference-aware concordance and specific robustness checks. Claude Code chose the specific implementations and wrote the pre-specified plan (`reports/phase2/analysis_plan.md`). |
| Debugging and dependency problems | AI | e.g. tool installation in WSL, a fastp–Salmon deadlock, a pandas index bug |
| Computational execution (all analyses, including Salmon re-quantification) | AI | Run by Claude Code on the team's machine |

### 3.5 Verification and robustness

| Activity | Who | Detail |
|---|---|---|
| Automated tests (pytest), validators, CI | AI | 33 unit tests plus metadata and methodology validators |
| Data-integrity, checksum and metadata consistency checks | AI | Recorded in reports and provenance files |
| Re-quantification pilot design and pre-registered decision rules | AI drafted | `reports/quantification_validation/pilot_protocol.md`; the team chose which samples to download |
| Sensitivity and robustness analyses; comparing statistical implementations | AI | `reports/phase2/robustness_report.md` |
| Errors found and corrected by AI | AI | e.g. a gene-ID check bug, a pathway-covariate bug, an invalid length-offset approach (recorded as a negative result) |
| **Independent scientific validation of results** | **Not yet documented** | Automated tests and AI-generated verification reports are **not** independent scientific validation. No record shows the team independently re-checking the Phase 2 results. **[TEAM TO CONFIRM]** what the team has reviewed. |

### 3.6 Figures and visualisation

| Activity | Who | Detail |
|---|---|---|
| Which figure types to produce | Team specified categories; AI chose specifics | The Phase 1 instructions requested histograms, boxplots, library-size, PCA and clustering figures. The Phase 2 figures (fig16–18) and diagnostics were chosen by Claude Code. |
| Plotting code, layout, colours, labels | AI | All figures were generated by AI-written code; a colour-vision-checked palette was used |
| What the figures establish | **Team responsibility** | Figures show all relevant data (e.g. contaminated samples are shown, not removed). **[TEAM TO CONFIRM]** that figure interpretations have been reviewed. |

### 3.7 Interpretation and scientific conclusions

| Activity | Who | Detail |
|---|---|---|
| Written interpretations in the reports and the final report draft | **AI drafted** | All interpretive text in `reports/` was drafted by Claude Code and labelled exploratory or draft |
| Draft conclusions (e.g. "PE mimics metabolic dysfunction is not supported") | **AI drafted** | `reports/phase2/interpretation.md`, `reports/final_report_draft.md` |
| Accepting interpretations and conclusions for submission | **Team** | **[TEAM TO CONFIRM]**. The final report draft states that it has not yet been reviewed by the team. |
| Claims about causality, statistical reliability or novelty | **Team** | Drafts deliberately avoid causal and novelty claims. The novelty check against the wider literature is outstanding. |

### 3.8 Writing, presentation and communication

| Activity | Who | Detail |
|---|---|---|
| README, reports, final report draft, this declaration | AI drafted | Written by Claude Code; the team is responsible for editing and approving them |
| Milestone slide outline (`milestone_progress.md`) | AI drafted | A 10-slide outline was proposed |
| Slide design and PowerPoint production | **Not recorded in this session** | The team states that presentation production was AI-assisted. **[TEAM TO CONFIRM]** the tools used and the extent. |

## 4. Decisions made by the team (recorded in the session)

These decisions are documented in the session and/or the repository. Other decisions may have been made outside these records.

| Date | Decision | Where recorded |
|---|---|---|
| 2026-11-02 | Project scope and requirements for Phase 1, the feasibility review and Phase 2 | Team instructions (session) |
| 2026-11-02 | Repository to be private, later made **public** | Session |
| 2026-11-02 | Merge PR #1 without a separate review ("bypass review") | Session; PR #1 |
| 2026-11-02 | Revisions to the methodology proposal (flowcell evidence, proxy naming, tximport offset, Cook's distance, power caveats, smoke-test scope, CI) | Team review of PR #2 (reported as done with ChatGPT); `reports/methodology_proposal.md` revision 2 |
| 2026-11-02 | Approve installing tools in WSL with micromamba | Session |
| 2026-11-02 | **Keep the pre-specified primary model** (no covariate) | `reports/phase2/analysis_plan.md` (Team decisions) |
| 2026-11-02 | Use the 5-sample pilot as verification; **proceed with the deposited count matrix** under assumptions A1–A3; no full re-quantification | `reports/quantification_validation/quantification_validation_report.md` §5.1 |
| 2026-11-02 | Choice of pilot samples (the team downloaded CD_1, CD_MP_3, MD_4, MD_MP_3, then CD_MP_5 on AI recommendation) | Validation report §4 |
| 2026-11-02 | Request a final report draft; keep it in the repository only | Session |

Not decided in the recorded material: **adopted hypotheses**, **acceptance of results and interpretations**, **final conclusions**. All are **[TEAM TO CONFIRM]**.

## 5. Statements for the team to confirm before submission

1. Which hypotheses the team formally adopts. The current H1 and H2 are AI suggestions.
2. Which results, figures and interpretations the team has independently reviewed and accepts.
3. Whether the preprocessing and methodological choices in the reports have been reviewed and approved by the team.
4. Which other AI tools were used (ChatGPT for the PR #2 review; any tools for slides, literature or the FASTQ download script) and for what.
5. Which literature the team assessed beyond the source paper.
6. Whether the BT3041 AI-use rules permit this usage. The guidelines were not checked when this was drafted.
7. Author names and contribution statements for the final report.

## 6. Short declaration for the final report (≈120 words)

> We used AI tools extensively in this project, mainly Claude Code (Anthropic, Claude Opus 5.5). Claude Code wrote and ran nearly all of the project's Python and R code. That covered data acquisition and auditing, preprocessing, DESeq2 factorial modelling, pathway, concordance and robustness analyses, automated tests, and all figures. It also drafted the project reports and this report's first draft. [TEAM TO CONFIRM: other tools, e.g. ChatGPT for reviewing the methodology proposal and AI assistance for slides.] AI-generated results, interpretations and candidate hypotheses were treated as drafts, not accepted findings. The research problem, project scope and key methodological decisions were set by the team, and the team takes responsibility for the hypotheses, interpretations and conclusions submitted. Full details and items still under review: `reports/ai_usage_declaration.md` in the project repository.

## 7. Presentation disclosure (last slide)

> AI tools (mainly Claude Code) were used extensively for code, analysis, figures and drafting; scientific decisions and conclusions are the team's. See the AI Assistance section of the repository README for the full disclosure.
