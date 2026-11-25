# AI usage record

**Project:** An Analysis of Microplastic Effects on the Liver (BT3041, IIT Madras; re-analysis of GEO GSE294180)
**Repository:** https://github.com/codingstarind/microplastic-liver-transcriptomics
**Prepared:** 2026-11-25, by Claude Code at the team's request, from the records listed in Section 1.

This file records which AI tools were used, which model versions the records show, what each tool was used for, and which software libraries and versions each task used. Where something cannot be verified from the records, it is marked **[TEAM TO CONFIRM]** instead of guessed. The earlier, more detailed declaration written during Phase 1 and Phase 2 is kept unchanged at [`reports/ai_usage_declaration.md`](../reports/ai_usage_declaration.md).

## 1. Records used

| Record | What it shows | Where |
|---|---|---|
| Claude Code session transcripts (JSONL, kept on the analysis machine) | Every prompt, every reply, and the model ID attached to each reply | Summarised in [`chat_histories/README.md`](chat_histories/README.md); prompts in [`prompts/`](prompts/); readable transcripts in [`chat_histories/`](chat_histories/) |
| Claude Code web link supplied by the team | The shared Claude session | https://claude.ai/code/cshk_01ED6AbAq8wfpf2wvCJVoiqv69jadvzqiRCTFmVueLHH45 |
| Git history | Commits from 2026-11-02 to 2026-11-25, made through Claude Code under the team's Git identity, carrying `Co-Authored-By: Claude Opus 5.5` lines | `git log` |
| Installed-environment records | Exact package versions, read from the environments on 2026-11-25 | [`reproducibility/environment_records/`](reproducibility/environment_records/) |
| Run-info files written by the analyses | Versions logged at run time (for example `"DESeq2": "1.46.0"`, `"sklearn": "1.9.1"`) | `results/**/run_info.json`, `research/**/run_info*.json` |

## 2. AI tools and model versions

| Tool | Version shown by the records | Used for | Evidence |
|---|---|---|---|
| **Claude Code** (Anthropic), desktop app, Code tab | App versions 2.1.309 to 2.1.330 (the `version` field of the transcripts) | Running all of the sessions listed below | Transcripts |
| **Claude Opus 5.5** (`claude-opus-5-5`) | Model ID on the replies of all 7 project sessions | Main model in every session: code, analysis runs, figures, literature summaries, drafting, repository management | Transcripts (per-session counts in [`chat_histories/README.md`](chat_histories/README.md)) |
| **Claude Opus 5** (`claude-opus-5`) | Model ID on 77 replies in the main literature-review conversation, and in sub-agents of three sessions | Literature review (2026-11-20) and sub-agent searches during research planning and published-gene replication | Transcripts of sessions `9228d9bf`, `5dbfd995` and `eb0c5c48` |
| **Claude Haiku 5.5** (`claude-haiku-5-5`) | Model ID in sub-agent transcripts | Short file-search and exploration sub-tasks started by Claude Code | Sub-agent transcripts of sessions `949232bd`, `5dbfd995`, `eb0c5c48` and `59dc8d15` |
| **ChatGPT** (OpenAI) | Not recorded | The team said it reviewed the Phase 2 methodology proposal (PR #2) with ChatGPT | Team message in session `949232bd`. **[TEAM TO CONFIRM]** model and extent |
| Tools used for the final LaTeX report's language editing, and for the slides | Not recorded in these sessions | The final report states that AI helped with "drafting and language editing, including grammar and spelling" | **[TEAM TO CONFIRM]** which tool and model were used for the final report text and for `presentation/` |
| GitHub Copilot reviewer | n/a | Requested on PR #1 and PR #2; produced no review (quota) | PR pages |

## 3. Use by task

The type of assistance is one of the following:
- **Code**: writing, running and debugging code.
- **Literature**: searching and summarising sources.
- **Analysis discussion**: proposing or explaining methods, checks and decisions.
- **Drafting**: writing reports or documentation.
- **Editing**: language and format changes to text the team wrote.

| # | Task (report section) | Dates | AI tool and model | Type of assistance | Libraries and installed versions | Session |
|---|---|---|---|---|---|---|
| 1 | Repository set-up; GEO data acquisition; count audit; sample metadata (§2.2) | 2026-11-02 | Claude Code, Opus 5.5 | Code; Drafting (audit and dataset reports) | Python 3.14.6, pandas 3.0.5, numpy 2.5.1, requests 2.32.3 | `949232bd` |
| 2 | Filtering, normalisation, EDA, PCA, hierarchical clustering (§2.2, §3.1, Fig. 1) | 2026-11-02, figure redrawn 2026-11-25 | Claude Code, Opus 5.5 | Code; Analysis discussion (preprocessing choices) | scikit-learn 1.9.1, scipy 1.18.1, matplotlib 3.11.1, seaborn 0.13.2 | `949232bd`, `59dc8d15` |
| 3 | Methodology feasibility proposal | 2026-11-02 | Claude Code, Opus 5.5; ChatGPT (team review) | Analysis discussion; Drafting | pydeseq2 0.5.4 (smoke test only) | `949232bd` |
| 4 | Re-quantification pilot (fastp → Salmon → tximport) | 2026-11-02 | Claude Code, Opus 5.5 | Code; Analysis discussion | salmon 1.10.3, fastp 0.23.4, tximport 1.34.0, R 4.4.3 | `949232bd`, `03c789f5` |
| 5 | DESeq2 factorial model, edgeR and limma-voom comparison (§3.1, §4.1) | 2026-11-02 | Claude Code, Opus 5.5 | Code | DESeq2 1.46.0, edgeR 4.4.0, limma 3.62.1, ashr 2.2.63, apeglm 1.28.0 | `03c789f5` |
| 6 | Pathway enrichment: directional fry and camera, plus fgsea (§3.1, Fig. 2) | 2026-11-02; audit 2026-11-24; figure 2026-11-25 | Claude Code, Opus 5.5 | Code | limma 3.62.1, fgsea 1.32.2, msigdbr 26.1.1 (MSigDB 2026.1.Mm) | `03c789f5`, `59dc8d15` |
| 7 | Robustness checks: leave-one-mouse-out; muscle, depth and clock covariates | 2026-11-02 | Claude Code, Opus 5.5 | Code | DESeq2 1.46.0, limma 3.62.1 | `03c789f5` |
| 8 | Literature review and research-gap audit; reference list (§1) | 2026-11-20 | Claude Code, Opus 5.5 and Opus 5 | Literature; Drafting (`research/literature/`) | none | `9228d9bf` |
| 9 | Replication of published microplastic liver genes | 2026-11-21 | Claude Code, Opus 5.5 and Opus 5 (sub-agents) | Code; Literature | DESeq2 1.46.0 | `eb0c5c48` |
| 10 | Research strategy and Phase 1 roadmap | 2026-11-21 | Claude Code, Opus 5.5 and Opus 5 (sub-agents) | Analysis discussion; Literature | none | `5dbfd995` |
| 11 | Reference datasets and signatures: MASH-1/2, FAST, CLOCK, PPARα; gene universe (§2.1) | 2026-11-22 | Claude Code, Opus 5.5 | Code | DESeq2 1.46.0, limma 3.62.1, pandas 3.0.5 | `59dc8d15` |
| 12 | Gene-level regression with F4, F8 and Blocks A/B/C: OLS, Ridge, Elastic Net (§3.4, Table 4) | 2026-11-22 | Claude Code, Opus 5.5 | Code | scikit-learn 1.9.1, joblib 1.6.0 | `59dc8d15` |
| 13 | CatBoost and HistGradientBoosting models (§3.5, Table 5) | 2026-11-22; follow-ups 2026-11-25 | Claude Code, Opus 5.5 | Code | catboost 1.2.10 (separate interpreter, `C:\Python314`), scikit-learn 1.9.1 | `59dc8d15` |
| 14 | Mouse-level eight-pathway scores and tests (§3.2, Table 2) | 2026-11-22 to 11-23 (signature scores); 2026-11-25 (eight-pathway tests) | Claude Code, Opus 5.5 | Code | singscore implemented in project code; scipy 1.18.1 | `59dc8d15` |
| 15 | Cross-study MASH-state classifier and pathway ablation (§3.3, Table 3) | 2026-11-23 (classifier); 2026-11-25 (ablation) | Claude Code, Opus 5.5 | Code | scikit-learn 1.9.1 | `59dc8d15` |
| 16 | Hong-PE transfer and external validation (§3.5) | 2026-11-23 | Claude Code, Opus 5.5 | Code | scikit-learn 1.9.1, DESeq2 1.46.0 | `59dc8d15` |
| 17 | Cholesterol and lipid gene prediction; nonlinear vs linear; interaction diagnostics (§4.2) | 2026-11-25 | Claude Code, Opus 5.5 | Code | scikit-learn 1.9.1, catboost 1.2.10 | `59dc8d15` |
| 18 | Paper evidence audit, figure recommendations, paper outline | 2026-11-23 | Claude Code, Opus 5.5 | Drafting (organising evidence for the report) | none | `59dc8d15` |
| 19 | Final report draft in the repository (`reports/final_report_draft.md`, Phase 1–2 only) | 2026-11-02 | Claude Code, Opus 5.5 | Drafting | none | `949232bd` |
| 20 | Final submitted report (LaTeX PDF) | not recorded | **[TEAM TO CONFIRM]** | Per the report: AI-assisted drafting and language editing | n/a | not in these sessions |
| 21 | Repository organisation, README, this supplement, report-number check | 2026-11-25 | Claude Code, Opus 5.5 | Code; Drafting | as in task 2 | `59dc8d15` |

Full package lists:
- Python: [`python_venv_pip_freeze.txt`](reproducibility/environment_records/python_venv_pip_freeze.txt) and [`python314_catboost_interpreter_packages.txt`](reproducibility/environment_records/python314_catboost_interpreter_packages.txt).
- R/Bioconductor and quantification tools: [`r_and_quant_tool_versions.txt`](reproducibility/environment_records/r_and_quant_tool_versions.txt) and the two `conda_env_*` files.

## 4. What the team did

As stated in the final report, the students chose the research questions, reviewed the analyses and made the final decisions about interpretation. The prompts in [`prompts/`](prompts/) show the team setting each phase's scope and constraints. For example, the team repeatedly told Claude not to add interpretations or hypotheses ("our team will do that"), not to retrain frozen models, and to exclude Phase 5C from the paper. Decisions recorded in the repository are listed in [`reports/ai_usage_declaration.md`](../reports/ai_usage_declaration.md) §4.

## 5. Items the records cannot confirm

1. Which tool and model were used for language editing of the final LaTeX report and for the slides.
2. The ChatGPT model version and the extent of its use.
3. Whether any other AI tool was used outside these Claude Code sessions.
4. The repository records show that Claude Code drafted research reports, the roadmap documents and an early report draft. The team should check that the report's declaration page agrees with this record.
