# Reproducibility: environments, versions and run order

## 1. Environments and versions

Versions were read from the installed environments on 2026-11-25. They agree with the versions written into the analyses' own `run_info.json` files.

| Environment | Used for | Key versions | Full record |
|---|---|---|---|
| Python virtual environment `.venv` (Windows 11) | All Python analyses except CatBoost fitting | Python 3.14.6; numpy 2.5.1; pandas 3.0.5; scipy 1.18.1; scikit-learn 1.9.1; matplotlib 3.11.1; seaborn 0.13.2; joblib 1.6.0; pytest 9.1.1 | [`python_venv_pip_freeze.txt`](environment_records/python_venv_pip_freeze.txt) |
| Second Python interpreter `C:\Python314\python.exe` | Steps that fit CatBoost (Phase 4 and 4B record "Python 3.14, scikit-learn 1.9.1, CatBoost 1.2.10"; interaction follow-up `ix01`). `.venv` has no CatBoost. | Python 3.14.6; catboost 1.2.10; same numpy, pandas, scipy and scikit-learn versions as `.venv` | [`python314_catboost_interpreter_packages.txt`](environment_records/python314_catboost_interpreter_packages.txt) |
| Conda env `mpliver-r` in WSL Ubuntu 24.04 ([`envs/r-analysis.yml`](../../envs/r-analysis.yml)) | DESeq2, edgeR, limma (fry, camera, voom), fgsea, tximport | R 4.4.3; DESeq2 1.46.0; limma 3.62.1; edgeR 4.4.0; fgsea 1.32.2; msigdbr 26.1.1; ashr 2.2.63; apeglm 1.28.0; tximport 1.34.0; BiocParallel 1.40.0 | [`r_and_quant_tool_versions.txt`](environment_records/r_and_quant_tool_versions.txt), [`conda_env_mpliver-r_packages.txt`](environment_records/conda_env_mpliver-r_packages.txt) |
| Conda env `mpliver-quant` in WSL ([`envs/quantification.yml`](../../envs/quantification.yml)) | Re-quantification pilot | salmon 1.10.3; fastp 0.23.4 | [`conda_env_mpliver-quant_packages.txt`](environment_records/conda_env_mpliver-quant_packages.txt) |

`requirements.txt` pins the same Python versions.

**Random seeds and settings** are fixed in each script and recorded in its folder's `README.md`, plan or `run_info*.json`. For example:
- Phase 4 outer splits use seeds 101–103.
- The follow-up analyses use seed 20261125.
- Locked plans and hashes are in `research/roadmap_phase4/results/protocol_lock.json`, `research/roadmap_phase4e/results/plan_lock.json` and `research/roadmap_phase5/results/plan_lock.json`.

## 2. Data

Raw and large prepared data are **not committed**. They are downloaded by the scripts:
- **GSE294180:** `scripts/01_download_data.py` downloads the files and verifies them against `metadata/data_provenance.json`.
- **Reference studies:** `research/roadmap_phase2/scripts/01_…08_*.py` download and prepare them in `D:/mpliver_work/phase2_refs/`. The location is set in `research/roadmap_phase2/scripts/p2_common.py` and can be overridden with the environment variable `MPLIVER_P2_DATA`.

NCBI sometimes blocks scripted downloads. If that happens, download the listed file in a browser and place it at the path the script prints.

## 3. Run order

Run every command from the repository root. Each folder's `README.md` gives the exact commands, inputs and outputs. All analysis steps after Phase 2 read earlier outputs read-only.

| Step | Folder | What it produces | Report section |
|---|---|---|---|
| 1 | `scripts/01–04` | Download, sample metadata, count audit, filtering and normalisation | §2.2 |
| 2 | `scripts/05–06` | EDA, PCA, clustering, muscle-gene list | §3.1, §4.1 |
| 3 | `analysis/deseq2/run_deseq2_factorial.R` (WSL) | DESeq2 factorial contrasts | §3.1, §4.1 |
| 4 | `analysis/pathways/pathway_analysis.R` (WSL) | fry, camera and fgsea for 50 Hallmark + 3 PPAR sets | §3.1, Fig. 2 |
| 5 | `analysis/concordance`, `analysis/interaction`, `analysis/robustness` | Concordance, interaction patterns, sensitivity variants | §2.2, §5.3 |
| 6 | `pipelines/quantification/` (WSL) | Re-quantification pilot (optional; 5 samples) | §2.2 |
| 7 | `research/roadmap_phase2/` | Reference datasets, gene universe | §2.1 |
| 8 | `research/roadmap_phase3/` (Python + WSL R) | Reference signatures: MASH, FAST, CLOCK, PPARα | §2.1, §3.4 |
| 9 | `research/roadmap_phase4/` | F4/F8 regression with OLS, Ridge, Elastic Net, CatBoost and HGB; module splits | §3.4, §3.5 |
| 10 | `research/roadmap_phase4b/`, `research/roadmap_phase4c_diagnostics/` | Blocks A/B/C, FULL (28 features), ablation and diagnostics | §3.4, Table 4 |
| 11 | `research/roadmap_phase4d/` | Per-mouse signature scores (singscore) | §3.2 |
| 12 | `research/roadmap_phase4e/` | Cross-study MASH-state classifier, frozen scores on Main PE | §3.3, Table 3 |
| 13 | `research/roadmap_phase5/`, `research/roadmap_phase5b/` | Hong-PE transfer and external validation | §3.5, §4.2 |
| 14 | `research/enrichment_audit/` | Gene-mapping audit and fry/camera sensitivity | §3.1 |
| 15 | `research/pathway_followup/` | Eight-pathway scores and tests; classifier pathway ablation | §3.2, Table 2, §4.2 |
| 16 | `research/lipid_prediction_followup/` | Lipid and cholesterol gene groups vs other genes | Table 4 |
| 17 | `research/nonlinear_chol_followup/` | CatBoost and HGB vs OLS for cholesterol-synthesis genes | Table 5 |
| 18 | `research/interaction_followup/` (`ix01` uses `C:\Python314\python.exe`) | Interaction diagnostics, pair-only HGB | §4.2 |
| 19 | `research/section3_1_figures/` | Report Figures 1 and 2 | Fig. 1–2 |
| 20 | `supplementary/reproducibility/verify_report_numbers.py` | Checks the report's numbers against steps 1–19 | all |

Not part of the report:
- `research/roadmap_phase5c/` (paused, nothing scored).
- `research/cross_study_replication/` (published-gene replication).
- `research/paper_evidence_audit/` (evidence index used while writing).

**Checking the report's numbers.** This script reads the saved outputs and checks the values printed in the final report. It does not re-run any model. On 2026-11-25 all 75 values matched; the results are in [`../tables/report_number_crosswalk.csv`](../tables/report_number_crosswalk.csv).

```powershell
.venv\Scripts\python.exe supplementary\reproducibility\verify_report_numbers.py
```

**Tests.**

```powershell
.venv\Scripts\python.exe -m pytest
```

## 4. Run times

Run times are given in each folder's README where they were recorded. The follow-up analyses take between 10 s and 11 min each.
