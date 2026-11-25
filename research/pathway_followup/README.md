# Pathway follow-up and R2 classifier ablation (exploratory, post hoc)

Computational outputs only. Earlier phase outputs were read, not modified. The Phase 4E code was imported unchanged.

## Deliverables
- **Concise results table:** `TABLE_results_concise.csv`.
- **Pathway-score figure:** `figures/fig1_pathway_scores_19_mice.{png,pdf,svg}`.
- **Classifier-ablation comparison:** `figures/fig2_classifier_ablation.{png,pdf,svg}`.
- **Evidence audit:** `EVIDENCE_AUDIT.md`.

## Scripts (in order)
1. `scripts/pf01_pathways.py`: 8 pathway scores, contrasts, gene-set overlap, correlations, per-mouse consistency, leave-one-out, Cook's distance, covariate/muscle sensitivity, external studies. Runs in about 1 min.
2. `scripts/pf02_classifier_ablation.py` with `pf02_lib.py`: Phase 4E reproduction plus 22 ablations, scaled effects and mouse bootstrap. Runs in about 1–2 min on 12 workers (env `PF_JOBS`).
3. `scripts/pf03_report.py`: builds the table and Figure 2.

Run each from the repository root with `.venv\Scripts\python.exe`.

## Inputs (read-only)
- **GSE294180:**
  - `data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz`
  - `results/deseq2/provisional_deposited/vst_matrix.csv.gz`
  - `metadata/sample_metadata.csv`
- **Earlier phases:**
  - `research/roadmap_phase3/results/gene_sets.csv`
  - `research/roadmap_phase4d/{signatures/FAST_main_n50.csv, results/per_mouse_scores.csv}`
  - `research/roadmap_phase5/{frozen/phase4e_R2_l2_logistic_params.json, per_sample_scores.csv}`
  - `research/roadmap_phase4e/scripts/e01_mash_state_models.py` (imported unchanged)
- **Other:**
  - MASH reference matrices in `D:/mpliver_work/phase2_refs`
  - Phase 5 external matrices in `D:/mpliver_work/phase5`
  - `research/enrichment_audit/tables/T2_*` (original fry q-values)
