# Phase 4: model comparison and validation

**Start with:** `phase4_decision_report.md`.

## Protocol documents (locked before any fit; hashes in `results/protocol_lock.json`; test-checked)
1. `research_question_and_prediction_targets.md`: task A (gene-effect prediction) vs task B (demonstration classifier); rows, labels, features, independent units; the feasibility gate.
2. `literature_model_justification.md`: why each model and validation choice was made, with original methodological citations.
3. `evaluation_protocol.md`: grouped nested CV, leakage rules, metrics, decision rules, variants and controls.
4. `model_search_spaces.csv`: bounded hyperparameter grids (≤ 16 configurations per model) with justifications.

## Scripts (`scripts/`)

| Script | Role |
|---|---|
| `p4_data.py` | Gene feature tables (read-only from Phases 2–3), targets, co-expression modules (Ward, 20 clusters, GSE117134) and chromosome groups, frozen outer splits, training-only τ² |
| `p4_models.py` | Model factory: pipelines and grids |
| `p4_run.py` | Parallel nested-CV runner (joblib/loky, 19 workers, single-threaded tasks) |
| `p4_alt_universe_refs.R` | Refits the references on the alternative muscle universe (variant V6) |
| `p4_depth_thinning.R` | Depth-only negative-control targets (binomial thinning) |
| `p4_classify.py` | Task B demonstration |
| `p4_analyze.py` | Rank tables, Nadeau–Bengio paired intervals, robustness counts, figures |

**Run order:**
```powershell
python p4_run.py --jobs 19 --only primary,positive
```
```powershell
python p4_run.py --jobs 19 --only neg,variants,loo
```
```powershell
python p4_classify.py
```
```powershell
python p4_classify.py hgb
```
```powershell
python p4_analyze.py
```

The two R scripts run first, via WSL as in Phase 3.

## Results
- **`results/taskA/`:**
  - `metrics_*.csv`: one row per (run, model, split), with held-out metrics and chosen hyperparameters;
  - `predictions/`: fold-by-fold held-out predictions;
  - `rank_table.csv`, `paired_comparisons.csv`, `robustness_counts.csv`, `variants_r2_oos.csv`, `negative_controls_r2_oos.csv`, `coefficient_stability.csv`, `importance_summary.csv`, `chosen_hyperparameters.csv`.
- **`results/taskB/`:** internal CV, external GSE246328, descriptive GSE294180 vehicle mice; plus the HGB fix files and the preserved original degenerate run.

## Tests
`tests/test_phase4.py`:
- protocol lock;
- no gene in both training and test;
- groups never straddle;
- each gene tested once per repeat;
- inner folds nested and group-disjoint;
- no forbidden or target-derived features;
- training-only τ² and scaler;
- single-axis models use one column;
- results completeness.
