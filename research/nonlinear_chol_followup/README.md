# Nonlinear vs linear cholesterol-synthesis prediction follow-up (exploratory; 2026-11-25)

This follow-up compares F8-OLS, F8-Ridge, F8-CatBoost and F8-HGB using the original Phase 4 held-out predictions, on the locked 20-module, 5-fold × 3-repeat splits. No model was refitted or tuned. Start with `REPORT.md`.

## Files

| Path | Content |
|---|---|
| `figures/fig_chol_linear_vs_nonlinear.{png,pdf,svg}` | Two panels (CD, MD): cholesterol-synthesis held-out R²_group with 95% module-bootstrap CI for the 4 models, matched-gene reference, ΔR² vs OLS with Holm p |
| `TABLE_results_concise.csv` | Cholesterol-synthesis metrics per diet × model, plus paired comparisons vs OLS and vs Ridge |
| `REPORT.md` | Results (numbers only), decision rules, controls, items not done |
| `results/group_metrics_by_model.csv` | R², RMSE, direction agreement and Spearman for cholesterol synthesis, lipogenesis, PPARα and other genes, with matched references |
| `results/paired_nonlinear_vs_linear.csv` | All gene-paired tests, Holm (primary) and BH (secondary) |
| `results/chol_dominance_controls.csv` | Leave-one-gene-out, leave-one-module-out and top-gene share |
| `results/chol_per_gene_predictions.csv` | Per-gene predictions and errors for all 4 models |
| `results/matching_balance.csv` | Balance of the matched reference genes |
| `results/reproduction_check_vs_phase4.csv` | Check against stored Phase 4 metrics |
| `results/run_info_nl01.json` | Run record |
| `logs/` | Run logs |
| `scripts/nl01_nonlinear_vs_linear.py` | Analysis (about 10 s) |
| `scripts/nl02_table_figure.py` | Table and figure |

## Reproduce

Run these from the repository root. Step 1 needs `research/lipid_prediction_followup/results/oof_predictions_with_null.csv.gz` (from lp01).

```powershell
.venv\Scripts\python.exe research\nonlinear_chol_followup\scripts\nl01_nonlinear_vs_linear.py
```
```powershell
.venv\Scripts\python.exe research\nonlinear_chol_followup\scripts\nl02_table_figure.py
```
