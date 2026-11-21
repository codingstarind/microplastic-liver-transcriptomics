# Cross-study gene replication (BT3041)

Do published microplastic liver gene findings show up in GSE294180? This directory checks
eight published genes against the existing DESeq2 factorial analysis. It does not change any
earlier result. Read [`replication_report.md`](replication_report.md) first.

| Panel | Genes | Source |
|---|---|---|
| Primary (pre-specified) | Fasn, Acly, Acaca, Ppara | Djouina et al. 2023, polyethylene, qPCR |
| Exploratory | Acot3, Abcc3, Nr1i3 | Li et al. 2025, polystyrene, RNA-seq + qPCR |
| Exploratory | Nr4a1 | Chiu et al. 2025, polystyrene (abstract only); supporting qPCR study |

## Files

| File | What it holds |
|---|---|
| `published_gene_evidence.csv` | Verified literature claims, with design context and how far each was verified |
| `gene_replication_results.csv` | Per gene and contrast: log2FC, 95% CI, p, genome-wide BH and panel Holm adjusted p, expression, direction agreement, category |
| `interaction_results.csv` | Diet × PE estimates (Wald, LRT, OLS on VST, bootstrap, muscle-adjusted), additive-vs-factorial comparison |
| `leave_one_out_results.csv` | 19 omitted mice × 8 genes × 3 contrasts |
| `leave_one_out_summary.csv` | How often each gene keeps its direction and significance; most influential mouse |
| `covariate_sensitivity.csv` | Primary fit vs the existing depth-, muscle-, clock-, yield-adjusted and drop-MD_1 fits |
| `bootstrap_summary.csv` | Exploratory stratified bootstrap intervals |
| `clustering_comparison.csv` | Gene clusters under four distance measures and two linkages |
| `panel_pca_summary.csv` | Variance explained, share attributable to diet / PE / individual mice, loadings |
| `gene_replication_forest.png`, `gene_expression_heatmap.png`, `gene_panel_pca.png` | Figures |
| `intermediate/` | Raw outputs of the R step (additive model, LRT, OLS, bootstrap draws, correlations) |
| `scripts/` | Code (below) |

## Reproduce

The R step needs the project's WSL environment (`mpliver-r`, DESeq2 1.46). It reads the
primary fit `results/deseq2/provisional_deposited/dds.rds`, so that fit must exist first
(`analysis/deseq2/run_deseq2_factorial.R`). The Python steps use the project `.venv`. The
leave-one-out and covariate-adjusted refits are the existing ones in `results/deseq2/`
(made by `analysis/robustness/run_sensitivity.sh`); they are read here, not refitted.

From the repository root, in Windows:

```powershell
.venv\Scripts\python research/cross_study_replication/scripts/00_write_published_evidence.py
```

Inside WSL (Ubuntu-24.04), from the repository root (about 25 minutes, mostly the 2000 bootstrap draws):

```bash
MAMBA_ROOT_PREFIX=~/mpliver-tools/mamba ~/mpliver-tools/bin/micromamba run -n mpliver-r \
  Rscript research/cross_study_replication/scripts/01_fit_panel_models.R "$(pwd)" 2000
```

Back in Windows:

```powershell
.venv\Scripts\python research/cross_study_replication/scripts/02_build_tables_and_figures.py
```

`00` is a curated record. Its values were checked by hand against the papers, and the
script only writes them out. The bootstrap uses a fixed seed (20261121).
