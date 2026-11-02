"""Step 7: build and execute notebooks/01_phase1_walkthrough.ipynb.

The notebook does no new computation: it reads the tables written by steps 1-6
and links to the figures, so the Phase 1 results can be reviewed in one place.
Rebuild it after re-running the pipeline.
"""

import sys

import nbformat
from nbclient import NotebookClient

from mpliver.config import REPO_ROOT

NOTEBOOK = REPO_ROOT / "notebooks" / "01_phase1_walkthrough.ipynb"

md = nbformat.v4.new_markdown_cell
code = nbformat.v4.new_code_cell

# Figures are linked rather than embedded so the notebook stays small in git.
EDA_FIGURES = [
    "fig01_library_sizes", "fig02_count_distributions",
    "fig03_sample_boxplots_before_after_normalisation", "fig04_expression_dominance",
    "fig05_mean_variance", "fig06_sample_qc_by_group", "fig07_selected_genes_by_group",
    "fig08_variance_share_by_factor", "fig09_tissue_marker_screen",
]
PCA_FIGURES = [
    "fig10_pca_primary", "fig11_pca_top_loadings", "fig12_pca_all_filtered_genes",
    "fig13_pca_muscle_genes_removed", "fig14_sample_correlation_heatmap",
    "fig15_dendrograms_sensitivity",
]


def figure_links(folder: str, names: list[str]) -> str:
    lines = ["Figures (PNG and PDF; each one is explained in `reports/eda_report.md`):\n"]
    lines += [f"- [{n}](../figures/{folder}/{n}.png)" for n in names]
    return "\n".join(lines)


CELLS = [
    md("# Phase 1 walkthrough: GSE294180 liver RNA-seq\n\n"
       "Read-only tour of the outputs of `scripts/01`–`06`. Interpretation lives in "
       "`reports/`; this notebook shows the underlying tables.\n\n"
       "Re-generate with `python scripts/07_build_walkthrough_notebook.py`."),
    code("import json\nfrom pathlib import Path\n\nimport pandas as pd\n\n"
         "ROOT = Path.cwd().parent if Path.cwd().name == 'notebooks' else Path.cwd()\n"
         "T = ROOT / 'results' / 'tables'\n"
         "pd.set_option('display.width', 160)"),
    md("## 1. Samples (from GEO records)"),
    code("meta = pd.read_csv(ROOT / 'metadata' / 'sample_metadata.csv')\n"
         "meta[['geo_accession', 'counts_column', 'group', 'diet', 'microplastic',\n"
         "      'library_name', 'geo_batch']]"),
    code("meta.groupby(['diet', 'microplastic']).size().unstack()"),
    md("## 2. Audit of the deposited count file"),
    code("audit = json.loads((T / 'audit_summary.json').read_text())\n"
         "pd.Series(audit['integrity'])"),
    code("pd.Series(audit['value_type_evidence'])"),
    code("pd.read_csv(T / 'counts_preview.csv').iloc[:5, :6]"),
    code("pd.read_csv(T / 'per_sample_audit.csv', index_col=0)[\n"
         "    ['group', 'total_counts', 'n_genes_ge_10', 'fraction_counts_chrM',\n"
         "     'top_gene_symbol']]"),
    code("pd.read_csv(T / 'sex_marker_counts.csv', index_col=0)"),
    md("## 3. Preprocessing"),
    code("json.loads((T / 'filter_summary.json').read_text())"),
    code("pd.read_csv(T / 'size_factors.csv', index_col=0).round(3)"),
    md("## 4. Descriptive EDA"),
    md(figure_links("eda", EDA_FIGURES)),
    code("pd.read_csv(T / 'variance_share_by_factor.csv')"),
    code("pd.read_csv(T / 'tissue_marker_cpm.csv', index_col=0).round(1)"),
    md("## 5. PCA and clustering"),
    md(figure_links("pca_clustering", PCA_FIGURES)),
    code("pd.read_csv(T / 'pca_variance_explained.csv').iloc[:, :7]"),
    code("assoc = pd.read_csv(T / 'pca_pc_associations.csv')\n"
         "assoc[assoc['PC'].isin(['PC1', 'PC2'])]"),
    code("agree = pd.read_csv(T / 'clustering_label_agreement.csv')\n"
         "agree.pivot_table(index=['clustering', 'k'], columns='label',\n"
         "                  values='adjusted_rand_index').round(2)"),
    md("## 6. Where to read the interpretation\n\n"
       "- `reports/dataset_audit_report.md`: data issues\n"
       "- `reports/eda_report.md`: every figure explained\n"
       "- `reports/preliminary_observations.md`: candidate questions for Phase 2"),
]


def main() -> int:
    nb = nbformat.v4.new_notebook(cells=CELLS)
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3",
                                 "language": "python"}
    NotebookClient(nb, timeout=300, kernel_name="python3",
                   resources={"metadata": {"path": str(NOTEBOOK.parent)}}).execute()
    NOTEBOOK.parent.mkdir(parents=True, exist_ok=True)
    nbformat.write(nb, NOTEBOOK)
    print(f"Wrote and executed {NOTEBOOK}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
