# data/processed

Derived matrices written by `scripts/04_preprocess_counts.py`. Not committed;
regenerate them from the raw file. File names state the transformation applied.
See `metadata/data_dictionary.md` for column definitions and
`reports/preprocessing_decisions.md` for the reasoning.

| File | Derived from | Transformation |
|---|---|---|
| `counts_original.csv.gz` | raw counts | none (sample columns put in metadata order) |
| `gene_filter.csv` | raw counts | per-gene filter statistics, `passes_filter` flag |
| `counts_filtered.csv.gz` | raw counts | low-expression genes removed; values unchanged |
| `normalized_counts_filtered.csv.gz` | filtered counts | divided by median-of-ratios size factor |
| `log2_normalized_filtered.csv.gz` | normalised counts | log2(x + 1) |
| `log2_cpm_filtered.csv.gz` | filtered counts | log2(CPM + 1); sensitivity check only |
