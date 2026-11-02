"""Step 4: derive analysis-ready matrices from the deposited counts.

The raw file is never overwritten. Every derived matrix is written to
data/processed/ with a name that says what was done to it:

    counts_original.csv.gz           deposited values, sample columns in metadata order
    gene_filter.csv                  per-gene filter statistics and the keep/drop flag
    size_factors.csv (results/tables) median-of-ratios size factors and library sizes
    counts_filtered.csv.gz           original values, low-expression genes removed
    normalized_counts_filtered.csv.gz  filtered counts / size factor
    log2_normalized_filtered.csv.gz  log2(normalized + 1), used for EDA, PCA, clustering
    log2_cpm_filtered.csv.gz         log2(CPM + 1), used only as a sensitivity check

Choices (see reports/preprocessing_decisions.md for the reasoning):
    * no rounding - estimated counts are kept as they are for exploration;
    * low-expression filter: >= 10 counts in >= 4 samples (smallest group size);
    * size factors from DESeq2's median-of-ratios method on the filtered genes;
    * log2(x + 1) for exploration only; no gene-wise standardisation is stored.
"""

import json
import sys

import pandas as pd

from mpliver.config import (
    COUNTS_FILE,
    GENE_ANNOTATION_FILE,
    PROCESSED_DIR,
    SAMPLE_METADATA_FILE,
    TABLES_DIR,
)
from mpliver.counts import read_counts_file, split_annotation_and_counts
from mpliver.normalization import (
    counts_per_million,
    filter_low_expression,
    library_sizes,
    log2_transform,
    median_of_ratios_size_factors,
    normalize_by_size_factors,
)

MIN_COUNT = 10
MIN_SAMPLES = 4


def main() -> int:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(SAMPLE_METADATA_FILE)
    _, counts = split_annotation_and_counts(read_counts_file(COUNTS_FILE))
    counts = counts[metadata["counts_column"]]  # fixed, documented sample order
    annotation = pd.read_csv(GENE_ANNOTATION_FILE, index_col="gene")

    counts.to_csv(PROCESSED_DIR / "counts_original.csv.gz")

    keep = filter_low_expression(counts, MIN_COUNT, MIN_SAMPLES)
    gene_filter = pd.DataFrame({
        "gene_name": annotation["gene_name"],
        "gene_type": annotation["gene_type"],
        "mean_count": counts.mean(axis=1),
        "n_samples_ge_min_count": (counts >= MIN_COUNT).sum(axis=1),
        "passes_filter": keep,
    })
    gene_filter.to_csv(PROCESSED_DIR / "gene_filter.csv")
    filtered = counts[keep]

    size_factors = median_of_ratios_size_factors(filtered)
    sf_table = pd.DataFrame({
        "library_size_all_genes": library_sizes(counts),
        "library_size_filtered_genes": library_sizes(filtered),
        "size_factor": size_factors,
    })
    sf_table["library_size_relative_to_mean"] = (
        sf_table["library_size_filtered_genes"] / sf_table["library_size_filtered_genes"].mean())
    sf_table = sf_table.join(metadata.set_index("counts_column")[["group"]])
    sf_table.index.name = "counts_column"
    sf_table.to_csv(TABLES_DIR / "size_factors.csv")

    normalized = normalize_by_size_factors(filtered, size_factors)
    log_normalized = log2_transform(normalized)
    log_cpm = log2_transform(counts_per_million(filtered))

    filtered.to_csv(PROCESSED_DIR / "counts_filtered.csv.gz")
    normalized.to_csv(PROCESSED_DIR / "normalized_counts_filtered.csv.gz")
    log_normalized.to_csv(PROCESSED_DIR / "log2_normalized_filtered.csv.gz")
    log_cpm.to_csv(PROCESSED_DIR / "log2_cpm_filtered.csv.gz")

    filter_summary = {
        "rule": f">= {MIN_COUNT} counts in >= {MIN_SAMPLES} samples",
        "n_genes_input": int(len(counts)),
        "n_genes_kept": int(keep.sum()),
        "n_genes_removed": int((~keep).sum()),
        "fraction_of_total_counts_kept": float(filtered.to_numpy().sum()
                                               / counts.to_numpy().sum()),
        "kept_by_gene_type": gene_filter[keep]["gene_type"].value_counts().head(8).to_dict(),
        "size_factor_range": [float(size_factors.min()), float(size_factors.max())],
    }
    (TABLES_DIR / "filter_summary.json").write_text(
        json.dumps(filter_summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(filter_summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
