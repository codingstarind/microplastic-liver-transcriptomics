"""Library-size normalisation, low-expression filtering and log transformation.

These functions return new objects and never alter the input matrix, so the
original deposited values stay available for comparison.
"""

import numpy as np
import pandas as pd


def library_sizes(counts: pd.DataFrame) -> pd.Series:
    return counts.sum(axis=0).rename("library_size")


def counts_per_million(counts: pd.DataFrame) -> pd.DataFrame:
    return counts.div(library_sizes(counts), axis=1) * 1e6


def median_of_ratios_size_factors(counts: pd.DataFrame) -> pd.Series:
    """DESeq2-style size factors (Anders & Huber 2010).

    For each gene, take the geometric mean across samples (genes with a zero in
    any sample are skipped). A sample's size factor is the median of its
    count / geometric-mean ratios. Robust to a few very highly expressed or
    strongly changing genes, unlike total-count scaling.
    """
    positive = counts[(counts > 0).all(axis=1)]
    if positive.empty:
        raise ValueError("no genes are non-zero in every sample")
    log_counts = np.log(positive)
    log_geo_means = log_counts.mean(axis=1)
    log_ratios = log_counts.sub(log_geo_means, axis=0)
    return np.exp(log_ratios.median(axis=0)).rename("size_factor")


def normalize_by_size_factors(counts: pd.DataFrame, size_factors: pd.Series) -> pd.DataFrame:
    return counts.div(size_factors[counts.columns], axis=1)


def log2_transform(values: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
    """log2(x + pseudocount). Compresses the dynamic range; the pseudocount avoids log(0)."""
    if (values < 0).any().any():
        raise ValueError("negative values cannot be log-transformed")
    return np.log2(values + pseudocount)


def filter_low_expression(counts: pd.DataFrame, min_count: float = 10,
                          min_samples: int = 4) -> pd.Series:
    """Boolean mask: genes with >= ``min_count`` in at least ``min_samples`` samples.

    ``min_samples`` defaults to the smallest group size (CD+Veh, n = 4) so that a
    gene expressed in only one group is still kept. This mirrors the filter
    suggested in the DESeq2 vignette.
    """
    return ((counts >= min_count).sum(axis=1) >= min_samples).rename("passes_filter")


def top_variable_genes(log_expr: pd.DataFrame, n: int = 500) -> pd.Index:
    """Genes with the highest variance across samples (on the log scale)."""
    return log_expr.var(axis=1).sort_values(ascending=False).index[:n]
