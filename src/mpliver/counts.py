"""Load the deposited count table and run data-integrity checks on it.

Nothing here modifies values: every function either reads the file or
summarises what is in it.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd

ENSEMBL_GENE_RE = re.compile(r"^ENSMUSG\d{11}(\.\d+)?(_PAR_Y)?$")
SAMPLE_COLUMN_RE = re.compile(r"^(CD|MD)(_MP)?_\d+_S\d+$")


def read_counts_file(path: Path) -> pd.DataFrame:
    """Read the deposited CSV exactly as provided (no index, no type coercion of IDs)."""
    return pd.read_csv(path, dtype={0: str})


def split_annotation_and_counts(raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Separate identifier/annotation columns from per-sample count columns.

    Sample columns are recognised by the naming scheme documented in GEO
    (e.g. ``CD_MP_1_S48``). The first column is used as the gene identifier.
    """
    sample_cols = [c for c in raw.columns if SAMPLE_COLUMN_RE.match(str(c))]
    annotation_cols = [c for c in raw.columns if c not in sample_cols]
    if not annotation_cols:
        raise ValueError("no gene identifier column found")
    id_col = annotation_cols[0]
    annotation = raw[annotation_cols].set_index(id_col)
    counts = raw[[id_col, *sample_cols]].set_index(id_col)
    counts.index.name = "gene"
    annotation.index.name = "gene"
    return annotation, counts


def classify_gene_ids(ids: pd.Index) -> pd.Series:
    """Label each identifier as versioned/unversioned Ensembl gene ID or 'symbol'."""
    def label(gene_id: str) -> str:
        if not ENSEMBL_GENE_RE.match(gene_id):
            return "symbol"
        return "ensembl_versioned" if "." in gene_id else "ensembl_unversioned"
    return pd.Series([label(str(i)) for i in ids], index=ids, name="id_format")


def strip_ensembl_version(ids: pd.Index) -> pd.Index:
    """Drop '.N' version suffixes from Ensembl IDs; symbols such as 'H2-M10.1' are kept."""
    return pd.Index([str(i).split(".")[0] if ENSEMBL_GENE_RE.match(str(i)) else str(i)
                     for i in ids], name=ids.name)


def non_integer_fraction(counts: pd.DataFrame, tol: float = 1e-9) -> float:
    """Fraction of non-missing values that are not whole numbers."""
    values = counts.to_numpy(dtype=float)
    values = values[~np.isnan(values)]
    return float(np.mean(np.abs(values - np.round(values)) > tol)) if values.size else 0.0


def check_sample_ids(count_columns, metadata_columns) -> dict:
    """Compare sample identifiers between the count matrix and the metadata table."""
    counts_set, meta_set = set(count_columns), set(metadata_columns)
    return {
        "n_counts_columns": len(counts_set),
        "n_metadata_samples": len(meta_set),
        "only_in_counts": sorted(counts_set - meta_set),
        "only_in_metadata": sorted(meta_set - counts_set),
        "all_match": counts_set == meta_set,
    }


def audit_counts(counts: pd.DataFrame) -> dict:
    """Integrity summary of a genes x samples matrix (values left untouched)."""
    values = counts.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    id_formats = classify_gene_ids(counts.index).value_counts().to_dict()
    unversioned = strip_ensembl_version(counts.index)
    row_sums = counts.sum(axis=1)

    return {
        "n_genes": int(counts.shape[0]),
        "n_samples": int(counts.shape[1]),
        "dtypes": sorted({str(t) for t in counts.dtypes}),
        "gene_id_formats": {k: int(v) for k, v in id_formats.items()},
        "n_duplicate_gene_ids": int(counts.index.duplicated().sum()),
        "n_duplicate_ids_after_version_strip": int(unversioned.duplicated().sum()),
        "n_duplicate_count_profiles": int(counts[row_sums > 0].duplicated().sum()),
        "n_missing_values": int(np.isnan(values).sum()),
        "n_infinite_values": int(np.isinf(values).sum()),
        "n_negative_values": int((finite < 0).sum()),
        "min_value": float(finite.min()),
        "max_value": float(finite.max()),
        "non_integer_fraction": non_integer_fraction(counts),
        "fraction_zero": float((values == 0).mean()),
        "n_genes_all_zero": int((row_sums == 0).sum()),
        "n_genes_zero_in_any_sample": int((counts == 0).any(axis=1).sum()),
    }


def value_type_evidence(counts: pd.DataFrame) -> dict:
    """Facts that distinguish raw/estimated counts from TPM or normalised values.

    - TPM/CPM: every column sums to ~1e6.
    - Integer read counts (featureCounts/HTSeq, or tximport counts rounded for
      DESeq2): no fractional values; smallest non-zero value is 1.
    - Salmon/tximport estimated counts (countsFromAbundance="no"): fractional
      values, many non-zero values below 1, column sums of the order of the
      number of assigned reads.
    - DESeq2 normalised counts: fractional, but within a sample every value is
      an integer multiple of 1/size_factor.
    """
    sums = counts.sum(axis=0)
    nonzero = counts.where(counts > 0)
    min_nonzero = nonzero.min(axis=0)
    # If a column were integer / size_factor, dividing by its smallest non-zero
    # value would give (near-)integers.
    ratios = nonzero.div(min_nonzero, axis=1)
    ratio_integer_fraction = float(
        (np.abs(ratios - np.round(ratios)) < 1e-6).sum().sum() / ratios.notna().sum().sum()
    )
    fractional = (counts - np.floor(counts)).to_numpy().ravel()
    fractional = fractional[fractional > 1e-9]
    return {
        "column_sum_min": float(sums.min()),
        "column_sum_max": float(sums.max()),
        "columns_sum_to_1e6": bool(np.allclose(sums, 1e6, rtol=1e-3)),
        "non_integer_fraction": non_integer_fraction(counts),
        "min_nonzero_value": float(min_nonzero.min()),
        "fraction_nonzero_below_1": float(((counts > 0) & (counts < 1)).sum().sum()
                                          / (counts > 0).sum().sum()),
        "fraction_values_integer_multiple_of_column_min": ratio_integer_fraction,
        "n_fractional_values": int(fractional.size),
        "fractional_part_quartiles": (np.quantile(fractional, [0.25, 0.5, 0.75]).tolist()
                                      if fractional.size else []),
    }


def per_sample_summary(counts: pd.DataFrame) -> pd.DataFrame:
    """Per-sample totals and detection counts on the original scale."""
    return pd.DataFrame({
        "total_counts": counts.sum(axis=0),
        "n_genes_detected": (counts > 0).sum(axis=0),
        "n_genes_ge_10": (counts >= 10).sum(axis=0),
        "fraction_zero": (counts == 0).mean(axis=0),
        "max_gene_count": counts.max(axis=0),
        "top_gene": counts.idxmax(axis=0),
        "non_integer_fraction": counts.apply(lambda col: non_integer_fraction(col.to_frame())),
    })
