import numpy as np
import pandas as pd

from mpliver.counts import (
    audit_counts,
    check_sample_ids,
    classify_gene_ids,
    non_integer_fraction,
    split_annotation_and_counts,
    strip_ensembl_version,
)


def make_raw():
    return pd.DataFrame({
        "gene_id": ["ENSMUSG00000000001.5", "ENSMUSG00000000003.16", "ENSMUSG00000000028.16"],
        "CD_1_S39": [10.0, 0.0, 3.5],
        "MD_MP_1_S53": [12.0, 0.0, 4.0],
    })


def test_split_annotation_and_counts():
    annotation, counts = split_annotation_and_counts(make_raw())
    assert list(counts.columns) == ["CD_1_S39", "MD_MP_1_S53"]
    assert counts.index.name == "gene"
    assert annotation.shape == (3, 0)


def test_classify_gene_ids():
    ids = pd.Index(["ENSMUSG00000000001.5", "ENSMUSG00000000001", "Actb"])
    assert classify_gene_ids(ids).tolist() == [
        "ensembl_versioned", "ensembl_unversioned", "symbol"]


def test_strip_version_leaves_dotted_symbols_alone():
    ids = pd.Index(["ENSMUSG00000000001.5", "H2-M10.1", "H2-M10.2"])
    assert strip_ensembl_version(ids).tolist() == ["ENSMUSG00000000001", "H2-M10.1", "H2-M10.2"]


def test_non_integer_fraction():
    assert non_integer_fraction(pd.DataFrame({"a": [1.0, 2.0, 2.5, np.nan]})) == 1 / 3


def test_audit_counts_detects_problems():
    counts = pd.DataFrame(
        {"s1": [0.0, 5.0, 5.0, np.nan], "s2": [0.0, 7.0, 7.0, 1.0]},
        index=pd.Index(["ENSMUSG00000000001.1", "ENSMUSG00000000002.1",
                        "ENSMUSG00000000002.2", "ENSMUSG00000000004.1"]),
    )
    audit = audit_counts(counts)
    assert audit["n_genes_all_zero"] == 1
    assert audit["n_missing_values"] == 1
    assert audit["n_duplicate_gene_ids"] == 0
    assert audit["n_duplicate_ids_after_version_strip"] == 1
    assert audit["n_duplicate_count_profiles"] == 1


def test_check_sample_ids():
    result = check_sample_ids(["a", "b"], ["b", "c"])
    assert not result["all_match"]
    assert result["only_in_counts"] == ["a"]
    assert result["only_in_metadata"] == ["c"]
