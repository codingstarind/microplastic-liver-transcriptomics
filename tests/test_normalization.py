import numpy as np
import pandas as pd
import pytest

from mpliver.normalization import (
    counts_per_million,
    filter_low_expression,
    log2_transform,
    median_of_ratios_size_factors,
    normalize_by_size_factors,
)


def test_size_factors_recover_known_depth_ratio():
    base = pd.Series([100.0, 50.0, 10.0, 7.0, 300.0])
    counts = pd.DataFrame({"a": base, "b": base * 2, "c": base * 0.5})
    sf = median_of_ratios_size_factors(counts)
    assert sf["b"] / sf["a"] == pytest.approx(2.0)
    assert sf["c"] / sf["a"] == pytest.approx(0.5)
    normalized = normalize_by_size_factors(counts, sf)
    assert np.allclose(normalized["a"], normalized["b"])


def test_size_factors_ignore_genes_with_zeros():
    counts = pd.DataFrame({"a": [10.0, 0.0, 20.0], "b": [20.0, 1000.0, 40.0]})
    sf = median_of_ratios_size_factors(counts)
    assert sf["b"] / sf["a"] == pytest.approx(2.0)


def test_size_factors_robust_to_one_dominant_gene():
    counts = pd.DataFrame({"a": [10.0, 20.0, 30.0, 40.0, 50.0],
                           "b": [10.0, 20.0, 30.0, 40.0, 50000.0]})
    sf = median_of_ratios_size_factors(counts)
    assert sf["b"] / sf["a"] == pytest.approx(1.0)


def test_cpm_columns_sum_to_one_million():
    counts = pd.DataFrame({"a": [1.0, 3.0], "b": [10.0, 30.0]})
    assert np.allclose(counts_per_million(counts).sum(), 1e6)


def test_log2_transform_does_not_modify_input():
    counts = pd.DataFrame({"a": [0.0, 1.0, 3.0]})
    out = log2_transform(counts)
    assert out["a"].tolist() == [0.0, 1.0, 2.0]
    assert counts["a"].tolist() == [0.0, 1.0, 3.0]


def test_log2_transform_rejects_negative():
    with pytest.raises(ValueError):
        log2_transform(pd.DataFrame({"a": [-1.0]}))


def test_filter_low_expression():
    counts = pd.DataFrame({
        "s1": [10, 0, 50], "s2": [12, 0, 1], "s3": [0, 9, 1], "s4": [11, 30, 1]})
    mask = filter_low_expression(counts, min_count=10, min_samples=3)
    assert mask.tolist() == [True, False, False]
