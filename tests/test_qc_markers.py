import pandas as pd

from mpliver.qc_markers import genes_tracking_score, marker_set_cpm


def test_marker_set_cpm_sums_set_genes_per_million():
    counts = pd.DataFrame({"s1": [900_000, 100_000, 0], "s2": [990_000, 0, 10_000]},
                          index=["Alb", "Myh1", "Acta1"])
    scores = marker_set_cpm(counts, {"muscle": ["Myh1", "Acta1", "NotAGene"]})
    assert scores.loc["s1", "muscle"] == 100_000
    assert scores.loc["s2", "muscle"] == 10_000


def test_genes_tracking_score_finds_correlated_genes():
    score = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0], index=list("abcde"))
    log_expr = pd.DataFrame({
        "a": [1, 5, 2], "b": [2, 4, 2], "c": [3, 3, 2], "d": [4, 2, 2], "e": [5, 1, 2.1]},
        index=["up", "down", "flat"], dtype=float)
    hits = genes_tracking_score(log_expr, score, min_abs_rho=0.9)
    assert set(hits) == {"up", "down"}
