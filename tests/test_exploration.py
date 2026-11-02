import numpy as np
import pandas as pd
import pytest

from mpliver.exploration import (
    cluster_label_agreement,
    correlation_distance,
    hierarchical_clustering,
    run_pca,
)


@pytest.fixture
def two_group_matrix():
    rng = np.random.default_rng(0)
    genes, per_group = 200, 4
    # Shared gene-level baseline (as in real expression data) + a group effect on 20 genes.
    baseline = rng.normal(scale=3, size=genes)[:, None]
    shift = np.zeros(genes)
    shift[:20] = 3
    a = baseline + rng.normal(scale=0.5, size=(genes, per_group))
    b = baseline + rng.normal(scale=0.5, size=(genes, per_group)) + shift[:, None]
    cols = [f"a{i}" for i in range(per_group)] + [f"b{i}" for i in range(per_group)]
    return pd.DataFrame(np.hstack([a, b]), columns=cols)


def test_pca_separates_groups_on_pc1(two_group_matrix):
    result = run_pca(two_group_matrix, n_components=3)
    pc1 = result["scores"]["PC1"]
    assert np.sign(pc1.filter(like="a")).nunique() == 1
    assert np.sign(pc1.filter(like="a").iloc[0]) != np.sign(pc1.filter(like="b").iloc[0])
    assert result["explained_variance_ratio"].is_monotonic_decreasing


def test_clustering_recovers_groups(two_group_matrix):
    dist = correlation_distance(two_group_matrix)
    assert np.allclose(np.diag(dist), 0)
    tree = hierarchical_clustering(dist, method="average")
    labels = pd.Series(["a"] * 4 + ["b"] * 4, name="group")
    agreement = cluster_label_agreement(tree, labels, k_values=[2])
    assert agreement.loc[0, "adjusted_rand_index"] == pytest.approx(1.0)
