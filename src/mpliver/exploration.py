"""PCA and hierarchical clustering helpers for exploratory analysis.

Samples are observations and genes are features throughout.
"""

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score


def run_pca(log_expr: pd.DataFrame, n_components: int = 10) -> dict:
    """PCA on a genes x samples log-expression matrix.

    Genes are mean-centred but not scaled to unit variance, so genes with larger
    log-scale variation contribute more (the same convention as DESeq2's plotPCA).
    """
    X = log_expr.T.to_numpy()
    n_components = min(n_components, X.shape[0], X.shape[1])
    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(X - X.mean(axis=0))
    pcs = [f"PC{i + 1}" for i in range(n_components)]
    return {
        "scores": pd.DataFrame(scores, index=log_expr.columns, columns=pcs),
        "explained_variance_ratio": pd.Series(pca.explained_variance_ratio_, index=pcs),
        "loadings": pd.DataFrame(pca.components_.T, index=log_expr.index, columns=pcs),
    }


def correlation_distance(log_expr: pd.DataFrame, method: str = "pearson") -> pd.DataFrame:
    """1 - correlation between samples."""
    return 1 - log_expr.corr(method=method)


def hierarchical_clustering(distance: pd.DataFrame, method: str = "average") -> np.ndarray:
    """Linkage matrix from a square sample distance matrix."""
    condensed = squareform(distance.to_numpy(), checks=False)
    return linkage(condensed, method=method)


def cluster_label_agreement(linkage_matrix: np.ndarray, labels: pd.Series,
                            k_values=range(2, 7)) -> pd.DataFrame:
    """Adjusted Rand index between tree cuts at several k and known sample labels.

    Purely descriptive: it shows how well unsupervised clusters line up with the
    design at each k, without choosing a "correct" k.
    """
    rows = []
    for k in k_values:
        clusters = fcluster(linkage_matrix, t=k, criterion="maxclust")
        rows.append({
            "k": k,
            "label": labels.name,
            "adjusted_rand_index": adjusted_rand_score(labels.to_numpy(), clusters),
            "cluster_sizes": ",".join(str(n) for n in np.bincount(clusters)[1:]),
        })
    return pd.DataFrame(rows)
