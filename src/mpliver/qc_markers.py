"""Marker-gene scores for tissue contamination and cell-type composition.

Liver dissections can carry neighbouring tissue (diaphragm muscle, fat, blood).
Summing reads from genes that are specific to that tissue, per million reads,
gives a simple per-sample contamination indicator.
"""

import pandas as pd

# Marker sets chosen to be specific to the named tissue/cell type and low in
# hepatocytes. Genes that the MASH diet induces in liver itself (Cidea, Car3 and
# Cfd/adipsin, which rises ~100x in every MD sample while Adipoq/Lep/Plin1 stay
# at ~0) are deliberately excluded from the adipocyte set.
MARKER_SETS = {
    "skeletal_muscle": ["Myh1", "Myh2", "Myh4", "Acta1", "Ckm", "Mb", "Ttn", "Tnnt3",
                        "Mylpf", "Atp2a1", "Actn3", "Tnni2"],
    "adipocyte": ["Adipoq", "Lep", "Retn", "Plin1"],
    "erythrocyte": ["Hba-a1", "Hbb-bs", "Hbb-bt", "Alas2"],
    "pancreas_exocrine": ["Pnlip", "Cpa1", "Prss2", "Cela2a", "Ctrb1"],
    "cholangiocyte": ["Krt19", "Krt7", "Epcam", "Spp1"],
    "kupffer_macrophage": ["Clec4f", "Vsig4", "Cd5l", "Marco"],
    "stellate_fibrogenic": ["Col1a1", "Col3a1", "Dcn", "Lrat"],
}


def marker_set_cpm(counts: pd.DataFrame, marker_sets: dict = MARKER_SETS) -> pd.DataFrame:
    """Samples x marker sets: summed counts of each set's genes per million reads."""
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
    return pd.DataFrame({
        name: cpm.loc[[g for g in genes if g in cpm.index]].sum(axis=0)
        for name, genes in marker_sets.items()
    })


def genes_tracking_score(log_expr: pd.DataFrame, score: pd.Series,
                         min_abs_rho: float = 0.7) -> pd.Index:
    """Genes whose Spearman correlation with ``score`` across samples is >= ``min_abs_rho``.

    Spearman is used because one heavily contaminated sample would dominate a
    Pearson correlation.
    """
    ranks = log_expr[score.index].rank(axis=1)
    score_rank = score.rank()
    rho = ranks.apply(lambda row: row.corr(score_rank), axis=1)
    return rho[rho.abs() >= min_abs_rho].index
