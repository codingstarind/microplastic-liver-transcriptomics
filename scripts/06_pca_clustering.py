"""Step 6: exploratory PCA and hierarchical clustering of samples.

Primary analysis
    features : 500 most variable genes of log2(median-of-ratios normalised + 1),
               after the low-expression filter
    PCA      : genes mean-centred, not scaled (as DESeq2::plotPCA)
    distance : 1 - Pearson correlation between samples, same 500 genes
    linkage  : average (UPGMA)

Sensitivity analyses (does the picture depend on these choices?)
    * top 1000 / 2000 / all filtered genes
    * log2 CPM instead of median-of-ratios
    * Euclidean distance + Ward linkage, Spearman correlation + average linkage
    * genes tracking the skeletal-muscle contamination score removed
    * the most muscle-contaminated sample (MD_1) removed

No number of clusters is imposed: tree cuts at k = 2..6 are compared with the
known labels descriptively (adjusted Rand index).
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import pdist, squareform

from mpliver.config import (
    FIGURES_DIR,
    GENE_ANNOTATION_FILE,
    GROUP_COLORS,
    PROCESSED_DIR,
    SAMPLE_METADATA_FILE,
    TABLES_DIR,
)
from mpliver.exploration import (
    cluster_label_agreement,
    correlation_distance,
    hierarchical_clustering,
    run_pca,
)
from mpliver.normalization import top_variable_genes
from mpliver.plotting import INK, MUTED, group_legend_handles, save, scatter_by_group, set_style
from mpliver.qc_markers import genes_tracking_score, marker_set_cpm

OUT = FIGURES_DIR / "pca_clustering"
N_TOP = 500
MUSCLE_RHO = 0.7  # |Spearman rho| with the muscle-marker score to call a gene muscle-tracking
GRID_LINE = "#bdbdbd"


def short_name(col: str) -> str:
    return col.rsplit("_", 1)[0]


def eta_squared(values: np.ndarray, labels: np.ndarray) -> float:
    """Share of the variance of ``values`` explained by grouping on ``labels``."""
    total = ((values - values.mean()) ** 2).sum()
    between = sum(((values[labels == g].mean() - values.mean()) ** 2) * (labels == g).sum()
                  for g in np.unique(labels))
    return between / total


def pc_associations(scores: pd.DataFrame, meta: pd.DataFrame, sf: pd.DataFrame,
                    muscle_score: pd.Series) -> pd.DataFrame:
    """Descriptive link between each PC and design / technical variables."""
    m = meta.set_index("counts_column").loc[scores.index]
    lib = sf.loc[scores.index, "library_size_all_genes"].to_numpy()
    lib_order = m["library_name"].str.lstrip("S").astype(int).to_numpy()
    muscle = muscle_score.loc[scores.index].to_numpy()
    rows = []
    for pc in scores.columns:
        v = scores[pc].to_numpy()
        rows.append({
            "PC": pc,
            "eta2_diet": eta_squared(v, m["diet"].to_numpy()),
            "eta2_microplastic": eta_squared(v, m["microplastic"].to_numpy()),
            "eta2_group": eta_squared(v, m["group"].to_numpy()),
            "spearman_r_library_size": pd.Series(v).corr(pd.Series(lib), method="spearman"),
            "spearman_r_library_index": pd.Series(v).corr(pd.Series(lib_order),
                                                          method="spearman"),
            "spearman_r_muscle_marker_cpm": pd.Series(v).corr(pd.Series(muscle),
                                                              method="spearman"),
        })
    return pd.DataFrame(rows)


def fig_pca(result, meta, title_suffix, path):
    scores, evr = result["scores"], result["explained_variance_ratio"]
    groups = meta.set_index("counts_column").loc[scores.index, "group"].to_numpy()
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.7),
                             gridspec_kw={"width_ratios": [1, 1.25, 1.25]})
    ax = axes[0]
    k = len(evr)
    ax.bar(range(1, k + 1), evr * 100, color="#2a78d6", width=0.7)
    ax.plot(range(1, k + 1), np.cumsum(evr) * 100, color=INK, marker="o", markersize=3,
            linewidth=1, label="cumulative")
    ax.set_xticks(range(1, k + 1))
    ax.set_xlabel("Principal component")
    ax.set_ylabel("% variance explained")
    ax.set_title("A  Scree plot")
    ax.legend(loc="center right")
    for ax, (a, b), letter in zip(axes[1:], [("PC1", "PC2"), ("PC3", "PC4")], "BC"):
        scatter_by_group(ax, scores[a].to_numpy(), scores[b].to_numpy(), groups)
        for name, (x, y) in scores[[a, b]].iterrows():
            ax.annotate(short_name(name), (x, y), xytext=(4, 3), textcoords="offset points",
                        fontsize=6, color=MUTED)
        ax.set_xlabel(f"{a} ({evr[a]:.1%} of variance)")
        ax.set_ylabel(f"{b} ({evr[b]:.1%} of variance)")
        ax.set_title(f"{letter}  {a} vs {b}")
        ax.axhline(0, color=GRID_LINE, linewidth=0.8)
        ax.axvline(0, color=GRID_LINE, linewidth=0.8)
    axes[1].legend(handles=group_legend_handles(), loc="best")
    fig.suptitle(f"PCA of liver samples - {title_suffix}", fontweight="bold", fontsize=10,
                 y=1.03)
    save(fig, path)


def fig_loadings(result, ann, path, n=12):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, pc in zip(axes, ["PC1", "PC2"]):
        load = result["loadings"][pc]
        top = pd.concat([load.nsmallest(n // 2), load.nlargest(n // 2)]).sort_values()
        labels = ann["gene_name"].reindex(top.index).fillna(pd.Series(top.index, index=top.index))
        colors = ["#eb6834" if v > 0 else "#2a78d6" for v in top]
        ax.barh(range(len(top)), top.to_numpy(), color=colors, height=0.7)
        ax.set_yticks(range(len(top)), labels, style="italic")
        ax.axvline(0, color=INK, linewidth=0.8)
        ax.set_xlabel(f"Loading on {pc}")
        ax.set_title(f"{'AB'[pc == 'PC2']}  Genes with the largest {pc} loadings")
        ax.grid(axis="y", visible=False)
    save(fig, path)


def fig_heatmap(log_expr, meta, distance, method, path, title):
    tree = hierarchical_clustering(distance, method=method)
    m = meta.set_index("counts_column").loc[distance.index]
    row_colors = pd.DataFrame({
        "Group": m["group"].map(GROUP_COLORS),
        "Diet": m["diet"].map({"CD": "#d9d9d9", "MD": "#525252"}),
        "PE": m["microplastic"].map({"Veh": "#d9d9d9", "PE": "#525252"}),
    })
    corr = 1 - distance
    labels = [short_name(c) for c in corr.index]
    grid = sns.clustermap(
        corr.set_axis(labels, axis=0).set_axis(labels, axis=1),
        row_linkage=tree, col_linkage=tree, cmap="Blues",
        row_colors=row_colors.set_axis(labels), col_colors=row_colors.set_axis(labels),
        figsize=(8.5, 8), cbar_kws={"label": "Pearson correlation\n(log2 normalised)"},
        dendrogram_ratio=0.15, linewidths=0.3, linecolor="white", xticklabels=True,
        yticklabels=True,
    )
    grid.ax_heatmap.tick_params(labelsize=7)
    grid.figure.suptitle(title, fontweight="bold", fontsize=10, y=1.02)
    grid.ax_heatmap.legend(handles=group_legend_handles(), loc="upper left",
                           bbox_to_anchor=(1.02, 1.25), fontsize=7)
    grid.figure.text(0.02, -0.01, "Side bars: group colour; diet and PE in grey (light = CD/Veh,"
                     " dark = MD/PE).", fontsize=7, color=MUTED)
    save(grid.figure, path)
    return tree


def fig_dendrograms(trees, meta, labels_order, path):
    groups = meta.set_index("counts_column")["group"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    for ax, (title, tree) in zip(axes.ravel(), trees.items()):
        d = dendrogram(tree, labels=[short_name(c) for c in labels_order], ax=ax,
                       color_threshold=0, above_threshold_color=INK, leaf_font_size=7,
                       leaf_rotation=90)
        for tick, idx in zip(ax.get_xticklabels(), d["leaves"]):
            tick.set_color(GROUP_COLORS[groups[labels_order[idx]]])
            tick.set_fontweight("bold")
        ax.set_title(title, fontsize=9)
        ax.grid(False)
        ax.set_ylabel("Merge height (distance)")
    fig.legend(handles=group_legend_handles(long_labels=True), loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.03))
    fig.suptitle("Sample dendrograms under different distance / linkage choices "
                 "(label colour = group)", fontweight="bold", fontsize=10)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save(fig, path)


def main() -> int:
    set_style()
    meta = pd.read_csv(SAMPLE_METADATA_FILE)
    order = meta["counts_column"].tolist()
    log_norm = pd.read_csv(PROCESSED_DIR / "log2_normalized_filtered.csv.gz", index_col=0)[order]
    log_cpm = pd.read_csv(PROCESSED_DIR / "log2_cpm_filtered.csv.gz", index_col=0)[order]
    sf = pd.read_csv(TABLES_DIR / "size_factors.csv", index_col=0)
    ann = pd.read_csv(GENE_ANNOTATION_FILE, index_col="gene")
    counts = pd.read_csv(PROCESSED_DIR / "counts_original.csv.gz", index_col=0)[order]

    # Skeletal-muscle contamination (see fig09): genes whose expression tracks the
    # muscle-marker score are removed in one sensitivity analysis.
    muscle_score = np.log2(marker_set_cpm(counts)["skeletal_muscle"] + 1)
    muscle_genes = genes_tracking_score(log_norm, muscle_score, MUSCLE_RHO)
    pd.DataFrame({"gene_name": ann["gene_name"].reindex(muscle_genes)}).to_csv(
        TABLES_DIR / "muscle_tracking_genes.csv")
    no_muscle = log_norm.drop(index=muscle_genes)
    without_md1 = log_norm.drop(columns=[c for c in order if c.startswith("MD_1_")])

    # --- PCA: primary + sensitivity ---------------------------------------
    feature_sets = {
        f"top{N_TOP}_log2norm": log_norm.loc[top_variable_genes(log_norm, N_TOP)],
        "top1000_log2norm": log_norm.loc[top_variable_genes(log_norm, 1000)],
        "top2000_log2norm": log_norm.loc[top_variable_genes(log_norm, 2000)],
        "allfiltered_log2norm": log_norm,
        f"top{N_TOP}_log2cpm": log_cpm.loc[top_variable_genes(log_cpm, N_TOP)],
        f"top{N_TOP}_log2norm_muscle_genes_removed":
            no_muscle.loc[top_variable_genes(no_muscle, N_TOP)],
        f"top{N_TOP}_log2norm_sample_MD_1_removed":
            without_md1.loc[top_variable_genes(without_md1, N_TOP)],
    }
    variance_rows, assoc_rows, score_tables = [], [], {}
    for name, features in feature_sets.items():
        result = run_pca(features, n_components=10)
        evr = result["explained_variance_ratio"]
        variance_rows.append({"feature_set": name, "n_genes": len(features),
                              **{f"{pc}_pct": v * 100 for pc, v in evr.items()}})
        assoc = pc_associations(result["scores"], meta, sf, muscle_score).head(5)
        assoc.insert(0, "feature_set", name)
        assoc_rows.append(assoc)
        score_tables[name] = result["scores"]
        if name == f"top{N_TOP}_log2norm":
            primary = result
            fig_pca(result, meta, f"{N_TOP} most variable genes, log2 normalised counts",
                    OUT / "fig10_pca_primary")
            fig_loadings(result, ann, OUT / "fig11_pca_top_loadings")
            loadings = result["loadings"][["PC1", "PC2", "PC3"]].join(
                ann[["gene_name", "gene_type"]])
            loadings.to_csv(TABLES_DIR / "pca_primary_loadings.csv")
        elif name == "allfiltered_log2norm":
            fig_pca(result, meta, f"all {len(features):,} filtered genes (sensitivity)",
                    OUT / "fig12_pca_all_filtered_genes")
        elif name.endswith("muscle_genes_removed"):
            fig_pca(result, meta, f"{N_TOP} most variable genes after removing "
                    f"{len(muscle_genes)} muscle-tracking genes (sensitivity)",
                    OUT / "fig13_pca_muscle_genes_removed")

    pd.DataFrame(variance_rows).round(2).to_csv(TABLES_DIR / "pca_variance_explained.csv",
                                                index=False)
    pd.concat(assoc_rows).round(3).to_csv(TABLES_DIR / "pca_pc_associations.csv", index=False)
    scores = primary["scores"].join(meta.set_index("counts_column")[["group"]])
    scores.round(4).to_csv(TABLES_DIR / "pca_primary_scores.csv")

    # Outlier screen: distance of each sample to its own group centroid in PC1-PC4.
    pcs = primary["scores"][["PC1", "PC2", "PC3", "PC4"]]
    m = meta.set_index("counts_column").loc[pcs.index]
    centroid = pcs.groupby(m["group"].to_numpy()).transform("mean")
    dist = np.sqrt(((pcs - centroid) ** 2).sum(axis=1))
    pooled = dist.groupby(m["group"].to_numpy()).transform("median")
    pd.DataFrame({"group": m["group"], "distance_to_group_centroid_PC1to4": dist,
                  "ratio_to_group_median_distance": dist / pooled}).round(3).sort_values(
        "ratio_to_group_median_distance", ascending=False).to_csv(
        TABLES_DIR / "pca_outlier_screen.csv")

    # --- Hierarchical clustering --------------------------------------------
    top = log_norm.loc[top_variable_genes(log_norm, N_TOP)]
    pearson = correlation_distance(top, "pearson")
    spearman = correlation_distance(top, "spearman")
    euclid = pd.DataFrame(squareform(pdist(top.T.to_numpy(), "euclidean")),
                          index=top.columns, columns=top.columns)
    tree_primary = fig_heatmap(top, meta, pearson, "average",
                               OUT / "fig14_sample_correlation_heatmap",
                               f"Sample-sample correlation, {N_TOP} most variable genes "
                               "(1 - Pearson, average linkage)")
    top_no_muscle = no_muscle.loc[top_variable_genes(no_muscle, N_TOP)]
    trees = {
        "1 - Pearson, average (primary)": tree_primary,
        "1 - Spearman, average": hierarchical_clustering(spearman, "average"),
        "Euclidean, Ward": hierarchical_clustering(euclid, "ward"),
        "1 - Pearson, average, muscle-tracking genes removed": hierarchical_clustering(
            correlation_distance(top_no_muscle, "pearson"), "average"),
    }
    fig_dendrograms(trees, meta, top.columns.tolist(), OUT / "fig15_dendrograms_sensitivity")

    m = meta.set_index("counts_column").loc[top.columns]
    agreement = []
    for tree_name, tree in trees.items():
        for label in ["group", "diet", "microplastic"]:
            df = cluster_label_agreement(tree, m[label].rename(label))
            df.insert(0, "clustering", tree_name)
            agreement.append(df)
    pd.concat(agreement).round(3).to_csv(TABLES_DIR / "clustering_label_agreement.csv",
                                         index=False)

    print(pd.DataFrame(variance_rows).round(1).iloc[:, :6].to_string(index=False))
    print(pd.concat(assoc_rows).round(2).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
