"""Step 5: descriptive statistics and exploratory figures.

No hypothesis tests are run here. Figures go to figures/eda/, tables to
results/tables/. Interpretation is in reports/eda_report.md.
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mpliver.config import (
    FIGURES_DIR,
    GENE_ANNOTATION_FILE,
    GROUP_COLORS,
    GROUP_ORDER,
    PROCESSED_DIR,
    SAMPLE_METADATA_FILE,
    TABLES_DIR,
)
from mpliver.plotting import INK, MUTED, group_legend_handles, save, scatter_by_group, set_style
from mpliver.qc_markers import MARKER_SETS, marker_set_cpm

OUT = FIGURES_DIR / "eda"

# Genes chosen a priori, before looking at the data: markers named in the
# original paper (Jung et al. 2026) plus textbook steatosis/inflammation/fibrosis
# markers. Shown descriptively only.
SELECTED_GENES = [
    ("Ppara", "paper: induced by PE"),
    ("Nr1i3", "paper: induced in MD+PE"),
    ("Anxa2", "paper: spatial PE marker"),
    ("Plin4", "paper: shared PE DEG"),
    ("Fasn", "paper: shared PE DEG"),
    ("Pnpla3", "paper: shared PE DEG"),
    ("Cd36", "lipid uptake / steatosis"),
    ("Cidea", "lipid droplet / steatosis"),
    ("Lcn2", "acute inflammation"),
    ("Col1a1", "fibrosis"),
    ("Cyp4a14", "PPARa target"),
    ("Cyp7a1", "bile acid synthesis"),
]


def load():
    meta = pd.read_csv(SAMPLE_METADATA_FILE)
    meta["group"] = pd.Categorical(meta["group"], GROUP_ORDER, ordered=True)
    counts = pd.read_csv(PROCESSED_DIR / "counts_original.csv.gz", index_col=0)
    filtered = pd.read_csv(PROCESSED_DIR / "counts_filtered.csv.gz", index_col=0)
    log_norm = pd.read_csv(PROCESSED_DIR / "log2_normalized_filtered.csv.gz", index_col=0)
    sf = pd.read_csv(TABLES_DIR / "size_factors.csv", index_col=0)
    ann = pd.read_csv(GENE_ANNOTATION_FILE, index_col="gene")
    audit = pd.read_csv(TABLES_DIR / "per_sample_audit.csv", index_col=0)
    return meta, counts, filtered, log_norm, sf, ann, audit


def short_name(col: str) -> str:
    """'MD_MP_3_S55' -> 'MD_MP_3' for axis labels."""
    return col.rsplit("_", 1)[0]


def fig_library_sizes(meta, sf):
    groups = meta.set_index("counts_column")["group"]
    order = meta["counts_column"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), gridspec_kw={"width_ratios": [2.2, 1]})
    ax = axes[0]
    totals = sf.loc[order, "library_size_all_genes"] / 1e6
    ax.bar(range(len(order)), totals, color=[GROUP_COLORS[g] for g in groups[order]],
           edgecolor="white", linewidth=1)
    ax.set_xticks(range(len(order)), [short_name(c) for c in order], rotation=60, ha="right")
    ax.set_ylabel("Total estimated counts (millions)")
    ax.set_title("A  Library size per sample")
    ax.axhline(totals.mean(), color=MUTED, linestyle="--", linewidth=1)
    ax.text(len(order) - 0.5, totals.mean(), " mean", va="bottom", ha="right", color=MUTED)
    ax.legend(handles=group_legend_handles(), ncol=4, loc="upper left",
              bbox_to_anchor=(0, -0.32))

    ax = axes[1]
    scatter_by_group(ax, sf.loc[order, "library_size_relative_to_mean"].to_numpy(),
                     sf.loc[order, "size_factor"].to_numpy(), groups[order].to_numpy())
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, color=MUTED, linewidth=1, linestyle=":")
    ax.set_xlabel("Library size / mean library size")
    ax.set_ylabel("Median-of-ratios size factor")
    ax.set_title("B  Size factor vs. library size")
    save(fig, OUT / "fig01_library_sizes")


def fig_count_distribution(meta, counts, gene_filter_threshold=10):
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    ax = axes[0]
    mean_counts = counts.mean(axis=1)
    bins = np.linspace(0, np.log10(mean_counts.max() + 1) + 0.1, 60)
    ax.hist(np.log10(mean_counts + 1), bins=bins, color="#2a78d6", edgecolor="white",
            linewidth=0.5)
    ax.axvline(np.log10(gene_filter_threshold + 1), color=INK, linestyle="--", linewidth=1)
    ax.text(np.log10(gene_filter_threshold + 1), ax.get_ylim()[1] * 0.92,
            f"  mean = {gene_filter_threshold}", color=INK)
    ax.set_xlabel("log10(mean count across 19 samples + 1)")
    ax.set_ylabel("Number of genes")
    ax.set_title(f"A  Gene-level mean counts (all {len(counts):,} genes)")

    ax = axes[1]
    groups = meta.set_index("counts_column")["group"]
    grid = np.linspace(0, np.log2(counts.to_numpy().max() + 1), 200)
    for col in meta["counts_column"]:
        values = np.log2(counts[col] + 1)
        hist, edges = np.histogram(values[values > 0], bins=grid, density=True)
        centers = (edges[:-1] + edges[1:]) / 2
        ax.plot(centers, hist, color=GROUP_COLORS[groups[col]], linewidth=1, alpha=0.8)
    ax.set_xlabel("log2(count + 1), non-zero genes only")
    ax.set_ylabel("Density")
    ax.set_title("B  Per-sample distribution (one line per sample)")
    ax.legend(handles=group_legend_handles(), loc="upper right")
    save(fig, OUT / "fig02_count_distributions")


def fig_boxplots_before_after(meta, filtered, log_norm):
    order = meta["counts_column"]
    groups = meta.set_index("counts_column")["group"]
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    panels = [
        (np.log2(filtered[order] + 1), "A  Before normalisation: log2(count + 1)"),
        (log_norm[order], "B  After median-of-ratios normalisation: log2(normalised count + 1)"),
    ]
    for ax, (data, title) in zip(axes, panels):
        parts = ax.boxplot([data[c] for c in order], patch_artist=True, widths=0.6,
                           showfliers=False, medianprops={"color": INK, "linewidth": 1.2})
        for patch, col in zip(parts["boxes"], order):
            patch.set_facecolor(GROUP_COLORS[groups[col]])
            patch.set_edgecolor(INK)
            patch.set_linewidth(0.6)
        ax.axhline(np.median(data.to_numpy()), color=MUTED, linestyle="--", linewidth=1)
        ax.set_ylabel("log2 expression")
        ax.set_title(title)
    axes[1].set_xticks(range(1, len(order) + 1), [short_name(c) for c in order],
                       rotation=60, ha="right")
    fig.legend(handles=group_legend_handles(long_labels=True), ncol=4, loc="lower center",
               bbox_to_anchor=(0.5, -0.07))
    fig.text(0.01, -0.1, f"{len(filtered):,} genes passing the low-expression filter. "
             "Whiskers: 1.5 x IQR; outlying genes not drawn.", color=MUTED, fontsize=8)
    save(fig, OUT / "fig03_sample_boxplots_before_after_normalisation")


def fig_expression_dominance(meta, counts, ann):
    groups = meta.set_index("counts_column")["group"]
    share = counts.div(counts.sum(), axis=1)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), gridspec_kw={"width_ratios": [1.2, 1]})
    ax = axes[0]
    ranks = np.arange(1, len(counts) + 1)
    for col in meta["counts_column"]:
        cumulative = np.cumsum(np.sort(share[col].to_numpy())[::-1])
        ax.plot(ranks, cumulative, color=GROUP_COLORS[groups[col]], linewidth=1, alpha=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("Gene rank within sample (log scale)")
    ax.set_ylabel("Cumulative fraction of library")
    ax.set_title("A  How concentrated are the reads?")
    for n in (10, 100, 1000):
        frac = np.mean([np.sort(share[c])[::-1][:n].sum() for c in share.columns])
        ax.annotate(f"top {n}: {frac:.0%}", (n, frac), xytext=(6, -12),
                    textcoords="offset points", fontsize=8, color=INK)
        ax.plot(n, frac, "o", color=INK, markersize=3)
    ax.legend(handles=group_legend_handles(), loc="lower right")

    ax = axes[1]
    top = share.mean(axis=1).sort_values(ascending=False).head(15)
    labels = ann["gene_name"].reindex(top.index).fillna(pd.Series(top.index, index=top.index))
    ax.barh(range(len(top))[::-1], top.to_numpy() * 100, color="#2a78d6", height=0.7)
    ax.set_yticks(range(len(top))[::-1], labels)
    ax.set_xlabel("Mean % of library (19 samples)")
    ax.set_title("B  Fifteen most abundant transcripts")
    ax.grid(axis="y", visible=False)
    save(fig, OUT / "fig04_expression_dominance")


def fig_mean_variance(filtered, log_norm, normalized):
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    ax = axes[0]
    mean, var = normalized.mean(axis=1), normalized.var(axis=1)
    ax.scatter(mean, var, s=3, alpha=0.25, color="#2a78d6", edgecolors="none")
    lims = np.array([mean[mean > 0].min(), mean.max()])
    ax.plot(lims, lims, color=INK, linewidth=1, linestyle="--", label="Poisson: var = mean")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Mean normalised count")
    ax.set_ylabel("Variance across 19 samples")
    ax.set_title("A  Count scale: variance grows with mean")
    ax.legend(loc="upper left")

    ax = axes[1]
    ax.scatter(log_norm.mean(axis=1), log_norm.std(axis=1), s=3, alpha=0.25,
               color="#2a78d6", edgecolors="none")
    ax.set_xlabel("Mean log2(normalised count + 1)")
    ax.set_ylabel("SD of log2 expression")
    ax.set_title("B  log2 scale: low-count genes get inflated SD")
    save(fig, OUT / "fig05_mean_variance")


def fig_sample_qc(meta, audit):
    df = audit.join(meta.set_index("counts_column")[["replicate"]])
    metrics = [
        ("n_genes_ge_10", "Genes with >= 10 counts"),
        ("fraction_counts_chrM", "Fraction of counts from chrM genes"),
        ("fraction_counts_top_gene", "Fraction of counts from top gene"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.4))
    rng = np.random.default_rng(1)
    for ax, (col, label) in zip(axes, metrics):
        for i, group in enumerate(GROUP_ORDER):
            vals = df.loc[df["group"] == group, col].to_numpy()
            x = i + rng.uniform(-0.12, 0.12, size=len(vals))
            scatter_by_group(ax, x, vals, np.array([group] * len(vals)), size=40)
            ax.hlines(np.median(vals), i - 0.25, i + 0.25, color=INK, linewidth=1.5)
        ax.set_xticks(range(4), GROUP_ORDER, rotation=30, ha="right")
        ax.set_title(label, fontsize=9)
        ax.grid(axis="x", visible=False)
    fig.suptitle("Per-sample quality indicators (bar = group median)", fontweight="bold",
                 fontsize=10, y=1.03)
    save(fig, OUT / "fig06_sample_qc_by_group")


def fig_selected_genes(meta, log_norm, ann):
    symbol_to_id = {s: i for i, s in ann["gene_name"].dropna().items() if i in log_norm.index}
    present = [(g, why) for g, why in SELECTED_GENES if g in symbol_to_id]
    ncols = 4
    nrows = int(np.ceil(len(present) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(10, 2.5 * nrows), squeeze=False)
    rng = np.random.default_rng(2)
    groups = meta.set_index("counts_column")["group"]
    for ax, (gene, why) in zip(axes.ravel(), present):
        values = log_norm.loc[symbol_to_id[gene], meta["counts_column"]]
        for i, group in enumerate(GROUP_ORDER):
            vals = values[groups[values.index] == group].to_numpy()
            x = i + rng.uniform(-0.12, 0.12, size=len(vals))
            scatter_by_group(ax, x, vals, np.array([group] * len(vals)), size=30)
            ax.hlines(np.mean(vals), i - 0.25, i + 0.25, color=INK, linewidth=1.5)
        ax.set_title(f"{gene}", fontsize=9, style="italic")
        ax.text(0.02, 0.97, why, transform=ax.transAxes, fontsize=7, color=MUTED, va="top")
        ax.set_xticks(range(4), GROUP_ORDER, rotation=35, ha="right", fontsize=7)
        ax.grid(axis="x", visible=False)
    for ax in axes.ravel()[len(present):]:
        ax.set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("log2(norm. count + 1)")
    fig.suptitle("Selected genes by group (bar = group mean; descriptive only, no tests)",
                 fontweight="bold", fontsize=10, y=1.01)
    fig.tight_layout()
    save(fig, OUT / "fig07_selected_genes_by_group")
    missing = [g for g, _ in SELECTED_GENES if g not in symbol_to_id]
    return missing


def variance_share(log_norm: pd.DataFrame, labels: pd.Series) -> pd.Series:
    """Per-gene eta^2: between-label sum of squares / total sum of squares."""
    centred = log_norm.sub(log_norm.mean(axis=1), axis=0)
    total = (centred ** 2).sum(axis=1)
    group_means = log_norm.T.groupby(labels.to_numpy()).transform("mean").T
    between = (group_means.sub(log_norm.mean(axis=1), axis=0) ** 2).sum(axis=1)
    return between / total


def fig_variance_partition(meta, log_norm):
    """Per-gene share of variance tied to each design factor, vs. shuffled labels."""
    log_norm = log_norm[meta["counts_column"]]
    expressed = log_norm[log_norm.mean(axis=1) > 3]  # avoid noisy low-count genes
    factors = {"Diet (CD vs MD)": meta["diet"], "Microplastic (Veh vs PE)": meta["microplastic"],
               "Group (4 levels)": meta["group"].astype(str)}
    rng = np.random.default_rng(3)
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.2), sharey=True)
    rows = []
    for ax, (name, labels) in zip(axes, factors.items()):
        observed = variance_share(expressed, labels)
        shuffled = variance_share(expressed, pd.Series(rng.permutation(labels.to_numpy())))
        bins = np.linspace(0, 1, 41)
        ax.hist(shuffled, bins=bins, color="#c9c9c9", label="labels shuffled once")
        ax.hist(observed, bins=bins, histtype="step", color="#2a78d6", linewidth=1.8,
                label="actual labels")
        ax.set_title(name, fontsize=9)
        ax.set_xlabel("Share of gene's variance (eta²)")
        rows.append({"factor": name, "median_eta2": observed.median(),
                     "median_eta2_shuffled": shuffled.median(),
                     "n_genes_eta2_gt_0.5": int((observed > 0.5).sum()),
                     "n_genes_eta2_gt_0.5_shuffled": int((shuffled > 0.5).sum()),
                     "n_genes_considered": len(observed)})
    axes[0].set_ylabel("Number of genes")
    axes[0].legend(loc="upper right")
    fig.suptitle(f"How much of each gene's variation lines up with the design? "
                 f"({len(expressed):,} genes, mean log2 > 3)", fontweight="bold",
                 fontsize=10, y=1.04)
    save(fig, OUT / "fig08_variance_share_by_factor")
    pd.DataFrame(rows).to_csv(TABLES_DIR / "variance_share_by_factor.csv", index=False)


def fig_tissue_markers(meta, counts):
    """Per-sample reads from tissue/cell-type marker genes (contamination screen)."""
    scores = marker_set_cpm(counts[meta["counts_column"]])
    table = scores.join(meta.set_index("counts_column")[["group"]])
    table.index.name = "counts_column"
    table.round(2).to_csv(TABLES_DIR / "tissue_marker_cpm.csv")

    order = meta["counts_column"]
    groups = meta.set_index("counts_column")["group"]
    names = list(MARKER_SETS)
    fig, axes = plt.subplots(2, 4, figsize=(11, 5.6), sharex=True)
    for ax, name in zip(axes.ravel(), names):
        values = scores.loc[order, name]
        ax.bar(range(len(order)), values + 1, color=[GROUP_COLORS[g] for g in groups[order]],
               edgecolor="white", linewidth=0.8)
        ax.set_yscale("log")
        ax.set_title(name.replace("_", " "), fontsize=9)
        ax.grid(axis="x", visible=False)
        ax.set_xticks(range(len(order)), [short_name(c) for c in order], rotation=90,
                      fontsize=6)
    for ax in axes[:, 0]:
        ax.set_ylabel("Marker reads per million + 1")
    legend_ax = axes.ravel()[-1]
    legend_ax.axis("off")
    legend_ax.legend(handles=group_legend_handles(long_labels=True), loc="center")
    fig.suptitle("Tissue / cell-type marker reads per sample (log scale): contamination screen",
                 fontweight="bold", fontsize=10, y=1.0)
    fig.tight_layout()
    save(fig, OUT / "fig09_tissue_marker_screen")


def summary_tables(meta, counts, filtered, log_norm, sf, audit):
    per_group = audit.groupby("group").agg(
        n=("total_counts", "size"),
        library_size_mean_millions=("total_counts", lambda s: s.mean() / 1e6),
        library_size_min_millions=("total_counts", lambda s: s.min() / 1e6),
        library_size_max_millions=("total_counts", lambda s: s.max() / 1e6),
        genes_ge_10_mean=("n_genes_ge_10", "mean"),
        chrM_fraction_mean=("fraction_counts_chrM", "mean"),
    ).reindex(GROUP_ORDER)
    per_group.to_csv(TABLES_DIR / "library_summary_by_group.csv")

    gene_level = pd.DataFrame({
        "matrix": ["original (all genes)", "original (filtered genes)",
                   "log2 normalised (filtered genes)"],
        "n_genes": [len(counts), len(filtered), len(log_norm)],
        "median_of_gene_means": [counts.mean(axis=1).median(), filtered.mean(axis=1).median(),
                                 log_norm.mean(axis=1).median()],
        "q90_of_gene_means": [counts.mean(axis=1).quantile(0.9),
                              filtered.mean(axis=1).quantile(0.9),
                              log_norm.mean(axis=1).quantile(0.9)],
        "max_gene_mean": [counts.mean(axis=1).max(), filtered.mean(axis=1).max(),
                          log_norm.mean(axis=1).max()],
    })
    gene_level.to_csv(TABLES_DIR / "gene_level_summary.csv", index=False)

    cv = filtered.std(axis=1) / filtered.mean(axis=1)
    corr_values = log_norm.corr().to_numpy().copy()
    np.fill_diagonal(corr_values, np.nan)
    sample_corr = pd.DataFrame(corr_values, index=log_norm.columns, columns=log_norm.columns)
    pd.DataFrame({
        "mean_correlation_to_other_samples": sample_corr.mean(),
        "min_correlation_to_other_samples": sample_corr.min(),
    }).join(meta.set_index("counts_column")[["group"]]).to_csv(
        TABLES_DIR / "sample_correlation_summary.csv")
    return cv


def main() -> int:
    set_style()
    meta, counts, filtered, log_norm, sf, ann, audit = load()
    normalized = pd.read_csv(PROCESSED_DIR / "normalized_counts_filtered.csv.gz", index_col=0)
    fig_library_sizes(meta, sf)
    fig_count_distribution(meta, counts)
    fig_boxplots_before_after(meta, filtered, log_norm)
    fig_expression_dominance(meta, counts, ann)
    fig_mean_variance(filtered, log_norm, normalized)
    fig_sample_qc(meta, audit)
    missing = fig_selected_genes(meta, log_norm, ann)
    fig_variance_partition(meta, log_norm)
    fig_tissue_markers(meta, counts)
    summary_tables(meta, counts, filtered, log_norm, sf, audit)
    if missing:
        print(f"Selected genes not plotted (filtered out or absent): {missing}")
    print(f"Figures written to {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
