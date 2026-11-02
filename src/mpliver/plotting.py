"""Shared figure style so every plot uses the same colours, fonts and legends."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from mpliver.config import GROUP_COLORS, GROUP_LABELS, GROUP_MARKERS, GROUP_ORDER  # noqa: E402

INK = "#222222"
MUTED = "#6b6b6b"
GRID = "#e6e6e6"


def set_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.labelsize": 9,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.axisbelow": True,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "legend.frameon": False,
        "legend.fontsize": 8,
    })


def group_legend_handles(groups=GROUP_ORDER, long_labels: bool = False) -> list:
    labels = GROUP_LABELS if long_labels else {g: g for g in groups}
    return [
        Line2D([], [], marker=GROUP_MARKERS[g], linestyle="", markersize=7,
               markerfacecolor=GROUP_COLORS[g], markeredgecolor=INK, markeredgewidth=0.6,
               label=labels[g])
        for g in groups
    ]


def scatter_by_group(ax, x, y, groups, size: float = 55) -> None:
    """Scatter points coloured and shaped by experimental group."""
    for group in GROUP_ORDER:
        mask = groups == group
        ax.scatter(x[mask], y[mask], s=size, c=GROUP_COLORS[group],
                   marker=GROUP_MARKERS[group], edgecolors=INK, linewidths=0.6,
                   label=group, zorder=3)


def save(fig, path: Path) -> None:
    """Save as PNG (for slides/README) and PDF (vector, for reports)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    # No timestamps in the files, so re-running unchanged code gives identical outputs.
    fig.savefig(path.with_suffix(".png"), metadata={"Software": None})
    fig.savefig(path.with_suffix(".pdf"), metadata={"CreationDate": None, "Producer": None})
    plt.close(fig)
