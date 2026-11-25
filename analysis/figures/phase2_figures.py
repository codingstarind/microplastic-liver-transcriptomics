"""Phase 2 figures (PE response by diet, PE-diet concordance, interaction contamination check).

Usage: python analysis/figures/phase2_figures.py <deseq2_label> [<covariate_label>]
"""

import itertools
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mpliver.config import REPO_ROOT, SAMPLE_METADATA_FILE, TABLES_DIR
from mpliver.plotting import INK, MUTED, save, set_style

BLUE, ORANGE, GREY = "#2a78d6", "#eb6834", "#9a9a9a"
CONTRAST_STYLE = {"PE_in_CD": ("PE effect, control diet", BLUE, "o"),
                  "PE_in_MD": ("PE effect, MASH diet", ORANGE, "^"),
                  "Interaction": ("Interaction (MD - CD)", INK, "s")}


def provisional_tag(label):
    return "  [PROVISIONAL: deposited counts, no length offset]" if "provisional" in label else ""


def fig_pathway_forest(label, out):
    pathway_dir = REPO_ROOT / "results" / "pathways" / label
    scores = pd.read_csv(pathway_dir / "pathway_score_effects.csv.gz")
    fry = pd.read_csv(pathway_dir / "limma_fry_camera.csv")
    focus = scores[(scores["theme"].notna()) & (scores["theme"] != "") &
                   scores["family"].isin(["Hallmark", "PPAR"])]
    order = (focus[focus["contrast"] == "PE_in_CD"].sort_values(["theme", "pathway"])
             ["pathway"].tolist())
    fig, ax = plt.subplots(figsize=(9, 0.42 * len(order) + 1.6))
    offsets = {"PE_in_CD": -0.22, "PE_in_MD": 0.0, "Interaction": 0.22}
    for contrast, (name, color, marker) in CONTRAST_STYLE.items():
        sub = focus[focus["contrast"] == contrast].set_index("pathway").loc[order]
        y = np.arange(len(order)) + offsets[contrast]
        ax.errorbar(sub["estimate"], y, xerr=[sub["estimate"] - sub["ci95_low"],
                                              sub["ci95_high"] - sub["estimate"]],
                    fmt=marker, color=color, ms=5, lw=1, capsize=0, label=name)
        f = fry[fry["contrast"] == contrast].set_index("pathway")
        for yi, pw in zip(y, order):
            if pw in f.index and f.loc[pw, "fry_fdr"] < 0.05:
                ax.text(sub.loc[pw, "ci95_high"] + 0.03, yi, "*", color=color, va="center")
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_yticks(range(len(order)),
                  [p.replace("HALLMARK_", "").replace("_", " ").lower() for p in order],
                  fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Effect on pathway score (mean z of set genes, mouse-level model), 95% CI")
    ax.set_title("PE response of focus pathways under each diet" + provisional_tag(label),
                 fontsize=9)
    ax.legend(loc="lower right", fontsize=7)
    fig.text(0.01, -0.02, "* fry (self-contained rotation test) FDR < 0.05 within contrast. "
             "Mean-z scores dilute sets whose genes move in opposite directions.",
             fontsize=7, color=MUTED)
    save(fig, out / "fig16_pathway_pe_response_forest")


def log2_diff(norm, a, b):
    return np.log2(norm[a] + 1).mean(axis=1) - np.log2(norm[b] + 1).mean(axis=1)


def fig_concordance(label, out):
    meta = pd.read_csv(SAMPLE_METADATA_FILE)
    norm = pd.read_csv(REPO_ROOT / "data" / "processed" / "normalized_counts_filtered.csv.gz",
                       index_col=0)[meta["counts_column"]]
    norm = norm[norm.mean(axis=1) > 50]
    g = meta.set_index("counts_column")["group"]
    m = {k: g.index[g == k].tolist() for k in g.unique()}
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, (pe_grp, ref_grp) in zip(axes, [("CD+PE", "CD+Veh"), ("MD+PE", "MD+Veh")]):
        shared = m["CD+Veh"] if ref_grp == "CD+Veh" else m["MD+Veh"]
        naive_pe = log2_diff(norm, m[pe_grp], m[ref_grp])
        naive_diet = log2_diff(norm, m["MD+Veh"], m["CD+Veh"])
        obs, null = [], []
        for part in itertools.combinations(shared, len(shared) // 2):
            a, b = list(part), [x for x in shared if x not in part]
            diet = (log2_diff(norm, m["MD+Veh"], b) if ref_grp == "CD+Veh"
                    else log2_diff(norm, b, m["CD+Veh"]))
            obs.append(log2_diff(norm, m[pe_grp], a).corr(diet, method="spearman"))
            pool = m[pe_grp] + a
            for ref in itertools.combinations(pool, len(a)):
                trt = [x for x in pool if x not in ref]
                null.append(log2_diff(norm, trt, list(ref)).corr(diet, method="spearman"))
        ax.hist(null, bins=30, color=GREY, label="relabelled PE contrast (null)")
        for i, o in enumerate(obs):
            ax.axvline(o, color=BLUE, lw=1.2, label="split-reference estimate" if i == 0 else None)
        ax.axvline(naive_pe.corr(naive_diet, method="spearman"), color=ORANGE, ls="--", lw=1.5,
                   label="naive (shared controls)")
        ax.set_title(f"{pe_grp} vs {ref_grp}  against  MD+Veh vs CD+Veh", fontsize=9)
        ax.set_xlabel("Spearman correlation of gene log2 fold changes")
        ax.set_ylabel("Number of relabellings")
        ax.legend(fontsize=7, loc="upper left")
    fig.suptitle("Does the PE response resemble the MASH-diet response?" + provisional_tag(label),
                 fontsize=10, fontweight="bold", y=1.03)
    save(fig, out / "fig17_pe_vs_diet_concordance")


def fig_interaction_scatter(labels, out):
    muscle = set(pd.read_csv(TABLES_DIR / "muscle_tracking_genes.csv").iloc[:, 0])
    fig, axes = plt.subplots(1, len(labels), figsize=(5 * len(labels), 4.6), squeeze=False)
    for ax, label in zip(axes[0], labels):
        cd = pd.read_csv(REPO_ROOT / "results" / "deseq2" / label / "results_PE_in_CD.csv.gz",
                         index_col="gene")
        md = pd.read_csv(REPO_ROOT / "results" / "deseq2" / label / "results_PE_in_MD.csv.gz",
                         index_col="gene")
        it = pd.read_csv(REPO_ROOT / "results" / "deseq2" / label / "results_Interaction.csv.gz",
                         index_col="gene")
        x, y = cd["log2FC_ashr"], md["log2FC_ashr"].reindex(cd.index)
        sig = it["padj"].reindex(cd.index) < 0.05
        mus = cd.index.isin(muscle)
        ax.scatter(x[~sig], y[~sig], s=3, color=GREY, alpha=0.4, edgecolors="none",
                   label="all genes")
        ax.scatter(x[sig & ~mus], y[sig & ~mus], s=14, color=BLUE, label="interaction FDR < 0.05")
        ax.scatter(x[sig & mus], y[sig & mus], s=18, color=ORANGE, marker="^",
                   label="... and muscle-tracking")
        lim = 3
        ax.plot([-lim, lim], [-lim, lim], color=MUTED, lw=0.8, ls=":")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("PE effect under control diet (shrunken log2FC)")
        ax.set_ylabel("PE effect under MASH diet (shrunken log2FC)")
        name = "with muscle covariate" if "cov_muscle" in label else "primary model"
        ax.set_title(f"{name}\n{int(sig.sum())} interaction genes, "
                     f"{int((sig & mus).sum())} muscle-tracking", fontsize=9)
        ax.legend(fontsize=7, loc="upper left")
    fig.suptitle("Diet x PE interaction genes and skeletal-muscle carry-over" +
                 provisional_tag(labels[0]), fontsize=10, fontweight="bold", y=1.02)
    save(fig, out / "fig18_interaction_vs_muscle")


def main() -> int:
    set_style()
    label = sys.argv[1]
    cov_label = sys.argv[2] if len(sys.argv) > 2 else None
    out = REPO_ROOT / "figures" / "phase2" / label
    fig_pathway_forest(label, out)
    fig_concordance(label, out)
    fig_interaction_scatter([label] + ([cov_label] if cov_label else []), out)
    print(f"figures written to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
