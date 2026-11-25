"""Section 3.1 figures: PCA (before/after muscle-gene exclusion) and fry enrichment.
Reuses the original Phase 1 PCA workflow (scripts/06_pca_clustering.py via mpliver) and the
original limma fry results (results/pathways/provisional_deposited/limma_fry_camera.csv).
No statistics are recomputed for enrichment; PCA is re-run with the identical functions and
checked against results/tables/pca_variance_explained.csv and pca_primary_scores.csv.
Run from repo root: .venv/Scripts/python.exe research/section3_1_figures/scripts/make_figures.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from mpliver.config import PROCESSED_DIR, SAMPLE_METADATA_FILE, TABLES_DIR  # noqa: E402
from mpliver.exploration import run_pca  # noqa: E402
from mpliver.normalization import top_variable_genes  # noqa: E402
from mpliver.qc_markers import genes_tracking_score, marker_set_cpm  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1]
FIG = OUT / "figures"
FIG.mkdir(exist_ok=True)
N_TOP, MUSCLE_RHO = 500, 0.7  # same constants as scripts/06_pca_clustering.py

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 300, "savefig.bbox": "tight",
})
INK = "#222222"
DIET_FILL = {"CD": "#1f4e8c", "MD": "#f4b183"}   # dark vs light: separable in greyscale
PE_MARK = {"Veh": "o", "PE": "^"}


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{name}.{ext}")
    plt.close(fig)


# ---------------------------------------------------------------- Figure 1: PCA
meta = pd.read_csv(ROOT / SAMPLE_METADATA_FILE.relative_to(ROOT)) if SAMPLE_METADATA_FILE.is_absolute() else pd.read_csv(SAMPLE_METADATA_FILE)
order = meta["counts_column"].tolist()
log_norm = pd.read_csv(PROCESSED_DIR / "log2_normalized_filtered.csv.gz", index_col=0)[order]
counts = pd.read_csv(PROCESSED_DIR / "counts_original.csv.gz", index_col=0)[order]
muscle_score = np.log2(marker_set_cpm(counts)["skeletal_muscle"] + 1)
muscle_genes = genes_tracking_score(log_norm, muscle_score, MUSCLE_RHO)
no_muscle = log_norm.drop(index=muscle_genes)
sets = {"before": log_norm.loc[top_variable_genes(log_norm, N_TOP)],
        "after": no_muscle.loc[top_variable_genes(no_muscle, N_TOP)]}
res = {k: run_pca(v, n_components=10) for k, v in sets.items()}

# verification against the original outputs
saved_ve = pd.read_csv(TABLES_DIR / "pca_variance_explained.csv").set_index("feature_set")
saved_sc = pd.read_csv(TABLES_DIR / "pca_primary_scores.csv", index_col=0)
saved_mg = pd.read_csv(TABLES_DIR / "muscle_tracking_genes.csv", index_col=0)
check = {
    "n_muscle_genes": len(muscle_genes),
    "muscle_gene_list_identical_to_saved": set(muscle_genes) == set(saved_mg.index),
    "max_abs_diff_scores_vs_saved_primary_PC1_PC2": float((res["before"]["scores"][["PC1", "PC2"]] - saved_sc[["PC1", "PC2"]]).abs().max().max()),
}
for k, name in (("before", "top500_log2norm"), ("after", "top500_log2norm_muscle_genes_removed")):
    ev = res[k]["explained_variance_ratio"] * 100
    check[f"{k}_PC1_pct"], check[f"{k}_PC2_pct"] = round(ev["PC1"], 2), round(ev["PC2"], 2)
    check[f"{k}_matches_saved_table"] = bool(abs(round(ev["PC1"], 2) - saved_ve.loc[name, "PC1_pct"]) < 0.006
                                             and abs(round(ev["PC2"], 2) - saved_ve.loc[name, "PC2_pct"]) < 0.006)

fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.0))
titles = {"before": "a  All genes (500 most variable)",
          "after": f"b  After excluding {len(muscle_genes)} muscle-associated genes"}
rows = []
for ax, k in zip(axes, ("before", "after")):
    sc, ev = res[k]["scores"], res[k]["explained_variance_ratio"] * 100
    for _, m in meta.iterrows():
        x, y = sc.loc[m.counts_column, "PC1"], sc.loc[m.counts_column, "PC2"]
        ax.scatter(x, y, s=38, marker=PE_MARK[m.microplastic], c=DIET_FILL[m.diet],
                   edgecolors=INK, linewidths=0.7, zorder=3)
        rows.append(dict(panel=k, sample=m.counts_column, diet=m.diet, treatment=m.microplastic, PC1=x, PC2=y,
                         PC1_pct=ev["PC1"], PC2_pct=ev["PC2"]))
    ax.set_xlabel(f"PC1 ({ev['PC1']:.1f}% of variance)")
    ax.set_ylabel(f"PC2 ({ev['PC2']:.1f}% of variance)")
    ax.set_title(titles[k], loc="left", fontweight="bold")
    ax.axhline(0, color="#bbbbbb", lw=0.5, zorder=1)
    ax.axvline(0, color="#bbbbbb", lw=0.5, zorder=1)
handles = [Line2D([], [], ls="", marker=PE_MARK[t], ms=6.5, mfc=DIET_FILL[d], mec=INK, mew=0.7,
                  label=f"{d} + {'vehicle' if t == 'Veh' else 'PE'} (n = {((meta.diet == d) & (meta.microplastic == t)).sum()})")
           for d in ("CD", "MD") for t in ("Veh", "PE")]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.06))
fig.tight_layout(rect=(0, 0.06, 1, 1))
save(fig, "fig_pca")
pd.DataFrame(rows).to_csv(OUT / "source_data_fig_pca.csv", index=False)

# ---------------------------------------------------------- Figure 2: fry enrichment
fr = pd.read_csv(ROOT / "results/pathways/provisional_deposited/limma_fry_camera.csv")
fr = fr[fr.contrast.isin(["PE_in_CD", "PE_in_MD"])].copy()
assert fr.groupby("contrast").size().eq(53).all() and fr.family.value_counts().to_dict() == {"Hallmark": 100, "PPAR": 6}
v = fr[(fr.pathway == "SANDERSON_PPARA_TARGETS") & (fr.contrast == "PE_in_CD")].fry_fdr.iloc[0]
check["PPARA_targets_CD_fry_FDR"] = round(float(v), 4)
assert round(v, 4) == 0.0406

SPECIAL = {"MTORC1": "mTORC1", "DNA": "DNA", "KRAS": "KRAS", "TGF": "TGF", "UV": "UV", "E2F": "E2F", "G2M": "G2M",
           "MYC": "MYC", "P53": "p53", "PI3K": "PI3K", "AKT": "AKT", "MTOR": "mTOR", "IL2": "IL2", "STAT5": "STAT5",
           "IL6": "IL6", "JAK": "JAK", "STAT3": "STAT3", "TNFA": "TNFα", "NFKB": "NF-κB", "WNT": "WNT", "V1": "v1",
           "V2": "v2", "DN": "(down)", "UP": "(up)", "BETA": "β", "ALPHA": "α", "GAMMA": "γ", "NOTCH": "NOTCH"}


def pretty(p):
    if p == "SANDERSON_PPARA_TARGETS":
        return "PPARα targets (Sanderson)"
    if p == "WP_PPAR_SIGNALING_PATHWAY":
        return "PPAR signaling pathway (WikiPathways)"
    if p == "REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA":
        return "Lipid metabolism regulation by PPARα (Reactome)"
    toks = p.replace("HALLMARK_", "").split("_")
    out = [SPECIAL.get(t, t.lower()) for t in toks]
    s = " ".join(out).replace("TGF β", "TGF-β").replace("WNT β", "WNT/β").replace("PI3K AKT mTOR", "PI3K/AKT/mTOR")
    s = s.replace("IL6 JAK STAT3", "IL6/JAK/STAT3").replace("IL2 STAT5", "IL2/STAT5")
    return s if toks[0] in SPECIAL else s[0].upper() + s[1:]


fr["label"] = fr.pathway.map(pretty)
fr["signed"] = np.where(fr.fry_direction == "Up", 1, -1) * -np.log10(fr.fry_fdr)
TH5, TH10 = -np.log10(0.05), -np.log10(0.10)


def draw(sub_paths, name, height, note):
    d = fr[fr.pathway.isin(sub_paths)]
    best = d.groupby("pathway").fry_fdr.min()
    fam = d.groupby("pathway").family.first()
    hall = best[fam == "Hallmark"].sort_values(ascending=False).index.tolist()
    ppar = best[fam == "PPAR"].sort_values(ascending=False).index.tolist()
    # bottom -> top: PPAR rows, PPAR header, empty gap row, Hallmark rows, Hallmark header.
    # Family headers sit on their own rows so they never share a line with a pathway.
    order_ = ppar + ["#PPAR"] + [None] + hall + ["#HALL"]
    ypos = {p: i for i, p in enumerate(order_) if p and not p.startswith("#")}
    headers = {order_.index("#HALL"): "Hallmark (50 sets; BH within family)",
               order_.index("#PPAR"): "PPAR-related (3 sets; BH within family)"}
    lim = max(2.0, np.ceil(d.signed.abs().max() * 2) / 2 + 0.25)
    fig, axes = plt.subplots(1, 2, figsize=(6.9, height), sharey=True)
    for ax, (c, ttl) in zip(axes, (("PE_in_CD", "a  PE vs vehicle,\n    control diet (CD)"),
                                   ("PE_in_MD", "b  PE vs vehicle,\n    metabolic-dysfunction diet (MD)"))):
        for xv, ls in ((TH5, "--"), (-TH5, "--"), (TH10, ":"), (-TH10, ":")):
            ax.axvline(xv, color="#555555", lw=0.7, ls=ls, zorder=1)
        ax.axvline(0, color=INK, lw=0.6, zorder=1)
        for _, r in d[d.contrast == c].iterrows():
            y = ypos[r.pathway]
            ax.plot([0, r.signed], [y, y], color="#9a9a9a", lw=0.8, zorder=2)
            fc = "#000000" if r.fry_fdr < 0.05 else ("#9a9a9a" if r.fry_fdr < 0.10 else "white")
            ax.scatter(r.signed, y, marker="^" if r.fry_direction == "Up" else "v", s=30, c=fc,
                       edgecolors=INK, linewidths=0.7, zorder=3)
        ax.set_xlim(-lim, lim)
        ax.set_title(ttl, loc="left", fontweight="bold")
        ax.set_xlabel("−log10(fry FDR), signed by direction\n← down in PE      up in PE →")
        ax.grid(axis="y", color="#eeeeee", lw=0.5)
        ax.axhline(order_.index(None), color="#bbbbbb", lw=0.6)
    axes[0].set_yticks(list(ypos.values()))
    axes[0].set_yticklabels([fr.loc[fr.pathway == p, "label"].iloc[0] for p in ypos])
    axes[0].tick_params(axis="y", length=0)
    axes[1].tick_params(axis="y", length=0)
    axes[0].set_ylim(-0.7, len(order_) - 0.3)
    for yy, txt in headers.items():
        axes[0].text(-0.02, yy, txt, transform=axes[0].get_yaxis_transform(), ha="right", va="center",
                     fontweight="bold", fontsize=7)
    h = [Line2D([], [], ls="", marker="^", ms=6, mfc="#000000", mec=INK, label="FDR < 0.05"),
         Line2D([], [], ls="", marker="^", ms=6, mfc="#9a9a9a", mec=INK, label="0.05 ≤ FDR < 0.10"),
         Line2D([], [], ls="", marker="^", ms=6, mfc="white", mec=INK, label="FDR ≥ 0.10"),
         Line2D([], [], ls="", marker="^", ms=6, mfc="white", mec="white", label="▲ up  ▼ down"),
         Line2D([], [], color="#555555", ls="--", lw=0.8, label="FDR = 0.05"),
         Line2D([], [], color="#555555", ls=":", lw=0.8, label="FDR = 0.10")]
    fig.legend(handles=h, loc="lower center", ncol=6, frameon=False, bbox_to_anchor=(0.5, -0.01),
               handletextpad=0.3, columnspacing=1.0)
    fig.text(0.5, -0.045, note, ha="center", fontsize=6.5, color="#555555")
    fig.tight_layout(rect=(0, 0.06 if height < 5 else 0.025, 1, 1))
    save(fig, name)
    return d


bestall = fr.groupby("pathway").fry_fdr.min()
main_paths = bestall[(bestall < 0.25) | fr.groupby("pathway").family.first().eq("PPAR")].index.tolist()
note = ("limma fry (directional, self-contained); BH-FDR within family (50 Hallmark; 3 PPAR) and contrast. "
        "camera results not shown.")
dm = draw(main_paths, "fig_enrichment", 3.6, note + " Shown: all PPAR sets and Hallmark sets with fry FDR < 0.25 in either diet.")
ds = draw(fr.pathway.unique().tolist(), "fig_enrichment_full_supplementary", 9.6, note + " All 53 sets.")
cols = ["family", "pathway", "label", "contrast", "n_genes", "fry_direction", "fry_pvalue", "fry_fdr", "signed"]
dm[cols].to_csv(OUT / "source_data_fig_enrichment.csv", index=False)
ds[cols].to_csv(OUT / "source_data_fig_enrichment_full.csv", index=False)
check["n_rows_main"] = int(dm.pathway.nunique())
check["n_rows_full"] = int(ds.pathway.nunique())
check["fry_FDR_lt_0.05"] = fr.loc[fr.fry_fdr < 0.05, ["pathway", "contrast", "fry_fdr"]].round(4).values.tolist()
check["fry_FDR_0.05_to_0.10"] = fr.loc[(fr.fry_fdr >= 0.05) & (fr.fry_fdr < 0.10), ["pathway", "contrast", "fry_fdr"]].round(4).values.tolist()
json.dump(check, open(OUT / "verification.json", "w"), indent=2, ensure_ascii=False)
print(json.dumps(check, indent=1, ensure_ascii=False))
