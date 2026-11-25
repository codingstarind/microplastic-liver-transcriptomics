"""Enrichment audit, step 3: one pathway comparison figure, PE vs Veh in CD and MD side by side (+ interaction).

x: mouse-level pathway effect (mean z of VST over set genes; lm ~ diet*microplastic; 95% CI) from
results/pathways/provisional_deposited/pathway_score_effects.csv.gz.
Fill: limma fry directional FDR < 0.05 (BH within gene-set family x contrast); hollow otherwise.
Square marker (right edge, outside the axes): camera FDR < 0.05 (same family), coloured by camera direction.
Rows: 50 Hallmark sets (ordered by mean of CD and MD effect) then the 3 PPAR sets.
Writes figures/pathway_PE_vs_Veh_CD_MD_interaction.{png,pdf,svg} under research/enrichment_audit/.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
EA = REPO / "research" / "enrichment_audit"
FIG = EA / "figures"
FIG.mkdir(exist_ok=True)
w = pd.read_csv(EA / "tables" / "T2_enrichment_wide_CD_MD_interaction.csv")

UP, DOWN, NS = "#e34948", "#2a78d6", "#8c8b87"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK, "svg.fonttype": "none",
                     "pdf.fonttype": 42})


def nice(name):
    n = name.replace("HALLMARK_", "")
    special = {"REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA": "Reactome: lipid metabolism by PPARα",
               "WP_PPAR_SIGNALING_PATHWAY": "WikiPathways: PPAR signaling",
               "SANDERSON_PPARA_TARGETS": "Sanderson: PPARα targets"}
    if name in special:
        return special[name]
    n = n.replace("_", " ").lower().capitalize()
    for a, b in [("Tnfa", "TNFα"), ("nfkb", "NF-κB"), ("Il6", "IL6"), ("jak", "JAK"), ("stat3", "STAT3"),
                 ("Il2", "IL2"), ("stat5", "STAT5"), ("Kras", "KRAS"), ("dn", "down"), ("Myc", "MYC"),
                 ("v1", "V1"), ("v2", "V2"), ("E2f", "E2F"), ("G2m", "G2M"), ("Dna", "DNA"), ("Uv", "UV"),
                 ("Mtorc1", "mTORC1"), ("Pi3k akt mtor", "PI3K–AKT–mTOR"), ("Tgf beta", "TGF-β"),
                 ("Wnt beta catenin", "WNT/β-catenin"), ("Il", "IL"), ("P53", "p53"),
                 ("Interferon alpha", "Interferon-α"), ("Interferon gamma", "Interferon-γ"),
                 ("Reactive oxygen species", "Reactive oxygen species"), ("Notch", "NOTCH")]:
        n = n.replace(a, b)
    return n


h = w[w.family == "Hallmark"].copy()
h["order"] = (h["PE_in_CD__meanz_est"] + h["PE_in_MD__meanz_est"]) / 2
h = h.sort_values("order")
p = w[w.family == "PPAR"].copy().sort_values("pathway", ascending=False)
gap = 1.2
y_h = np.arange(len(h), dtype=float) + len(p) + gap
y_p = np.arange(len(p), dtype=float)
rows = pd.concat([p, h])
ys = np.concatenate([y_p, y_h])

panels = [("PE_in_CD", "PE vs Veh · control diet (CD)"),
          ("PE_in_MD", "PE vs Veh · MASH diet (MD)"),
          ("Interaction", "Diet × PE interaction (MD − CD)")]
lo = min(rows[f"{k}__meanz_ci95_low"].min() for k, _ in panels)
hi = max(rows[f"{k}__meanz_ci95_high"].max() for k, _ in panels)
pad = 0.05 * (hi - lo)
fig, axes = plt.subplots(1, 3, figsize=(10.2, 10.6), sharey=True, gridspec_kw={"wspace": 0.16})
fig.patch.set_facecolor(SURF)
for ax, (k, title) in zip(axes, panels):
    ax.set_facecolor(SURF)
    est, l, u = rows[f"{k}__meanz_est"].values, rows[f"{k}__meanz_ci95_low"].values, rows[f"{k}__meanz_ci95_high"].values
    q, d = rows[f"{k}__fry_fdr"].values, rows[f"{k}__fry_direction"].values
    qc = rows[f"{k}__camera_fdr"].values
    sig = q < 0.05
    col = np.where(sig, np.where(d == "Up", UP, DOWN), NS)
    ax.axvline(0, color=INK2, lw=0.8, zorder=1)
    for yy in ys:
        ax.axhline(yy, color=GRID, lw=0.4, zorder=0)
    ax.hlines(ys, l, u, color=col, lw=1.2, zorder=2)
    ax.scatter(est[~sig], ys[~sig], s=22, facecolor=SURF, edgecolor=NS, lw=1.1, zorder=3)
    ax.scatter(est[sig], ys[sig], s=30, facecolor=col[sig], edgecolor=SURF, lw=1.0, zorder=4)
    cam = qc < 0.05
    cd_ = rows[f"{k}__camera_direction"].values
    ccol = np.where(cd_ == "Up", UP, DOWN)
    tr = matplotlib.transforms.blended_transform_factory(ax.transAxes, ax.transData)
    ax.scatter(np.full(cam.sum(), 1.035), ys[cam], marker="s", s=12, color=ccol[cam], transform=tr,
               clip_on=False, zorder=3)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_title(title, fontsize=8.5, color=INK, loc="left", pad=6)
    ax.set_xlabel("Pathway score effect (mean z, 95% CI)", color=INK2)
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.axhline(len(p) - 1 + gap / 2 + 0.5, color=INK2, lw=0.6, ls=(0, (2, 2)))
axes[0].set_yticks(ys)
axes[0].set_yticklabels([nice(x) for x in rows.pathway])
axes[0].set_ylim(-0.8, ys.max() + 0.8)

from matplotlib.lines import Line2D
handles = [Line2D([], [], marker="o", ls="", mfc=UP, mec=SURF, ms=6, label="fry FDR < 0.05, up"),
           Line2D([], [], marker="o", ls="", mfc=DOWN, mec=SURF, ms=6, label="fry FDR < 0.05, down"),
           Line2D([], [], marker="o", ls="", mfc=SURF, mec=NS, ms=5.5, label="fry FDR ≥ 0.05"),
           Line2D([], [], marker="s", ls="", color=UP, ms=3.5, label="camera FDR < 0.05, up"),
           Line2D([], [], marker="s", ls="", color=DOWN, ms=3.5, label="camera FDR < 0.05, down")]
fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=7.5, bbox_to_anchor=(0.55, 0.005))
fig.text(0.01, 0.985, "Hallmark (50) and PPAR (3) gene sets, GSE294180 (n = 4/5/5/5 mice)", fontsize=9.5,
         color=INK, va="top", weight="bold")
fig.text(0.01, 0.967, "Points: per-mouse mean-z score effect from lm(score ~ diet × PE). Fill: limma fry directional test, "
         "BH within gene-set family × contrast. Squares at right edge: camera FDR < 0.05 (same family), coloured by camera direction. Hallmark rows ordered by mean of CD and MD effects; "
         "PPAR sets below the dashed line.", fontsize=6.6, color=INK2, va="top", wrap=True)
fig.subplots_adjust(left=0.25, right=0.975, top=0.925, bottom=0.07)
for ext in ["png", "pdf", "svg"]:
    fig.savefig(FIG / f"pathway_PE_vs_Veh_CD_MD_interaction.{ext}", dpi=300, facecolor=SURF)
print("saved")
