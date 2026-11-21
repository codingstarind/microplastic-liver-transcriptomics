# ruff: noqa: E501  (long captions and curated citation strings)
"""Cross-study gene replication: tables and figures.

Run from the repository root after 01_fit_panel_models.R:

    python research/cross_study_replication/scripts/02_build_tables_and_figures.py

Inputs (all read-only):
  results/deseq2/provisional_deposited/            primary factorial fit (existing)
  results/deseq2/loo/provisional_deposited__without_<mouse>/   19 leave-one-out refits (existing)
  results/deseq2/provisional_deposited__cov_<x>/   covariate-adjusted refits (existing)
  results/deseq2/covariates.csv                    technical covariates (existing)
  data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz
  research/cross_study_replication/published_gene_evidence.csv  (curated, hand-verified)
  research/cross_study_replication/intermediate/*   outputs of 01_fit_panel_models.R

Outputs go to research/cross_study_replication/ (see README.md there).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from scipy.cluster.hierarchy import cophenet, dendrogram, fcluster, linkage  # noqa: E402
from scipy.spatial.distance import pdist, squareform  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "research" / "cross_study_replication"
INTER = OUT / "intermediate"
PRIMARY = REPO / "results" / "deseq2" / "provisional_deposited"
LOO_DIR = REPO / "results" / "deseq2" / "loo"

DJOUINA = ["Fasn", "Acly", "Acaca", "Ppara"]
PS_PANEL = ["Acot3", "Abcc3", "Nr1i3", "Nr4a1"]
PANEL = DJOUINA + PS_PANEL
PANEL_NAME = {g: "Djouina 2023 (primary)" for g in DJOUINA} | {
    g: "Polystyrene studies (exploratory)" for g in PS_PANEL
}

CONTRASTS = {
    "PE_in_CD": "CD+PE vs CD+Veh (PE effect under control diet)",
    "PE_in_MD": "MD+PE vs MD+Veh (PE effect under MASH diet)",
    "Diet_in_Veh": "MD+Veh vs CD+Veh (diet effect without PE)",
    "Interaction": "(MD+PE - MD+Veh) - (CD+PE - CD+Veh), diet x PE",
}
# Published claims are about PE (or PS) versus unexposed mice on a non-MASH diet, so they map
# onto the PE contrasts only. Diet and interaction have no published counterpart.
COMPARABLE = {
    "PE_in_CD": "closest analogue (no metabolic-stress diet)",
    "PE_in_MD": "partial analogue (our MASH diet has no counterpart in the published study)",
    "PE_main_additive": "average PE effect over both diets (additive model)",
}
GROUP_ORDER = ["CD+Veh", "CD+PE", "MD+Veh", "MD+PE"]
ALPHA = 0.05
Z = 1.959963984540054

# Reference palette (dataviz skill, light mode).
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
RED, NEUTRAL = "#e34948", "#f0efec"
GROUP_COLOR = {"CD+Veh": BLUE, "CD+PE": ORANGE, "MD+Veh": AQUA, "MD+PE": YELLOW}
GROUP_MARKER = {"CD+Veh": "o", "CD+PE": "o", "MD+Veh": "s", "MD+PE": "s"}

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.edgecolor": AXIS,
        "axes.labelcolor": INK2,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "text.color": INK,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "figure.facecolor": "white",
        "axes.facecolor": SURFACE,
        "savefig.dpi": 220,
    }
)


# ---------------------------------------------------------------------------------------
def read_results(run_dir: Path, contrast: str) -> pd.DataFrame:
    return pd.read_csv(run_dir / f"results_{contrast}.csv.gz").set_index("gene")


def holm(p: pd.Series) -> pd.Series:
    """Holm step-down adjustment; NaN p-values are left out of the family."""
    out = pd.Series(np.nan, index=p.index)
    q = p.dropna().sort_values()
    m = len(q)
    adj = np.maximum.accumulate([(m - i) * v for i, v in enumerate(q.values)])
    out[q.index] = np.minimum(adj, 1.0)
    return out


def classify(lfc, padj, pub_dir, reason_na=None) -> str:
    if reason_na:
        return f"Not assessable ({reason_na})"
    if pub_dir not in ("up", "down"):
        return "Not assessable (no published direction for this contrast)"
    if pd.isna(lfc):
        return "Not assessable (no estimate)"
    if pd.isna(padj):
        return (
            "Not assessable (no adjusted p-value: removed by independent filtering or outlier flag)"
        )
    same = (lfc > 0) == (pub_dir == "up")
    supported = padj < ALPHA
    return ("Same direction, " if same else "Opposite direction, ") + (
        "statistically supported" if supported else "statistically inconclusive"
    )


def meta_table() -> pd.DataFrame:
    meta = pd.read_csv(REPO / "metadata" / "sample_metadata.csv")
    meta["mouse"] = meta["counts_column"].str.replace(r"_S\d+$", "", regex=True)
    return meta.set_index("counts_column")


# ---------------------------------------------------------------------------------------
def build_replication(meta: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    infl = pd.read_csv(PRIMARY / "sample_influence.csv").set_index("sample")
    counts = pd.read_csv(
        REPO / "data" / "raw" / "GSE294180_MP_raw_counts_gene_level.csv.gz", index_col=0
    )[infl.index]
    norm = counts.loc[PANEL] / infl["size_factor"]
    grp_means = norm.T.groupby(meta.loc[norm.columns, "group"]).mean().T[GROUP_ORDER]

    pub = evidence[evidence.use_for_comparison == "yes"].set_index("gene")
    assert pub.index.is_unique and set(pub.index) == set(PANEL)
    rows = []
    for contrast, meaning in CONTRASTS.items():
        res = read_results(PRIMARY, contrast).reindex(PANEL)
        p_holm = pd.concat([holm(res.loc[DJOUINA, "pvalue"]), holm(res.loc[PS_PANEL, "pvalue"])])
        for g in PANEL:
            r = res.loc[g]
            pub_dir = pub.loc[g, "published_direction"] if contrast.startswith("PE_") else None
            rows.append(
                {
                    "panel": PANEL_NAME[g],
                    "gene": g,
                    "contrast": contrast,
                    "contrast_meaning": meaning,
                    "published_direction": pub_dir if pub_dir else "no published claim",
                    "published_source": pub.loc[g, "source_short"],
                    "context_comparability": COMPARABLE.get(
                        contrast, "not comparable: no published estimate"
                    ),
                    "baseMean": r.baseMean,
                    **{f"mean_norm_count_{k}": grp_means.loc[g, k] for k in GROUP_ORDER},
                    "log2FC": r.log2FC,
                    "lfcSE": r.lfcSE,
                    "ci95_low": r.ci95_low,
                    "ci95_high": r.ci95_high,
                    "pvalue": r.pvalue,
                    "padj_genomewide_BH": r.padj,
                    "padj_panel_holm": p_holm[g],
                    "log2FC_ashr_shrunk": r.log2FC_ashr,
                    "fold_change": 2**r.log2FC,
                    "observed_direction": "up" if r.log2FC > 0 else "down",
                    "ci_excludes_zero": bool(r.ci95_low > 0 or r.ci95_high < 0),
                    "direction_agreement": (
                        "same" if (r.log2FC > 0) == (pub_dir == "up") else "opposite"
                    )
                    if pub_dir in ("up", "down")
                    else "not assessable",
                    "category_primary_genomewide_BH": classify(r.log2FC, r.padj, pub_dir),
                    "category_secondary_panel_holm": classify(r.log2FC, p_holm[g], pub_dir),
                }
            )
    # Additive-model average PE effect (fitted in 01_fit_panel_models.R).
    add = pd.read_csv(INTER / "additive_model_panel.csv").set_index("gene")
    p_holm = pd.concat([holm(add.loc[DJOUINA, "pvalue"]), holm(add.loc[PS_PANEL, "pvalue"])])
    for g in PANEL:
        r = add.loc[g]
        pub_dir = pub.loc[g, "published_direction"]
        rows.append(
            {
                "panel": PANEL_NAME[g],
                "gene": g,
                "contrast": "PE_main_additive",
                "contrast_meaning": "PE effect averaged over diets, ~ diet + PE (assumes no interaction)",
                "published_direction": pub_dir,
                "published_source": pub.loc[g, "source_short"],
                "context_comparability": COMPARABLE["PE_main_additive"],
                "baseMean": rows[PANEL.index(g)]["baseMean"],
                **{f"mean_norm_count_{k}": grp_means.loc[g, k] for k in GROUP_ORDER},
                "log2FC": r.log2FC,
                "lfcSE": r.lfcSE,
                "ci95_low": r.ci95_low,
                "ci95_high": r.ci95_high,
                "pvalue": r.pvalue,
                "padj_genomewide_BH": r.padj_genomewide,
                "padj_panel_holm": p_holm[g],
                "log2FC_ashr_shrunk": np.nan,
                "fold_change": 2**r.log2FC,
                "observed_direction": "up" if r.log2FC > 0 else "down",
                "ci_excludes_zero": bool(r.ci95_low > 0 or r.ci95_high < 0),
                "direction_agreement": "same"
                if (r.log2FC > 0) == (pub_dir == "up")
                else "opposite",
                "category_primary_genomewide_BH": classify(r.log2FC, r.padj_genomewide, pub_dir),
                "category_secondary_panel_holm": classify(r.log2FC, p_holm[g], pub_dir),
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------
def build_interaction() -> pd.DataFrame:
    res = read_results(PRIMARY, "Interaction").reindex(PANEL)
    cd = read_results(PRIMARY, "PE_in_CD").reindex(PANEL)
    md = read_results(PRIMARY, "PE_in_MD").reindex(PANEL)
    mus = read_results(
        REPO / "results/deseq2/provisional_deposited__cov_muscle", "Interaction"
    ).reindex(PANEL)
    lrt = pd.read_csv(INTER / "interaction_lrt_panel.csv").set_index("gene")
    ols = pd.read_csv(INTER / "ols_vst_panel.csv").set_index("gene")
    add = pd.read_csv(INTER / "additive_model_panel.csv").set_index("gene")
    boot = pd.read_csv(INTER / "bootstrap_panel.csv.gz")
    bq = boot.groupby("gene")["Interaction"].quantile([0.025, 0.975]).unstack()
    out = pd.DataFrame(
        {
            "panel": [PANEL_NAME[g] for g in PANEL],
            "gene": PANEL,
            "PE_in_CD_log2FC": cd.log2FC.values,
            "PE_in_MD_log2FC": md.log2FC.values,
            "interaction_log2FC": res.log2FC.values,
            "interaction_lfcSE": res.lfcSE.values,
            "interaction_ci95_low": res.ci95_low.values,
            "interaction_ci95_high": res.ci95_high.values,
            "interaction_wald_pvalue": res.pvalue.values,
            "interaction_padj_genomewide_BH": res.padj.values,
            "interaction_padj_panel_holm": pd.concat(
                [holm(res.loc[DJOUINA, "pvalue"]), holm(res.loc[PS_PANEL, "pvalue"])]
            ).values,
            "interaction_log2FC_ashr_shrunk": res.log2FC_ashr.values,
            "lrt_factorial_vs_additive_stat": lrt.loc[PANEL, "lrt_stat"].values,
            "lrt_pvalue": lrt.loc[PANEL, "lrt_pvalue"].values,
            "lrt_padj_genomewide_BH": lrt.loc[PANEL, "lrt_padj_genomewide"].values,
            "additive_PE_main_log2FC": add.loc[PANEL, "log2FC"].values,
            "additive_PE_main_ci95_low": add.loc[PANEL, "ci95_low"].values,
            "additive_PE_main_ci95_high": add.loc[PANEL, "ci95_high"].values,
            "ols_vst_interaction_estimate": ols.loc[PANEL, "interaction_estimate_vst"].values,
            "ols_vst_interaction_ci95_low": ols.loc[PANEL, "interaction_ci95_low_vst"].values,
            "ols_vst_interaction_ci95_high": ols.loc[PANEL, "interaction_ci95_high_vst"].values,
            "ols_interaction_F_pvalue": ols.loc[PANEL, "interaction_F_pvalue"].values,
            "ols_r2_additive": ols.loc[PANEL, "r2_additive"].values,
            "ols_r2_factorial": ols.loc[PANEL, "r2_factorial"].values,
            "ols_adj_r2_additive": ols.loc[PANEL, "adj_r2_additive"].values,
            "ols_adj_r2_factorial": ols.loc[PANEL, "adj_r2_factorial"].values,
            "ols_aic_additive": ols.loc[PANEL, "aic_additive"].values,
            "ols_aic_factorial": ols.loc[PANEL, "aic_factorial"].values,
            "ols_residual_shapiro_pvalue": ols.loc[PANEL, "shapiro_pvalue_residuals"].values,
            "bootstrap_interaction_p2.5": bq.loc[PANEL, 0.025].values,
            "bootstrap_interaction_p97.5": bq.loc[PANEL, 0.975].values,
            "muscle_adjusted_interaction_log2FC": mus.log2FC.values,
            "muscle_adjusted_interaction_padj": mus.padj.values,
        }
    )
    out["interaction_ci_excludes_zero"] = (out.interaction_ci95_low > 0) | (
        out.interaction_ci95_high < 0
    )
    out["interpretation"] = np.where(
        out.interaction_padj_genomewide_BH < ALPHA,
        "PE effect differs between diets (genome-wide FDR < 0.05)",
        np.where(
            out.interaction_wald_pvalue < ALPHA,
            "nominal evidence only (raw p < 0.05, not FDR-significant)",
            "no evidence that the PE effect differs between diets; interval is wide",
        ),
    )
    return out


# ---------------------------------------------------------------------------------------
def build_loo(meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    full = {
        c: read_results(PRIMARY, c).reindex(PANEL) for c in ("PE_in_CD", "PE_in_MD", "Interaction")
    }
    rows = []
    for d in sorted(LOO_DIR.glob("provisional_deposited__without_*")):
        mouse = d.name.split("__without_")[1]
        grp = meta.loc[meta.mouse == mouse, "group"].iloc[0]
        for c, f in full.items():
            r = read_results(d, c).reindex(PANEL)
            for g in PANEL:
                rows.append(
                    {
                        "panel": PANEL_NAME[g],
                        "gene": g,
                        "contrast": c,
                        "omitted_mouse": mouse,
                        "omitted_group": grp,
                        "full_log2FC": f.loc[g, "log2FC"],
                        "full_padj": f.loc[g, "padj"],
                        "loo_log2FC": r.loc[g, "log2FC"],
                        "loo_lfcSE": r.loc[g, "lfcSE"],
                        "loo_ci95_low": r.loc[g, "ci95_low"],
                        "loo_ci95_high": r.loc[g, "ci95_high"],
                        "loo_pvalue": r.loc[g, "pvalue"],
                        "loo_padj": r.loc[g, "padj"],
                        "delta_log2FC": r.loc[g, "log2FC"] - f.loc[g, "log2FC"],
                        "same_direction_as_full": bool(
                            np.sign(r.loc[g, "log2FC"]) == np.sign(f.loc[g, "log2FC"])
                        ),
                        "loo_ci_excludes_zero": bool(
                            r.loc[g, "ci95_low"] > 0 or r.loc[g, "ci95_high"] < 0
                        ),
                        "loo_padj_lt_0.05": bool(r.loc[g, "padj"] < ALPHA)
                        if pd.notna(r.loc[g, "padj"])
                        else False,
                    }
                )
    loo = pd.DataFrame(rows)
    assert loo.omitted_mouse.nunique() == 19, "expected 19 leave-one-out refits"

    def summarise(s: pd.DataFrame) -> pd.Series:
        i = s.delta_log2FC.abs().idxmax()
        return pd.Series(
            {
                "full_log2FC": s.full_log2FC.iloc[0],
                "full_padj": s.full_padj.iloc[0],
                "n_refits": len(s),
                "n_same_direction": int(s.same_direction_as_full.sum()),
                "n_ci_excludes_zero": int(s.loo_ci_excludes_zero.sum()),
                "n_padj_lt_0.05": int(s["loo_padj_lt_0.05"].sum()),
                "loo_log2FC_min": s.loo_log2FC.min(),
                "loo_log2FC_max": s.loo_log2FC.max(),
                "most_influential_mouse": s.loc[i, "omitted_mouse"],
                "max_abs_delta_log2FC": abs(s.loc[i, "delta_log2FC"]),
                "max_abs_delta_over_full_SE": abs(s.loc[i, "delta_log2FC"]) / s.full_SE.iloc[0],
            }
        )

    se = pd.concat({c: f["lfcSE"] for c, f in full.items()}, names=["contrast", "gene"]).rename(
        "full_SE"
    )
    loo2 = loo.join(se, on=["contrast", "gene"])
    summ = loo2.groupby(["panel", "gene", "contrast"], sort=False).apply(
        summarise, include_groups=False
    )
    summ = summ.reset_index()
    summ["single_mouse_drives_conclusion"] = np.where(
        (summ.full_padj < ALPHA) & (summ["n_padj_lt_0.05"] < summ.n_refits),
        "significance lost in some refit",
        np.where(
            (summ.full_padj >= ALPHA) & (summ["n_padj_lt_0.05"] > 0),
            "significance gained in some refit",
            np.where(summ.n_same_direction < summ.n_refits, "direction flips in some refit", "no"),
        ),
    )
    return loo, summ


# ---------------------------------------------------------------------------------------
def build_covariate_sensitivity(meta: pd.DataFrame) -> pd.DataFrame:
    variants = {
        "primary (no covariate)": PRIMARY,
        "+ sequencing depth": REPO / "results/deseq2/provisional_deposited__cov_depth",
        "+ skeletal-muscle score": REPO / "results/deseq2/provisional_deposited__cov_muscle",
        "+ clock-phase score": REPO / "results/deseq2/provisional_deposited__cov_clock_phase",
        "+ yield (quality) proxy": REPO / "results/deseq2/provisional_deposited__cov_quality_proxy",
        "drop MD_1 (muscle-contaminated)": REPO / "results/deseq2/provisional_deposited__drop_MD_1",
    }
    cov = pd.read_csv(REPO / "results/deseq2/covariates.csv").set_index("counts_column")
    # How much of each covariate's variance the 2x2 design explains (R^2 of covariate ~ group).
    grp = meta.loc[cov.index, "group"]
    r2 = {}
    for c in ["depth", "muscle", "clock_phase", "quality_proxy"]:
        y = cov[c]
        r2[c] = (
            1 - ((y - y.groupby(grp).transform("mean")) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        )
        pe = meta.loc[cov.index, "microplastic"] == "PE"
        r2[c + "_pe_only"] = np.corrcoef(y, pe.astype(float))[0, 1] ** 2
    r2_of = {
        "+ sequencing depth": "depth",
        "+ skeletal-muscle score": "muscle",
        "+ clock-phase score": "clock_phase",
        "+ yield (quality) proxy": "quality_proxy",
    }
    rows = []
    for v, d in variants.items():
        for c in ("PE_in_CD", "PE_in_MD", "Interaction"):
            r = read_results(d, c).reindex(PANEL)
            for g in PANEL:
                k = r2_of.get(v)
                rows.append(
                    {
                        "panel": PANEL_NAME[g],
                        "gene": g,
                        "contrast": c,
                        "variant": v,
                        "log2FC": r.loc[g, "log2FC"],
                        "lfcSE": r.loc[g, "lfcSE"],
                        "ci95_low": r.loc[g, "ci95_low"],
                        "ci95_high": r.loc[g, "ci95_high"],
                        "pvalue": r.loc[g, "pvalue"],
                        "padj": r.loc[g, "padj"],
                        "covariate_R2_explained_by_group": r2[k] if k else np.nan,
                        "covariate_R2_with_PE_indicator": r2[k + "_pe_only"] if k else np.nan,
                    }
                )
    out = pd.DataFrame(rows)
    base = out[out.variant == "primary (no covariate)"].set_index(["gene", "contrast"])
    out["primary_log2FC"] = [base.loc[(g, c), "log2FC"] for g, c in zip(out.gene, out.contrast)]
    out["same_direction_as_primary"] = np.sign(out.log2FC) == np.sign(out.primary_log2FC)
    out["se_inflation_vs_primary"] = out.lfcSE / [
        base.loc[(g, c), "lfcSE"] for g, c in zip(out.gene, out.contrast)
    ]
    return out


def build_bootstrap_summary() -> pd.DataFrame:
    boot = pd.read_csv(INTER / "bootstrap_panel.csv.gz")
    ref = pd.read_csv(INTER / "bootstrap_reference_fit.csv").set_index("gene")
    rows = []
    for g in PANEL:
        b = boot[boot.gene == g]
        for c in CONTRASTS:
            est = ref.loc[g, c]
            rows.append(
                {
                    "panel": PANEL_NAME[g],
                    "gene": g,
                    "contrast": c,
                    "reference_log2FC": est,
                    "n_bootstrap": len(b),
                    "boot_median": b[c].median(),
                    "boot_p2.5": b[c].quantile(0.025),
                    "boot_p97.5": b[c].quantile(0.975),
                    "boot_sd": b[c].std(),
                    "fraction_same_sign_as_reference": float(
                        (np.sign(b[c]) == np.sign(est)).mean()
                    ),
                }
            )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------
def forest_plot(rep: pd.DataFrame) -> None:
    panels = [
        ("PE_in_CD", "PE effect, control diet"),
        ("PE_in_MD", "PE effect, MASH diet"),
        ("Interaction", "Diet × PE interaction"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(11, 5.0), sharey=True)
    ypos = {g: i for i, g in enumerate(reversed(PANEL))}
    pub = rep.drop_duplicates("gene").set_index("gene")
    lim = (
        max(
            rep.loc[rep.contrast.isin([p for p, _ in panels]), "ci95_high"].abs().max(),
            rep.loc[rep.contrast.isin([p for p, _ in panels]), "ci95_low"].abs().max(),
        )
        + 0.2
    )
    for ax, (c, title) in zip(axes, panels):
        d = rep[rep.contrast == c].set_index("gene")
        ax.axvline(0, color=AXIS, lw=1, zorder=1)
        ax.axhspan(-0.5, len(PS_PANEL) - 0.5, color=NEUTRAL, zorder=0)
        for g in PANEL:
            r = d.loc[g]
            y = ypos[g]
            sig = pd.notna(r.padj_genomewide_BH) and r.padj_genomewide_BH < ALPHA
            ax.plot(
                [r.ci95_low, r.ci95_high],
                [y, y],
                color=BLUE,
                lw=2,
                solid_capstyle="round",
                zorder=2,
            )
            ax.scatter(
                r.log2FC,
                y,
                s=55,
                zorder=3,
                color=BLUE if sig else "white",
                edgecolors=BLUE,
                linewidths=2,
            )
        ax.set_title(title, loc="left")
        ax.set_xlim(-lim, lim)
        ax.grid(axis="x", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.set_xlabel("log2 fold change (95% Wald CI)")
    labels = []
    for g in reversed(PANEL):
        p = pub.loc[g, "published_direction"]
        arrow = {"up": "↑", "down": "↓"}.get(p, "?")
        labels.append(f"{g}  (published {arrow})")
    axes[0].set_yticks(range(len(PANEL)), labels)
    axes[0].text(
        -0.78,
        (len(PANEL) - 1 + len(PS_PANEL)) / 2,
        "Djouina 2023\n(primary)",
        transform=axes[0].get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=8,
        color=MUTED,
        rotation=90,
    )
    axes[0].text(
        -0.78,
        (len(PS_PANEL) - 1) / 2,
        "PS studies\n(exploratory)",
        transform=axes[0].get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=8,
        color=MUTED,
        rotation=90,
    )
    h1 = axes[2].scatter(
        [], [], s=55, color=BLUE, edgecolors=BLUE, linewidths=2, label="genome-wide FDR < 0.05"
    )
    h2 = axes[2].scatter(
        [], [], s=55, color="white", edgecolors=BLUE, linewidths=2, label="FDR ≥ 0.05 or not tested"
    )
    fig.legend(
        handles=[h1, h2],
        loc="upper left",
        frameon=False,
        ncol=2,
        bbox_to_anchor=(0.005, 0.955),
        fontsize=8,
    )
    fig.suptitle(
        "Published lipid/xenobiotic genes in GSE294180 (primary DESeq2 factorial model)",
        x=0.01,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    fig.text(
        0.01,
        0.005,
        "Positive = higher with PE (or larger PE effect under MASH diet). Arrows give the published "
        "direction for PE/PS vs unexposed; the interaction has no published counterpart.\n"
        "Deposited count matrix (provisional). Intervals are unadjusted 95% Wald CIs. n = 4/5/5/5 mice.",
        fontsize=7.5,
        color=INK2,
    )
    fig.tight_layout(rect=(0.05, 0.06, 1, 0.9))
    fig.savefig(OUT / "gene_replication_forest.png")
    plt.close(fig)


def standardized_panel(meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    vst = pd.read_csv(PRIMARY / "vst_matrix.csv.gz").set_index("gene").loc[PANEL]
    order = sorted(vst.columns, key=lambda s: (GROUP_ORDER.index(meta.loc[s, "group"]), s))
    vst = vst[order]
    zs = vst.sub(vst.mean(axis=1), axis=0).div(vst.std(axis=1, ddof=1), axis=0)
    return zs, meta.loc[order, "group"]


def heatmap_and_clustering(meta: pd.DataFrame) -> pd.DataFrame:
    zs, grp = standardized_panel(meta)
    X = zs.values
    # Three similarity measures for clustering genes.
    d_pear = pdist(X, metric="correlation")  # 1 - Pearson r
    rho = spearmanr(X, axis=1).statistic
    d_spear = squareform(1 - rho, checks=False)  # 1 - Spearman rho
    d_euc = pdist(X, metric="euclidean")  # on z-scores: monotone in 1 - r
    # Within-group (residual) correlation: removes the shared treatment means.
    resid = zs.T.groupby(grp.values).transform(lambda s: s - s.mean()).T
    d_resid = pdist(resid.values, metric="correlation")
    rows = []
    trees = {}
    for name, d in [
        ("1 - Pearson r (primary)", d_pear),
        ("1 - Spearman rho", d_spear),
        ("Euclidean on z-scores", d_euc),
        ("1 - within-group Pearson r (treatment means removed)", d_resid),
    ]:
        for method in ("average", "complete"):
            Z_ = linkage(d, method=method)
            trees[(name, method)] = Z_
            coph = cophenet(Z_, d)[0]
            for k in (2, 3):
                lab = fcluster(Z_, k, criterion="maxclust")
                clusters = {}
                for g, lb in zip(PANEL, lab):
                    clusters.setdefault(lb, []).append(g)
                rows.append(
                    {
                        "distance": name,
                        "linkage": method,
                        "k": k,
                        "cophenetic_correlation": coph,
                        "clusters": " | ".join(sorted(",".join(v) for v in clusters.values())),
                    }
                )
    clus = pd.DataFrame(rows)
    corr_total = pd.DataFrame(np.corrcoef(X), index=PANEL, columns=PANEL)
    corr_within = pd.DataFrame(np.corrcoef(resid.values), index=PANEL, columns=PANEL)
    corr_total.to_csv(INTER / "gene_correlation_total.csv")
    corr_within.to_csv(INTER / "gene_correlation_within_group.csv")

    # Figure: z-score heatmap, genes clustered by 1 - Pearson r, average linkage.
    Zp = trees[("1 - Pearson r (primary)", "average")]
    fig = plt.figure(figsize=(10, 4.6))
    gs = fig.add_gridspec(
        3,
        4,
        width_ratios=[1.3, 9, 0.9, 0.25],
        height_ratios=[0.35, 0.35, 8],
        wspace=0.03,
        hspace=0.06,
    )
    axd = fig.add_subplot(gs[2, 0])
    dn = dendrogram(
        Zp,
        orientation="left",
        no_labels=True,
        ax=axd,
        color_threshold=0,
        above_threshold_color=INK2,
    )
    axd.axis("off")
    leaf = [PANEL[i] for i in dn["leaves"]]  # bottom -> top
    axh = fig.add_subplot(gs[2, 1])
    cmap = LinearSegmentedColormap.from_list("div", [BLUE, NEUTRAL, RED])
    M = zs.loc[leaf].values
    im = axh.imshow(
        M, aspect="auto", cmap=cmap, vmin=-2.5, vmax=2.5, origin="lower", interpolation="nearest"
    )
    axh.set_yticks(range(len(leaf)), leaf)
    axh.yaxis.tick_right()
    axh.set_xticks(
        range(zs.shape[1]), [meta.loc[s, "mouse"] for s in zs.columns], rotation=90, fontsize=7
    )
    for i in range(1, 4):
        axh.axvline((grp.values == GROUP_ORDER[i]).argmax() - 0.5, color="white", lw=2)
    for s in axh.spines.values():
        s.set_visible(False)
    diet_c = {"CD": "#c3c2b7", "MD": INK2}
    pe_c = {"Veh": "#e1e0d9", "PE": ORANGE}
    for row, (lab, f, cmap_) in enumerate([("Diet", "diet", diet_c), ("PE", "microplastic", pe_c)]):
        ax = fig.add_subplot(gs[row, 1], sharex=axh)
        vals = [meta.loc[s, f] for s in zs.columns]
        for j, v in enumerate(vals):
            ax.add_patch(plt.Rectangle((j - 0.5, 0), 1, 1, color=cmap_[v], lw=0))
        ax.set_ylim(0, 1)
        ax.set_yticks([0.5], [lab])
        ax.tick_params(axis="x", bottom=False, labelbottom=False)
        ax.tick_params(axis="y", length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        for j0, v in enumerate(vals):
            if j0 == 0 or vals[j0 - 1] != v:
                ax.text(
                    j0 - 0.4,
                    0.5,
                    v,
                    va="center",
                    fontsize=7,
                    color="white" if cmap_[v] in (INK2, ORANGE) else INK,
                )
    cax = fig.add_subplot(gs[2, 3])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("z-score of VST expression", color=INK2)
    cb.outline.set_visible(False)
    fig.suptitle(
        "Published candidate genes: gene-standardised VST expression, genes clustered by 1 − Pearson r "
        "(average linkage)",
        x=0.01,
        ha="left",
        fontsize=10,
        fontweight="bold",
    )
    fig.text(
        0.01,
        0.01,
        "Mice ordered by group (CD+Veh, CD+PE, MD+Veh, MD+PE), not clustered. With 8 genes and "
        "19 mice the dendrogram is descriptive only.",
        fontsize=7.5,
        color=INK2,
    )
    fig.subplots_adjust(left=0.03, right=0.92, top=0.88, bottom=0.2)
    fig.savefig(OUT / "gene_expression_heatmap.png")
    plt.close(fig)
    return clus


def panel_pca(meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    zs, grp = standardized_panel(meta)
    X = zs.T.values  # mice x genes, columns centred & scaled
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    var = S**2 / (S**2).sum()
    scores = pd.DataFrame(U * S, index=zs.columns, columns=[f"PC{i + 1}" for i in range(len(S))])
    loadings = pd.DataFrame(Vt.T, index=PANEL, columns=scores.columns)
    # Orient each PC so that Fasn loads positively (sign of a PC is arbitrary).
    for pc in scores.columns:
        if loadings.loc["Fasn", pc] < 0:
            loadings[pc] *= -1
            scores[pc] *= -1
    # Share of each PC's variance attributable to diet, PE, diet x PE and residual (individual
    # mice), from sequential (type I) sums of squares in the order diet, PE, interaction.
    diet = (meta.loc[zs.columns, "diet"] == "MD").astype(float).values
    pe = (meta.loc[zs.columns, "microplastic"] == "PE").astype(float).values
    designs = [np.ones((len(diet), 1))]
    designs.append(np.column_stack([designs[-1], diet]))
    designs.append(np.column_stack([designs[-1], pe]))
    designs.append(np.column_stack([designs[-1], diet * pe]))
    rows = []
    for i, pc in enumerate(scores.columns[:4]):
        y = scores[pc].values
        rss = [((y - D @ np.linalg.lstsq(D, y, rcond=None)[0]) ** 2).sum() for D in designs]
        tss = rss[0]
        rows.append(
            {
                "PC": pc,
                "variance_explained": var[i],
                "share_diet": (rss[0] - rss[1]) / tss,
                "share_PE": (rss[1] - rss[2]) / tss,
                "share_diet_x_PE": (rss[2] - rss[3]) / tss,
                "share_residual_individual": rss[3] / tss,
                **{f"loading_{g}": loadings.loc[g, pc] for g in PANEL},
            }
        )
    pcs = pd.DataFrame(rows)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.6), gridspec_kw={"width_ratios": [1.25, 1]})
    for gname in GROUP_ORDER:
        m = grp.values == gname
        a1.scatter(
            scores.loc[m, "PC1"],
            scores.loc[m, "PC2"],
            s=60,
            color=GROUP_COLOR[gname],
            marker=GROUP_MARKER[gname],
            edgecolors="white",
            linewidths=1.5,
            label=gname,
            zorder=3,
        )
    for s in scores.index:
        a1.annotate(
            meta.loc[s, "mouse"],
            (scores.loc[s, "PC1"], scores.loc[s, "PC2"]),
            fontsize=6,
            color=MUTED,
            xytext=(4, 3),
            textcoords="offset points",
        )
    a1.set_xlabel(f"PC1 ({var[0] * 100:.0f}% of panel variance)")
    a1.set_ylabel(f"PC2 ({var[1] * 100:.0f}%)")
    a1.axhline(0, color=GRID, lw=0.8, zorder=1)
    a1.axvline(0, color=GRID, lw=0.8, zorder=1)
    a1.legend(frameon=False, fontsize=8, loc="best")
    a1.set_title("Mouse scores", loc="left")
    for sp in ("top", "right"):
        a1.spines[sp].set_visible(False)
    yy = np.arange(len(PANEL))[::-1]
    a2.barh(yy + 0.18, loadings["PC1"], height=0.34, color=BLUE, label="PC1")
    a2.barh(yy - 0.18, loadings["PC2"], height=0.34, color=ORANGE, label="PC2")
    a2.set_yticks(yy, PANEL)
    a2.axvline(0, color=AXIS, lw=1)
    a2.set_xlabel("loading")
    a2.set_title("Gene loadings", loc="left")
    a2.legend(frameon=False, fontsize=8, loc="lower right")
    a2.grid(axis="x", color=GRID, lw=0.6)
    a2.set_axisbelow(True)
    for sp in ("top", "right", "left"):
        a2.spines[sp].set_visible(False)
    a2.tick_params(axis="y", length=0)
    r = pcs.iloc[0]
    r2 = pcs.iloc[1]
    fig.suptitle(
        "Exploratory PCA restricted to the 8 published genes (descriptive only)",
        x=0.01,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    fig.text(
        0.01,
        0.01,
        f"Share of PC1 variance: diet {r.share_diet:.0%}, PE {r.share_PE:.0%}, diet×PE {r.share_diet_x_PE:.0%}, "
        f"residual (individual mice) {r.share_residual_individual:.0%}.  PC2: diet {r2.share_diet:.0%}, "
        f"PE {r2.share_PE:.0%}, diet×PE {r2.share_diet_x_PE:.0%}, residual {r2.share_residual_individual:.0%}. "
        "Type I sums of squares (unbalanced 4/5/5/5).\nInput: VST (blind = FALSE) values, each gene z-scored. "
        "Eight genes chosen from the literature: not a validation or classification test.",
        fontsize=7.5,
        color=INK2,
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    fig.savefig(OUT / "gene_panel_pca.png")
    plt.close(fig)
    scores.round(4).to_csv(INTER / "panel_pca_scores.csv")
    return pcs, loadings


def main() -> None:
    meta = meta_table()
    evidence = pd.read_csv(OUT / "published_gene_evidence.csv")
    rep = build_replication(meta, evidence)
    rep.to_csv(OUT / "gene_replication_results.csv", index=False, float_format="%.6g")
    build_interaction().to_csv(OUT / "interaction_results.csv", index=False, float_format="%.6g")
    loo, loo_summary = build_loo(meta)
    loo.to_csv(OUT / "leave_one_out_results.csv", index=False, float_format="%.6g")
    loo_summary.to_csv(OUT / "leave_one_out_summary.csv", index=False, float_format="%.6g")
    build_covariate_sensitivity(meta).to_csv(
        OUT / "covariate_sensitivity.csv", index=False, float_format="%.6g"
    )
    build_bootstrap_summary().to_csv(
        OUT / "bootstrap_summary.csv", index=False, float_format="%.6g"
    )
    forest_plot(rep)
    heatmap_and_clustering(meta).to_csv(
        OUT / "clustering_comparison.csv", index=False, float_format="%.4g"
    )
    pcs, _ = panel_pca(meta)
    pcs.to_csv(OUT / "panel_pca_summary.csv", index=False, float_format="%.4g")
    print("done")


if __name__ == "__main__":
    main()
