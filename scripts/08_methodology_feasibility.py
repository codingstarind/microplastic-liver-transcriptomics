"""Step 8: feasibility checks for the Phase 2 methodology proposal.

No hypothesis tests are run and no p-values are computed for genes. This script
only describes the design and the noise level, so we can judge which
comparisons the data can support:

    1. sequencing depth vs gene-assigned counts (where the library-size gap arises)
    2. confounding / collinearity between design factors and technical variables
    3. gene-wise dispersion (biological variability) from within-group variances
    4. minimum detectable log2 fold change for each candidate comparison

Outputs go to results/methodology/ (existing Phase 1 results are not touched).
"""

import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

from mpliver.config import METADATA_DIR, PROCESSED_DIR, REPO_ROOT, SAMPLE_METADATA_FILE
from mpliver.qc_markers import marker_set_cpm

OUT = REPO_ROOT / "results" / "methodology"
CELLS = ["CD+Veh", "CD+PE", "MD+Veh", "MD+PE"]


def sequencing_depth_table(meta, counts):
    runs = pd.read_csv(METADATA_DIR / "sra_run_info.csv")
    df = meta.merge(runs, on="geo_accession")
    df["gene_assigned_counts"] = counts.sum().reindex(df["counts_column"]).to_numpy()
    df["assigned_fraction"] = df["gene_assigned_counts"] / df["read_pairs"]
    return df[["counts_column", "group", "library_name", "read_pairs", "gene_assigned_counts",
               "assigned_fraction", "flowcell", "run_number"]]


def design_matrix(meta):
    diet = (meta["diet"] == "MD").astype(float).to_numpy()
    pe = (meta["microplastic"] == "PE").astype(float).to_numpy()
    return pd.DataFrame({"intercept": 1.0, "diet_MD": diet, "PE": pe,
                         "diet_MD:PE": diet * pe}, index=meta["counts_column"])


def variance_inflation(X: pd.DataFrame, extra: pd.Series) -> dict:
    """R^2 of an extra covariate regressed on the design, and its VIF."""
    beta, *_ = np.linalg.lstsq(X.to_numpy(), extra.to_numpy(), rcond=None)
    resid = extra.to_numpy() - X.to_numpy() @ beta
    r2 = 1 - resid.var() / extra.to_numpy().var()
    return {"r2_explained_by_design": float(r2), "vif": float(1 / (1 - r2)) if r2 < 1 else None}


def method_of_moments_dispersion(norm_counts, groups):
    """Pooled within-group NB dispersion: (var - mean) / mean^2, per gene."""
    rows = []
    for g in CELLS:
        sub = norm_counts.loc[:, groups == g]
        rows.append(((sub.var(axis=1) - sub.mean(axis=1)) / sub.mean(axis=1) ** 2)
                    * (sub.shape[1] - 1))
    df_total = sum((groups == g).sum() - 1 for g in CELLS)
    return (sum(rows) / df_total).clip(lower=1e-4)


def minimum_detectable_lfc(mean_count, dispersion, n_cells, contrast, alpha, power=0.8):
    """Smallest |log2FC| detectable with a Wald test on the log-mean scale.

    Var(log mu_hat) for a cell with n mice is approximately (1/mu + phi) / n.
    The contrast's variance is sum(c_i^2 * Var_i). This is an approximation
    (it ignores dispersion-estimation error and shrinkage) meant for planning only.
    """
    var_cell = (1 / mean_count + dispersion) / np.asarray(n_cells, dtype=float)
    se = np.sqrt(np.sum(np.asarray(contrast) ** 2 * var_cell))
    return (norm.ppf(1 - alpha / 2) + norm.ppf(power)) * se / np.log(2)


CONTRASTS = {
    "PE within CD (CD+PE vs CD+Veh)": [-1, 1, 0, 0],
    "PE within MD (MD+PE vs MD+Veh)": [0, 0, -1, 1],
    "Diet within Veh (MD+Veh vs CD+Veh)": [-1, 0, 1, 0],
    "PE main effect (average over diets)": [-0.5, 0.5, -0.5, 0.5],
    "Diet x PE interaction (difference of PE effects)": [1, -1, -1, 1],
}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(SAMPLE_METADATA_FILE)
    order = meta["counts_column"]
    counts = pd.read_csv(PROCESSED_DIR / "counts_original.csv.gz", index_col=0)[order]
    norm_counts = pd.read_csv(PROCESSED_DIR / "normalized_counts_filtered.csv.gz",
                              index_col=0)[order]
    groups = meta.set_index("counts_column")["group"].reindex(order).to_numpy()
    summary = {}

    # 1. Sequencing depth vs assignment ----------------------------------------
    depth = sequencing_depth_table(meta, counts)
    depth.to_csv(OUT / "sequencing_depth_by_sample.csv", index=False)
    by_group = depth.groupby("group")[["read_pairs", "gene_assigned_counts",
                                       "assigned_fraction"]].agg(["mean", "min", "max"])
    by_group.reindex(CELLS).round(3).to_csv(OUT / "sequencing_depth_by_group.csv")
    summary["flowcells_seen"] = sorted(depth["flowcell"].dropna().unique().tolist())
    summary["n_runs_with_header"] = int(depth["flowcell"].notna().sum())

    # 2. Confounding and collinearity -------------------------------------------
    X = design_matrix(meta)
    summary["design_rank"] = int(np.linalg.matrix_rank(X.to_numpy()))
    summary["design_columns"] = list(X.columns)
    muscle = np.log2(marker_set_cpm(counts)["skeletal_muscle"] + 1).reindex(order)
    technical = pd.DataFrame({
        "log2_read_pairs": np.log2(depth.set_index("counts_column")["read_pairs"]),
        "assigned_fraction": depth.set_index("counts_column")["assigned_fraction"],
        "log2_muscle_cpm": muscle,
        "library_index": meta.set_index("counts_column")["library_name"]
                             .str.lstrip("S").astype(int),
    }).reindex(order)
    technical["group"] = groups
    technical.round(4).to_csv(OUT / "technical_covariates_by_sample.csv")
    summary["technical_vs_design"] = {
        col: variance_inflation(X, technical[col]) for col in
        ["log2_read_pairs", "assigned_fraction", "log2_muscle_cpm", "library_index"]
    }
    pe = (meta["microplastic"] == "PE").to_numpy()
    reads = technical["log2_read_pairs"].to_numpy()
    summary["read_pairs_overlap_PE_vs_Veh"] = {
        "max_PE": float(2 ** reads[pe].max()), "min_Veh": float(2 ** reads[~pe].min()),
        "perfect_separation": bool(reads[pe].max() < reads[~pe].min()),
    }
    contaminated = muscle > np.log2(100 + 1)
    summary["muscle_cpm_gt_100_by_group"] = (
        pd.Series(contaminated.to_numpy(), index=groups).groupby(level=0).sum()
        .reindex(CELLS).astype(int).to_dict())

    # 3. Dispersion ---------------------------------------------------------------
    disp = method_of_moments_dispersion(norm_counts, pd.Series(groups, index=order))
    disp_no_md1 = method_of_moments_dispersion(
        norm_counts.drop(columns=[c for c in order if c.startswith("MD_1_")]),
        pd.Series(groups, index=order).drop([c for c in order if c.startswith("MD_1_")]))
    mean_expr = norm_counts.mean(axis=1)
    bins = pd.cut(mean_expr, [10, 50, 100, 500, 1000, 5000, 1e9],
                  labels=["10-50", "50-100", "100-500", "500-1k", "1k-5k", ">5k"])
    disp_table = pd.DataFrame({
        "n_genes": bins.value_counts().sort_index(),
        "median_dispersion": disp.groupby(bins, observed=True).median(),
        "median_dispersion_without_MD_1": disp_no_md1.groupby(bins, observed=True).median(),
    })
    disp_table["median_biological_CV"] = np.sqrt(disp_table["median_dispersion"])
    disp_table.index.name = "mean_normalised_count"
    disp_table.round(4).to_csv(OUT / "dispersion_by_expression.csv")
    summary["overall_median_dispersion"] = float(disp.median())
    summary["overall_median_dispersion_without_MD_1"] = float(disp_no_md1.median())

    # 4. Minimum detectable effects --------------------------------------------
    n_all = [4, 5, 5, 5]
    n_drop_md1 = [4, 5, 4, 5]
    rows = []
    for label, c in CONTRASTS.items():
        for mean_count in [50, 200, 1000]:
            phi = float(disp[(mean_expr > mean_count / 2) & (mean_expr < mean_count * 2)]
                        .median())
            for alpha, alpha_label in [(0.05, "0.05 (single gene)"),
                                       (0.001, "0.001 (~genome-wide FDR scale)")]:
                rows.append({
                    "contrast": label, "mean_count": mean_count,
                    "typical_dispersion": round(phi, 4), "alpha": alpha_label,
                    "mde_log2fc_all_19": minimum_detectable_lfc(mean_count, phi, n_all, c, alpha),
                    "mde_log2fc_without_MD_1": minimum_detectable_lfc(
                        mean_count, phi, n_drop_md1, c, alpha),
                })
    pd.DataFrame(rows).round(2).to_csv(OUT / "minimum_detectable_effects.csv", index=False)

    (OUT / "feasibility_summary.json").write_text(json.dumps(summary, indent=2) + "\n",
                                                   encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(by_group.reindex(CELLS).round(3).to_string())
    print(disp_table.round(3).to_string())
    print(pd.DataFrame(rows).round(2).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
