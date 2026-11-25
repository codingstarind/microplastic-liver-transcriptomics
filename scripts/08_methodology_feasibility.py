"""Step 8: feasibility checks for the Phase 2 methodology proposal.

No hypothesis tests are run and no p-values are computed for genes. This script
describes the design and the noise level, so we can judge which comparisons the
data can support:

    1. sequencing depth vs summed gene counts (where the library-size gap arises)
    2. confounding / collinearity between design factors and technical variables
    3. rough gene-wise dispersion, with leave-one-mouse-out uncertainty
    4. rough minimum detectable log2 fold change for each candidate comparison

All dispersion and power numbers are planning approximations; see the
docstrings in src/mpliver/feasibility.py and the proposal for their limits.
Outputs go to results/methodology/ (Phase 1 results are not touched).
"""

import json
import sys

import numpy as np
import pandas as pd

from mpliver.config import (
    METADATA_DIR,
    PROCESSED_DIR,
    REPO_ROOT,
    SAMPLE_METADATA_FILE,
    TABLES_DIR,
)
from mpliver.feasibility import (
    CELLS,
    CONTRASTS,
    design_matrix,
    jackknife_interval,
    leave_one_mouse_out,
    minimum_detectable_lfc,
    moment_dispersion,
    r2_on_design,
)
from mpliver.qc_markers import marker_set_cpm

OUT = REPO_ROOT / "results" / "methodology"
EXPRESSION_BINS = [10, 50, 100, 500, 1000, 5000, 1e9]
BIN_LABELS = ["10-50", "50-100", "100-500", "500-1k", "1k-5k", ">5k"]
# Not a real proportion: Salmon counts fragments, but reads removed by fastp,
# multi-mapping redistribution, rounding and non-gene targets all enter the ratio.
PROXY = "gene_counts_per_read_pair_proxy"


def sequencing_depth_table(meta, counts):
    runs = pd.read_csv(METADATA_DIR / "sra_run_info.csv")
    df = meta.merge(runs, on="geo_accession")
    df["summed_gene_counts"] = counts.sum().reindex(df["counts_column"]).to_numpy()
    df[PROXY] = df["summed_gene_counts"] / df["read_pairs"]
    return df[["counts_column", "group", "library_name", "read_pairs", "summed_gene_counts",
               PROXY, "sequencer_id", "run_number", "flowcell", "lane"]]


def dispersion_by_bin(disp, mean_expr):
    bins = pd.cut(mean_expr, EXPRESSION_BINS, labels=BIN_LABELS)
    return disp.groupby(bins, observed=True).quantile([0.25, 0.5, 0.75]).unstack()


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(SAMPLE_METADATA_FILE)
    order = meta["counts_column"]
    counts = pd.read_csv(PROCESSED_DIR / "counts_original.csv.gz", index_col=0)[order]
    norm_counts = pd.read_csv(PROCESSED_DIR / "normalized_counts_filtered.csv.gz",
                              index_col=0)[order]
    size_factors = pd.read_csv(TABLES_DIR / "size_factors.csv", index_col=0)["size_factor"]
    groups = meta.set_index("counts_column")["group"].reindex(order)
    summary = {}

    # 1. Sequencing depth -------------------------------------------------------
    depth = sequencing_depth_table(meta, counts)
    depth.to_csv(OUT / "sequencing_depth_by_sample.csv", index=False)
    by_group = (depth.groupby("group")[["read_pairs", "summed_gene_counts", PROXY]]
                .agg(["mean", "min", "max"]).reindex(CELLS))
    by_group.round(3).to_csv(OUT / "sequencing_depth_by_group.csv")
    verified = depth.dropna(subset=["flowcell"])
    summary["flowcell_evidence"] = {
        "n_runs_total": int(len(depth)),
        "n_runs_with_header_read": int(len(verified)),
        "flowcells_seen": sorted(verified["flowcell"].unique().tolist()),
        "run_numbers_seen": sorted(verified["run_number"].astype(str).unique().tolist()),
        "groups_with_header_read": verified["group"].value_counts().reindex(CELLS)
                                                    .fillna(0).astype(int).to_dict(),
        "runs_without_header": depth.loc[depth["flowcell"].isna(), "counts_column"].tolist(),
    }

    # 2. Confounding and collinearity -------------------------------------------
    X = design_matrix(meta["diet"].eq("MD"), meta["microplastic"].eq("PE"), order)
    summary["design_rank"] = int(np.linalg.matrix_rank(X.to_numpy()))
    muscle = np.log2(marker_set_cpm(counts)["skeletal_muscle"] + 1).reindex(order)
    technical = pd.DataFrame({
        "log2_read_pairs": np.log2(depth.set_index("counts_column")["read_pairs"]),
        PROXY: depth.set_index("counts_column")[PROXY],
        "log2_muscle_cpm": muscle,
        "library_index": meta.set_index("counts_column")["library_name"]
                             .str.lstrip("S").astype(int),
    }).reindex(order)
    technical["group"] = groups
    technical.round(4).to_csv(OUT / "technical_covariates_by_sample.csv")
    summary["technical_vs_design"] = {
        col: r2_on_design(X, technical[col])
        for col in ["log2_read_pairs", PROXY, "log2_muscle_cpm", "library_index"]
    }
    pe = meta["microplastic"].eq("PE").to_numpy()
    reads = depth.set_index("counts_column").loc[order, "read_pairs"].to_numpy()
    summary["read_pairs_PE_vs_Veh"] = {
        "max_PE": int(reads[pe].max()), "min_Veh": int(reads[~pe].min()),
        "perfect_separation": bool(reads[pe].max() < reads[~pe].min()),
    }
    summary["muscle_cpm_gt_100_by_group"] = (
        pd.Series((muscle > np.log2(101)).to_numpy(), index=groups.to_numpy())
        .groupby(level=0).sum().reindex(CELLS).astype(int).to_dict())

    # 3. Dispersion, with leave-one-mouse-out uncertainty ---------------------------
    disp = moment_dispersion(norm_counts, groups, size_factors)
    mean_expr = norm_counts.mean(axis=1)
    point = dispersion_by_bin(disp, mean_expr)
    loo = {}
    for mouse, rest in leave_one_mouse_out(groups):
        cols = rest.index
        d = moment_dispersion(norm_counts[cols], rest, size_factors[cols])
        loo[mouse] = dispersion_by_bin(d, mean_expr)[0.5]
    loo = pd.DataFrame(loo)
    loo.round(5).to_csv(OUT / "dispersion_median_leave_one_mouse_out.csv")
    jk = pd.DataFrame({b: jackknife_interval(loo.loc[b].to_numpy(), point.loc[b, 0.5])
                       for b in point.index}, index=["low", "high"]).T
    disp_table = pd.DataFrame({
        "n_genes": pd.cut(mean_expr, EXPRESSION_BINS, labels=BIN_LABELS)
                     .value_counts().reindex(BIN_LABELS),
        "dispersion_q25": point[0.25], "dispersion_median": point[0.5],
        "dispersion_q75": point[0.75],
        "median_jackknife_low": jk["low"].clip(lower=1e-4),
        "median_jackknife_high": jk["high"],
        "median_loo_min": loo.min(axis=1), "median_loo_max": loo.max(axis=1),
        "mouse_moving_median_most": (loo.sub(point[0.5], axis=0)).abs().idxmax(axis=1),
    })
    disp_table["median_biological_CV"] = np.sqrt(disp_table["dispersion_median"])
    disp_table.index.name = "mean_normalised_count"
    disp_table.round(4).to_csv(OUT / "dispersion_by_expression.csv")
    summary["overall_median_dispersion"] = float(disp.median())

    # 4. Minimum detectable effects, with dispersion uncertainty -------------------
    bin_of = {50: "10-50", 200: "100-500", 1000: "500-1k"}
    rows = []
    for label, contrast in CONTRASTS.items():
        for mean_count, b in bin_of.items():
            for alpha, alpha_label in [(0.05, "0.05 (single gene)"),
                                       (0.001, "0.001 (illustrative genome-wide)")]:
                def mde(phi, n=(4, 5, 5, 5), mean_count=mean_count, contrast=contrast,
                        alpha=alpha):
                    return minimum_detectable_lfc(mean_count, phi, list(n), contrast, alpha)
                rows.append({
                    "contrast": label, "mean_count": mean_count, "alpha": alpha_label,
                    "mde_median_dispersion": mde(disp_table.loc[b, "dispersion_median"]),
                    "mde_q25_dispersion": mde(disp_table.loc[b, "dispersion_q25"]),
                    "mde_q75_dispersion": mde(disp_table.loc[b, "dispersion_q75"]),
                    "mde_jackknife_low": mde(disp_table.loc[b, "median_jackknife_low"]),
                    "mde_jackknife_high": mde(disp_table.loc[b, "median_jackknife_high"]),
                    "mde_median_dispersion_without_MD_1":
                        mde(disp_table.loc[b, "dispersion_median"], n=(4, 5, 4, 5)),
                })
    mde_table = pd.DataFrame(rows)
    mde_table.round(2).to_csv(OUT / "minimum_detectable_effects.csv", index=False)

    (OUT / "feasibility_summary.json").write_text(json.dumps(summary, indent=2) + "\n",
                                                   encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(by_group.round(3).to_string())
    print(disp_table.round(4).to_string())
    print(mde_table[mde_table["mean_count"] == 200].round(2).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
