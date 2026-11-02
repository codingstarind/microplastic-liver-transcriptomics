"""Consistency checks on the committed methodology-feasibility outputs (runs in CI).

Needs no downloaded data: it checks that the committed tables are complete,
internally consistent and use the agreed (qualified) column names.
"""

import json
import sys

import pandas as pd

from mpliver.config import METADATA_DIR, REPO_ROOT, SAMPLE_METADATA_FILE

OUT = REPO_ROOT / "results" / "methodology"
EXPECTED = {
    "sequencing_depth_by_sample.csv": ["counts_column", "group", "read_pairs",
                                       "summed_gene_counts", "gene_counts_per_read_pair_proxy",
                                       "flowcell"],
    "technical_covariates_by_sample.csv": ["log2_read_pairs", "gene_counts_per_read_pair_proxy",
                                           "log2_muscle_cpm", "library_index", "group"],
    "dispersion_by_expression.csv": ["mean_normalised_count", "dispersion_median",
                                     "median_jackknife_low", "median_jackknife_high"],
    "minimum_detectable_effects.csv": ["contrast", "mean_count", "alpha",
                                       "mde_median_dispersion", "mde_jackknife_low",
                                       "mde_jackknife_high"],
}


def main() -> int:
    errors = []
    meta = pd.read_csv(SAMPLE_METADATA_FILE)

    runs = pd.read_csv(METADATA_DIR / "sra_run_info.csv")
    if len(runs) != 19 or runs["run_accession"].duplicated().any():
        errors.append("sra_run_info.csv should list 19 unique runs")
    if set(runs["geo_accession"]) != set(meta["geo_accession"]):
        errors.append("sra_run_info.csv GSM accessions do not match sample_metadata.csv")

    tables = {}
    for name, columns in EXPECTED.items():
        path = OUT / name
        if not path.exists():
            errors.append(f"missing {name}")
            continue
        tables[name] = df = pd.read_csv(path)
        missing = [c for c in columns if c not in df.columns]
        if missing:
            errors.append(f"{name} lacks columns {missing}")
        if any("assigned_fraction" in c for c in df.columns):
            errors.append(f"{name} uses the unqualified name 'assigned_fraction'")

    depth = tables.get("sequencing_depth_by_sample.csv")
    if depth is not None and (depth["summed_gene_counts"] > depth["read_pairs"]).any():
        errors.append("summed gene counts exceed read pairs for some sample")

    mde = tables.get("minimum_detectable_effects.csv")
    if mde is not None:
        if mde[["mde_median_dispersion", "mde_jackknife_low", "mde_jackknife_high"]] \
                .isna().any().any():
            errors.append("minimum_detectable_effects.csv has missing values")
        if not (mde["mde_jackknife_low"] <= mde["mde_jackknife_high"]).all():
            errors.append("jackknife interval bounds are reversed")
        key = ["mean_count", "alpha"]
        pivot = mde.pivot_table(index=key, columns="contrast", values="mde_median_dispersion")
        main_eff = pivot["PE main effect (average over diets)"]
        within = pivot["PE within MD (MD+PE vs MD+Veh)"]
        inter = pivot["Diet x PE interaction (difference of PE effects)"]
        if not ((main_eff < within) & (within < inter)).all():
            errors.append("expected main effect < simple effect < interaction MDE ordering")

    summary = json.loads((OUT / "feasibility_summary.json").read_text(encoding="utf-8"))
    if summary["flowcell_evidence"]["n_runs_total"] != 19:
        errors.append("feasibility_summary.json should cover 19 runs")

    for error in errors:
        print(f"[error] {error}")
    if not errors:
        print("Methodology output checks passed.")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
