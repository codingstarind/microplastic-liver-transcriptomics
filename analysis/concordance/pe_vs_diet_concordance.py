"""Direction A: does the PE response resemble the diet (MASH) response?

Shared-reference problem: PE_in_CD = CD+PE - CD+Veh and Diet_in_Veh = MD+Veh - CD+Veh both
subtract the CD+Veh mean, so their noise is positively correlated and naive concordance is
biased upward. The primary analysis therefore uses a **split reference**: the 4 CD+Veh mice
are divided into two pairs, one pair serving as the PE reference and the other as the diet
reference, so the two contrasts share no mice. All 6 ordered splits are evaluated.

Effect estimates are log2 differences of group means of size-factor-normalised counts
(pseudocount 1, genes with mean normalised count > 50), which track the DESeq2 MLEs closely
for expressed genes. Reference distribution: within each split, the PE contrast's 7 mice
(5 CD+PE + 2 reference) are exhaustively relabelled (21 ways); this removes any PE effect
while keeping gene-gene correlation and the independence from the diet contrast.

The same is done for PE_in_MD vs diet, where MD+Veh is shared: MD+Veh (n = 5) is split 2/3.

Usage: python analysis/concordance/pe_vs_diet_concordance.py <normalized_counts.csv.gz> <label>
       [--drop=SAMPLE ...]   (sensitivity analyses; e.g. --drop=MD_1)
"""

import itertools
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mpliver.config import REPO_ROOT, SAMPLE_METADATA_FILE

TOP_N = 500
MIN_MEAN = 50


def log2_diff(norm, treated, reference):
    return (np.log2(norm[treated] + 1).mean(axis=1) - np.log2(norm[reference] + 1).mean(axis=1))


def hallmark_means(lfc: pd.Series, sets: dict) -> pd.Series:
    return pd.Series({s: lfc.reindex(g).dropna().mean() for s, g in sets.items()
                      if lfc.reindex(g).notna().sum() >= 10})


def concordance(pe: pd.Series, diet: pd.Series, diet_selector: pd.Series, sets: dict) -> dict:
    """Agreement statistics; genes for the top-N check are chosen with diet_selector only."""
    rho = spearmanr(pe, diet).statistic
    top = diet_selector.abs().sort_values(ascending=False).index[:TOP_N]
    sign_agree = float((np.sign(pe[top]) == np.sign(diet[top])).mean())
    hm_pe, hm_diet = hallmark_means(pe, sets), hallmark_means(diet, sets)
    common = hm_pe.index.intersection(hm_diet.index)
    return {"spearman_all_genes": float(rho),
            f"sign_agreement_top{TOP_N}_diet_genes": sign_agree,
            "hallmark_mean_lfc_spearman":
                float(spearmanr(hm_pe[common], hm_diet[common]).statistic)}


def run(norm, groups, sets, pe_treated, shared_name, shared, diet_other, pe_is_reference_side):
    """Evaluate all splits of the shared group; return per-split stats and null quantiles."""
    rows = []
    k = len(shared) // 2
    sizes = [k] if len(shared) % 2 == 0 else [k, k + 1]
    for size in sizes:
        for part in itertools.combinations(shared, size):
            a = list(part)                                  # used in the PE contrast
            b = [m for m in shared if m not in a]           # used in the diet contrast
            if pe_is_reference_side:                        # shared group = PE reference
                pe = log2_diff(norm, pe_treated, a)
                diet = log2_diff(norm, b, diet_other) if shared_name == "MD+Veh" else \
                    log2_diff(norm, diet_other, b)
            else:
                raise ValueError("unsupported layout")
            obs = concordance(pe, diet, diet, sets)
            null = []
            pool = pe_treated + a
            for ref in itertools.combinations(pool, len(a)):
                trt = [m for m in pool if m not in ref]
                null.append(concordance(log2_diff(norm, trt, list(ref)), diet, diet, sets))
            null = pd.DataFrame(null)
            row = {"shared_group": shared_name, "pe_reference": "|".join(a),
                   "diet_side": "|".join(b), **obs}
            for col in null.columns:
                row[f"{col}_null_median"] = float(null[col].median())
                row[f"{col}_null_q95"] = float(null[col].quantile(0.95))
                row[f"{col}_perm_p_one_sided"] = float((null[col] >= obs[col]).mean())
            rows.append(row)
    return pd.DataFrame(rows)


def main() -> int:
    norm_path, label = sys.argv[1], sys.argv[2]
    drop = [a.split("=", 1)[1] for a in sys.argv[3:] if a.startswith("--drop=")]
    out = REPO_ROOT / "results" / "concordance" / label
    out.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(SAMPLE_METADATA_FILE)
    meta = meta[~meta["counts_column"].str.rsplit("_", n=1).str[0].isin(drop)]
    norm = pd.read_csv(norm_path, index_col=0)[meta["counts_column"]]
    norm = norm[norm.mean(axis=1) > MIN_MEAN]
    groups = meta.set_index("counts_column")["group"]
    members = {g: groups.index[groups == g].tolist() for g in groups.unique()}
    hallmark = pd.read_csv(REPO_ROOT / "results" / "concordance" / "hallmark_sets_2026.1.Mm.csv")
    sets = hallmark.groupby("gs_name")["gene_symbol"].apply(list).to_dict()

    # Naive (shared-reference) values, for comparison with the split estimates.
    naive_cd = concordance(log2_diff(norm, members["CD+PE"], members["CD+Veh"]),
                           log2_diff(norm, members["MD+Veh"], members["CD+Veh"]),
                           log2_diff(norm, members["MD+Veh"], members["CD+Veh"]), sets)
    naive_md = concordance(log2_diff(norm, members["MD+PE"], members["MD+Veh"]),
                           log2_diff(norm, members["MD+Veh"], members["CD+Veh"]),
                           log2_diff(norm, members["MD+Veh"], members["CD+Veh"]), sets)

    # PE_in_CD vs diet: CD+Veh shared (PE reference vs diet reference).
    cd = run(norm, groups, sets, members["CD+PE"], "CD+Veh", members["CD+Veh"],
             members["MD+Veh"], pe_is_reference_side=True)
    # PE_in_MD vs diet: MD+Veh shared (PE reference vs diet treated group).
    md = run(norm, groups, sets, members["MD+PE"], "MD+Veh", members["MD+Veh"],
             members["CD+Veh"], pe_is_reference_side=True)
    cd.insert(0, "comparison", "PE_in_CD vs Diet_in_Veh")
    md.insert(0, "comparison", "PE_in_MD vs Diet_in_Veh")
    splits = pd.concat([cd, md], ignore_index=True)
    splits.round(4).to_csv(out / "split_reference_concordance.csv", index=False)

    stats = ["spearman_all_genes", f"sign_agreement_top{TOP_N}_diet_genes",
             "hallmark_mean_lfc_spearman"]
    summary = {
        "label": label, "dropped": drop, "n_genes": int(len(norm)),
        "min_mean_normalised_count": MIN_MEAN,
        "naive_shared_reference": {"PE_in_CD": naive_cd, "PE_in_MD": naive_md},
        "split_reference": {
            comp: {s: {"min": float(g[s].min()), "median": float(g[s].median()),
                       "max": float(g[s].max()),
                       "null_median_of_medians": float(g[f"{s}_null_median"].median()),
                       "perm_p_range": [float(g[f"{s}_perm_p_one_sided"].min()),
                                        float(g[f"{s}_perm_p_one_sided"].max())]}
                   for s in stats}
            for comp, g in splits.groupby("comparison")},
    }
    (out / "concordance_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
