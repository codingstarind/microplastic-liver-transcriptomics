"""Stage A: audit the GSE316023 FPKM-derived effective-length offset.

Checks (analysis_plan.md section 2):
 (i)   count and FPKM matrices: identical samples, genes, zero pattern
 (ii)  implied length L = count*1e9/(FPKM*N) in bp range; correlation with GENCODE vM35
       median transcript length
 (iii) identifiability: between-sample variation beyond a per-sample constant
 (iv)  instability from 4-decimal FPKM rounding (relative length error bound 5e-5/FPKM)
 (v)   DESeq2-style normalisation factors (L / rowGeoMean(L)) finite, positive, centred
 (vi)  rounding of estimated counts
Outputs: results/offset_audit/audit.json, offset_gene_table.csv.gz, figures/offset_audit.png
"""
import gzip
import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

P3 = Path(__file__).resolve().parents[1]
REPO = P3.parents[1]
P2 = REPO / "research" / "roadmap_phase2"
D = Path("D:/mpliver_work/phase2_refs/GSE316023")
OUT = P3 / "results" / "offset_audit"
OUT.mkdir(parents=True, exist_ok=True)


def gencode_tx_length() -> pd.Series:
    tx_len, tx_gene = {}, {}
    with gzip.open(REPO / "data" / "raw" / "gencode.vM35.annotation.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "exon":
                continue
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            g = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
            tx_len[t] = tx_len.get(t, 0) + int(f[4]) - int(f[3]) + 1
            tx_gene[t] = g
    s = pd.Series(tx_len)
    return s.groupby(pd.Series(tx_gene)).median()


def main():
    c = pd.read_csv(D / "GSE316023_COUNT_NCOMMS-25-40823-T.csv.gz", index_col=0)
    f = pd.read_csv(D / "GSE316023_FPKM_NCOMMS-25-40823-T.csv.gz", index_col=0)
    meta = pd.read_csv(P2 / "datasets" / "GSE316023" / "sample_metadata.csv")
    case = meta.loc[meta.contrast_role == "MASH_primary:case", "sample_id"].tolist()
    ctrl = meta.loc[meta.contrast_role == "MASH_primary:control", "sample_id"].tolist()
    s16 = case + ctrl
    uni = pd.read_csv(P2 / "metadata" / "phase3_gene_universe.csv", index_col=0)
    U = uni.index[uni.universe_main]
    mp = pd.read_csv(P2 / "datasets" / "GSE316023" / "gene_mapping.csv.gz")
    kept = mp[mp.kept_in_ready].set_index("original_id")["vm35_symbol"]
    res = {}

    # (i)
    res["i_same_samples_order"] = bool(c.columns.equals(f.columns))
    res["i_same_gene_ids_order"] = bool(c.index.equals(f.index))
    zc, zf = (c.values == 0), (f.values == 0)
    res["i_zero_pattern_mismatch_cells"] = int((zc != zf).sum())
    res["i_zero_mismatch_count_pos_fpkm0"] = int((~zc & zf).sum())
    res["i_zero_mismatch_count0_fpkm_pos"] = int((zc & ~zf).sum())
    # where count>0 but FPKM rounds to 0 (precision floor)
    N = c.sum()
    L = (c * 1e9 / (f * N)).where((c > 0) & (f > 0))
    # restrict to unique original ids kept in ready matrix, universe genes, 16 samples
    uniq = ~c.index.duplicated(keep=False)
    uniq_ids = set(c.index[uniq])
    keep_ids = [i for i in kept.index if i in uniq_ids]
    sym = kept.loc[keep_ids].values

    def to_u(m):
        x = m.loc[keep_ids].copy()
        x.index = sym
        return x.reindex(U)[s16]

    Lu, cu, fu = to_u(L), to_u(c), to_u(f)
    res["universe_genes"] = int(len(U))
    res["universe_genes_found_in_GSE316023"] = int(cu.notna().all(axis=1).sum())

    # (ii)
    medL = Lu.median(axis=1)
    res["ii_implied_length_quantiles_bp"] = {q: float(medL.quantile(q)) for q in (0.01, 0.25, 0.5, 0.75, 0.99)}
    gl = gencode_tx_length()
    both = medL.dropna().index.intersection(gl.index)
    rho, _ = spearmanr(medL[both], gl[both])
    res["ii_spearman_implied_vs_gencode_median_tx_length"] = float(rho)
    res["ii_n_genes_compared"] = int(len(both))
    ratio = (medL[both] / gl[both])
    res["ii_median_ratio_implied_to_gencode"] = float(ratio.median())

    # (iii) per-sample constant: log L_ij - row mean = col effect + residual
    lg = np.log(Lu.dropna())
    lg = lg[(cu.loc[lg.index] >= 10).all(axis=1)]
    rc = lg.sub(lg.mean(axis=1), axis=0)
    col = rc.median(axis=0)
    resid = rc.sub(col, axis=1)
    res["iii_genes_all_counts_ge10"] = int(len(lg))
    res["iii_sample_scale_factors_exp_range"] = [float(np.exp(col.min())), float(np.exp(col.max()))]
    res["iii_var_share_column_constant"] = float((rc.sub(resid)).values.var() / rc.values.var())
    res["iii_residual_sd_log"] = float(resid.values.std())
    gm = lg[case].mean(axis=1) - lg[ctrl].mean(axis=1) - (col[case].mean() - col[ctrl].mean())
    res["iii_genes_abs_log2_group_length_ratio_gt_0.263"] = int((gm.abs() / np.log(2) > 0.263).sum())
    res["iii_genes_abs_log2_group_length_ratio_gt_0.1"] = int((gm.abs() / np.log(2) > 0.1).sum())

    # (iv) rounding instability
    relerr = 5e-5 / fu.where(fu > 0)
    unstable = (relerr > 0.01) & (cu > 0)
    res["iv_cells_rel_len_err_gt_1pct"] = int(unstable.sum().sum())
    res["iv_genes_any_cell_err_gt_1pct"] = int(unstable.any(axis=1).sum())
    res["iv_genes_any_cell_err_gt_1pct_and_detected"] = int((unstable.any(axis=1) & ((cu >= 10).sum(axis=1) >= 4)).sum())
    res["iv_cells_count_pos_fpkm_zero"] = int(((cu > 0) & (fu == 0)).sum().sum())
    res["iv_genes_no_length_any_sample"] = int(Lu.isna().any(axis=1).sum())

    # (v) normalisation factors as DESeq2 builds them from avgTxLength
    Lf = Lu.apply(lambda r: r.fillna(r.mean()), axis=1).dropna()
    nm = Lf.div(np.exp(np.log(Lf).mean(axis=1)), axis=0)
    res["v_normMatrix_finite"] = bool(np.isfinite(nm.values).all())
    res["v_normMatrix_positive"] = bool((nm.values > 0).all())
    res["v_normMatrix_row_geomean_max_abs_log"] = float(np.abs(np.log(nm).mean(axis=1)).max())
    res["v_normMatrix_range"] = [float(nm.values.min()), float(nm.values.max())]
    res["v_normMatrix_q01_q99"] = [float(np.quantile(nm.values, 0.01)), float(np.quantile(nm.values, 0.99))]

    # (vi) rounding
    s = c[s16]
    res["vi_noninteger_fraction"] = float(((s % 1) != 0).values.mean())
    res["vi_max_rel_library_change"] = float(((np.rint(s) - s).abs().sum() / s.sum()).max())

    # documentation
    res["documentation"] = ("GEO data_processing: 'FPKM was generated for following process'; mm10; quantification "
                            "software, length definition and FPKM denominator NOT STATED. Paper (PMC13392390) methods: "
                            "no RNA-seq quantification details.")
    # decision per plan section 2
    checks = dict(
        i=res["i_same_samples_order"] and res["i_same_gene_ids_order"],
        ii=0.5 <= res["ii_spearman_implied_vs_gencode_median_tx_length"],
        iii=res["iii_var_share_column_constant"] < 0.9,
        iv=res["iv_genes_any_cell_err_gt_1pct_and_detected"] < 0.01 * res["universe_genes"],
        v=res["v_normMatrix_finite"] and res["v_normMatrix_positive"] and res["v_normMatrix_row_geomean_max_abs_log"] < 1e-8,
        vi=res["vi_max_rel_library_change"] < 1e-4,
    )
    res["checks_pass"] = checks
    res["provenance_documented"] = False
    res["decision"] = ("no_offset_primary__offset_approximate_sensitivity" if all(checks.values())
                       else "no_offset_primary__offset_not_used_as_primary_checks_failed")
    (OUT / "audit.json").write_text(json.dumps(res, indent=2))
    tab = pd.DataFrame({"median_implied_length": medL, "gencode_median_tx_length": gl.reindex(medL.index),
                        "log2_GAN_vs_chow_length_ratio": (gm / np.log(2)).reindex(medL.index),
                        "any_cell_len_err_gt_1pct": unstable.any(axis=1)})
    tab.to_csv(OUT / "offset_gene_table.csv.gz")
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].loglog(gl[both], medL[both], ".", ms=2, alpha=.3)
    ax[0].set_xlabel("GENCODE vM35 median transcript length (bp)")
    ax[0].set_ylabel("implied effective length (bp), median of 16")
    ax[0].set_title(f"(ii) Spearman {rho:.2f}", fontsize=9)
    mean_ct = np.log10(cu.loc[gm.index].mean(axis=1))
    ax[1].scatter(mean_ct, gm / np.log(2), s=2, alpha=.3)
    ax[1].axhline(0.263, ls="--", c="k", lw=.7)
    ax[1].axhline(-0.263, ls="--", c="k", lw=.7)
    ax[1].set_xlabel("log10 mean count (16 samples)")
    ax[1].set_ylabel("log2 GAN/chow implied length (sample-scale removed)")
    ax[1].set_title("(iii) group difference in implied length", fontsize=9)
    fig.tight_layout()
    fig.savefig(P3 / "figures" / "offset_audit.png", dpi=110)
    # export the length matrix used for the sensitivity model, universe x 16 samples
    Lf.to_csv(OUT / "gse316023_length_matrix_universe16.csv.gz")
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
