"""POST-HOC VERIFICATION (review addendum, 2026-11-22). Does not change the locked decision.

Exhaustively enumerates all choose(16,8) = 12,870 allocations of the 16 GSE316023 primary-contrast mice
into 8 'case' vs 8 'control', using EXACTLY the estimator of the original Monte Carlo null
(p3lib.relabel_null / c01_positive_control.py):
  statistic = Pearson correlation of average ranks of
              (mean log2 normalised count in 'case' - mean in 'control')  and  GSE294180 Diet_in_Veh log2FC,
  genes = main universe present in the MASH1 log2norm matrix with a non-NA target log2FC.
Outputs: results/posthoc_audit/exact_permutation.json, exact_permutation_by_overlap.csv
The original Monte Carlo result (results/positive_control/positive_control.json) is left untouched.
"""
import itertools
import json

import numpy as np
import pandas as pd
from p3lib import P3, read_gz_csv, target, universe
from scipy.stats import rankdata

OUT = P3 / "results" / "posthoc_audit"
OUT.mkdir(parents=True, exist_ok=True)
LABEL = "MASH1_GSE316023_GANvsChow_TA"


def main():
    U = universe("main")
    tgt = target("Diet_in_Veh")
    ln = read_gz_csv(P3 / "results" / "mash" / f"{LABEL}_log2norm.csv.gz", index_col=0)
    ln = ln.loc[ln.index.intersection(list(U))]
    g = ln.index.intersection(tgt["log2FC"].dropna().index)
    X = ln.loc[g].values
    cols = list(ln.columns)
    true_case = np.array([c.startswith("SD_") for c in cols])  # SD_ = GAN diet, NC_ = normal chow
    assert true_case.sum() == 8 and len(cols) == 16
    t = pd.Series(tgt.loc[g, "log2FC"].values).rank().values   # average ranks, as in relabel_null
    tc = (t - t.mean()) / np.sqrt(((t - t.mean()) ** 2).sum())

    def stat_for(masks):  # masks: (16, k) boolean
        W = masks / 8.0 - (~masks) / 8.0
        D = X @ W
        R = rankdata(D, axis=0)
        R = R - R.mean(0)
        R = R / np.sqrt((R ** 2).sum(0))
        return tc @ R

    obs = float(stat_for(true_case[:, None])[0])
    combos = np.array(list(itertools.combinations(range(16), 8)))
    stats = np.empty(len(combos))
    overlap = np.empty(len(combos), int)
    for s in range(0, len(combos), 1000):
        cb = combos[s:s + 1000]
        M = np.zeros((16, len(cb)), bool)
        M[cb.T, np.arange(len(cb))] = True
        stats[s:s + 1000] = stat_for(M)
        overlap[s:s + 1000] = (M & true_case[:, None]).sum(0)
    n_ge = int((stats >= obs - 1e-12).sum())
    is_true = overlap == 8
    res = dict(
        n_allocations=int(len(stats)), n_genes=int(len(g)),
        case_columns=[c for c, k in zip(cols, true_case) if k], control_columns=[c for c, k in zip(cols, true_case) if not k],
        observed_stat=obs,
        observed_matches_original_montecarlo=None,
        observed_allocation_stat_from_enumeration=float(stats[is_true][0]),
        exact_one_sided_p=n_ge / len(stats), n_allocations_ge_observed=n_ge,
        exact_q99=float(np.quantile(stats, .99)), exact_q95=float(np.quantile(stats, .95)),
        exact_max=float(stats.max()), exact_median=float(np.median(stats)),
        balanced_4v4_n=int((overlap == 4).sum()), balanced_4v4_max=float(stats[overlap == 4].max()),
        note="Post-hoc verification; locked decision (FAIL on c1) unchanged.",
    )
    orig = json.loads((P3 / "results" / "positive_control" / "positive_control.json").read_text())["M5_relabel_null"]
    res["original_montecarlo"] = orig
    res["observed_matches_original_montecarlo"] = bool(abs(orig["observed_rho_meandiff"] - obs) < 1e-9)
    by = pd.DataFrame({"overlap_true_GAN_in_case": overlap, "stat": stats}).groupby("overlap_true_GAN_in_case").stat.agg(
        n="size", mean="mean", min="min", max="max", n_ge_observed=lambda s: int((s >= obs - 1e-12).sum()))
    by.to_csv(OUT / "exact_permutation_by_overlap.csv")
    np.save(OUT / "exact_permutation_stats.npy", stats)
    (OUT / "exact_permutation.json").write_text(json.dumps(res, indent=2))
    print(json.dumps({k: v for k, v in res.items() if k not in ("case_columns", "control_columns")}, indent=1))
    print(by.round(4).to_string())


if __name__ == "__main__":
    main()
