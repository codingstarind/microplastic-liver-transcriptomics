"""Training-label permutation null for external AUC (R2, logistic C=1), 200 per primary split.
Deviation D1 (technical): with permuted labels the |diff|>=1 rule can select <10 genes (empty set -> NaN score);
fallback = top 50 by t. Real training sets always had 50+50 genes passing (verified), so real fits are unaffected."""
import json
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import e01_mash_state_models as E
from scipy import stats

_orig = E.welch_signature


def welch_fallback(X, y, n=50, min_diff=1.0):
    up, dn = _orig(X, y, n, min_diff)
    if len(up) < 10 or len(dn) < 10:
        a, b = X[y == 1], X[y == 0]
        t = np.nan_to_num(stats.ttest_ind(a, b, equal_var=False).statistic)
        o = np.argsort(t)
        up, dn = (o[::-1][:n] if len(up) < 10 else up), (o[:n] if len(dn) < 10 else dn)
    return up, dn


def job(i, sp):
    E.welch_signature = welch_fallback
    if "mash" not in E.CTX:
        E.build_ctx("mash")
    return E.perm_one(i, sp)


if __name__ == "__main__":
    splits = ["S1_train316TA_test246", "S2_train246_test316TA"]
    pj = [(i, sp) for sp in splits for i in range(200)]
    pr = Parallel(n_jobs=E.N_JOBS, verbose=2, batch_size=10)(delayed(job)(i, sp) for i, sp in pj)
    perms = {sp: [r for (i, s_), r in zip(pj, pr) if s_ == sp] for sp in splits}
    pd.DataFrame(perms).to_csv(E.RES / "permutation_null_auc_R2_logistic.csv", index=False)
    m = pd.read_csv(E.RES / "external_metrics.csv")
    obs = m[(m.rep == "R2") & (m.model == "logistic_default")].set_index("split").auc
    out = {sp: dict(observed=float(obs[sp]), null_mean=float(np.nanmean(x)), null_q95=float(np.nanquantile(x, .95)),
                    p=float((1 + np.sum(np.array(x) >= obs[sp])) / (1 + len(x)))) for sp, x in perms.items()}
    json.dump(out, open(E.RES / "permutation_summary.json", "w"), indent=2)
    print(out)
