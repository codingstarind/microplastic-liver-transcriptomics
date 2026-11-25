"""Step 4 (Q1): frozen gene-level predictions vs independent external plastic log2FC vectors (locked_plan.md section 4).
No model is refitted. Per (contrast x frozen model): weighted R2 vs zero / training-mean / own-mean baselines, unweighted R2,
MAE, RMSE, Spearman + module-block bootstrap CI, sign agreement, calibration slope (descriptive), noise ceiling, and a
mouse-label permutation null (all relabellings) on a log2CPM mean-difference proxy."""
import itertools

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from p5_common import FRZ, N_JOBS, P5, RES, RS, SEED, load_prepared
from scipy import stats

DE = RES / "external_de"
P = pd.read_csv(FRZ / "gene_level_frozen_predictions.csv.gz", index_col=0)
CO = pd.read_csv(FRZ / "gene_level_frozen_coefficients.csv")
TRAIN_MEAN = CO[CO.feature == "(training weighted mean y)"].set_index("model").coef_per_SD
mods = pd.read_csv(RS / "roadmap_phase4/results/gene_groups.csv", index_col=0)["module"]

CONTRASTS = [  # (study, file tag, exposed groups, role)
    ("HONG", "PE_vs_CON", ["PE"], "PRIMARY"),
    ("HONG", "PE_vs_CON__6mice_only", ["PE"], "sensitivity (6-mouse model)"),
    ("HONG", "PS_vs_CON", ["PS"], "exploratory specificity (same lab, shared controls)"),
    ("HONG", "PP_vs_CON", ["PP"], "exploratory specificity (same lab, shared controls)"),
    ("HONG", "PB_vs_CON", ["PB"], "exploratory specificity (same lab, shared controls)"),
    ("GSE245069", "PSpooled_vs_CON", ["PS_80nm", "PS_0.5um", "PS_5um"], "SECONDARY cross-polymer"),
    ("GSE245069", "PSpooled_vs_CON__limmavoom_unrounded", ["PS_80nm", "PS_0.5um", "PS_5um"], "sensitivity (limma-voom, unrounded)"),
    ("GSE245069", "PS80nm_vs_CON", ["PS_80nm"], "exploratory size (shared control)"),
    ("GSE245069", "PS05um_vs_CON", ["PS_0.5um"], "exploratory size (shared control)"),
    ("GSE245069", "PS5um_vs_CON", ["PS_5um"], "exploratory size (shared control)"),
    ("GSE277772", "PS_vs_CON", ["PS"], "exploratory (normalised values; not validation)"),
]
MODELS = ["F8_ols", "F8C_ridge", "FAST_only_ols", "MASH_only_ols", "observed_GSE294180_ashr", "zero"]


def dl_tau2(y, s):
    w = 1 / s**2
    yb = (w * y).sum() / w.sum()
    q = (w * (y - yb) ** 2).sum()
    c = w.sum() - (w**2).sum() / w.sum()
    return max(0.0, (q - (len(y) - 1)) / c)


def load_target(study, tag):
    t = pd.read_csv(DE / f"{study}__{tag}.csv.gz").set_index("gene")
    t = t[np.isfinite(t.log2FC) & np.isfinite(t.lfcSE) & (t.lfcSE > 0)]
    return t


def block_boot_rho(y, p, grp, n=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    u = np.unique(grp)
    idx = {g: np.where(grp == g)[0] for g in u}
    out = []
    for _ in range(n):
        ii = np.concatenate([idx[g] for g in rng.choice(u, len(u))])
        out.append(stats.spearmanr(y[ii], p[ii]).statistic)
    return np.quantile(out, [.025, .975])


def proxy_null(study, exposed, genes, pred_rank):
    """All relabellings of (CON + exposed) mice; Spearman(pred, mean log2CPM difference). Returns (obs, null array)."""
    X, S = load_prepared(study)
    keep_m = S.group.isin(["CON"] + exposed)
    X, S = X.loc[:, keep_m.values], S[keep_m]
    L = np.log2(X / X.sum() * 1e6 + 1).reindex(genes).values  # genes x mice (genes already prefiltered by DE step)
    y = (S.group != "CON").values
    n, k = len(y), int((~y).sum())
    obs_diff = L[:, y].mean(1) - L[:, ~y].mean(1)
    rho = lambda d: stats.spearmanr(pred_rank, d).statistic  # noqa: E731
    null = []
    for ctrl in itertools.combinations(range(n), k):
        m = np.ones(n, bool)
        m[list(ctrl)] = False
        null.append(rho(L[:, m].mean(1) - L[:, ~m].mean(1)))
    return rho(obs_diff), np.array(null)


def evaluate(study, tag, exposed, role, target_name):
    t = load_target(study, tag)
    genes = P.index.intersection(t.index)
    y, s = t.loc[genes, "log2FC"].values, t.loc[genes, "lfcSE"].values
    tau2 = dl_tau2(y, s)
    w = 1 / (s**2 + tau2)
    ybar_w = np.average(y, weights=w)
    signal_share = 1 - np.mean(s**2) / np.var(y)
    grp = mods.reindex(genes).fillna(-1).astype(int).values
    rel = s <= np.quantile(s, 1 / 3)
    rows = []
    for m in MODELS:
        key = f"{m}__{target_name}"
        p = np.zeros(len(genes)) if m == "zero" else P.loc[genes, key].values
        ok = np.isfinite(p)
        yy, pp, ww, gg = y[ok], p[ok], w[ok], grp[ok]
        sse_w = np.sum(ww * (yy - pp) ** 2)
        tm = TRAIN_MEAN.get(f"F8_ols__{target_name}", 0.0)
        r = dict(study=study, contrast=tag, role=role, training_target=target_name, model=m, n_genes=int(ok.sum()),
                 tau2_ext=tau2, signal_share=signal_share,
                 R2w_vs_zero=1 - sse_w / np.sum(ww * yy**2),
                 R2w_vs_training_mean=1 - sse_w / np.sum(ww * (yy - tm) ** 2),
                 R2w_vs_own_mean=1 - sse_w / np.sum(ww * (yy - np.average(yy, weights=ww)) ** 2),
                 R2_unweighted_vs_own_mean=1 - np.sum((yy - pp) ** 2) / np.sum((yy - yy.mean()) ** 2),
                 MAE=np.mean(np.abs(yy - pp)), RMSE=np.sqrt(np.mean((yy - pp) ** 2)), pred_sd=pp.std(), obs_sd=yy.std())
        if m != "zero":
            r["spearman"] = stats.spearmanr(yy, pp).statistic
            r["spearman_ci_low"], r["spearman_ci_high"] = block_boot_rho(yy, pp, gg)
            r["sign_agree_all"] = float(np.mean(np.sign(yy) == np.sign(pp)))
            top = np.abs(pp) >= np.quantile(np.abs(pp), .9)
            r["sign_agree_top10pct_pred"] = float(np.mean(np.sign(yy[top]) == np.sign(pp[top])))
            r["calibration_slope_w"] = float(np.cov(pp, yy, aweights=ww)[0, 1] / np.cov(pp, aweights=ww))
            r["spearman_reliable_tercile"] = stats.spearmanr(yy[rel[ok]], pp[rel[ok]]).statistic
            r["R2w_vs_own_mean_reliable_tercile"] = 1 - np.sum(ww[rel[ok]] * (yy[rel[ok]] - pp[rel[ok]]) ** 2) / np.sum(
                ww[rel[ok]] * (yy[rel[ok]] - np.average(yy[rel[ok]], weights=ww[rel[ok]])) ** 2)
            if "limmavoom" not in tag and "6mice" not in tag:
                obs, null = proxy_null(study, exposed, genes[ok], stats.rankdata(pp))
                r.update(perm_rho_proxy_obs=obs, perm_n_labellings=len(null), perm_p_one_sided=float(np.mean(null >= obs - 1e-12)),
                         perm_null_q95=float(np.quantile(null, .95)))
        rows.append(r)
    return rows


jobs = [(st, tag, ex, role, tgt) for st, tag, ex, role in CONTRASTS for tgt in ("PE_in_CD", "PE_in_MD")]
out = Parallel(n_jobs=N_JOBS, backend="loky", verbose=5)(delayed(evaluate)(*j) for j in jobs)
R = pd.DataFrame([r for o in out for r in o])


def verdict(r):
    if r.model == "zero":
        return ""
    if r.signal_share <= 0.2:
        return "UNINTERPRETABLE (target mostly noise)"
    fail = (r.R2w_vs_own_mean <= 0) or (r.spearman_ci_low <= 0 <= r.spearman_ci_high) or (r.spearman_ci_high < 0)
    if fail:
        return "FAIL"
    return "PASS-direction" if (r.spearman > 0 and r.get("perm_p_one_sided", 1) <= 0.10) else "positive but permutation-unsupported"


R["verdict_locked_rule"] = R.apply(verdict, axis=1)
R.to_csv(P5 / "gene_level_external_results.csv", index=False)
# external contrasts table
cr = []
for st, tag, ex, role in CONTRASTS:
    t = load_target(st, tag)
    X, S = load_prepared(st)
    cr.append(dict(study=st, contrast=tag, role=role, exposed_groups=";".join(ex), n_control=int((S.group == "CON").sum()),
                   n_exposed=int(S.group.isin(ex).sum()), method="limma-voom (unrounded)" if "limmavoom" in tag else (
                       "limma-trend on log2(norm+1)" if st == "GSE277772" else "DESeq2 Wald MLE"),
                   genes_tested=len(t), genes_in_frozen_feature_table=len(P.index.intersection(t.index)),
                   n_padj_lt_0_05=int((t.padj < .05).sum()), n_padj_lt_0_10=int((t.padj < .10).sum()),
                   median_lfcSE=float(t.lfcSE.median()),
                   signal_share=float(1 - np.mean(t.lfcSE**2) / np.var(t.log2FC))))
pd.DataFrame(cr).to_csv(P5 / "external_contrasts.csv", index=False)
# per-gene table for the primary contrast (figures)
t = load_target("HONG", "PE_vs_CON")
g = P.index.intersection(t.index)
pd.concat([t.loc[g, ["log2FC", "lfcSE", "padj"]], P.loc[g]], axis=1).to_csv(RES / "gene_level_primary_merged.csv.gz")
cols = ["study", "contrast", "training_target", "model", "n_genes", "signal_share", "R2w_vs_zero", "R2w_vs_training_mean", "R2w_vs_own_mean",
        "spearman", "spearman_ci_low", "spearman_ci_high", "sign_agree_all", "perm_p_one_sided", "perm_n_labellings", "verdict_locked_rule"]
with pd.option_context("display.width", 300, "display.max_rows", 200):
    print(pd.DataFrame(cr).round(3).to_string(index=False))
    print(R[cols].round(3).to_string(index=False))
