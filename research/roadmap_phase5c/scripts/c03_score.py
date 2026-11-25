"""Phase 5C step 3: LOCKED evaluation of frozen Phase 5 / 4D / 4E artefacts on PRJNA916285 (Z23).
Implements 03_locked_protocol.md exactly. Nothing is fitted, tuned or reoriented on Z23.

Unit of inference = mouse. Gene-level: Spearman(frozen prediction, Z23 DESeq2 MLE log2FC) with an exhaustive
mouse-relabelling null on the log2 group-mean-CPM-ratio proxy (as Phase 5/5B) and a within-arm mouse bootstrap.
Mouse-level: frozen singscores / logistic log-odds, exact permutation over all labellings, LOO, Hedges g.
Usage: python c03_score.py [--dry-run]   (dry run = labels replaced by a fixed non-biological split; prints nothing numeric)
"""
import itertools
import json
import sys

import numpy as np
import pandas as pd
from c00_common import (FRZ, P5, P5B, PAR, RES, RS, SEED, TAB, WORK, holm, load_counts, log2cpm, score_4d, score_4e, verify_frozen,
                        welch_ci)
from scipy import stats

DRY = "--dry-run" in sys.argv
B = 1000
OUT = RES / ("dryrun" if DRY else "")
OUT.mkdir(parents=True, exist_ok=True)
TOUT = OUT if DRY else TAB
frozen_ok = verify_frozen()
X, S = load_counts()
if DRY:  # label-blind functional test: alternate assignment unrelated to treatment
    S = S.copy()
    S["group"] = np.where(np.arange(len(S)) % 2 == 0, S.group.str.replace("MPs", "").replace({"EGCG": "EGCG", "": "PBS"}), S.group)

DE = RES / "de"
P = pd.read_csv(FRZ / "gene_level_frozen_predictions.csv.gz", index_col=0)
CAND5B = pd.read_csv(P5B / "results/candidate_predictions_all_genes.csv.gz", index_col=0)
mods = pd.read_csv(RS / "roadmap_phase4/results/gene_groups.csv", index_col=0)["module"]

# ---------------- prespecified model sets ----------------
PRIMARY_GENE = "F8_ols__PE_in_CD"
SECONDARY_GENE = ["F8_ols__PE_in_MD", "F8C_ridge__PE_in_CD", "FAST_only_ols__PE_in_CD", "MASH_only_ols__PE_in_CD", "observed_GSE294180_ashr__PE_in_CD"]
EXPLORATORY_GENE = [c for c in CAND5B.columns if c.startswith(("I2b", "I2c"))]  # Phase 5B post-hoc models
CONTRASTS = [  # name, control group(s), exposed group(s), stratified?, role
    ("MPs_vs_PBS", ["PBS"], ["MPs"], False, "PRIMARY"),
    ("EGCGMPs_vs_EGCG", ["EGCG"], ["EGCGMPs"], False, "SECONDARY (MP on EGCG background)"),
    ("MPmain_additive", ["PBS", "EGCG"], ["MPs", "EGCGMPs"], True, "SECONDARY (16 mice, stratified by EGCG)"),
    ("EGCG_vs_PBS", ["PBS"], ["EGCG"], False, "SPECIFICITY COMPARATOR (non-plastic perturbation)"),
]
SCORES_PRIMARY = "FAST_singscore"
SCORES_SECONDARY = ["MASH_singscore", "PPARA_singscore", "R0_logistic_logodds", "R2_logistic_logodds"]


def mice_for(ctrl, expo):
    m = S[S.group.isin(ctrl + expo)]
    return m, m.group.isin(expo).values, (m.group.isin(["EGCG", "EGCGMPs"]).values)


def labellings(y, strata, stratified):
    """All relabellings preserving group sizes (within EGCG stratum if stratified). Yields boolean exposed vectors."""
    n = len(y)
    if not stratified:
        k = int(y.sum())
        for ex in itertools.combinations(range(n), k):
            v = np.zeros(n, bool)
            v[list(ex)] = True
            yield v
        return
    i0, i1 = np.where(~strata)[0], np.where(strata)[0]
    k0, k1 = int(y[i0].sum()), int(y[i1].sum())
    for a in itertools.combinations(i0, k0):
        for b in itertools.combinations(i1, k1):
            v = np.zeros(n, bool)
            v[list(a) + list(b)] = True
            yield v


def proxy_diff(L, y, strata, stratified):
    if not stratified:
        return L[:, y].mean(1) - L[:, ~y].mean(1)
    return np.mean([L[:, y & (strata == s)].mean(1) - L[:, ~y & (strata == s)].mean(1) for s in (False, True)], axis=0)


def arith_proxy(Xc, y, strata, stratified):
    """log2 ratio of group-mean CPM (+1): the fast MLE surrogate used in Phase 5B."""
    C = (Xc / Xc.sum(0) * 1e6)
    if not stratified:
        return np.log2(C[:, y].mean(1) + 1) - np.log2(C[:, ~y].mean(1) + 1)
    return np.mean([np.log2(C[:, y & (strata == s)].mean(1) + 1) - np.log2(C[:, ~y & (strata == s)].mean(1) + 1) for s in (False, True)], axis=0)


def dl_tau2(y, s):
    w = 1 / s**2
    yb = (w * y).sum() / w.sum()
    q = (w * (y - yb) ** 2).sum()
    c = w.sum() - (w**2).sum() / w.sum()
    return max(0.0, (q - (len(y) - 1)) / c)


# ---------------- gene level ----------------
gene_rows, null_store, boot_store = [], {}, {}
rng = np.random.default_rng(SEED)
for name, ctrl, expo, strat, role in CONTRASTS:
    f = DE / f"Z23__{name}.csv.gz"
    if not f.exists() or DRY:
        if not DRY:
            gene_rows.append(dict(contrast=name, role=role, status="SKIPPED (contrast not estimable: mice unavailable)"))
            continue
    m, y, strata = mice_for(ctrl, expo)
    if y.sum() < 3 or (~y).sum() < 3 or (strat and len(m) < 16):
        gene_rows.append(dict(contrast=name, role=role, status="SKIPPED (mice unavailable)"))
        continue
    if DRY:
        t = pd.DataFrame({"log2FC": rng.normal(size=len(P)), "lfcSE": 0.3}, index=P.index)
    else:
        t = pd.read_csv(f).set_index("gene")
        t = t[np.isfinite(t.log2FC) & np.isfinite(t.lfcSE) & (t.lfcSE > 0)]
    genes = P.index.intersection(t.index).intersection(X.index)
    yv, sv = t.loc[genes, "log2FC"].values, t.loc[genes, "lfcSE"].values
    tau2 = dl_tau2(yv, sv)
    w = 1 / (sv**2 + tau2)
    Xc = X.loc[genes, m.index].values.astype(float)
    Lab = list(labellings(y, strata, strat))
    true_idx = next(i for i, v in enumerate(Lab) if np.array_equal(v, y))
    diffs = np.array([arith_proxy(Xc, v, strata, strat) for v in Lab], dtype=np.float32)  # labellings x genes
    # mouse bootstrap (within arm [and stratum]) of the proxy target
    groups = [np.where((y == a) & (strata == s))[0] for a in (False, True) for s in (False, True) if ((y == a) & (strata == s)).any()]
    bidx = [np.concatenate([rng.choice(g, len(g), replace=True) for g in groups]) for _ in range(B)]
    bdiffs = [arith_proxy(Xc[:, ii], y[ii], strata[ii], strat) for ii in bidx]
    models = [PRIMARY_GENE] + SECONDARY_GENE + EXPLORATORY_GENE
    for mod in models:
        p = (CAND5B if mod in EXPLORATORY_GENE else P).reindex(genes)[mod].values
        ok = np.isfinite(p)
        pr = stats.rankdata(p[ok])
        rho = stats.spearmanr(yv[ok], p[ok]).statistic
        null = np.array([stats.spearmanr(pr, d[ok]).statistic for d in diffs])
        obs_proxy = null[true_idx]
        rank_true = int((null > obs_proxy + 1e-12).sum() + 1)
        boot = np.array([stats.spearmanr(pr, d[ok]).statistic for d in bdiffs]) if mod in (PRIMARY_GENE, *SECONDARY_GENE[:1]) else np.array([np.nan])
        sse = np.sum(w[ok] * (yv[ok] - p[ok]) ** 2)
        r = dict(contrast=name, role=role, model=mod,
                 model_status="PRIMARY frozen" if mod == PRIMARY_GENE else ("secondary frozen" if mod in SECONDARY_GENE else "EXPLORATORY post-hoc (Phase 5B)"),
                 n_genes=int(ok.sum()), n_ctrl=int((~y).sum()), n_exp=int(y.sum()), signal_share=float(1 - np.mean(sv**2) / np.var(yv)),
                 spearman_mle=rho, spearman_proxy=obs_proxy, n_labellings=len(Lab), rank_true_of_labellings=rank_true,
                 p_perm_upper=rank_true / len(Lab), p_perm_lower=float(np.mean(null <= obs_proxy + 1e-12)),
                 null_q95=float(np.quantile(null, .95)), null_median=float(np.median(null)),
                 boot_lo=float(np.nanquantile(boot, .025)) if np.isfinite(boot).any() else np.nan,
                 boot_hi=float(np.nanquantile(boot, .975)) if np.isfinite(boot).any() else np.nan,
                 R2w_vs_own_mean=1 - sse / np.sum(w[ok] * (yv[ok] - np.average(yv[ok], weights=w[ok])) ** 2),
                 R2w_vs_zero=1 - sse / np.sum(w[ok] * yv[ok] ** 2),
                 sign_agree=float(np.mean(np.sign(yv[ok]) == np.sign(p[ok]))),
                 sign_agree_top10pct_pred=float(np.mean((np.sign(yv[ok]) == np.sign(p[ok]))[np.abs(p[ok]) >= np.quantile(np.abs(p[ok]), .9)])),
                 calibration_slope_w=float(np.cov(p[ok], yv[ok], aweights=w[ok])[0, 1] / np.cov(p[ok], aweights=w[ok])),
                 MAE=float(np.mean(np.abs(yv[ok] - p[ok]))), pred_sd=float(p[ok].std()), obs_sd=float(yv[ok].std()))
        gene_rows.append(r)
        if mod == PRIMARY_GENE:
            null_store[name] = dict(null=null.tolist(), true_idx=true_idx, boot=boot.tolist())
    # mouse-state diagnostic: how much of the relabelling rho is explained by mouse PCs (as Phase 5B)
    if name in ("MPs_vs_PBS", "EGCGMPs_vs_EGCG", "EGCG_vs_PBS"):
        L = log2cpm(X.loc[genes, m.index]).values
        Lc = L - L.mean(1, keepdims=True)
        U, s_, Vt = np.linalg.svd(Lc, full_matrices=False)
        pcs = Vt[:3].T * s_[:3]
        D = np.array([pcs[v].mean(0) - pcs[~v].mean(0) for v in Lab])
        nl = np.array(null_store[name]["null"])
        Z = np.column_stack([np.ones(len(D)), D])
        beta, *_ = np.linalg.lstsq(Z, nl, rcond=None)
        r2 = 1 - np.sum((nl - Z @ beta) ** 2) / np.sum((nl - nl.mean()) ** 2)
        null_store[name].update(pc_var_explained=(s_[:3] ** 2 / np.sum(s_**2)).tolist(), r2_rho_on_pc_diffs=float(r2),
                                pc_scores={mid: pcs[i].tolist() for i, mid in enumerate(m.index)})
G = pd.DataFrame(gene_rows)


def gene_verdict(r):
    if r.get("model") != PRIMARY_GENE or not isinstance(r.get("spearman_mle"), float):
        return ""
    if r.spearman_proxy > 0 and r.p_perm_upper <= 0.05 and r.spearman_mle > 0:
        return "SUPPORTIVE"
    if r.spearman_proxy < 0 and r.p_perm_lower <= 0.05:
        return "CONTRADICTORY"
    return "INCONCLUSIVE"


if len(G) and "model" in G:
    G["verdict_locked"] = G.apply(gene_verdict, axis=1)
    G["p_perm_upper_holm_secondary"] = np.nan
    for c, d in G.groupby("contrast"):
        i = d.index[d.model.isin(SECONDARY_GENE)]
        if len(i):
            G.loc[i, "p_perm_upper_holm_secondary"] = holm(G.loc[i, "p_perm_upper"].values)
G.to_csv(TOUT / "z23_gene_level_results.csv", index=False)
(OUT / "z23_gene_level_nulls.json").write_text(json.dumps(null_store))

# ---------------- mouse level ----------------
d4, cov4, nb4 = score_4d(X, "counts")
SC = d4.rename(columns=lambda c: f"{c}_singscore")
feats = None
cov = [dict(score=f"Phase4D_{k}", background=nb4, **v) for k, v in cov4.items()]
for rep in ("R0", "R2"):
    lo, F, cv, nb = score_4e(X, rep)
    SC[f"{rep}_logistic_logodds"] = lo
    if rep == "R2":
        feats = F
        cov.append(dict(score="Phase4E_R0R2", background=nb, background_frozen=len(PAR[rep]["background_genes"]), **cv))
SC = SC.join(S[["group", "mp_exposed", "egcg"]])
SC.index.name = "mouse_id"
SC.to_csv(OUT / "z23_per_mouse_scores.csv")
pd.DataFrame(cov).to_csv(OUT / "z23_score_gene_coverage.csv", index=False)
p5s = pd.read_csv(P5 / "sample_level_external_results.csv")
ref = p5s[(p5s.study == "HONG") & (p5s.contrast == "PE_vs_CON")].set_index("score")
mrows = []
for name, ctrl, expo, strat, role in CONTRASTS:
    m, y, strata = mice_for(ctrl, expo)
    if y.sum() < 3 or (~y).sum() < 3 or (strat and len(m) < 16):
        mrows.append(dict(contrast=name, role=role, status="SKIPPED (mice unavailable)"))
        continue
    Lab = list(labellings(y, strata, strat))
    for sc in [SCORES_PRIMARY] + SCORES_SECONDARY:
        v = SC.loc[m.index, sc].values
        stat = lambda yy: (np.mean([v[yy & (strata == s)].mean() - v[~yy & (strata == s)].mean() for s in (False, True)])  # noqa: E731
                           if strat else v[yy].mean() - v[~yy].mean())
        obs = stat(y)
        null = np.array([stat(l_) for l_ in Lab])
        w = welch_ci(v[y], v[~y]) if not strat else dict(estimate=obs)
        loo = []
        for i in range(len(v)):
            k = np.arange(len(v)) != i
            vv, yy, ss = v[k], y[k], strata[k]
            loo.append(np.sign(np.mean([vv[yy & (ss == s)].mean() - vv[~yy & (ss == s)].mean() for s in (False, True)]) if strat
                               else vv[yy].mean() - vv[~yy].mean()))
        bidx = [np.concatenate([rng.choice(np.where((y == a) & (strata == s))[0], int(((y == a) & (strata == s)).sum()))
                                for a in (False, True) for s in (False, True) if ((y == a) & (strata == s)).any()]) for _ in range(B)]
        bt = np.array([(np.mean([v[ii][y[ii] & (strata[ii] == s)].mean() - v[ii][~y[ii] & (strata[ii] == s)].mean() for s in (False, True)])
                        if strat else v[ii][y[ii]].mean() - v[ii][~y[ii]].mean()) for ii in bidx])
        g294 = float(ref.loc[sc, "gse294180_PE_in_CD"]) if sc in ref.index else np.nan
        mrows.append(dict(contrast=name, role=role, score=sc, score_status="PRIMARY" if sc == SCORES_PRIMARY else "secondary",
                          n_ctrl=int((~y).sum()), n_exp=int(y.sum()), **{k_: w.get(k_) for k_ in ("estimate", "se", "ci_low", "ci_high", "hedges_g")},
                          boot_lo=float(np.quantile(bt, .025)), boot_hi=float(np.quantile(bt, .975)), n_labellings=len(Lab),
                          p_perm_two_sided=float(np.mean(np.abs(null) >= abs(obs) - 1e-12)), p_perm_min=float(np.mean(np.abs(null) >= np.abs(null).max() - 1e-12)),
                          loo_same_sign=int(np.sum(np.array(loo) == np.sign(obs))), loo_n=len(loo),
                          gse294180_PE_in_CD=g294, hong_PE_vs_CON=float(ref.loc[sc, "estimate"]) if sc in ref.index else np.nan,
                          same_sign_as_gse294180=bool(np.sign(obs) == np.sign(g294)) if np.isfinite(g294) else None))
M = pd.DataFrame(mrows)


def mouse_verdict(r):
    if r.get("score") != SCORES_PRIMARY or r.contrast == "EGCG_vs_PBS":
        return ""
    if r.same_sign_as_gse294180 and r.p_perm_two_sided <= 0.05 and r.loo_same_sign == r.loo_n:
        return "SUPPORTIVE"
    if (not r.same_sign_as_gse294180) and r.p_perm_two_sided <= 0.05:
        return "CONTRADICTORY"
    return "INCONCLUSIVE"


if "score" in M:
    M["verdict_locked"] = M.apply(mouse_verdict, axis=1)
    M["p_perm_holm_secondary4"] = np.nan
    for c, d in M.groupby("contrast"):
        i = d.index[d.score.isin(SCORES_SECONDARY)]
        if len(i):
            M.loc[i, "p_perm_holm_secondary4"] = holm(M.loc[i, "p_perm_two_sided"].values)
M.to_csv(TOUT / "z23_mouse_level_results.csv", index=False)

# R2 decomposition (descriptive), same convention as Phase 5
coef_sd = pd.Series(np.array(PAR["R2"]["coef"]) / np.array(PAR["R2"]["scaler_scale"]), index=PAR["R2"]["feature_names"])
dec = []
for name, ctrl, expo, strat, role in CONTRASTS[:2] + CONTRASTS[3:]:
    m, y, _ = mice_for(ctrl, expo)
    if y.sum() < 3 or (~y).sum() < 3:
        continue
    F = feats.loc[m.index]
    diff = F[y].mean() - F[~y].mean()
    for f_, c_ in (coef_sd * diff).items():
        dec.append(dict(contrast=name, feature=f_, feature_diff=diff[f_], contribution_logodds=c_))
pd.DataFrame(dec).to_csv(OUT / "z23_r2_decomposition.csv", index=False)
json.dump(dict(frozen_hashes_verified=all(frozen_ok.values()), dry_run=DRY, B=B, seed=SEED, n_mice=len(S)), open(OUT / "c03_run_info.json", "w"), indent=2)
if DRY:
    print("dry run OK:", len(G), "gene rows,", len(M), "mouse rows")
else:
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(G[["contrast", "model", "n_genes", "spearman_mle", "spearman_proxy", "rank_true_of_labellings", "n_labellings", "p_perm_upper",
                 "boot_lo", "boot_hi", "R2w_vs_own_mean", "verdict_locked"]].round(3).to_string(index=False))
        print(M[["contrast", "score", "estimate", "ci_low", "ci_high", "hedges_g", "p_perm_two_sided", "p_perm_min", "loo_same_sign", "loo_n",
                 "gse294180_PE_in_CD", "hong_PE_vs_CON", "verdict_locked"]].round(4).to_string(index=False))
