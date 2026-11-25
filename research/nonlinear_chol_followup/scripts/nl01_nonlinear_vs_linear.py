"""Nonlinear vs linear held-out prediction for lipid gene groups (exploratory; 2026-11-25).

No model is fitted or tuned. Reads the ORIGINAL Phase 4 task-A held-out predictions
(research/roadmap_phase4/results/taskA/predictions/primary__PE_in_{CD,MD}.csv.gz):
  F8-OLS = 'ols_all', F8-Ridge = 'ridge', F8-CatBoost = 'catboost', F8-HGB = 'hgb'
  (F8 = 4 reference effects + 4 SEs: MASH, FAST, PPARA, CIRC; no pathway-membership features).
Splits: 20 co-expression modules, 5 outer folds x 3 repeats (seeds 101-103); tuning by inner GroupKFold(4) on training
genes only (Phase 4 protocol). Null baseline per split and module labels: from research/lipid_prediction_followup
(results/oof_predictions_with_null.csv.gz; identical splits).

Decision rules fixed before any group-level result was computed:
  Primary family (8 tests): CHOL group, gene-paired gain (MSE_linear - MSE_nonlinear) for
    {CatBoost, HGB} x {OLS, Ridge} x {CD, MD}; p = module-cluster t (one-sample mean of per-gene gain,
    CR1-type cluster SE, t with G-1 df); Holm within family. Sensitivity p: exact module sign-flip.
  Secondary family: same pairs for LIPO, PPARA, OTHER and CHOL-minus-OTHER / CHOL-minus-matched
    difference-in-gain; BH within family.
  Feature-combination follow-up runs only if a primary test has Holm p < 0.05 with gain > 0.
Group metrics: R2_group = 1 - SSE/SSE_null, RMSE, direction agreement (share of gene x repeat with
  sign(pred) == sign(y)), Spearman(y, repeat-mean pred). CIs: module block bootstrap (B = 2000).
Matched reference: 10 nearest OTHER genes per group gene on repeat-mean null squared error (as lp03).
Run: .venv/Scripts/python.exe research/nonlinear_chol_followup/scripts/nl01_nonlinear_vs_linear.py
"""
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

RS = Path(__file__).resolve().parents[2]
OUT = RS / "nonlinear_chol_followup"
RES = OUT / "results"
RES.mkdir(parents=True, exist_ok=True)
SEED, B, K = 20261125, 2000, 10
MODELS = {"ols_all": "F8-OLS", "ridge": "F8-Ridge", "catboost": "F8-CatBoost", "hgb": "F8-HGB"}
LIN, NL = ["F8-OLS", "F8-Ridge"], ["F8-CatBoost", "F8-HGB"]
SETS = {"CHOL": "CHOL_SYNTHESIS_CORE", "LIPO": "LIPOGENESIS_CORE", "PPARA": "SANDERSON_PPARA_TARGETS"}
TARGETS = ("PE_in_CD", "PE_in_MD")

# ---------------- load original predictions + null + modules
pr = []
for t in TARGETS:
    d = pd.read_csv(RS / f"roadmap_phase4/results/taskA/predictions/primary__{t}.csv.gz")
    d = d[d.model.isin(MODELS)].copy()
    d["target"], d["model_id"] = t, d.model.map(MODELS)
    pr.append(d)
P = pd.concat(pr, ignore_index=True)
N = pd.read_csv(RS / "lipid_prediction_followup/results/oof_predictions_with_null.csv.gz",
                usecols=["gene", "target", "repeat", "fold", "y", "null", "module", "model_id"])
N = N[N.model_id == "F8-OLS"].drop(columns="model_id").rename(columns={"y": "y_lp"})
P = P.merge(N, on=["gene", "target", "repeat", "fold"], how="left", validate="many_to_one")
assert P.null.notna().all(), "fold assignment mismatch between Phase 4 predictions and lipid follow-up splits"
assert np.allclose(P.y, P.y_lp)
assert (P.groupby(["model_id", "target", "repeat"]).gene.nunique() == 11758).all()
# F8-OLS here must equal F8-OLS used in the lipid follow-up (Phase 4C file)
lpo = pd.read_csv(RS / "lipid_prediction_followup/results/oof_predictions_with_null.csv.gz")
lpo = lpo[lpo.model_id == "F8-OLS"].set_index(["gene", "target", "repeat"]).pred
p4o = P[P.model_id == "F8-OLS"].set_index(["gene", "target", "repeat"]).pred
max_diff_ols_vs_lp = float((p4o - lpo.loc[p4o.index]).abs().max())

# ---------------- reproduction check vs stored Phase 4 per-split metrics
M = pd.read_csv(RS / "roadmap_phase4/results/taskA/metrics_primary_positive.csv")
chk = []
for (mo, t, r, k), g in P.groupby(["model", "target", "repeat", "fold"]):
    m = M[(M.run == f"primary__{t}") & (M.model == mo) & (M.repeat == r) & (M.fold == k)]
    sse, sse0 = ((g.y - g.pred) ** 2).sum(), ((g.y - g.null) ** 2).sum()
    chk.append(dict(model=MODELS[mo], target=t, repeat=r, fold=k, rmse_stored=m.rmse.iloc[0],
                    rmse_recomputed=float(np.sqrt(((g.y - g.pred) ** 2).mean())), r2_stored=m.r2_oos.iloc[0],
                    r2_recomputed=float(1 - sse / sse0), sign_stored=m.sign_agree.iloc[0],
                    sign_recomputed=float((np.sign(g.pred) == np.sign(g.y)).mean())))
C = pd.DataFrame(chk)
for a in ("rmse", "r2", "sign"):
    C[f"absdiff_{a}"] = (C[f"{a}_stored"] - C[f"{a}_recomputed"]).abs()
C.to_csv(RES / "reproduction_check_vs_phase4.csv", index=False)

# ---------------- per-gene summaries
P["se"], P["se0"] = (P.y - P.pred) ** 2, (P.y - P.null) ** 2
P["dir"] = (np.sign(P.pred) == np.sign(P.y)).astype(float)
G = P.groupby(["model_id", "target", "gene"]).agg(y=("y", "first"), module=("module", "first"), mse=("se", "mean"),
                                                   mse0=("se0", "mean"), pred=("pred", "mean"), dir=("dir", "mean")).reset_index()
gs = pd.read_csv(RS / "roadmap_phase3/results/gene_sets.csv")
genes = np.array(sorted(G.gene.unique()))
mods = np.array(sorted(G.module.unique()))
GROUPS = {k: np.isin(genes, gs.gene[gs.set == v]) for k, v in SETS.items()}
GROUPS["OTHER"] = ~(GROUPS["CHOL"] | GROUPS["LIPO"] | GROUPS["PPARA"])
W = {(m, t): G[(G.model_id == m) & (G.target == t)].set_index("gene").loc[genes] for m in MODELS.values() for t in TARGETS}
modv = W[("F8-OLS", "PE_in_CD")].module.values

# matched OTHER (null squared error, model-independent)
MATCH = {}
bal = []
for t in TARGETS:
    m0 = W[("F8-OLS", t)].mse0.values
    for g in ("CHOL", "LIPO", "PPARA"):
        gm, avail, ref = GROUPS[g], GROUPS["OTHER"].copy(), np.zeros(len(genes), bool)
        for i in np.argsort(-m0 * gm)[: gm.sum()]:
            cand = np.where(avail)[0]
            pick = cand[np.argsort(np.abs(m0[cand] - m0[i]))[:K]]
            ref[pick], avail[pick] = True, False
        MATCH[(t, g)] = ref
        bal.append(dict(diet=t, group=g, n_group=int(gm.sum()), n_matched=int(ref.sum()), median_mse0_group=float(np.median(m0[gm])),
                        median_mse0_matched=float(np.median(m0[ref])), n_modules_matched=int(len(set(modv[ref])))))
pd.DataFrame(bal).to_csv(RES / "matching_balance.csv", index=False)

rng = np.random.default_rng(SEED)
BOOTW = rng.multinomial(len(mods), np.ones(len(mods)) / len(mods), size=B)



def group_stats(df, mask):
    x = df[mask]
    S = x.groupby("module").agg(n=("mse", "size"), sse=("mse", "sum"), sse0=("mse0", "sum"), dirs=("dir", "sum"),
                                 ss=("mse", "sum")).reindex(mods, fill_value=0).values
    bs = BOOTW @ S
    with np.errstate(invalid="ignore", divide="ignore"):
        br2, brmse, bdir = 1 - bs[:, 1] / bs[:, 2], np.sqrt(bs[:, 1] / bs[:, 0]), bs[:, 3] / bs[:, 0]
    q = lambda a: (np.nanquantile(a, .025), np.nanquantile(a, .975))  # noqa: E731
    return dict(n_genes=int(mask.sum()), n_modules=int(x.module.nunique()), R2_group=1 - x.mse.sum() / x.mse0.sum(),
                R2_ci=q(br2), rmse=float(np.sqrt(x.mse.mean())), rmse_ci=q(brmse), rmse_null=float(np.sqrt(x.mse0.mean())),
                direction_agreement=float(x.dir.mean()), dir_ci=q(bdir),
                spearman=float(stats.spearmanr(x.y, x.pred).statistic))


def cluster_mean_test(v, m):
    """One-sample mean of per-gene values v, clustered by module m: CR1 SE, t(G-1); exact/MC sign-flip at module level."""
    u = np.unique(m)
    Gc, n = len(u), len(v)
    mean = v.mean()
    sums = np.array([(v[m == g] - mean).sum() for g in u])
    se = np.sqrt((Gc / (Gc - 1)) * (sums**2).sum()) / n
    t = mean / se if se > 0 else np.nan
    p_t = float(2 * stats.t.sf(abs(t), Gc - 1)) if Gc > 1 else np.nan
    tot = np.array([v[m == g].sum() for g in u])
    if Gc <= 16:
        flips = np.array(list(itertools.product([1, -1], repeat=Gc)))
    else:
        flips = np.random.default_rng(SEED).choice([1, -1], size=(20000, Gc))
    null = flips @ tot
    p_sf = float(np.mean(np.abs(null) >= abs(tot.sum()) - 1e-12))
    qt = stats.t.ppf(.975, Gc - 1) if Gc > 1 else np.nan
    return mean, mean - qt * se, mean + qt * se, p_t, p_sf, Gc


def diff_test(v, a_mask, b_mask, m):
    """Difference in mean of v between masks a and b with module FE (CR1 cluster SE by module, t(G-1))."""
    keep = a_mask | b_mask
    y, x, mm = v[keep], a_mask[keep].astype(float), m[keep]
    D = pd.get_dummies(mm).values.astype(float)
    X = np.column_stack([x, D])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    e = y - X @ b
    XtXi = np.linalg.pinv(X.T @ X)
    u = np.unique(mm)
    meat = sum(np.outer(X[mm == g].T @ e[mm == g], X[mm == g].T @ e[mm == g]) for g in u)
    n, k, Gc = len(y), X.shape[1], len(u)
    V = XtXi @ meat @ XtXi * (Gc / (Gc - 1)) * ((n - 1) / (n - k))
    se = np.sqrt(V[0, 0])
    qt = stats.t.ppf(.975, Gc - 1)
    return b[0], b[0] - qt * se, b[0] + qt * se, float(2 * stats.t.sf(abs(b[0] / se), Gc - 1))


# ---------------- (A) group metrics for every model
rows = []
for t in TARGETS:
    for mid in MODELS.values():
        df = W[(mid, t)]
        for g in ("CHOL", "LIPO", "PPARA", "OTHER"):
            s = group_stats(df, GROUPS[g])
            ms = group_stats(df, MATCH[(t, g)]) if g != "OTHER" else None
            rows.append(dict(diet=t, model=mid, group=g, n_genes=s["n_genes"], n_modules=s["n_modules"],
                             R2_group=s["R2_group"], R2_ci_low=s["R2_ci"][0], R2_ci_high=s["R2_ci"][1],
                             rmse=s["rmse"], rmse_ci_low=s["rmse_ci"][0], rmse_ci_high=s["rmse_ci"][1], rmse_null=s["rmse_null"],
                             direction_agreement=s["direction_agreement"], dir_ci_low=s["dir_ci"][0], dir_ci_high=s["dir_ci"][1],
                             spearman=s["spearman"],
                             matched_R2_group=None if ms is None else ms["R2_group"],
                             matched_R2_ci_low=None if ms is None else ms["R2_ci"][0],
                             matched_R2_ci_high=None if ms is None else ms["R2_ci"][1],
                             matched_direction_agreement=None if ms is None else ms["direction_agreement"],
                             R2_minus_matched=None if ms is None else s["R2_group"] - ms["R2_group"]))
GM = pd.DataFrame(rows)
GM.to_csv(RES / "group_metrics_by_model.csv", index=False)

# ---------------- (B) gene-paired nonlinear vs linear
pair_rows = []
for t in TARGETS:
    for nl in NL:
        for li in LIN:
            gain = W[(li, t)].mse.values - W[(nl, t)].mse.values  # >0: nonlinear lower held-out squared error
            ddir = W[(nl, t)].dir.values - W[(li, t)].dir.values
            for g in ("CHOL", "LIPO", "PPARA", "OTHER"):
                mk = GROUPS[g]
                mean, lo, hi, p_t, p_sf, Gc = cluster_mean_test(gain[mk], modv[mk])
                r2nl = 1 - W[(nl, t)].mse.values[mk].sum() / W[(nl, t)].mse0.values[mk].sum()
                r2li = 1 - W[(li, t)].mse.values[mk].sum() / W[(li, t)].mse0.values[mk].sum()
                S1 = pd.DataFrame({"m": modv[mk], "a": W[(nl, t)].mse.values[mk], "b": W[(li, t)].mse.values[mk],
                                   "z": W[(li, t)].mse0.values[mk]}).groupby("m").sum().reindex(mods, fill_value=0).values
                bs = BOOTW @ S1
                with np.errstate(invalid="ignore", divide="ignore"):
                    bd = (bs[:, 1] - bs[:, 0]) / bs[:, 2]
                row = dict(diet=t, nonlinear=nl, linear=li, group=g, n_genes=int(mk.sum()), n_modules=Gc,
                           R2_nonlinear=r2nl, R2_linear=r2li, dR2=r2nl - r2li, dR2_ci_low=np.nanquantile(bd, .025),
                           dR2_ci_high=np.nanquantile(bd, .975), mean_gain_mse=mean, gain_ci_low=lo, gain_ci_high=hi,
                           p_cluster_t=p_t, p_signflip=p_sf, signflip_min_p=2 / 2**Gc if Gc <= 16 else 1 / 20000,
                           n_genes_nonlinear_better=int((gain[mk] > 0).sum()),
                           d_direction_agreement=float(ddir[mk].mean()))
                if g != "OTHER":
                    for ref, rmask in (("OTHER", GROUPS["OTHER"]), ("matched", MATCH[(t, g)])):
                        b_, l_, h_, p_ = diff_test(gain, mk, rmask, modv)
                        row.update({f"gain_minus_{ref}": b_, f"gain_minus_{ref}_ci_low": l_, f"gain_minus_{ref}_ci_high": h_,
                                    f"gain_minus_{ref}_p": p_})
                pair_rows.append(row)
PR = pd.DataFrame(pair_rows)


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj, run = np.empty(len(p)), 0
    for r, j in enumerate(o):
        run = max(run, (len(p) - r) * p[j])
        adj[j] = min(1, run)
    return adj


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


prim = PR.group == "CHOL"
PR.loc[prim, "p_holm_primary8"] = holm(PR.loc[prim, "p_cluster_t"])
PR.loc[prim, "p_signflip_holm_primary8"] = holm(PR.loc[prim, "p_signflip"])
sec = []
for i, r in PR.iterrows():
    if r.group != "CHOL":
        sec.append((i, "p_cluster_t", r.p_cluster_t))
    if r.group != "OTHER":
        sec.append((i, "gain_minus_OTHER_p", r.gain_minus_OTHER_p))
        sec.append((i, "gain_minus_matched_p", r.gain_minus_matched_p))
qs = bh([s[2] for s in sec])
for (i, col, _), q in zip(sec, qs):
    PR.loc[i, f"{col}_BH_secondary"] = q
PR.to_csv(RES / "paired_nonlinear_vs_linear.csv", index=False)

# ---------------- (C) dominance controls for CHOL: leave-one-gene-out / leave-one-module-out, top-gene share
dom = []
for t in TARGETS:
    for nl in NL:
        for li in LIN:
            mk = GROUPS["CHOL"]
            a, b, z = W[(nl, t)].mse.values[mk], W[(li, t)].mse.values[mk], W[(li, t)].mse0.values[mk]
            gn, mm = genes[mk], modv[mk]
            full = (b.sum() - a.sum()) / z.sum()
            gain = b - a
            loo_g = [((np.delete(b, i).sum() - np.delete(a, i).sum()) / np.delete(z, i).sum()) for i in range(len(a))]
            loo_m = {m_: ((b[mm != m_].sum() - a[mm != m_].sum()) / z[mm != m_].sum()) for m_ in np.unique(mm)}
            srt = np.sort(np.abs(gain))[::-1]
            dom.append(dict(diet=t, nonlinear=nl, linear=li, dR2_full=full, loo_gene_min=min(loo_g), loo_gene_max=max(loo_g),
                            loo_gene_most_influential=gn[int(np.argmax(np.abs(np.array(loo_g) - full)))],
                            loo_gene_n_same_sign=int(np.sum(np.sign(loo_g) == np.sign(full))),
                            loo_module_min=min(loo_m.values()), loo_module_max=max(loo_m.values()),
                            loo_module_n_same_sign=int(sum(np.sign(v) == np.sign(full) for v in loo_m.values())),
                            n_modules=len(loo_m), share_abs_gain_top3_genes=float(srt[:3].sum() / srt.sum()),
                            genes_by_module=json.dumps(pd.Series(mm).value_counts().sort_index().to_dict(), default=int)))
pd.DataFrame(dom).to_csv(RES / "chol_dominance_controls.csv", index=False)

# per-gene CHOL table
pg = []
for t in TARGETS:
    for gname in genes[GROUPS["CHOL"]]:
        r = dict(diet=t, gene=gname, module=int(W[("F8-OLS", t)].loc[gname, "module"]), y=W[("F8-OLS", t)].loc[gname, "y"],
                 mse_null=W[("F8-OLS", t)].loc[gname, "mse0"])
        for mid in MODELS.values():
            r[f"pred_{mid}"] = W[(mid, t)].loc[gname, "pred"]
            r[f"mse_{mid}"] = W[(mid, t)].loc[gname, "mse"]
            r[f"dir_{mid}"] = W[(mid, t)].loc[gname, "dir"]
        pg.append(r)
pd.DataFrame(pg).to_csv(RES / "chol_per_gene_predictions.csv", index=False)

info = dict(seed=SEED, B=B, k_matched=K, max_abs_diff_F8OLS_vs_lipid_followup=max_diff_ols_vs_lp,
            reproduction_max_absdiff=dict(rmse=float(C.absdiff_rmse.max()), r2=float(C.absdiff_r2.max()),
                                          sign=float(C.absdiff_sign.max())),
            group_sizes={k: int(v.sum()) for k, v in GROUPS.items()},
            primary_family="CHOL x {CatBoost,HGB} x {OLS,Ridge} x {CD,MD}; module-cluster t; Holm(8)",
            followup_trigger=bool(((PR.p_holm_primary8 < .05) & (PR.mean_gain_mse > 0)).any()))
json.dump(info, open(RES / "run_info_nl01.json", "w"), indent=2)
pd.set_option("display.width", 250)
print(json.dumps(info, indent=1))
print(GM[GM.group == "CHOL"].round(4).to_string(index=False))
print(PR[PR.group == "CHOL"].round(4).to_string(index=False))
print(pd.DataFrame(dom).round(4).to_string(index=False))
