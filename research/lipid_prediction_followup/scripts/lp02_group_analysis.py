"""Lipid-gene prediction follow-up, step 2: held-out prediction quality for lipid-metabolism gene sets vs other genes.
Reads results/oof_predictions_with_null.csv.gz (lp01). No model is fitted, tuned or selected here.

Models (fixed before any group metric was computed; designations from earlier phases, not from these results):
  F8-OLS     Phase 4 primary gene model
  F8+C-Ridge Phase 4B recommended / Phase 5 frozen secondary gene model
  FULL-OLS   Phase 4B primary FULL comparison row (28 features; includes block-B pathway-membership indicators
             for exactly LIPOGENESIS_CORE, CHOL_SYNTHESIS_CORE and SANDERSON_PPARA_TARGETS)
  secondary estimators also reported: F8+C-OLS, FULL-Ridge.
Gene sets: research/roadmap_phase3/results/gene_sets.csv (the frozen sets used for block B).
Per gene: squared error averaged over the 3 outer repeats (each gene is held out exactly once per repeat).
Null baseline: weighted training-gene mean of the split (Phase 4/4B definition).
Group metrics: RMSE_model, RMSE_null, dMSE = MSE_null - MSE_model, R2_group = 1 - SSE_model/SSE_null.
Uncertainty / dependence (20 co-expression modules = outer-split grouping unit):
  (a) module block bootstrap, B = 2000 (resample modules with replacement; genes travel with their module);
  (b) within-module label permutation, 2000 permutations (group labels shuffled within each module);
  (c) per-gene dMSE ~ group + module fixed effects, CR1 module-cluster-robust SE, t with G-1 df.
Overlap: exclusive-membership categories and a joint regression with all three indicators + module FE.
Run: .venv/Scripts/python.exe research/lipid_prediction_followup/scripts/lp02_group_analysis.py
"""
import json
from pathlib import Path

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

RS = Path(__file__).resolve().parents[2]
OUT = RS / "lipid_prediction_followup"
RES, FIG = OUT / "results", OUT / "figures"
SEED, B, NPERM = 20261125, 2000, 2000
PRIMARY = ["F8-OLS", "F8+C-Ridge", "FULL-OLS"]
SECONDARY = ["F8+C-OLS", "FULL-Ridge"]
SETS = {"LIPO": "LIPOGENESIS_CORE", "CHOL": "CHOL_SYNTHESIS_CORE", "PPARA": "SANDERSON_PPARA_TARGETS"}
B_SETS = ["LIPOGENESIS_CORE", "CHOL_SYNTHESIS_CORE", "SANDERSON_PPARA_TARGETS", "HALLMARK_FATTY_ACID_METABOLISM",
          "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
          "HALLMARK_XENOBIOTIC_METABOLISM", "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_PEROXISOME",
          "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_MTORC1_SIGNALING"]
PAIRS = [("F8+C-Ridge", "F8-OLS"), ("FULL-OLS", "F8-OLS"), ("FULL-OLS", "F8+C-Ridge")]

P = pd.read_csv(RES / "oof_predictions_with_null.csv.gz")
P["se"], P["se0"] = (P.y - P.pred) ** 2, (P.y - P.null) ** 2
G = P.groupby(["model_id", "target", "gene"]).agg(y=("y", "first"), module=("module", "first"), mse=("se", "mean"),
                                                   mse0=("se0", "mean"), pred=("pred", "mean")).reset_index()
G["d"] = G.mse0 - G.mse
gs = pd.read_csv(RS / "roadmap_phase3/results/gene_sets.csv")
genes = G.gene.unique()
mem = pd.DataFrame({k: np.isin(genes, gs.gene[gs.set == v]) for k, v in SETS.items()}, index=genes)
mem["any_B12"] = np.isin(genes, gs.gene[gs.set.isin(B_SETS)])
nset = mem[list(SETS)].sum(1)
GROUPS = {"LIPO": mem.LIPO, "CHOL": mem.CHOL, "PPARA": mem.PPARA, "LIPID_any": nset > 0,
          "LIPO_only": mem.LIPO & (nset == 1), "CHOL_only": mem.CHOL & (nset == 1), "PPARA_only": mem.PPARA & (nset == 1),
          "LIPID_multi": nset >= 2, "OTHER": nset == 0, "OTHER_noB12": ~mem.any_B12}
GROUPS = {k: v.values for k, v in GROUPS.items()}
mods = np.array(sorted(G.module.unique()))
rng = np.random.default_rng(SEED)
BOOTW = rng.multinomial(len(mods), np.ones(len(mods)) / len(mods), size=B)  # B x 20 module counts

# ---- overlap table
ov = []
for a in SETS:
    for b in SETS:
        ov.append(dict(set_a=a, set_b=b, n_shared=int((mem[a] & mem[b]).sum())))
ovT = pd.DataFrame(ov).pivot(index="set_a", columns="set_b", values="n_shared")
ovT.to_csv(RES / "gene_set_overlap_in_prediction_universe.csv")


def eff_modules(m):
    c = pd.Series(m).value_counts().values
    p = c / c.sum()
    return float(1 / (p**2).sum())


def module_sums(df, mask):
    """Per-module sums for a gene mask: n, SSE, SSE0, sum d."""
    x = df[mask]
    s = x.groupby("module").agg(n=("d", "size"), sse=("mse", "sum"), sse0=("mse0", "sum"), sd=("d", "sum")).reindex(mods, fill_value=0)
    return s.values  # 20 x 4


def stats_from_sums(S):  # S: (..., 4) summed
    n, sse, sse0, sd = S[..., 0], S[..., 1], S[..., 2], S[..., 3]
    with np.errstate(invalid="ignore", divide="ignore"):
        return dict(r2=1 - sse / sse0, dmse=sd / n, rmse=np.sqrt(sse / n), rmse0=np.sqrt(sse0 / n))


def boot(S):  # S: 20 x 4 -> B x 4
    return BOOTW @ S


def cr1_fe(yv, xv, module, extra=None):
    """OLS of yv on xv (+extra columns) with module fixed effects; CR1 cluster-robust SE by module."""
    D = pd.get_dummies(module).values.astype(float)
    X = np.column_stack([xv] + ([] if extra is None else [extra]) + [D])
    b, *_ = np.linalg.lstsq(X, yv, rcond=None)
    e = yv - X @ b
    XtXi = np.linalg.pinv(X.T @ X)
    Gc, n, k = len(np.unique(module)), len(yv), X.shape[1]
    meat = np.zeros((k, k))
    for g in np.unique(module):
        i = module == g
        u = X[i].T @ e[i]
        meat += np.outer(u, u)
    V = XtXi @ meat @ XtXi * (Gc / (Gc - 1)) * ((n - 1) / (n - k))
    q = 1 if extra is None else 1 + extra.shape[1]
    se = np.sqrt(np.diag(V)[:q])
    t = b[:q] / se
    return b[:q], se, 2 * stats.t.sf(np.abs(t), Gc - 1), stats.t.ppf(.975, Gc - 1)


modv_all = G[(G.model_id == PRIMARY[0]) & (G.target == "PE_in_CD")].set_index("gene").loc[genes].module.values
PERM_LABELS = {}
for ref in ("OTHER", "OTHER_noB12"):
    for g in ("LIPO", "CHOL", "PPARA", "LIPID_any", "LIPO_only", "CHOL_only", "PPARA_only", "LIPID_multi"):
        keep = GROUPS[g] | GROUPS[ref]
        mk, lab = modv_all[keep], GROUPS[g][keep]
        order = np.argsort(mk, kind="stable")
        Lm = np.empty((NPERM, keep.sum()), bool)
        for i in range(NPERM):
            pi = np.lexsort((rng.random(len(mk)), mk))
            Lm[i, order] = lab[pi]
        PERM_LABELS[(g, ref)] = Lm
rows, cmp_rows, joint_rows, pair_rows, rep_rows = [], [], [], [], []
for tgt in ("PE_in_CD", "PE_in_MD"):
    for mid in PRIMARY + SECONDARY:
        df = G[(G.model_id == mid) & (G.target == tgt)].set_index("gene").loc[genes].reset_index()
        modv = df.module.values
        sums = {g: module_sums(df, m) for g, m in GROUPS.items()}
        bst = {g: stats_from_sums(boot(S)) for g, S in sums.items()}
        for g, m in GROUPS.items():
            pt = stats_from_sums(sums[g].sum(0))
            x = df[m]
            rows.append(dict(diet=tgt, model=mid, primary=mid in PRIMARY, group=g, n_genes=int(m.sum()),
                             n_modules=int(x.module.nunique()), eff_n_modules=eff_modules(x.module), rmse_model=pt["rmse"],
                             rmse_null=pt["rmse0"], dMSE=pt["dmse"], dMSE_ci_low=np.nanquantile(bst[g]["dmse"], .025),
                             dMSE_ci_high=np.nanquantile(bst[g]["dmse"], .975), R2_group=pt["r2"],
                             R2_ci_low=np.nanquantile(bst[g]["r2"], .025), R2_ci_high=np.nanquantile(bst[g]["r2"], .975),
                             frac_genes_improved=float((x.d > 0).mean()), spearman_y_pred=stats.spearmanr(x.y, x.pred).statistic,
                             mean_abs_y=float(x.y.abs().mean())))
        # per-repeat spread of R2_group
        pr = P[(P.model_id == mid) & (P.target == tgt)]
        for r_, q in pr.groupby("repeat"):
            q = q.set_index("gene").loc[genes]
            for g in ("LIPO", "CHOL", "PPARA", "LIPID_any", "OTHER"):
                qq = q[GROUPS[g]]
                rep_rows.append(dict(diet=tgt, model=mid, repeat=r_, group=g, R2_group=1 - qq.se.sum() / qq.se0.sum(),
                                     dMSE=float((qq.se0 - qq.se).mean())))
        # group vs OTHER (and vs OTHER_noB12)
        for ref in ("OTHER", "OTHER_noB12"):
            for g in ("LIPO", "CHOL", "PPARA", "LIPID_any", "LIPO_only", "CHOL_only", "PPARA_only", "LIPID_multi"):
                a, o = stats_from_sums(sums[g].sum(0)), stats_from_sums(sums[ref].sum(0))
                dr2, ddm = a["r2"] - o["r2"], a["dmse"] - o["dmse"]
                bd_r2 = bst[g]["r2"] - bst[ref]["r2"]
                bd_dm = bst[g]["dmse"] - bst[ref]["dmse"]
                keep = GROUPS[g] | GROUPS[ref]
                L = PERM_LABELS[(g, ref)].astype(float)  # NPERM x n_keep, labels permuted within module among kept genes
                Lc = 1.0 - L
                mse_, mse0_, d_ = df.mse.values[keep], df.mse0.values[keep], df.d.values[keep]
                n1 = L.sum(1)
                nr2 = (1 - (L @ mse_) / (L @ mse0_)) - (1 - (Lc @ mse_) / (Lc @ mse0_))
                ndm = (L @ d_) / n1 - (Lc @ d_) / (len(d_) - n1)
                sub = df[keep]
                b_, se_, p_, qt = cr1_fe(sub.d.values, GROUPS[g][keep].astype(float), sub.module.values)
                cmp_rows.append(dict(diet=tgt, model=mid, primary=mid in PRIMARY, group=g, reference=ref,
                                     n_group=int(GROUPS[g].sum()), n_ref=int(GROUPS[ref].sum()),
                                     n_modules_group=int(df[GROUPS[g]].module.nunique()),
                                     diff_R2=dr2, diff_R2_ci_low=np.nanquantile(bd_r2, .025), diff_R2_ci_high=np.nanquantile(bd_r2, .975),
                                     diff_R2_perm_p=float((np.sum(np.abs(nr2) >= abs(dr2) - 1e-15) + 1) / (NPERM + 1)),
                                     diff_dMSE=ddm, diff_dMSE_ci_low=np.nanquantile(bd_dm, .025), diff_dMSE_ci_high=np.nanquantile(bd_dm, .975),
                                     diff_dMSE_perm_p=float((np.sum(np.abs(ndm) >= abs(ddm) - 1e-15) + 1) / (NPERM + 1)),
                                     moduleFE_diff_dMSE=b_[0], moduleFE_ci_low=b_[0] - qt * se_[0], moduleFE_ci_high=b_[0] + qt * se_[0],
                                     moduleFE_p_CR1=p_[0]))
        # joint regression: three indicators + module FE, among LIPID_any U OTHER
        keep = GROUPS["LIPID_any"] | GROUPS["OTHER"]
        sub = df[keep]
        Xi = np.column_stack([GROUPS[k][keep].astype(float) for k in SETS])
        b_, se_, p_, qt = cr1_fe(sub.d.values, Xi[:, 0], sub.module.values, extra=Xi[:, 1:])
        for j, k in enumerate(SETS):
            joint_rows.append(dict(diet=tgt, model=mid, term=k, coef_dMSE=b_[j], ci_low=b_[j] - qt * se_[j],
                                   ci_high=b_[j] + qt * se_[j], p_CR1=p_[j], n_modules=len(mods)))

    # model-vs-model difference-in-difference (gene-paired): improvement of A over B, lipid minus other
    for a_, b0 in PAIRS:
        da = G[(G.model_id == a_) & (G.target == tgt)].set_index("gene").loc[genes]
        db = G[(G.model_id == b0) & (G.target == tgt)].set_index("gene").loc[genes]
        gain = (db.mse - da.mse).values  # >0: A has lower held-out error than B
        dd = pd.DataFrame({"module": da.module.values, "gain": gain})
        for g in ("LIPO", "CHOL", "PPARA", "LIPID_any"):
            for ref in ("OTHER",):
                S1 = dd[GROUPS[g]].groupby("module").gain.agg(["size", "sum"]).reindex(mods, fill_value=0).values
                S0 = dd[GROUPS[ref]].groupby("module").gain.agg(["size", "sum"]).reindex(mods, fill_value=0).values
                pt1, pt0 = S1[:, 1].sum() / S1[:, 0].sum(), S0[:, 1].sum() / S0[:, 0].sum()
                b1, b0_ = BOOTW @ S1, BOOTW @ S0
                with np.errstate(invalid="ignore", divide="ignore"):
                    bdiff = b1[:, 1] / b1[:, 0] - b0_[:, 1] / b0_[:, 0]
                    bg1 = b1[:, 1] / b1[:, 0]
                keep = GROUPS[g] | GROUPS[ref]
                bb, ss, pp, qt = cr1_fe(gain[keep], GROUPS[g][keep].astype(float), da.module.values[keep])
                pair_rows.append(dict(diet=tgt, model_A=a_, model_B=b0, group=g, reference=ref,
                                      mean_gain_group=pt1, gain_group_ci_low=np.nanquantile(bg1, .025), gain_group_ci_high=np.nanquantile(bg1, .975),
                                      mean_gain_ref=pt0, diff_gain=pt1 - pt0, diff_ci_low=np.nanquantile(bdiff, .025),
                                      diff_ci_high=np.nanquantile(bdiff, .975), moduleFE_diff=bb[0], moduleFE_ci_low=bb[0] - qt * ss[0],
                                      moduleFE_ci_high=bb[0] + qt * ss[0], moduleFE_p_CR1=pp[0]))

R = pd.DataFrame(rows)
CMP = pd.DataFrame(cmp_rows)
R.to_csv(RES / "group_metrics_all.csv", index=False)
CMP.to_csv(RES / "group_vs_reference_comparisons.csv", index=False)
pd.DataFrame(joint_rows).to_csv(RES / "joint_overlap_regression_moduleFE.csv", index=False)
pd.DataFrame(pair_rows).to_csv(RES / "model_pair_gain_lipid_vs_other.csv", index=False)
pd.DataFrame(rep_rows).to_csv(RES / "group_metrics_per_repeat.csv", index=False)

# ---- concise table: primary models, main groups, vs OTHER
main = ["LIPO", "CHOL", "PPARA", "LIPID_any", "OTHER"]
T = R[R.primary & R.group.isin(main)].copy()
c = CMP[(CMP.reference == "OTHER") & CMP.primary][["diet", "model", "group", "diff_R2", "diff_R2_ci_low", "diff_R2_ci_high",
                                                    "diff_R2_perm_p", "diff_dMSE", "diff_dMSE_ci_low", "diff_dMSE_ci_high",
                                                    "diff_dMSE_perm_p", "moduleFE_diff_dMSE", "moduleFE_ci_low",
                                                    "moduleFE_ci_high", "moduleFE_p_CR1"]]
T = T.merge(c, on=["diet", "model", "group"], how="left")
T["model"] = pd.Categorical(T.model, PRIMARY)
T["group"] = pd.Categorical(T.group, main)
T = T.sort_values(["diet", "model", "group"])
keepcols = ["diet", "model", "group", "n_genes", "n_modules", "eff_n_modules", "rmse_model", "rmse_null", "R2_group", "R2_ci_low",
            "R2_ci_high", "dMSE", "dMSE_ci_low", "dMSE_ci_high", "frac_genes_improved", "diff_R2", "diff_R2_ci_low",
            "diff_R2_ci_high", "diff_R2_perm_p", "diff_dMSE", "diff_dMSE_ci_low", "diff_dMSE_ci_high", "diff_dMSE_perm_p",
            "moduleFE_diff_dMSE", "moduleFE_ci_low", "moduleFE_ci_high", "moduleFE_p_CR1"]
T[keepcols].to_csv(OUT / "TABLE_results_concise.csv", index=False)

# figure: see lp03_matched_and_figure.py
json.dump({"seed": SEED, "B_bootstrap": B, "n_permutations": NPERM, "n_genes": int(len(genes)), "n_modules": int(len(mods)),
           "primary_models": PRIMARY, "secondary_models": SECONDARY,
           "group_sizes": {k: int(v.sum()) for k, v in GROUPS.items()}}, open(RES / "run_info_lp02.json", "w"), indent=2)
print(T[keepcols].round(4).to_string(index=False))
print(ovT)
