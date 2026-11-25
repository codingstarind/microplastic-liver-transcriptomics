"""Tasks 2-4: prediction-error diagnostics, target quality/identifiability, leakage/dependence audit.
Inputs: Phase 4C fix-experiment OOF predictions (F8-OLS, F8+C-Ridge, F8+B+C-Ridge), Phase 4 saved OOF predictions
(null, CatBoost-F8), GSE294180 DESeq2 results + 19 LOO refits, Phase 4 depth-thinning targets, Phase 3 gene lengths."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from c00_common import (  # noqa: E402
    B_COLS,
    REPO,
    RS,
    F,
    T,
    _csv,
    block_boot,
    dl_tau2,
    gene_groups,
    loo_names,
    raw_table,
    target,
)
from scipy.stats import spearmanr  # noqa: E402
from sklearn.linear_model import LinearRegression  # noqa: E402
from sklearn.model_selection import KFold  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

TG = ["PE_in_CD", "PE_in_MD"]
grp = gene_groups()
Xraw = raw_table()
res = {}

# ---------------- OOF predictions (mean over 3 repeats per gene) ----------------
fx = pd.read_csv(T / "fix_experiment_oof_predictions.csv.gz")
fx["label"] = fx.rep.replace({"FULL-A": "F8+B+C"}) + "_" + fx.model
p4 = pd.concat([_csv(RS / "roadmap_phase4" / "results" / "taskA" / "predictions" / f"primary__{t}.csv.gz").assign(target=t) for t in TG])
p4["model"] = p4["model"].fillna("null")
p4 = p4[p4.model.isin(["null", "catboost"])].assign(label=lambda d: d.model.replace({"catboost": "F8_catboost"}))
oof = pd.concat([fx[["gene", "y", "pred", "target", "label", "repeat"]], p4[["gene", "y", "pred", "target", "label", "repeat"]]])
O = oof.groupby(["target", "label", "gene"]).agg(y=("y", "first"), pred=("pred", "mean")).reset_index()
info = {t: target(t).join(_csv(REPO / "results/deseq2/provisional_deposited" / f"results_{t}.csv.gz").set_index("gene")[["baseMean", "padj"]]) for t in TG}
labs = ["F8_ols", "F8+C_ridge", "F8+B+C_ridge", "F8_catboost"]
rows = []
for t in TG:
    nul = O[(O.target == t) & (O.label == "null")].set_index("gene").pred
    for lab in labs:
        d = O[(O.target == t) & (O.label == lab)].set_index("gene")
        d = d.join(info[t][["s", "baseMean"]]).join(grp.module).assign(null=nul)
        d["resid"] = d.y - d.pred
        slope = LinearRegression().fit(d[["pred"]], d.y).coef_[0]
        cal = block_boot(d.reset_index(), lambda x: LinearRegression().fit(x[["pred"]], x.y).coef_[0], n=300)
        hs = block_boot(d.reset_index(), lambda x: spearmanr(x.resid.abs(), x.s).statistic, n=300)
        hb = block_boot(d.reset_index(), lambda x: spearmanr(x.resid.abs(), x.baseMean).statistic, n=300)
        rows.append(dict(target=t, model=lab, n=len(d), obs_sd=d.y.std(), pred_sd=d.pred.std(),
                         pred_to_obs_sd=d.pred.std() / d.y.std(), calibration_slope=slope, calib_ci=(cal[1], cal[2]),
                         rho_absresid_targetSE=hs[0], rho_absresid_SE_ci=(hs[1], hs[2]),
                         rho_absresid_baseMean=hb[0], rho_absresid_baseMean_ci=(hb[1], hb[2]),
                         r2_oos_gene_mean=1 - ((d.y - d.pred) ** 2).sum() / ((d.y - d.null) ** 2).sum()))
E = pd.DataFrame(rows)
E.to_csv(T / "error_overview.csv", index=False)

# subgroup errors (F8_ols primary, F8+B+C_ridge for comparison)
sub = []
gs = Xraw[B_COLS]
for t in TG:
    nul = O[(O.target == t) & (O.label == "null")].set_index("gene").pred
    for lab in ["F8_ols", "F8+B+C_ridge"]:
        d = O[(O.target == t) & (O.label == lab)].set_index("gene").join(info[t][["s", "baseMean", "padj"]]).join(grp.module)
        d["null"] = nul
        d["se_tercile"] = pd.qcut(d.s, 3, labels=["low_SE_reliable", "mid_SE", "high_SE_noisy"])
        d["basemean_quintile"] = pd.qcut(np.log10(d.baseMean + 1), 5, labels=[f"bm_Q{i}" for i in range(1, 6)])
        d["absobs_quintile_POSTHOC"] = pd.qcut(d.y.abs(), 5, labels=[f"absy_Q{i}" for i in range(1, 6)])
        d["pathway"] = "no_set"
        for c in B_COLS[::-1]:
            d.loc[gs[c].reindex(d.index) == 1, "pathway"] = c[2:]
        for gcol, kind in [("pathway", "predefined_pathway (first-listed set if several)"), ("se_tercile", "target_SE_tercile"),
                           ("basemean_quintile", "baseMean_quintile"), ("module", "co-expression module"),
                           ("absobs_quintile_POSTHOC", "DIAGNOSTIC-ONLY post-hoc |observed| bin")]:
            for g, dd in d.groupby(gcol, observed=True):
                sse, sse0 = ((dd.y - dd.pred) ** 2).sum(), ((dd.y - dd.null) ** 2).sum()
                sub.append(dict(target=t, model=lab, grouping=kind, group=str(g), n=len(dd),
                                rmse=float(np.sqrt(((dd.y - dd.pred) ** 2).mean())), r2_oos=float(1 - sse / sse0),
                                sign_acc=float((np.sign(dd.y) == np.sign(dd.pred)).mean()),
                                mean_resid=float(dd.resid.mean()) if "resid" in dd else float((dd.y - dd.pred).mean())))
S = pd.DataFrame(sub)
S.to_csv(T / "subgroup_errors.csv", index=False)

# ---------------- target quality ----------------
tq, infl = [], []
for t in TG:
    full = target(t)
    L = pd.DataFrame({d: target(f"LOO:{d}:{t}").y for d in loo_names()})
    g = full.index.intersection(L.index).intersection(Xraw.index)
    L, y, s = L.loc[g], full.loc[g, "y"], full.loc[g, "s"]
    n = L.shape[1]
    jk = np.sqrt((n - 1) / n * ((L.sub(L.mean(axis=1), axis=0)) ** 2).sum(axis=1))
    rel_deseq = 1 - (s**2).mean() / y.var()
    rel_jk = 1 - (jk**2).mean() / y.var()
    contrast_mice = [d for d in loo_names() if (t == "PE_in_CD" and d.startswith("CD")) or (t == "PE_in_MD" and d.startswith("MD"))]
    tq.append(dict(target=t, n_genes=len(g), var_log2fc=y.var(), mean_deseq_se2=(s**2).mean(), mean_jackknife_se2=(jk**2).mean(),
                   median_jk_over_deseq_se=float((jk / s).median()), reliability_deseq=rel_deseq, reliability_jackknife=rel_jk,
                   n_contrast_mice=len(contrast_mice),
                   note="jackknife over all 19 overlapping refits; only contrast mice change the estimate materially; approximate"))
    for d in loo_names():
        dlt = L[d] - y
        infl.append(dict(target=t, dropped_mouse=d, in_contrast=d in contrast_mice, rms_change=float(np.sqrt((dlt**2).mean())),
                         max_abs_change=float(dlt.abs().max()), spearman_with_full=float(spearmanr(L[d], y, nan_policy="omit").statistic), n_na=int(L[d].isna().sum())))
TQ = pd.DataFrame(tq)
TQ.to_csv(T / "target_quality.csv", index=False)
I = pd.DataFrame(infl)
I.to_csv(T / "single_mouse_influence_on_target.csv", index=False)
# LOO influence on the model fit (Phase 4 saved LOO metrics, ols_all; Phase 4B sensitivity LOO for F8+B+C ridge)
m4 = pd.read_csv(RS / "roadmap_phase4" / "results" / "taskA" / "metrics_neg_variants_loo.csv")
m4["model"] = m4["model"].fillna("null")
lo = m4[m4.run.str.startswith("loo_") & (m4.model == "ols_all")].groupby("run").r2_oos.mean().reset_index()
lo["dropped"] = lo.run.str.extract(r"loo_(.*)__")[0]
lo["target"] = lo.run.str.split("__").str[1]
lo.to_csv(T / "loo_model_r2_F8_ols.csv", index=False)

# confound probes (descriptive)
lengths = _csv(RS / "roadmap_phase3" / "results" / "offset_audit" / "offset_gene_table.csv.gz", index_col=0).gencode_median_tx_length
probe = []
for t, dep in [("PE_in_CD", "DEPTH_CD"), ("PE_in_MD", "DEPTH_MD"), ("Diet_in_Veh", None)]:
    y = target(t).y.reindex(Xraw.index)
    row = dict(target=t, rho_with_tx_length=float(spearmanr(y, lengths.reindex(Xraw.index), nan_policy="omit").statistic),
               rho_with_baseMean=float(spearmanr(y, info.get(t, target(t)).get("baseMean", pd.Series(dtype=float)).reindex(Xraw.index)
                                                 if t in info else _csv(REPO / "results/deseq2/provisional_deposited/results_Diet_in_Veh.csv.gz").set_index("gene").baseMean.reindex(Xraw.index),
                                                 nan_policy="omit").statistic))
    if dep:
        row["rho_with_depth_thinning_effect"] = float(spearmanr(y, target(dep).y.reindex(Xraw.index), nan_policy="omit").statistic)
    probe.append(row)
pd.DataFrame(probe).to_csv(T / "confound_probes.csv", index=False)

# ---------------- leakage / dependence ----------------
def icc(v, g):
    df = pd.DataFrame({"v": v, "g": g}).dropna()
    k = df.groupby("g").size()
    ms_b = (k * (df.groupby("g").v.mean() - df.v.mean()) ** 2).sum() / (len(k) - 1)
    ms_w = ((df.v - df.groupby("g").v.transform("mean")) ** 2).sum() / (len(df) - len(k))
    n0 = (len(df) - (k**2).sum() / len(df)) / (len(k) - 1)
    return float((ms_b - ms_w) / (ms_b + (n0 - 1) * ms_w))


dep = []
for t in TG:
    d = O[(O.target == t) & (O.label == "F8_ols")].set_index("gene").join(grp.module)
    dep.append(dict(target=t, icc_target_by_module=icc(d.y, d.module), icc_oof_resid_by_module=icc(d.y - d.pred, d.module)))
    # diagnostic-only random-gene CV (ignores modules) to quantify optimism
    j = Xraw.join(target(t), how="inner").dropna(subset=["y", "s"])
    X, y, s = j[[c for c in Xraw.columns if c.split("_")[0] in ("MASH", "FAST", "PPARA", "CIRC")][:8]].values, j.y.values, j.s.values
    r2 = []
    for tr, te in KFold(5, shuffle=True, random_state=101).split(X):
        w = 1 / (s**2 + dl_tau2(y[tr], s[tr]))
        sc = StandardScaler().fit(X[tr])
        mdl = LinearRegression().fit(sc.transform(X[tr]), y[tr], sample_weight=w[tr])
        p = mdl.predict(sc.transform(X[te]))
        null = np.average(y[tr], weights=w[tr])
        r2.append(1 - ((y[te] - p) ** 2).sum() / ((y[te] - null) ** 2).sum())
    dep[-1]["r2_random_gene_cv_DIAGNOSTIC"] = float(np.mean(r2))
    dep[-1]["r2_module_grouped_cv"] = float(pd.read_csv(T / "fix_experiment_metrics_by_fold.csv").query("target==@t and rep=='F8' and model=='ols'").r2_oos.mean())
pd.DataFrame(dep).to_csv(T / "dependence_audit.csv", index=False)

# ---------------- figures ----------------
fig, ax = plt.subplots(2, 4, figsize=(18, 8))
for i, t in enumerate(TG):
    d = O[(O.target == t) & (O.label == "F8_ols")].set_index("gene").join(info[t][["s", "baseMean"]])
    d["resid"] = d.y - d.pred
    ax[i, 0].scatter(d.pred, d.y, s=2, alpha=.2)
    lim = np.quantile(d.y.abs(), .995)
    ax[i, 0].plot([-lim, lim], [-lim, lim], "k--", lw=.6)
    ax[i, 0].set_xlim(-lim, lim)
    ax[i, 0].set_ylim(-lim, lim)
    ax[i, 0].set_xlabel("OOF predicted (F8 OLS)")
    ax[i, 0].set_ylabel(f"observed {t}")
    q = pd.qcut(d.pred, 10)
    b = d.groupby(q, observed=True).agg(p=("pred", "mean"), r=("resid", "mean"))
    ax[i, 1].scatter(d.pred, d.resid, s=2, alpha=.15)
    ax[i, 1].plot(b.p, b.r, "r-o", ms=3)
    ax[i, 1].axhline(0, c="k", lw=.6)
    ax[i, 1].set_xlabel("fitted")
    ax[i, 1].set_ylabel("residual (red: decile means)")
    ax[i, 2].scatter(d.s, d.resid, s=2, alpha=.15)
    ax[i, 2].set_xscale("log")
    ax[i, 2].set_xlabel("target lfcSE (log)")
    ax[i, 2].set_ylabel("residual")
    ax[i, 3].scatter(d.baseMean, d.resid, s=2, alpha=.15)
    ax[i, 3].set_xscale("log")
    ax[i, 3].set_xlabel("GSE294180 baseMean (log)")
    ax[i, 3].set_ylabel("residual")
    for a in ax[i]:
        a.set_title(t, fontsize=8)
fig.tight_layout()
fig.savefig(F / "fig2_residual_diagnostics_F8_ols.png", dpi=100)

fig, ax = plt.subplots(1, 3, figsize=(16, 4.3))
for i, t in enumerate(TG):
    s_ = S[(S.target == t) & (S.model == "F8_ols") & (S.grouping.isin(["target_SE_tercile", "baseMean_quintile"]))]
    ax[0].plot(s_.group, s_.r2_oos, "o-", label=t)
ax[0].axhline(0, c="k", lw=.6)
ax[0].tick_params(axis="x", rotation=45, labelsize=7)
ax[0].set_title("R2_oos by target-SE tercile and baseMean quintile (F8 OLS)", fontsize=8)
ax[0].legend(fontsize=7)
for t in TG:
    ii = I[I.target == t]
    ax[1].scatter(ii.dropped_mouse, ii.rms_change, label=t, s=20)
ax[1].tick_params(axis="x", rotation=90, labelsize=6)
ax[1].set_ylabel("RMS change of target log2FC when mouse dropped")
ax[1].legend(fontsize=7)
ax[1].set_title("Single-mouse influence on the PE target", fontsize=8)
for t in TG:
    ll = lo[lo.target == t]
    ax[2].scatter(ll.dropped, ll.r2_oos, label=t, s=20)
ax[2].tick_params(axis="x", rotation=90, labelsize=6)
ax[2].set_ylabel("grouped-CV R2_oos (F8 OLS, Phase 4 LOO runs)")
ax[2].legend(fontsize=7)
ax[2].set_title("Single-mouse influence on model fit", fontsize=8)
fig.tight_layout()
fig.savefig(F / "fig3_subgroups_and_influence.png", dpi=110)

with pd.option_context("display.width", 230):
    print(E.round(3).to_string(index=False))
    print(S[S.grouping.isin(["target_SE_tercile", "baseMean_quintile", "DIAGNOSTIC-ONLY post-hoc |observed| bin"]) & (S.model == "F8_ols")].round(3).to_string(index=False))
    print(S[(S.grouping.str.startswith("predefined")) & (S.model == "F8_ols")].round(3).to_string(index=False))
    print(TQ.round(4).to_string(index=False))
    print(I.sort_values("rms_change", ascending=False).head(8).round(4).to_string(index=False))
    print(pd.read_csv(T / "confound_probes.csv").round(3).to_string(index=False))
    print(pd.read_csv(T / "dependence_audit.csv").round(4).to_string(index=False))
    print(lo.round(4).sort_values("r2_oos").head(6).to_string(index=False))
