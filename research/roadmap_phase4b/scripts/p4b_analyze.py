"""Phase 4B analysis: comparisons vs F8 / F4 / OLS-F8, block ablations, controls, reproducibility, 3 figures.
Rules from comparison_protocol.md section 3."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

P4B = Path(__file__).resolve().parents[1]
R = P4B / "results"
FIG = P4B / "figures"
REPS = ["F4", "F8", "A4", "F8+A", "F8+B", "F8+C", "FULL", "FULL-A", "FULL-B", "FULL-C"]
MODELS = ["ols", "ridge", "elasticnet", "hgb", "catboost"]
TGTS = ["PE_in_CD", "PE_in_MD"]


def nb(d, ntest, ntrain):
    J = len(d)
    m, s2 = d.mean(), d.var(ddof=1)
    se = np.sqrt((1 / J + np.mean(ntest) / np.mean(ntrain)) * s2)
    t = stats.t.ppf(0.975, J - 1)
    return m, m - t * se, m + t * se


def paired(m, tgt, a, b):
    """a, b = (rep, model). Returns stats for metric(a) - metric(b)."""
    A = m[(m.target == tgt) & (m.rep == a[0]) & (m.model == a[1])].sort_values(["repeat", "fold"])
    B = m[(m.target == tgt) & (m.rep == b[0]) & (m.model == b[1])].sort_values(["repeat", "fold"])
    d = A.rmse.values - B.rmse.values
    dr2 = A.r2_oos.values - B.r2_oos.values
    mean, lo, hi = nb(d, A.n_test.values, A.n_train.values)
    rule = (hi < 0) and ((d < 0).sum() >= 12) and (mean <= -0.005) and (dr2.mean() >= 0.01)
    return dict(target=tgt, rep=a[0], model=a[1], vs_rep=b[0], vs_model=b[1], d_rmse=mean, ci_low=lo, ci_high=hi,
                n_better=int((d < 0).sum()), n_splits=len(d), d_r2_oos=float(dr2.mean()),
                d_mae=float((A.mae.values - B.mae.values).mean()), meets_rule=bool(rule),
                stat_only=bool(hi < 0 and (d < 0).sum() >= 12))


def main():
    m = pd.read_csv(R / "metrics_main.csv")
    s = m.groupby(["target", "rep", "model"])[["rmse", "mae", "r2_oos", "wr2_oos", "spearman"]].agg(["mean", "std"])
    s.columns = [f"{a}_{b}" for a, b in s.columns]
    s = s.reset_index()
    rows = []
    for tgt in TGTS:
        for rep in REPS:
            for mo in MODELS:
                if rep != "F8":
                    rows.append({**paired(m, tgt, (rep, mo), ("F8", mo)), "comparison": "vs_F8_same_model"})
                if rep != "F4":
                    rows.append({**paired(m, tgt, (rep, mo), ("F4", mo)), "comparison": "vs_F4_same_model"})
                if not (rep == "F8" and mo == "ols"):
                    rows.append({**paired(m, tgt, (rep, mo), ("F8", "ols")), "comparison": "vs_OLS_F8"})
    C = pd.DataFrame(rows)
    comp = s.merge(C[C.comparison == "vs_F8_same_model"][["target", "rep", "model", "d_rmse", "ci_low", "ci_high", "n_better",
                                                           "d_r2_oos", "meets_rule", "stat_only"]],
                   on=["target", "rep", "model"], how="left")
    comp2 = C[C.comparison == "vs_OLS_F8"][["target", "rep", "model", "d_rmse", "ci_low", "ci_high", "n_better", "meets_rule"]]
    comp2.columns = ["target", "rep", "model"] + [f"vsOLSF8_{c}" for c in comp2.columns[3:]]
    comp = comp.merge(comp2, on=["target", "rep", "model"], how="left")
    comp.to_csv(P4B / "model_comparison_results.csv", index=False)
    C.to_csv(R / "all_paired_comparisons.csv", index=False)

    ab = []
    for tgt in TGTS:
        for mo in MODELS:
            for blk in ("A", "B", "C"):
                p = paired(m, tgt, (f"FULL-{blk}", mo), ("FULL", mo))  # >0 means removing block hurts
                ab.append(dict(target=tgt, model=mo, block_removed=blk, d_rmse_removed_minus_full=p["d_rmse"],
                               ci_low=p["ci_low"], ci_high=p["ci_high"], n_splits_full_better=15 - p["n_better"],
                               d_r2_oos_removed_minus_full=p["d_r2_oos"], block_useful=bool(p["ci_low"] > 0)))
            for blk, rep in (("A", "F8+A"), ("B", "F8+B"), ("C", "F8+C")):
                p = paired(m, tgt, (rep, mo), ("F8", mo))
                ab.append(dict(target=tgt, model=mo, block_removed=f"add_{blk}_to_F8", d_rmse_removed_minus_full=-p["d_rmse"],
                               ci_low=-p["ci_high"], ci_high=-p["ci_low"], n_splits_full_better=p["n_better"],
                               d_r2_oos_removed_minus_full=-p["d_r2_oos"], block_useful=bool(p["ci_high"] < 0)))
    A = pd.DataFrame(ab)
    A.to_csv(P4B / "feature_ablation_results.csv", index=False)

    # reproducibility vs Phase 4
    p4 = pd.read_csv(P4B.parent / "roadmap_phase4" / "results" / "taskA" / "metrics_primary_positive.csv")
    p4["model"] = p4["model"].fillna("null").replace({"ols_all": "ols"})
    p4 = p4[p4.run.str.startswith("primary__")].assign(target=lambda d: d.run.str.split("__").str[1])
    j = m[m.rep == "F8"].merge(p4, on=["target", "model", "repeat", "fold"], suffixes=("", "_p4"))
    repro = dict(n=len(j), max_abs_rmse_diff=float((j.rmse - j.rmse_p4).abs().max()))
    # controls
    ctl = {}
    if (R / "metrics_controls.csv").exists():
        c = pd.read_csv(R / "metrics_controls.csv")
        ctl = c.groupby(["stage", "target", "model"]).r2_oos.mean().round(4).reset_index().to_dict("records")
        # compare shuffled-B with real F8+B (repeat 0 only)
        real = m[(m.rep == "F8+B") & (m.repeat == 0) & m.model.isin(["ridge", "hgb"])].groupby(["target", "model"]).r2_oos.mean()
        shuf = c[c.stage == "ctrl_shuffled_B"].groupby(["target", "model"]).r2_oos.mean()
        f8 = m[(m.rep == "F8") & (m.repeat == 0) & m.model.isin(["ridge", "hgb"])].groupby(["target", "model"]).r2_oos.mean()
        ctl_b = pd.DataFrame({"F8": f8, "F8+B_real": real, "F8+B_shuffled": shuf}).round(4)
        ctl_b.to_csv(R / "control_shuffled_B.csv")
    coef = pd.read_csv(R / "coefficients_main.csv")
    cs = coef[coef.model == "ols"].groupby(["target", "rep", "feature"]).coef_per_SD.agg(["mean", "std", "min", "max"]).reset_index()
    cs["sign_stable"] = np.sign(cs["min"]) == np.sign(cs["max"])
    cs.to_csv(R / "ols_coefficient_stability.csv", index=False)
    shortlist = comp[comp.meets_rule == True][["target", "rep", "model", "d_rmse", "d_r2_oos"]]  # noqa: E712
    summary = dict(reproducibility_vs_phase4=repro, shortlist=shortlist.to_dict("records"), controls=ctl)
    (R / "analysis_summary.json").write_text(json.dumps(summary, indent=2, default=float))

    # ---------------- figures (3) ----------------
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
    for i, tgt in enumerate(TGTS):
        p = s[s.target == tgt].pivot(index="rep", columns="model", values="r2_oos_mean").reindex(REPS)[MODELS]
        im = ax[i].imshow(p.values, cmap="viridis", aspect="auto")
        for a in range(p.shape[0]):
            for b in range(p.shape[1]):
                c = C[(C.target == tgt) & (C.rep == p.index[a]) & (C.model == MODELS[b]) & (C.comparison == "vs_F8_same_model")]
                mark = ""
                if len(c):
                    mark = "*" if c.meets_rule.iloc[0] else ("+" if c.stat_only.iloc[0] else "")
                ax[i].text(b, a, f"{p.values[a, b]:.3f}{mark}", ha="center", va="center", fontsize=7, color="w")
        ax[i].set_xticks(range(len(MODELS)))
        ax[i].set_xticklabels(MODELS)
        ax[i].set_yticks(range(len(REPS)))
        ax[i].set_yticklabels(REPS)
        ax[i].set_title(f"{tgt}: mean R2_oos (15 grouped splits)\n* meets locked rule vs F8 (same model); + statistical only", fontsize=8)
        fig.colorbar(im, ax=ax[i])
    fig.tight_layout()
    fig.savefig(FIG / "fig1_representation_x_model_r2.png", dpi=110)

    fig, ax = plt.subplots(1, 2, figsize=(14, 4.6))
    for i, tgt in enumerate(TGTS):
        d = A[(A.target == tgt) & A.block_removed.isin(["A", "B", "C"])]
        x = np.arange(len(MODELS))
        for k, blk in enumerate(["A", "B", "C"]):
            dd = d[d.block_removed == blk].set_index("model").reindex(MODELS)
            ax[i].bar(x + (k - 1) * .27, dd.d_rmse_removed_minus_full, .27, label=f"remove block {blk}")
            ax[i].errorbar(x + (k - 1) * .27, dd.d_rmse_removed_minus_full,
                           yerr=[dd.d_rmse_removed_minus_full - dd.ci_low, dd.ci_high - dd.d_rmse_removed_minus_full], fmt="none", c="k", lw=.7)
        ax[i].axhline(0, c="k", lw=.6)
        ax[i].set_xticks(x)
        ax[i].set_xticklabels(MODELS)
        ax[i].set_ylabel("RMSE(FULL minus block) - RMSE(FULL), log2FC\n(>0: block helps; NB 95% CI)")
        ax[i].set_title(f"{tgt}: feature-block ablation", fontsize=9)
        ax[i].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_feature_block_ablation.png", dpi=110)

    fig, ax = plt.subplots(1, 2, figsize=(14, 4.6))
    for i, tgt in enumerate(TGTS):
        data, labels = [], []
        for rep in ["F8+A", "F8+B", "F8+C", "FULL"]:
            for mo in ["ols", "hgb"]:
                A_ = m[(m.target == tgt) & (m.rep == rep) & (m.model == mo)].sort_values(["repeat", "fold"]).rmse.values
                B_ = m[(m.target == tgt) & (m.rep == "F8") & (m.model == mo)].sort_values(["repeat", "fold"]).rmse.values
                data.append(A_ - B_)
                labels.append(f"{rep}\n{mo}")
        ax[i].boxplot(data, tick_labels=labels)
        ax[i].axhline(0, c="k", lw=.6)
        ax[i].axhline(-0.005, c="r", ls="--", lw=.6)
        ax[i].tick_params(axis="x", labelsize=7)
        ax[i].set_ylabel("per-split RMSE difference vs F8 (same model)")
        ax[i].set_title(f"{tgt}: stability over 15 splits (red: practical threshold -0.005)", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "fig3_stability_per_split.png", dpi=110)
    with pd.option_context("display.width", 250):
        print(comp[["target", "rep", "model", "r2_oos_mean", "rmse_mean", "d_rmse", "ci_low", "ci_high", "n_better", "d_r2_oos",
                    "meets_rule", "vsOLSF8_d_rmse", "vsOLSF8_meets_rule"]].round(4).to_string())
        print(A.round(4).to_string())
        print(json.dumps(summary, indent=1, default=float)[:3000])


if __name__ == "__main__":
    main()
