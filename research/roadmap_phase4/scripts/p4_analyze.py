"""Aggregate task A/B results: rank tables, Nadeau-Bengio paired intervals, robustness counts, controls, figures.
Decision rules exactly as evaluation_protocol.md section 6 (with deviation D1: 5-split variants/LOO)."""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from p4_data import P4, RES  # noqa: E402
from p4_models import BASELINES, COMPLEXITY, TUNED  # noqa: E402
from scipy import stats  # noqa: E402

A = RES / "taskA"
FIG = P4 / "figures"
FIG.mkdir(exist_ok=True)
SINGLE = ["ols_MASH", "ols_FAST", "ols_PPARA", "ols_CIRC"]
ALL = BASELINES + TUNED


def nb_interval(d, n_test, n_train):
    """Nadeau-Bengio corrected resampled t: var = (1/J + n_test/n_train) * s^2."""
    J = len(d)
    m, s2 = d.mean(), d.var(ddof=1)
    se = np.sqrt((1 / J + np.mean(n_test) / np.mean(n_train)) * s2)
    t = stats.t.ppf(0.975, J - 1)
    p = 2 * stats.t.sf(abs(m / se), J - 1) if se > 0 else np.nan
    return m, m - t * se, m + t * se, p


def paired(m, run, model, comp, metric="rmse"):
    a = m[(m.run == run) & (m.model == model)].set_index(["repeat", "fold"]).sort_index()
    b = m[(m.run == run) & (m.model == comp)].set_index(["repeat", "fold"]).sort_index()
    d = (a[metric] - b[metric]).values
    mean, lo, hi, p = nb_interval(d, a.n_test.values, a.n_train.values)
    better = (hi < 0) and ((d < 0).sum() >= int(np.ceil(0.8 * len(d))))
    return dict(run=run, model=model, comparator=comp, metric=metric, mean_diff=mean, ci_low=lo, ci_high=hi, p_nb=p,
                n_splits=len(d), n_splits_better=int((d < 0).sum()), better_by_rule=bool(better))


def main():
    m1 = pd.read_csv(A / "metrics_primary_positive.csv")
    m2 = pd.read_csv(A / "metrics_neg_variants_loo.csv")
    # the model name "null" is parsed as NA by pandas; restore it
    for d in (m1, m2):
        d["model"] = d["model"].fillna("null")
    m = pd.concat([m1, m2], ignore_index=True)
    keys = ["rmse", "mae", "r2_oos", "r2_heldout_mean", "wr2_oos", "spearman", "sign_agree", "sign_agree_deg"]
    summ = m.groupby(["run", "model"])[keys].agg(["mean", "std"])
    summ.columns = [f"{a}_{b}" for a, b in summ.columns]
    summ = summ.reset_index()
    summ["complexity"] = summ.model.map(COMPLEXITY)
    summ.to_csv(A / "summary_by_run_model.csv", index=False)

    out, ranks = {}, []
    comps = []
    for run in ["primary__PE_in_CD", "primary__PE_in_MD", "positive__Diet_in_Veh"]:
        s = summ[summ.run == run].set_index("model")
        best_single = s.loc[SINGLE, "rmse_mean"].idxmin()
        out[run] = {"best_single_axis": best_single}
        r = s[["rmse_mean", "rmse_std", "mae_mean", "r2_oos_mean", "r2_oos_std", "spearman_mean", "sign_agree_deg_mean"]].copy()
        r["rank_rmse"] = r.rmse_mean.rank()
        r.insert(0, "run", run)
        ranks.append(r.reset_index())
        for model in ALL:
            for comp in ("null", best_single, "ols_all"):
                if model != comp:
                    comps.append(paired(m, run, model, comp))
    rank = pd.concat(ranks)
    rank.to_csv(A / "rank_table.csv", index=False)
    C = pd.DataFrame(comps)
    C.to_csv(A / "paired_comparisons.csv", index=False)

    # robustness: variants (5 splits) and LOO (5 splits): mean dRMSE vs ols_all and vs null
    var = summ[summ.run.str.match(r"V\d")].copy()
    var["target"] = var.run.str.split("__").str[1]
    rob = []
    for tgt in ("PE_in_CD", "PE_in_MD"):
        for model in ["ols_all", "ridge", "rf", "hgb"] + SINGLE:
            row = dict(target=tgt, model=model)
            v = var[(var.target == tgt)]
            vm = v[v.model == model].set_index("run").rmse_mean
            vo = v[v.model == "ols_all"].set_index("run").rmse_mean
            vn = v[v.model == "null"].set_index("run").rmse_mean
            row["variants_better_than_ols_all"] = int((vm < vo.reindex(vm.index)).sum()) if model != "ols_all" else None
            row["variants_better_than_null"] = int((vm < vn.reindex(vm.index)).sum())
            row["n_variants"] = int(len(vm))
            lo = summ[summ.run.str.startswith("loo_") & summ.run.str.endswith(tgt)]
            lm = lo[lo.model == model].set_index("run").rmse_mean
            lo_o = lo[lo.model == "ols_all"].set_index("run").rmse_mean
            lo_n = lo[lo.model == "null"].set_index("run").rmse_mean
            if len(lm):
                row["loo_better_than_ols_all"] = int((lm < lo_o.reindex(lm.index)).sum()) if model != "ols_all" else None
                row["loo_better_than_null"] = int((lm < lo_n.reindex(lm.index)).sum())
                row["n_loo"] = int(len(lm))
                row["loo_r2_oos_min"] = float(lo[lo.model == model].r2_oos_mean.min())
                row["loo_r2_oos_max"] = float(lo[lo.model == model].r2_oos_mean.max())
            rob.append(row)
    R = pd.DataFrame(rob)
    R.to_csv(A / "robustness_counts.csv", index=False)
    vt = var.pivot_table(index=["target", "run"], columns="model", values="r2_oos_mean")
    vt.to_csv(A / "variants_r2_oos.csv")

    # controls
    ctrl = summ[summ.run.str.startswith("neg_")].pivot_table(index="run", columns="model", values="r2_oos_mean")
    ctrl.to_csv(A / "negative_controls_r2_oos.csv")

    # coefficient & importance stability (primary runs)
    coef = pd.read_csv(A / "coefficients_primary_positive.csv")
    cs = coef.groupby(["run", "model", "feature"]).coef_per_SD.agg(["mean", "std", "min", "max"]).reset_index()
    cs["sign_stable"] = (np.sign(cs["min"]) == np.sign(cs["max"]))
    cs.to_csv(A / "coefficient_stability.csv", index=False)
    imp = pd.read_csv(A / "importance_primary_positive.csv")
    imp.groupby(["run", "model", "feature"]).perm_importance_mse.agg(["mean", "std"]).reset_index().to_csv(
        A / "importance_summary.csv", index=False)
    # chosen hyperparameters
    hp = m1[m1.model.isin(TUNED)].groupby(["run", "model"]).best_params.agg(lambda s: s.value_counts().to_dict())
    hp.reset_index().to_csv(A / "chosen_hyperparameters.csv", index=False)

    # ---------------- figures ----------------
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6), sharey=False)
    for i, run in enumerate(["primary__PE_in_CD", "primary__PE_in_MD", "positive__Diet_in_Veh"]):
        s = summ[summ.run == run].set_index("model").reindex(ALL)
        nb = [nb_interval((m[(m.run == run) & (m.model == mo)].sort_values(["repeat", "fold"]).r2_oos.values),
                          m[(m.run == run) & (m.model == mo)].n_test.values, m[(m.run == run) & (m.model == mo)].n_train.values)
              for mo in ALL]
        mean = np.array([x[0] for x in nb])
        lo = np.array([x[1] for x in nb])
        hi = np.array([x[2] for x in nb])
        col = ["grey"] + ["tab:blue"] * 4 + ["tab:green"] + ["tab:orange"] * 6
        ax[i].bar(range(len(ALL)), mean, color=col)
        ax[i].errorbar(range(len(ALL)), mean, yerr=[mean - lo, hi - mean], fmt="none", c="k", lw=.8)
        ax[i].axhline(0, c="k", lw=.6)
        ax[i].set_xticks(range(len(ALL)))
        ax[i].set_xticklabels(ALL, rotation=60, fontsize=7)
        ax[i].set_title(f"{run}: out-of-sample R2 vs training mean\n(15 grouped outer splits; NB-corrected 95% CI)", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "taskA_performance_r2.png", dpi=110)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
    for i, run in enumerate(["primary__PE_in_CD", "primary__PE_in_MD"]):
        s = summ[summ.run == run].set_index("model")
        bs = out[run]["best_single_axis"]
        for mo in ALL:
            ax[i].scatter(COMPLEXITY[mo] + np.random.default_rng(hash(mo) % 99).uniform(-.15, .15),
                          s.loc[mo, "r2_oos_mean"] - s.loc[bs, "r2_oos_mean"], s=30)
            ax[i].annotate(mo, (COMPLEXITY[mo], s.loc[mo, "r2_oos_mean"] - s.loc[bs, "r2_oos_mean"]), fontsize=7)
        ax[i].axhline(0, c="k", lw=.6)
        ax[i].set_xlabel("model complexity (0 null ... 5 boosting)")
        ax[i].set_ylabel(f"mean R2_oos gain vs best single axis ({bs})")
        ax[i].set_title(run, fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "taskA_complexity_vs_gain.png", dpi=110)

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.3))
    for i, tgt in enumerate(["PE_in_CD", "PE_in_MD"]):
        lo = summ[summ.run.str.startswith("loo_") & summ.run.str.endswith(tgt)]
        data = [lo[lo.model == mo].r2_oos_mean.values for mo in BASELINES + ["ridge", "hgb"]]
        ax[i].boxplot(data, tick_labels=BASELINES + ["ridge", "hgb"])
        ax[i].axhline(0, c="k", lw=.6)
        ax[i].tick_params(axis="x", rotation=45, labelsize=7)
        ax[i].set_title(f"{tgt}: R2_oos across 19 leave-one-mouse-out targets", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "taskA_loo_stability.png", dpi=110)

    fig, ax = plt.subplots(figsize=(9, 6))
    vt2 = vt[["null", "ols_MASH", "ols_FAST", "ols_PPARA", "ols_CIRC", "ols_all", "ridge", "rf", "hgb"]]
    im = ax.imshow(vt2.values, cmap="RdBu_r", vmin=-0.15, vmax=0.15, aspect="auto")
    ax.set_yticks(range(len(vt2)))
    ax.set_yticklabels([f"{a} | {b.split('__')[0]}" for a, b in vt2.index], fontsize=7)
    ax.set_xticks(range(vt2.shape[1]))
    ax.set_xticklabels(vt2.columns, rotation=45, fontsize=8)
    for i in range(vt2.shape[0]):
        for j2 in range(vt2.shape[1]):
            ax.text(j2, i, f"{vt2.values[i, j2]:.3f}", ha="center", va="center", fontsize=6)
    fig.colorbar(im, label="mean R2_oos (5 grouped splits)")
    ax.set_title("Sensitivity variants", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "taskA_variant_stability.png", dpi=110)

    json.dump(out, open(A / "best_single_axis.json", "w"), indent=2)
    print(rank.round(4).to_string())
    print(C[C.comparator.isin(["null"]) | C.model.isin(["ols_all", "ridge", "rf", "hgb", "catboost", "elasticnet", "tree"])]
          .round(4).to_string())
    print(R.to_string())
    print(ctrl.round(3).to_string())
    print(cs[cs.model == "ols_all"].round(3).to_string())


if __name__ == "__main__":
    main()
