"""Concise table + two-panel figure (CD, MD) for cholesterol-synthesis genes. Reads nl01 outputs only.
Run: .venv/Scripts/python.exe research/nonlinear_chol_followup/scripts/nl02_table_figure.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OUT = Path(__file__).resolve().parents[1]
RES, FIG = OUT / "results", OUT / "figures"
FIG.mkdir(exist_ok=True)
ORDER = ["F8-OLS", "F8-Ridge", "F8-CatBoost", "F8-HGB"]
GM = pd.read_csv(RES / "group_metrics_by_model.csv")
PR = pd.read_csv(RES / "paired_nonlinear_vs_linear.csv")

# ---- concise table (CHOL)
T = GM[GM.group == "CHOL"][["diet", "model", "n_genes", "n_modules", "R2_group", "R2_ci_low", "R2_ci_high", "rmse", "rmse_ci_low",
                            "rmse_ci_high", "rmse_null", "direction_agreement", "spearman", "matched_R2_group",
                            "matched_R2_ci_low", "matched_R2_ci_high"]].copy()
for li in ("F8-OLS", "F8-Ridge"):
    s = PR[(PR.group == "CHOL") & (PR.linear == li)][["diet", "nonlinear", "dR2", "dR2_ci_low", "dR2_ci_high", "mean_gain_mse",
                                                       "gain_ci_low", "gain_ci_high", "p_cluster_t", "p_holm_primary8",
                                                       "p_signflip", "n_genes_nonlinear_better", "gain_minus_matched",
                                                       "gain_minus_matched_ci_low", "gain_minus_matched_ci_high",
                                                       "gain_minus_matched_p_BH_secondary"]]
    tag = "vsOLS" if li == "F8-OLS" else "vsRidge"
    s = s.rename(columns={"nonlinear": "model", **{c: f"{tag}_{c}" for c in s.columns if c not in ("diet", "nonlinear")}})
    T = T.merge(s, on=["diet", "model"], how="left")
T["model"] = pd.Categorical(T.model, ORDER)
T = T.sort_values(["diet", "model"])
T.to_csv(OUT / "TABLE_results_concise.csv", index=False)

# ---- figure
COL = {"F8-OLS": "#2a6ebb", "F8-Ridge": "#7fa7d9", "F8-CatBoost": "#d0782a", "F8-HGB": "#9c4f12"}
fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), sharey=True)
x = np.arange(len(ORDER))
for ax, tgt, ttl in ((axes[0], "PE_in_CD", "PE vs Veh within CD"), (axes[1], "PE_in_MD", "PE vs Veh within MD")):
    d = GM[(GM.diet == tgt) & (GM.group == "CHOL")].set_index("model").loc[ORDER]
    ax.axhline(0, color="#4a4a4a", lw=.8)
    ax.axvspan(1.5, 3.5, color="#f3efe9", zorder=0, lw=0)
    for i, m in enumerate(ORDER):
        ax.vlines(i, d.loc[m, "R2_ci_low"], d.loc[m, "R2_ci_high"], color=COL[m], lw=2)
        ax.scatter(i, d.loc[m, "R2_group"], s=48, color=COL[m], zorder=3, label="Cholesterol-synthesis genes (21)" if i == 0 else None)
        ax.vlines(i + .22, d.loc[m, "matched_R2_ci_low"], d.loc[m, "matched_R2_ci_high"], color="#9a9a9a", lw=1.2)
        ax.scatter(i + .22, d.loc[m, "matched_R2_group"], s=26, color="white", edgecolor="#6b6b6b", lw=1.2, zorder=3,
                   label="Matched other genes (210; same null error)" if i == 0 else None)
        ax.text(i, -0.17, f"R² {d.loc[m, 'R2_group']:.3f}\nRMSE {d.loc[m, 'rmse']:.3f}\ndir {d.loc[m, 'direction_agreement']:.2f}",
                ha="center", va="top", fontsize=6.6, color="#3a3a3a")
    for i, m in ((2, "F8-CatBoost"), (3, "F8-HGB")):
        r = PR[(PR.diet == tgt) & (PR.group == "CHOL") & (PR.nonlinear == m) & (PR.linear == "F8-OLS")].iloc[0]
        ax.text(i, 1.02, f"ΔR² vs OLS {r.dR2:+.3f}\n[{r.dR2_ci_low:+.3f}, {r.dR2_ci_high:+.3f}]\nHolm p {r.p_holm_primary8:.3g}",
                ha="center", va="bottom", fontsize=6.6, color="#3a3a3a")
    ax.set_xticks(x)
    ax.set_xticklabels(["F8-OLS\n(linear)", "F8-Ridge\n(linear)", "F8-CatBoost\n(nonlinear)", "F8-HGB\n(nonlinear)"], fontsize=7.5)
    ax.set_xlim(-.5, 3.6)
    ax.set_ylim(-.42, 1.22)
    ax.set_title(ttl, fontsize=9.5, loc="left")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.tick_params(axis="y", labelsize=7.5)
axes[0].set_ylabel("Held-out R² within gene group\n1 − SSE(model) / SSE(null: training mean)", fontsize=8)
h, l_ = axes[0].get_legend_handles_labels()
fig.legend(h, l_, loc="upper left", bbox_to_anchor=(.01, .93), ncol=2, fontsize=7.2, frameon=False)
fig.suptitle("Original Phase 4 held-out predictions (F8 reference features; 20 co-expression modules, 5 folds × 3 repeats). "
             "Bars: 95% module-block bootstrap CI (B = 2000). dir = sign agreement. Holm over 8 primary tests.",
             fontsize=7.6, x=.01, ha="left")
fig.tight_layout(rect=(0, 0, 1, .89))
for e in ("png", "pdf", "svg"):
    fig.savefig(FIG / f"fig_chol_linear_vs_nonlinear.{e}", dpi=300)
print(T.round(4).to_string(index=False))
