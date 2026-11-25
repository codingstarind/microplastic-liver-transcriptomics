"""Figures 2-4 for Phase 5B (fig1 is written by b03)."""
import joblib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from b00_common import FIG, P5, RES, RS, TAB  # noqa: E402
import c00_common  # noqa: E402,F401  (adds Phase 4B paths for unpickling RepBuilder)

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
                     "savefig.facecolor": "#fcfcfb"})
INK, GREY = "#0b0b0b", "#8a8984"
C3 = {"GSE294180": "#8a8984", "HONG": "#2a78d6", "GSE245069": "#eb6834"}

# ---- fig2: domain shift + uncertainty ----
fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
s4d = pd.read_csv(RS / "roadmap_phase4d/results/per_mouse_scores.csv", index_col=0)
s5 = pd.read_csv(P5 / "per_sample_scores.csv", index_col=0)
groups = [("GSE294180 CD+Veh", s4d[s4d.group == "CD+Veh"].FAST, "GSE294180"), ("GSE294180 CD+PE", s4d[s4d.group == "CD+PE"].FAST, "GSE294180"),
          ("Hong CON", s5[(s5.study == "HONG") & (s5.group == "CON")].FAST_singscore, "HONG"), ("Hong PE", s5[(s5.study == "HONG") & (s5.group == "PE")].FAST_singscore, "HONG"),
          ("GSE245069 CON", s5[(s5.study == "GSE245069") & (s5.group == "CON")].FAST_singscore, "GSE245069"),
          ("GSE245069 PS", s5[(s5.study == "GSE245069") & (s5.group != "CON")].FAST_singscore, "GSE245069")]
for k, (lab, v, st) in enumerate(groups):
    ax[0].scatter(np.full(len(v), k) + np.linspace(-.12, .12, len(v)), v, s=28, color=C3[st])
    ax[0].hlines(v.mean(), k - .3, k + .3, color=INK)
ax[0].set_xticks(range(len(groups)))
ax[0].set_xticklabels([g[0] for g in groups], rotation=35, ha="right", fontsize=8)
ax[0].set_ylabel("frozen FAST singscore")
ax[0].set_title("Baseline state differs between studies\n(controls: -0.11 vs +0.10 vs +0.04)")
rr = pd.read_csv(TAB / "gap_hong_target_resampling.csv")
for k, (col, lab) in enumerate((("transfer_spearman_F8_CD", "F8-OLS (PE_in_CD)"), ("transfer_spearman_F8_MD", "F8-OLS (PE_in_MD)"))):
    v = rr[rr.fit != "full"][col]
    ax[1].scatter(np.full(len(v), k) + np.linspace(-.15, .15, len(v)), v, s=24, color=["#2a78d6", "#eb6834"][k])
    ax[1].scatter(k, rr[rr.fit == "full"][col], marker="D", s=60, color=INK, zorder=3)
ax[1].set_xticks([0, 1])
ax[1].set_xticklabels(["F8-OLS\n(PE_in_CD)", "F8-OLS\n(PE_in_MD)"])
ax[1].set_ylabel("Hong transfer Spearman")
ax[1].set_title("15 mouse-resampled DESeq2 targets\n(diamond = all 6 mice)")
u = pd.read_csv(TAB / "count_unit_per_sample.csv")
u = u[u.study != "GSE277772"]
for k, st in enumerate(("GSE294180", "GSE245069", "HONG")):
    v = u[u.study == st].counts_per_read_pair
    ax[2].scatter(np.full(len(v), k) + np.linspace(-.15, .15, len(v)), v, s=24, color=C3[st])
ax[2].axhspan(.6, 1.0, color="#e6e5e0", zorder=0)
ax[2].set_xticks(range(3))
ax[2].set_xticklabels(["GSE294180", "GSE245069", "Hong"])
ax[2].set_ylabel("gene counts per sequenced read pair")
ax[2].set_title("Count unit: Hong ~2x (coverage-derived)")
fig.tight_layout()
fig.savefig(FIG / "fig2_domain_shift_uncertainty.png", dpi=130)
plt.close(fig)

# ---- fig3: before/after on the same task ----
V = pd.read_csv(TAB / "intervention_verdicts.csv").set_index("candidate")
order = ["B0_F8_OLS", "B1_F8C_Ridge", "B2_F8_HGB", "I2a_robustCD_Ridge", "I2b_MLEavg_Ridge", "I2c_robustavg_Ridge", "I3_Huber_normalscore"]
fig, ax = plt.subplots(1, 2, figsize=(13, 4.2))
x = np.arange(len(order))
ax[0].plot(x - .2, V.loc[order, "train_spearman"], "s", color=GREY, label="train (in-fold)")
ax[0].plot(x, V.loc[order, "cv_spearman"], "o", color=INK, label="outer grouped CV (held-out genes)")
ax[0].errorbar(x + .2, V.loc[order, "rho_hong_pe"], yerr=[V.loc[order, "rho_hong_pe"] - V.loc[order, "rho_hong_lo"], V.loc[order, "rho_hong_hi"] - V.loc[order, "rho_hong_pe"]],
               fmt="D", color="#2a78d6", label="Hong PE (exploratory; bars = mouse bootstrap 95%)")
ax[0].axhline(V.loc["B0_F8_OLS", "rho_hong_pe"], color="#2a78d6", lw=.6, ls=":")
ax[0].set_xticks(x)
ax[0].set_xticklabels([o.split("_")[0] for o in order], fontsize=8)
ax[0].set_ylabel("Spearman vs PE log2FC")
ax[0].set_title("Same gene-effect task: train / CV / external")
ax[0].legend(fontsize=7, frameon=False, loc="lower left")
ax[1].errorbar(x, V.loc[order, "rho_controlfree_median"], yerr=[V.loc[order, "rho_controlfree_median"] - V.loc[order, "rho_controlfree_lo"],
                                                                V.loc[order, "rho_controlfree_hi"] - V.loc[order, "rho_controlfree_median"]], fmt="o", color="#eb6834")
for i, o in enumerate(order):
    ax[1].text(i, V.loc[o, "rho_controlfree_hi"] + .01, f"rank {int(V.loc[o, 'rank_of_true'])}/20", ha="center", fontsize=7, color=INK)
ax[1].axhline(0, color=GREY, lw=.6)
ax[1].set_xticks(x)
ax[1].set_xticklabels([o.split("_")[0] for o in order], fontsize=8)
ax[1].set_ylabel("Spearman vs PE minus other polymers (no shared control)")
ax[1].set_title("PE-specificity checks (bootstrap 95%; text = true-labelling rank)")
fig.tight_layout()
fig.savefig(FIG / "fig3_before_after_interventions.png", dpi=130)
plt.close(fig)

# ---- fig4: interpretability: standardized coefficients ----
coef = {}
for cid in ("B0_F8_OLS", "I2a_robustCD_Ridge", "I2b_MLEavg_Ridge", "I2c_robustavg_Ridge"):
    m = joblib.load(RES / "models" / f"{cid}.joblib")
    coef[cid] = pd.Series(m.named_steps["model"].coef_, index=m.named_steps["rep"].get_feature_names_out())
CF = pd.DataFrame(coef)
CF.to_csv(TAB / "candidate_coefficients_per_SD.csv")
fig, a = plt.subplots(figsize=(9, 4))
y = np.arange(len(CF))
for k, (c, col) in enumerate(zip(CF.columns, ["#0b0b0b", "#1baf7a", "#2a78d6", "#eb6834"])):
    a.barh(y + (k - 1.5) * .2, CF[c], height=.19, color=col, label=c)
a.set_yticks(y)
a.set_yticklabels(CF.index)
a.invert_yaxis()
a.axvline(0, color=INK, lw=.6)
a.set_xlabel("coefficient per SD of feature (log2FC units)")
a.legend(fontsize=7, frameon=False)
a.set_title("What the gene-level models use: FAST (fasting) and MASH reference effects dominate")
fig.tight_layout()
fig.savefig(FIG / "fig4_coefficients.png", dpi=130)
plt.close(fig)
print(CF.round(4).to_string())
