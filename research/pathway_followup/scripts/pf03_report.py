"""Pathway follow-up, step 3: concise results table + classifier-ablation figure (reads pf01/pf02 outputs only).
Run: .venv/Scripts/python.exe research/pathway_followup/scripts/pf03_report.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OUT = Path(__file__).resolve().parents[1]
RES, FIG = OUT / "results", OUT / "figures"

# ---------------- Table 1: pathway contrasts (primary singscore) + robustness counts
T = pd.read_csv(RES / "pathway_contrasts.csv")
T = T[(T.score_type == "singscore_primary") & T.contrast.isin(["PE_in_CD", "PE_in_MD", "Interaction"])]
mz = pd.read_csv(RES / "pathway_contrasts.csv")
mz = mz[mz.score_type == "meanz_VST_secondary"].set_index(["pathway", "contrast"])
loo = pd.read_csv(RES / "loo_summary.csv").set_index(["pathway", "contrast"])
sens = pd.read_csv(RES / "covariate_and_muscle_sensitivity.csv")
sc = sens.groupby(["pathway", "contrast"]).agg(n_sens=("same_sign_as_primary", "size"), n_sens_same_sign=("same_sign_as_primary", "sum"),
                                               sens_max_p=("p", "max")).reset_index().set_index(["pathway", "contrast"])
uq = pd.read_csv(RES / "pathway_contrasts_unique_genes.csv").set_index(["pathway", "contrast"])
ext = pd.read_csv(RES / "external_studies_pathway_contrasts.csv")
hpe = ext[(ext.study == "HONG") & (ext.exposure == "PE")].set_index("pathway")
rows = []
for _, r in T.iterrows():
    k = (r.pathway, r.contrast)
    rows.append(dict(pathway=r.pathway, contrast=r.contrast, estimate=r.estimate, ci95=f"[{r.ci_low:.4f}, {r.ci_high:.4f}]",
                     std_effect=r.std_effect, p=r.p, p_exact_perm=r.p_exact_perm, q_BH_8scores=r.p_BH_8,
                     q_lt_0_05=r.p_BH_8 < .05, q_lt_0_10=r.p_BH_8 < .10,
                     orig_fry_direction=r.orig_fry_direction, orig_fry_q=r.orig_fry_q_family, orig_fry_family=r.orig_fry_family,
                     meanz_VST_estimate=mz.loc[k, "estimate"], meanz_VST_p=mz.loc[k, "p"],
                     unique_genes_estimate=uq.loc[k, "estimate"], unique_genes_p=uq.loc[k, "p"], n_unique_genes=uq.loc[k, "n_genes_unique"],
                     loo_n_same_sign_of_19=loo.loc[k, "n_same_sign"], loo_max_p=loo.loc[k, "loo_max_p"], loo_most_influential=loo.loc[k, "most_influential"],
                     sensitivity_n_same_sign=f"{int(sc.loc[k, 'n_sens_same_sign'])}/{int(sc.loc[k, 'n_sens'])}", sensitivity_max_p=sc.loc[k, "sens_max_p"],
                     hong_PE_vs_CON_estimate=hpe.loc[r.pathway, "estimate"], hong_PE_p_exact=hpe.loc[r.pathway, "p_exact_perm"]))
R = pd.DataFrame(rows)
R.to_csv(OUT / "TABLE_results_concise.csv", index=False)

# ---------------- Figure 2: ablation comparison on a common scale
B = pd.read_csv(RES / "ablation_scaled_effects_bootstrap.csv")
C = pd.read_csv(RES / "ablation_gse294180_contrasts.csv")
M = pd.read_csv(RES / "ablation_external_metrics.csv")
order = ["full_R2", "drop_R0", "drop_CHOL_SYN", "drop_CHOL_HOM", "drop_CHOL_both", "drop_MTORC1", "drop_CHOL_both+MTORC1",
         "drop_PPARA", "drop_PPARA_group(PPARA+FAM+PEROX)", "drop_LIPO", "drop_sterol_lipogenic(CHOL_both+LIPO)",
         "drop_5_metabolic(CHOL_both+LIPO+PPARA+MTORC1)", "drop_all_lipid(CHOL_both+LIPO+PPARA+FAM+PEROX+BILE)", "drop_all_lipid+MTORC1",
         "keep_only_R0+INFL+EMT", "drop_FAM", "drop_PEROX", "drop_BILE", "drop_INFL", "drop_EMT", "drop_XENO", "drop_OXPHOS"]
y = np.arange(len(order))[::-1]
fig, axes = plt.subplots(1, 3, figsize=(13.5, 7.6), gridspec_kw={"width_ratios": [1, 1, .75], "wspace": .42}, sharey=True)
for ax, c, title in ((axes[0], "PE_in_CD", "PE vs Veh · CD"), (axes[1], "PE_in_MD", "PE vs Veh · MD")):
    d = B[B.contrast == c].set_index("ablation").loc[order]
    h = C[C.contrast == c].set_index("ablation").loc[order]
    sig = h.p_holm.values < .05
    col = np.where(sig, "#1c5cab", "#8c8b87")
    ax.axvline(0, color="#52514e", lw=.8)
    ax.axvline(d.loc["full_R2", "scaled_effect"], color="#1c5cab", lw=.6, ls=(0, (2, 2)))
    ax.hlines(y, d.scaled_ci_low, d.scaled_ci_high, color=col, lw=1.3)
    ax.scatter(d.scaled_effect, y, s=26, color=np.where(sig, "#1c5cab", "white"), edgecolor=col, lw=1.1, zorder=3)
    for yy, (_, r) in zip(y, h.iterrows()):
        ax.text(1.01, yy, f"{r.estimate:+.2f} lo\npH {r.p_holm:.2g}", transform=ax.get_yaxis_transform(), ha="left", clip_on=False, linespacing=.9, va="center", fontsize=5.8, color="#52514e")
    ax.set_title(title, fontsize=8.5, loc="left")
    ax.set_xlabel("PE contrast / model's training MASH−control separation\n(mouse bootstrap 95% CI)", fontsize=7)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
axes[0].set_yticks(y)
axes[0].set_yticklabels(order, fontsize=6.8)
m = M.pivot_table(index="ablation", columns="split", values=["bal_acc", "brier"]).loc[order]
ax = axes[2]
for j, (sp, lab, mk) in enumerate((("S1_train316TA_test246", "S1", "o"), ("S2_train246_test316TA", "S2", "s"))):
    ax.scatter(m[("bal_acc", sp)], y, marker=mk, s=18, color=["#1c5cab", "#c0392b"][j], label=f"{lab} balanced accuracy", zorder=3)
for yy, (_, r) in zip(y, m.iterrows()):
    ax.text(1.02, yy, f"Brier {r[('brier', 'S1_train316TA_test246')]:.3f}/{r[('brier', 'S2_train246_test316TA')]:.3f}",
            transform=ax.get_yaxis_transform(), fontsize=5.8, va="center", color="#52514e")
ax.set_xlim(.7, 1.02)
ax.set_title("Study-held-out validity (AUC = 1.0 in all)", fontsize=8.5, loc="left")
ax.set_xlabel("balanced accuracy", fontsize=7)
ax.legend(fontsize=6, frameon=False, loc="lower left")
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)
fig.suptitle("R2 L2-logistic MASH-state classifier: leave-feature(s)-out ablation, frozen models applied to GSE294180. "
             "Filled = Holm p < 0.05 (3 PE contrasts). Text: raw log-odds, Holm p.", fontsize=8, x=.01, ha="left")
fig.subplots_adjust(left=.24, right=.92, top=.93, bottom=.1)
for e in ("png", "pdf", "svg"):
    fig.savefig(FIG / f"fig2_classifier_ablation.{e}", dpi=300)
print(R.round(4).to_string(index=False))
