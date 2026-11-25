"""Task B summary figure + merged table (original run with the degenerate HGB replaced by the D2 fix)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from p4_data import P4, RES  # noqa: E402

B = RES / "taskB"
f = pd.read_csv(B / "internal_cv_folds.csv")
f = pd.concat([f[f.model != "hgb"], pd.read_csv(B / "internal_cv_folds_hgb_fix.csv")])
e = pd.read_csv(B / "external_GSE246328.csv")
e = pd.concat([e[e.model != "hgb"], pd.read_csv(B / "external_GSE246328_hgb_fix.csv")])
d = pd.read_csv(B / "descriptive_GSE294180_vehicle.csv")
d = pd.concat([d[d.model != "hgb"], pd.read_csv(B / "descriptive_GSE294180_vehicle_hgb_fix.csv")])
order = ["dummy", "logistic_l2", "gaussian_nb", "tree_gini", "tree_entropy", "rf", "hgb"]
tab = f.groupby("model")[["balanced_accuracy", "roc_auc", "brier"]].agg(["mean", "std"]).reindex(order).round(3)
tab.columns = [f"internal_{a}_{b}" for a, b in tab.columns]
tab = tab.join(e.set_index("model")[["balanced_accuracy", "roc_auc", "sensitivity", "specificity", "brier", "tn", "fp", "fn", "tp"]]
               .add_prefix("external_").round(3))
tab = tab.join(d.groupby(["model", "diet"]).p_GAN_state.mean().unstack().add_prefix("GSE294180_mean_pGAN_").round(2))
tab.to_csv(B / "taskB_summary_merged.csv")
fig, ax = plt.subplots(1, 2, figsize=(12, 4))
x = range(len(order))
ax[0].bar([i - .2 for i in x], tab["internal_balanced_accuracy_mean"], .4, yerr=tab["internal_balanced_accuracy_std"], label="internal CV (GSE316023, n=32)")
ax[0].bar([i + .2 for i in x], tab["external_balanced_accuracy"], .4, label="external GSE246328 (n=15)")
ax[0].set_xticks(list(x))
ax[0].set_xticklabels(order, rotation=40, fontsize=8)
ax[0].set_ylim(0, 1.05)
ax[0].legend(fontsize=7)
ax[0].set_title("Task B (demonstration): GAN vs chow; ceiling effect, not a model ranking", fontsize=9)
for i, mo in enumerate(order):
    s = d[d.model == mo]
    ax[1].scatter([i - .1] * (s.diet == "CD").sum(), s[s.diet == "CD"].p_GAN_state, c="tab:blue", s=15, label="CD+Veh" if i == 0 else None)
    ax[1].scatter([i + .1] * (s.diet == "MD").sum(), s[s.diet == "MD"].p_GAN_state, c="tab:red", s=15, label="MD+Veh" if i == 0 else None)
ax[1].set_xticks(list(x))
ax[1].set_xticklabels(order, rotation=40, fontsize=8)
ax[1].set_ylabel("predicted P(GAN state)")
ax[1].legend(fontsize=7)
ax[1].set_title("Descriptive: GSE294180 vehicle mice (purified control diet; not comparable)", fontsize=9)
fig.tight_layout()
fig.savefig(P4 / "figures" / "taskB_classification_demo.png", dpi=110)
print(tab.to_string())
