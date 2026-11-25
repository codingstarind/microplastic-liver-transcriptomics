"""Diagnostics (no new model selection): which R2 features drive the frozen MASH-state score in GSE294180; per-feature
PE contrasts on the 19 mice; cross-study performance, before/after and calibration figures."""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from e01_mash_state_models import CTX, FIG, MODELS, RES, SEED, SETS, Pipeline, StandardScaler, build_ctx, clone  # noqa: E402
from e03_apply_to_gse294180 import BothStudySignature, contr  # noqa: E402  (re-runs e03 deterministically; outputs unchanged)
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold  # noqa: E402

T = CTX["mash"]["T"]
train = ((T.study == "GSE246328") | ((T.study == "GSE316023") & (T.contrast == "MASH_primary"))).values
tgt = np.where((T.study == "GSE294180").values)[0]
base, grid = MODELS["l2_logistic"]
pipe = Pipeline([("feat", BothStudySignature("R2", "mash")), ("sc", StandardScaler()), ("clf", clone(base))])
gs = GridSearchCV(pipe, grid, cv=RepeatedStratifiedKFold(n_splits=4, n_repeats=3, random_state=SEED), scoring="neg_brier_score")
tr = np.where(train)[0]
gs.fit(tr.reshape(-1, 1), T.y.values[tr])
fz = gs.best_estimator_
names = ["R0_signature"] + [f"R1_{s}" for s in SETS]
coef = pd.DataFrame({"feature": names, "coef_per_SD": fz.named_steps["clf"].coef_[0]})
F = fz.named_steps["feat"].transform(tgt.reshape(-1, 1))
Z = fz.named_steps["sc"].transform(F)
sm = pd.read_csv(RES / "gse294180_frozen_model_scores.csv", index_col=0)
MD, PE = (sm.diet == "MD").astype(int).values, (sm.microplastic == "PE").astype(int).values
rows = []
for j, n in enumerate(names):
    contrib = Z[:, j] * coef.coef_per_SD.iloc[j]  # this feature's additive contribution to the log-odds
    c = contr(contrib, MD, PE)
    rows.append(dict(feature=n, coef_per_SD=coef.coef_per_SD.iloc[j], contrib_PE_in_CD=c["PE_in_CD"]["estimate"],
                     contrib_PE_in_MD=c["PE_in_MD"]["estimate"], contrib_Diet_in_Veh=c["Diet_in_Veh"]["estimate"]))
D = pd.DataFrame(rows).sort_values("contrib_PE_in_MD", key=abs, ascending=False)
D.to_csv(RES / "r2_feature_contributions_gse294180.csv", index=False)

M = pd.read_csv(RES / "external_metrics.csv")
prim = M[M.split.isin(["S1_train316TA_test246", "S2_train246_test316TA"])]
fig, ax = plt.subplots(1, 2, figsize=(14, 4.5))
for i, met in enumerate(["auc", "brier"]):
    p = prim.pivot_table(index=["rep", "model"], columns="split", values=met)
    p = p.drop(index=[k for k in p.index if k[1] == "majority"])
    y = np.arange(len(p))
    ax[i].scatter(p.iloc[:, 0], y, label="S1 train 316023 -> test 246328")
    ax[i].scatter(p.iloc[:, 1], y, marker="x", label="S2 train 246328 -> test 316023")
    ax[i].set_yticks(y)
    ax[i].set_yticklabels([f"{a}|{b}" for a, b in p.index], fontsize=6)
    ax[i].set_xlabel(met + (" (higher better)" if met == "auc" else " (lower better)"))
    ax[i].legend(fontsize=7)
ax[0].set_title("Study-held-out AUC: ceiling for every non-PCA representation", fontsize=9)
ax[1].set_title("Brier: where representations differ (calibration across studies)", fontsize=9)
fig.tight_layout()
fig.savefig(FIG / "fig1_cross_study_performance.png", dpi=110)

P = pd.read_csv(RES / "external_predictions.csv")
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
for a, (rep, mo) in zip(ax, [("R0", "R0_rule"), ("R1", "l2_logistic"), ("PCA", "l2_logistic")]):
    for sp, mk in [("S1_train316TA_test246", "o"), ("S2_train246_test316TA", "x")]:
        d = P[(P.rep == rep) & (P.model == mo) & (P.split == sp)]
        a.scatter(d.y + np.random.default_rng(1).uniform(-.08, .08, len(d)), d.p, marker=mk, label=sp.split("_")[0])
    a.axhline(.5, ls="--", c="k", lw=.6)
    a.set_xticks([0, 1])
    a.set_xticklabels(["chow", "GAN"])
    a.set_ylabel("predicted P(GAN)")
    a.set_title(f"{rep} {mo}: external calibration", fontsize=9)
    a.legend(fontsize=7)
fig.tight_layout()
fig.savefig(FIG / "fig2_before_after_calibration.png", dpi=110)

perm = pd.read_csv(RES / "permutation_null_auc_R2_logistic.csv")
fig, ax = plt.subplots(figsize=(6, 4))
for c in perm.columns:
    ax.hist(perm[c].dropna(), bins=20, alpha=.5, label=c.split("_")[0])
ax.axvline(1.0, c="r")
ax.set_xlabel("external AUC with permuted TRAINING labels (R2, logistic C=1)")
ax.legend()
ax.set_title("Permutation null; red = observed (1.0)", fontsize=9)
fig.tight_layout()
fig.savefig(FIG / "fig3_permutation_null.png", dpi=110)
print(D.round(3).to_string(index=False))
json.dump({"r2_final_best_params": gs.best_params_}, open(RES / "r2_final_model_info.json", "w"))
