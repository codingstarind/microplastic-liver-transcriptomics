"""C: apply frozen MASH-state models to the 19 GSE294180 mice ONLY IF they passed results/validity.csv.
Final models are retrained on both MASH studies (31 mice; no GSE294180 data). The R0 signature uses genes passing the
training rule in BOTH studies (same sign, |diff|>=1), ranked by the smaller |t|. Output = decision function (log-odds),
analysed with score ~ diet*PE (Holm over the 3 PE contrasts) + 19 leave-one-out refits + technical associations."""
import json

import numpy as np
import pandas as pd
from e01_mash_state_models import (CTX, MODELS, RES, FIG, SEED, FeatureBuilder, Pipeline, StandardScaler, build_ctx, clone,
                                   singscore_signed)
from scipy import stats
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold

APPLY = [("R0", "l2_logistic"), ("R2", "l2_logistic")]  # predeclared clarification (05_ doc), fixed before external results
v = pd.read_csv(RES / "validity.csv")
ok = [(r, m) for r, m in APPLY if bool(v[(v.rep == r) & (v.model == m)].valid.iloc[0])]
info = {"requested": APPLY, "valid_and_applied": ok}
if not ok:
    info["decision"] = "NOT APPLIED: no predeclared model passed the study-held-out validity criteria"
    json.dump(info, open(RES / "application_decision.json", "w"), indent=2)
    print(info)
    raise SystemExit(0)

build_ctx("mash")
T = CTX["mash"]["T"]
X, R = CTX["mash"]["X"], CTX["mash"]["R"]
train = ((T.study == "GSE246328") | ((T.study == "GSE316023") & (T.contrast == "MASH_primary"))).values
tgt = (T.study == "GSE294180").values


class BothStudySignature(FeatureBuilder):
    """R0 built from genes concordant in both MASH studies (training data only)."""

    def fit(self, idx, y):
        idx = np.asarray(idx).ravel().astype(int)
        st = T.study.values[idx]
        ts, ds = [], []
        for s_ in ("GSE316023", "GSE246328"):
            j = idx[st == s_]
            yy = np.asarray(y)[st == s_]
            a, b = X[j][yy == 1], X[j][yy == 0]
            ts.append(np.nan_to_num(stats.ttest_ind(a, b, equal_var=False).statistic))
            ds.append(a.mean(0) - b.mean(0))
        t1, t2 = ts
        d1, d2 = ds
        upc = np.where((d1 >= 1) & (d2 >= 1))[0]
        dnc = np.where((d1 <= -1) & (d2 <= -1))[0]
        self.up_ = upc[np.argsort(-np.minimum(t1[upc], t2[upc]))][:50]
        self.dn_ = dnc[np.argsort(np.maximum(t1[dnc], t2[dnc]))][:50]
        return self


rows, scores, loo = [], {}, []
sm = pd.read_csv(RS_ := (RES.parents[1] / "roadmap_phase4d/results/per_mouse_scores.csv"), index_col=0)
for rep, model in ok:
    base, grid = MODELS[model]
    pipe = Pipeline([("feat", BothStudySignature(rep, "mash")), ("sc", StandardScaler()), ("clf", clone(base))])
    tr = np.where(train)[0]
    gs = GridSearchCV(pipe, grid, cv=RepeatedStratifiedKFold(n_splits=4, n_repeats=3, random_state=SEED), scoring="neg_brier_score")
    gs.fit(tr.reshape(-1, 1), T.y.values[tr])
    fz = gs.best_estimator_
    te = np.where(tgt)[0]
    s = fz.decision_function(te.reshape(-1, 1))
    key = f"{rep}_{model}"
    scores[key] = pd.Series(s, index=T.sample_id.values[te])
    rows.append(dict(model=key, best_params=json.dumps(gs.best_params_), n_up=len(fz.named_steps["feat"].up_), n_dn=len(fz.named_steps["feat"].dn_)))
S = pd.DataFrame(scores).join(sm[["group", "diet", "microplastic", "depth", "quality_proxy", "muscle", "clock_phase", "MASH", "FAST"]])
S.to_csv(RES / "gse294180_frozen_model_scores.csv")
MD, PE = (S.diet == "MD").astype(int).values, (S.microplastic == "PE").astype(int).values


def contr(y, md, pe):
    Xd = np.column_stack([np.ones(len(y)), md, pe, md * pe])
    b, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    r = y - Xd @ b
    df = len(y) - 4
    cov = (r @ r / df) * np.linalg.inv(Xd.T @ Xd)
    out = {}
    for k, l in {"PE_in_CD": [0, 0, 1, 0], "PE_in_MD": [0, 0, 1, 1], "Interaction": [0, 0, 0, 1], "Diet_in_Veh": [0, 1, 0, 0]}.items():
        l = np.array(l, float)
        e, se = l @ b, np.sqrt(l @ cov @ l)
        q = stats.t.ppf(.975, df)
        out[k] = dict(estimate=e, ci_low=e - q * se, ci_high=e + q * se, p=2 * stats.t.sf(abs(e / se), df), std_effect=e / np.sqrt(r @ r / df))
    return out


res = []
for key in list(scores) + ["MASH"]:
    y = S[key].values
    for k, d in contr(y, MD, PE).items():
        res.append(dict(score=key if key != "MASH" else "Phase4D_MASH_singscore", contrast=k, **d))
    for i in range(19):
        m = np.arange(19) != i
        for k, d in contr(y[m], MD[m], PE[m]).items():
            loo.append(dict(score=key, dropped=S.index[i], contrast=k, estimate=d["estimate"]))
Rz = pd.DataFrame(res)
Rz["p_holm"] = np.nan
for sc_ in Rz.score.unique():
    i = (Rz.score == sc_) & Rz.contrast.isin(["PE_in_CD", "PE_in_MD", "Interaction"])
    p = Rz.loc[i, "p"].values
    o = np.argsort(p)
    adj, run = np.empty(3), 0
    for r_, j in enumerate(o):
        run = max(run, (3 - r_) * p[j])
        adj[j] = min(1, run)
    Rz.loc[i, "p_holm"] = adj
Rz.to_csv(RES / "gse294180_frozen_model_contrasts.csv", index=False)
L = pd.DataFrame(loo)
Ls = L.merge(Rz[["score", "contrast", "estimate"]].replace({"Phase4D_MASH_singscore": "MASH"}), on=["score", "contrast"], suffixes=("", "_full"))
Ls = Ls.groupby(["score", "contrast"]).apply(lambda d: pd.Series(dict(min=d.estimate.min(), max=d.estimate.max(),
                                                                         n_sign_same=int((np.sign(d.estimate) == np.sign(d.estimate_full)).sum()))),
                                             include_groups=False).reset_index()
Ls.to_csv(RES / "gse294180_frozen_model_loo.csv", index=False)
tech = pd.DataFrame([dict(score=k, covariate=c, spearman=stats.spearmanr(S[k], S[c]).statistic) for k in list(scores) + ["MASH"]
                     for c in ("depth", "quality_proxy", "muscle", "clock_phase", "MASH", "FAST")])
tech.to_csv(RES / "gse294180_frozen_model_technical.csv", index=False)
info.update(models=rows, decision="APPLIED")
json.dump(info, open(RES / "application_decision.json", "w"), indent=2, default=str)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

order = ["CD+Veh", "CD+PE", "MD+Veh", "MD+PE"]
keys = list(scores) + ["MASH"]
fig, ax = plt.subplots(1, len(keys), figsize=(5 * len(keys), 4))
for a, k in zip(np.atleast_1d(ax), keys):
    for j, g in enumerate(order):
        vv = S[S.group == g][k]
        a.scatter(np.full(len(vv), j), vv, s=25)
        a.hlines(vv.mean(), j - .25, j + .25, color="k")
    a.set_xticks(range(4))
    a.set_xticklabels(order)
    a.set_title(k if k != "MASH" else "Phase 4D MASH singscore", fontsize=9)
fig.tight_layout()
fig.savefig(FIG / "fig4_gse294180_2x2_frozen_scores.png", dpi=110)
with pd.option_context("display.width", 200):
    print(pd.DataFrame(rows).to_string(index=False))
    print(Rz.round(4).to_string(index=False))
    print(Ls.round(4).to_string(index=False))
    print(tech.round(3).pivot(index="score", columns="covariate", values="spearman").to_string())
