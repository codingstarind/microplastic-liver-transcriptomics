"""Library for pf02 (importable by joblib workers). Wraps the unchanged Phase 4E e01 module."""
import itertools
import json

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.base import clone
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import e01_mash_state_models as e01

FEATS = ["R0_signature"] + [f"R1_{s}" for s in e01.SETS]
G = {"CHOL_SYN": "R1_CHOL_SYNTHESIS_CORE", "CHOL_HOM": "R1_HALLMARK_CHOLESTEROL_HOMEOSTASIS", "LIPO": "R1_LIPOGENESIS_CORE",
     "PPARA": "R1_SANDERSON_PPARA_TARGETS", "FAM": "R1_HALLMARK_FATTY_ACID_METABOLISM", "PEROX": "R1_HALLMARK_PEROXISOME",
     "BILE": "R1_HALLMARK_BILE_ACID_METABOLISM", "MTORC1": "R1_HALLMARK_MTORC1_SIGNALING", "INFL": "R1_HALLMARK_INFLAMMATORY_RESPONSE",
     "EMT": "R1_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "XENO": "R1_HALLMARK_XENOBIOTIC_METABOLISM",
     "OXPHOS": "R1_HALLMARK_OXIDATIVE_PHOSPHORYLATION", "R0": "R0_signature"}
ABLATIONS = {"full_R2": []}
ABLATIONS.update({f"drop_{k}": [v] for k, v in G.items()})
ABLATIONS.update({
    "drop_CHOL_both": [G["CHOL_SYN"], G["CHOL_HOM"]],
    "drop_CHOL_both+MTORC1": [G["CHOL_SYN"], G["CHOL_HOM"], G["MTORC1"]],
    "drop_PPARA_group(PPARA+FAM+PEROX)": [G["PPARA"], G["FAM"], G["PEROX"]],
    "drop_sterol_lipogenic(CHOL_both+LIPO)": [G["CHOL_SYN"], G["CHOL_HOM"], G["LIPO"]],
    "drop_5_metabolic(CHOL_both+LIPO+PPARA+MTORC1)": [G["CHOL_SYN"], G["CHOL_HOM"], G["LIPO"], G["PPARA"], G["MTORC1"]],
    "drop_all_lipid(CHOL_both+LIPO+PPARA+FAM+PEROX+BILE)": [G[k] for k in ("CHOL_SYN", "CHOL_HOM", "LIPO", "PPARA", "FAM", "PEROX", "BILE")],
    "drop_all_lipid+MTORC1": [G[k] for k in ("CHOL_SYN", "CHOL_HOM", "LIPO", "PPARA", "FAM", "PEROX", "BILE", "MTORC1")],
    "keep_only_R0+INFL+EMT": [f for f in FEATS if f not in (G["R0"], G["INFL"], G["EMT"])],
})
SPLITS = ["S1_train316TA_test246", "S2_train246_test316TA"]


class AblFB(e01.FeatureBuilder):
    def __init__(self, rep="R2", task="mash", drop=()):
        super().__init__(rep, task)
        self.drop = drop

    def transform(self, idx):
        F = super().transform(idx)
        keep = [i for i, f in enumerate(FEATS) if f not in self.drop]
        return F[:, keep]


class AblBoth(AblFB):
    """e03 final-model R0: genes concordant (|diff|>=1, same sign) in both MASH studies, ranked by the smaller |t|."""

    def fit(self, idx, y):
        idx = np.asarray(idx).ravel().astype(int)
        T, X = e01.CTX["mash"]["T"], e01.CTX["mash"]["X"]
        st = T.study.values[idx]
        ts, ds = [], []
        for s_ in ("GSE316023", "GSE246328"):
            j = idx[st == s_]
            yy = np.asarray(y)[st == s_]
            a, b = X[j][yy == 1], X[j][yy == 0]
            ts.append(np.nan_to_num(stats.ttest_ind(a, b, equal_var=False).statistic))
            ds.append(a.mean(0) - b.mean(0))
        (t1, t2), (d1, d2) = ts, ds
        upc = np.where((d1 >= 1) & (d2 >= 1))[0]
        dnc = np.where((d1 <= -1) & (d2 <= -1))[0]
        self.up_ = upc[np.argsort(-np.minimum(t1[upc], t2[upc]))][:50]
        self.dn_ = dnc[np.argsort(np.maximum(t1[dnc], t2[dnc]))][:50]
        return self


def _pipe(fb):
    base, grid = e01.MODELS["l2_logistic"]
    return Pipeline([("feat", fb), ("sc", StandardScaler()), ("clf", clone(base))]), grid


def _fit(fb, tr, y):
    pipe, grid = _pipe(fb)
    gs = GridSearchCV(pipe, grid, cv=RepeatedStratifiedKFold(n_splits=4, n_repeats=3, random_state=e01.SEED),
                      scoring="neg_brier_score", n_jobs=1, error_score="raise")
    gs.fit(tr.reshape(-1, 1), y)
    return gs.best_estimator_, gs.best_params_


def calib(y, p):
    lp = np.log(np.clip(p, 1e-6, 1 - 1e-6) / (1 - np.clip(p, 1e-6, 1 - 1e-6)))
    if len(set(y)) < 2:
        return np.nan, np.nan
    m = LogisticRegression(C=1e6, max_iter=10000).fit(lp[:, None], y)
    return float(m.intercept_[0]), float(m.coef_[0, 0])


def contr(y, md, pe, cov=None):
    cols = [np.ones(len(y)), md, pe, md * pe] + ([] if cov is None else [cov])
    Xd = np.column_stack(cols)
    b, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    r = y - Xd @ b
    df = len(y) - Xd.shape[1]
    s2 = r @ r / df
    cov_ = s2 * np.linalg.inv(Xd.T @ Xd)
    out = {}
    for k, l in {"PE_in_CD": [0, 0, 1, 0], "PE_in_MD": [0, 0, 1, 1], "Interaction": [0, 0, 0, 1], "Diet_in_Veh": [0, 1, 0, 0]}.items():
        l = np.array(l + [0] * (Xd.shape[1] - 4), float)
        e, se = l @ b, np.sqrt(l @ cov_ @ l)
        q = stats.t.ppf(.975, df)
        out[k] = dict(estimate=e, ci_low=e - q * se, ci_high=e + q * se, p=2 * stats.t.sf(abs(e / se), df), std_effect=e / np.sqrt(s2))
    return out


def holm3(p):
    p = np.asarray(p, float); o = np.argsort(p); adj = np.empty(3); run = 0
    for r_, j in enumerate(o):
        run = max(run, (3 - r_) * p[j]); adj[j] = min(1, run)
    return adj


def run_ablation(name, drop):
    if "mash" not in e01.CTX:
        e01.build_ctx("mash")
    T = e01.CTX["mash"]["T"]
    met = []
    for sp in SPLITS:
        trm, tem = e01.split_masks(T, sp)
        tr, te = np.where(trm.values)[0], np.where(tem.values)[0]
        fz, best = _fit(AblFB("R2", "mash", tuple(drop)), tr, T.y.values[tr])
        p = fz.predict_proba(te.reshape(-1, 1))[:, 1]
        s = fz.decision_function(te.reshape(-1, 1))
        y = T.y.values[te]
        m = e01.metrics(y, p, s)
        ci, cs = calib(y, p)
        met.append(dict(ablation=name, dropped=";".join(drop), n_features=len(FEATS) - len(drop), split=sp, best_params=json.dumps(best),
                        **m, calib_intercept=ci, calib_slope=cs, test_logodds_separation=s[y == 1].mean() - s[y == 0].mean()))
    # final model: both MASH studies (as e03), applied to GSE294180
    train = ((T.study == "GSE246328") | ((T.study == "GSE316023") & (T.contrast == "MASH_primary"))).values
    tgt = (T.study == "GSE294180").values
    tr, te = np.where(train)[0], np.where(tgt)[0]
    fz, best = _fit(AblBoth("R2", "mash", tuple(drop)), tr, T.y.values[tr])
    s_tr = fz.decision_function(tr.reshape(-1, 1))
    sep = float(s_tr[T.y.values[tr] == 1].mean() - s_tr[T.y.values[tr] == 0].mean())
    s = pd.Series(fz.decision_function(te.reshape(-1, 1)), index=T.sample_id.values[te], name=name)
    sm = pd.read_csv(e01.REPO / "metadata/sample_metadata.csv").set_index("counts_column").loc[s.index]
    p4d = pd.read_csv(e01.RS / "roadmap_phase4d/results/per_mouse_scores.csv", index_col=0).loc[s.index]
    md, pe = (sm.diet == "MD").astype(int).values, (sm.microplastic == "PE").astype(int).values
    cc = contr(s.values, md, pe)
    hol = holm3([cc[k]["p"] for k in ("PE_in_CD", "PE_in_MD", "Interaction")])
    rows = []
    for j, (k, d) in enumerate(cc.items()):
        row = dict(ablation=name, contrast=k, **d, p_holm=hol[j] if k != "Diet_in_Veh" else np.nan,
                   train_separation=sep, scaled_by_train_separation=d["estimate"] / sep,
                   scaled_by_diet_in_veh=d["estimate"] / cc["Diet_in_Veh"]["estimate"])
        if k in ("PE_in_CD", "PE_in_MD"):
            msk = (sm.diet == k[-2:]).values
            yy, pp = s.values[msk], pe[msk]
            obs = yy[pp == 1].mean() - yy[pp == 0].mean()
            dd = np.array([yy[list(c)].mean() - np.delete(yy, list(c)).mean() for c in itertools.combinations(range(len(yy)), int(pp.sum()))])
            row["p_exact_perm"] = float(np.mean(np.abs(dd) >= abs(obs) - 1e-12))
        rows.append(row)
    loo = []
    for i in range(len(s)):
        m_ = np.arange(len(s)) != i
        o = contr(s.values[m_], md[m_], pe[m_])
        for k in ("PE_in_CD", "PE_in_MD", "Interaction"):
            loo.append(dict(ablation=name, dropped_mouse=s.index[i], contrast=k, estimate=o[k]["estimate"], p=o[k]["p"]))
    cv = []
    for c in ("depth", "quality_proxy", "muscle", "clock_phase", "FAST"):
        z = p4d[c].values
        o = contr(s.values, md, pe, cov=(z - z.mean()) / z.std())
        for k in ("PE_in_CD", "PE_in_MD", "Interaction"):
            cv.append(dict(ablation=name, adjustment=f"cov_{c}", contrast=k, **o[k]))
    clf = fz.named_steps["clf"]
    info = dict(best_params=best, train_separation=sep, coef_per_SD=dict(zip([f for f in FEATS if f not in drop], clf.coef_[0].tolist())),
                intercept=float(clf.intercept_[0]), n_up=len(fz.named_steps["feat"].up_), n_dn=len(fz.named_steps["feat"].dn_))
    return dict(name=name, metrics=pd.DataFrame(met), contrasts=pd.DataFrame(rows), scores=s.to_frame(), loo=pd.DataFrame(loo),
                cov=pd.DataFrame(cv), info=info)
