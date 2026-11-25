"""Phase 4E: study-held-out MASH-state classifier (B1), fasting demonstration (B2) and, if valid, frozen application
to the 19 GSE294180 mice (C). Plan: 03_predeclared_model_evaluation_plan.md. All data-dependent steps (R0 signature
genes, PCA loadings, scaling) are fitted INSIDE each training fold via FeatureBuilder, which receives sample indices."""
import json
import os
import sys
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "1"
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy import stats  # noqa: E402
from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin, clone  # noqa: E402
from sklearn.dummy import DummyClassifier  # noqa: E402
from sklearn.ensemble import RandomForestClassifier  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import balanced_accuracy_score, brier_score_loss, log_loss, roc_auc_score  # noqa: E402
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold  # noqa: E402
from sklearn.naive_bayes import GaussianNB  # noqa: E402
from sklearn.pipeline import Pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402
from sklearn.tree import DecisionTreeClassifier  # noqa: E402

P4E = Path(__file__).resolve().parents[1]
RS, REPO = P4E.parent, P4E.parents[1]
RES, TAB, FIG = P4E / "results", P4E / "tables", P4E / "figures"
DROOT = Path("D:/mpliver_work/phase2_refs")
SEED = 20261122
N_JOBS = int(os.environ.get("P4E_JOBS", "12"))
SETS = ["LIPOGENESIS_CORE", "CHOL_SYNTHESIS_CORE", "SANDERSON_PPARA_TARGETS", "HALLMARK_FATTY_ACID_METABOLISM",
        "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_XENOBIOTIC_METABOLISM", "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_PEROXISOME",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_MTORC1_SIGNALING"]

# ---------------- data (cached as module globals; read-only in workers) ----------------
el = pd.read_csv(RS / "roadmap_phase2/metadata/phase3_eligible_samples.csv", dtype={"sample_id": str})
U = pd.read_csv(RS / "roadmap_phase2/metadata/phase3_gene_universe.csv", index_col=0)
UMAIN = set(U.index[U.universe_main])
gsets = pd.read_csv(RS / "roadmap_phase3/results/gene_sets.csv")


def ready(acc):
    x = pd.read_csv(DROOT / acc / "prepared/expression_ready.csv.gz", index_col=0)
    x.columns = x.columns.astype(str)
    return x


def eligible(cid):
    e = el[el.contrast_id == cid]
    return list(e.sample_id), (e.arm == "case").astype(int).values


raw = {"GSE316023": ready("GSE316023"), "GSE246328": ready("GSE246328"), "GSE107787": ready("GSE107787"), "GSE225695": ready("GSE225695")}
g294 = pd.read_csv(REPO / "data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz", index_col=0)
sm294 = pd.read_csv(REPO / "metadata/sample_metadata.csv").set_index("counts_column").loc[list(g294.columns)]
raw["GSE294180"] = g294
GENES = {"mash": sorted(UMAIN & set(raw["GSE316023"].index) & set(raw["GSE246328"].index) & set(g294.index)),
         "fast": sorted(UMAIN & set(raw["GSE107787"].index) & set(raw["GSE225695"].index))}


def logexpr(acc, cols, genes, kind):
    x = raw[acc].loc[genes, cols].astype(float)
    return np.log2(x / x.sum() * 1e6 + 1) if kind == "counts" else np.log2(x + 1)


def make_samples(task):
    """Returns (sample table, dict sample_key -> (log-expression vector, within-sample ranks)) for a task."""
    genes = GENES[task]
    rows, L = [], {}
    spec = {"mash": [("GSE316023", "MASH_primary", "counts"), ("GSE316023", "MASH_sens_TN", "counts"),
                     ("GSE246328", "MASH_secondary", "counts")],
            "fast": [("GSE107787", "FASTING_primary", "rpkm"), ("GSE225695", "FASTING_secondary", "counts")]}[task]
    for acc, cid, kind in spec:
        ids, y = eligible(cid)
        e = el[el.contrast_id == cid].set_index("sample_id")
        M = logexpr(acc, ids, genes, kind)
        for s, lab in zip(ids, y):
            k = f"{acc}|{s}"
            rows.append(dict(key=k, study=acc, contrast=cid, sample_id=s, gsm=e.loc[s, "gsm"], y=int(lab),
                             clean=bool(e.loc[s, "in_clean_sensitivity"]), libsize=float(raw[acc].loc[genes, s].sum())))
            L[k] = M[s].values
    if task == "mash":
        M = logexpr("GSE294180", list(g294.columns), genes, "counts")
        for s in g294.columns:
            k = f"GSE294180|{s}"
            rows.append(dict(key=k, study="GSE294180", contrast="target_19", sample_id=s, gsm=sm294.loc[s, "geo_accession"],
                             y=-1, clean=True, libsize=float(g294.loc[genes, s].sum())))
            L[k] = M[s].values
    T = pd.DataFrame(rows).drop_duplicates("key").reset_index(drop=True)
    X = np.stack([L[k] for k in T.key])  # samples x genes (log expression)
    R = np.apply_along_axis(lambda v: stats.rankdata(v), 1, X)  # within-sample ranks
    return T, X, R, genes


def singscore_up(R, idx):
    N, n = R.shape[1], len(idx)
    smin, smax = (n + 1) / 2, (2 * N - n + 1) / 2
    return (R[:, idx].mean(1) - smin) / (smax - smin) - 0.5


def singscore_signed(R, up, dn):
    N = R.shape[1]
    return singscore_up(R, up) + singscore_up(N + 1 - R, dn)


def welch_signature(X, y, n=50, min_diff=1.0):
    a, b = X[y == 1], X[y == 0]
    t = stats.ttest_ind(a, b, equal_var=False).statistic
    d = a.mean(0) - b.mean(0)
    t = np.nan_to_num(t)
    up = np.where((d >= min_diff))[0]
    dn = np.where((d <= -min_diff))[0]
    up = up[np.argsort(-t[up])][:n]
    dn = dn[np.argsort(t[dn])][:n]
    return up, dn


_MEMO = {}


def _memo(kind, task, idx, y, fn):
    """Per-process cache: identical training folds (re-used by every hyperparameter setting) are computed once."""
    k = (kind, task, np.asarray(idx).tobytes(), np.asarray(y).tobytes())
    if k not in _MEMO:
        _MEMO[k] = fn()
    return _MEMO[k]


class FeatureBuilder(BaseEstimator, TransformerMixin):
    """X = integer sample indices into module-level CTX (X log-expr, R ranks). Fits signature/PCA on training rows only."""

    def __init__(self, rep="R2", task="mash"):
        self.rep, self.task = rep, task

    def fit(self, idx, y):
        idx = np.asarray(idx).ravel().astype(int)
        Xg, Rg, gi = CTX[self.task]["X"], CTX[self.task]["R"], CTX[self.task]["set_idx"]
        if self.rep in ("R0", "R2"):
            self.up_, self.dn_ = _memo("sig", self.task, idx, y, lambda: welch_signature(Xg[idx], np.asarray(y)))
        if self.rep == "PCA":
            def _pca():
                Z = Xg[idx]
                Z = (Z - Z.mean(1, keepdims=True)) / Z.std(1, keepdims=True)
                mu = Z.mean(0)
                _, _, vt = np.linalg.svd(Z - mu, full_matrices=False)
                return mu, vt[:3]
            self.mu_, self.load_ = _memo("pca", self.task, idx, y, _pca)
        return self

    def transform(self, idx):
        idx = np.asarray(idx).ravel().astype(int)
        Xg, Rg, gi = CTX[self.task]["X"], CTX[self.task]["R"], CTX[self.task]["set_idx"]
        R = Rg[idx]
        f = []
        if self.rep in ("R0", "R2"):
            f.append(singscore_signed(R, self.up_, self.dn_)[:, None])
        if self.rep in ("R1", "R2"):
            f.append(CTX[self.task]["R1"][idx])
        if self.rep == "PCA":
            Z = Xg[idx]
            Z = (Z - Z.mean(1, keepdims=True)) / Z.std(1, keepdims=True)
            f.append((Z - self.mu_) @ self.load_.T)
        return np.hstack(f)


class R0Rule(BaseEstimator, ClassifierMixin):
    """Baseline: class = 1 if score > midpoint of training class means (score = single R0 feature)."""

    def fit(self, X, y):
        X = np.asarray(X)[:, 0]
        self.thr_ = (X[y == 1].mean() + X[y == 0].mean()) / 2
        self.classes_ = np.array([0, 1])
        return self

    def decision_function(self, X):
        return np.asarray(X)[:, 0] - self.thr_

    def predict_proba(self, X):
        p = (self.decision_function(X) > 0).astype(float)
        return np.column_stack([1 - p, p])

    def predict(self, X):
        return (self.decision_function(X) > 0).astype(int)


MODELS = {
    "majority": (DummyClassifier(strategy="prior"), {}),
    "R0_rule": (R0Rule(), {}),
    "logistic_weak": (LogisticRegression(C=100, max_iter=5000), {}),
    "logistic_default": (LogisticRegression(C=1, max_iter=5000), {}),
    "l2_logistic": (LogisticRegression(max_iter=5000), {"clf__C": [0.01, 0.1, 1, 10]}),
    "gaussian_nb": (GaussianNB(), {"clf__var_smoothing": [1e-9, 1e-3]}),
    "tree": (DecisionTreeClassifier(random_state=SEED), {"clf__max_depth": [1, 2], "clf__criterion": ["gini", "entropy"]}),
    "random_forest": (RandomForestClassifier(n_estimators=300, random_state=SEED, n_jobs=1),
                      {"clf__max_depth": [1, 2], "clf__min_samples_leaf": [1, 3]}),
}
CTX = {}


def build_ctx(task):
    T, X, R, genes = make_samples(task)
    gi = {g: i for i, g in enumerate(genes)}
    set_idx = {s: np.array([gi[g] for g in gsets.gene[gsets.set == s] if g in gi]) for s in SETS}
    CTX[task] = dict(T=T, X=X, R=R, genes=genes, set_idx=set_idx,
                     R1=np.column_stack([singscore_up(R, set_idx[s]) for s in SETS]))  # data-independent: computed once


def fit_eval(task, rep, model, train_mask, test_mask, label, y_override=None):
    T = CTX[task]["T"]
    tr, te = np.where(train_mask)[0], np.where(test_mask)[0]
    ytr = T.y.values[tr] if y_override is None else y_override
    base, grid = MODELS[model]
    pipe = Pipeline([("feat", FeatureBuilder(rep, task)), ("sc", StandardScaler()), ("clf", clone(base))])
    best = {}
    if grid:
        cv = RepeatedStratifiedKFold(n_splits=4, n_repeats=3, random_state=SEED)
        gs = GridSearchCV(pipe, grid, cv=cv, scoring="neg_brier_score", n_jobs=1, error_score="raise")
        gs.fit(tr.reshape(-1, 1), ytr)
        pipe, best = gs.best_estimator_, gs.best_params_
    else:
        pipe.fit(tr.reshape(-1, 1), ytr)
    p = pipe.predict_proba(te.reshape(-1, 1))[:, 1]
    s = pipe.decision_function(te.reshape(-1, 1)) if hasattr(pipe, "decision_function") and model != "majority" else p
    return pipe, best, p, np.asarray(s, float)


def metrics(y, p, s, boot=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    both = len(set(y)) == 2
    out = dict(n=len(y), n_pos=int(y.sum()), auc=roc_auc_score(y, s) if both and np.std(s) > 0 else 0.5,
               bal_acc=balanced_accuracy_score(y, (p >= .5).astype(int)), sens=float(((p >= .5) & (y == 1)).sum() / max(1, y.sum())),
               spec=float(((p < .5) & (y == 0)).sum() / max(1, (y == 0).sum())), brier=brier_score_loss(y, p),
               logloss=log_loss(y, np.clip(p, 1e-6, 1 - 1e-6), labels=[0, 1]))
    if both:
        pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
        aucs = []
        for _ in range(boot):
            i = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
            aucs.append(roc_auc_score(y[i], s[i]) if np.std(s[i]) > 0 else 0.5)
        out["auc_ci_low"], out["auc_ci_high"] = np.quantile(aucs, [.025, .975])
    return out


def split_masks(T, name):
    tr316 = (T.study == "GSE316023") & (T.contrast == "MASH_primary")
    tr316b = (T.study == "GSE316023")
    s246 = T.study == "GSE246328"
    return {"S1_train316TA_test246": (tr316, s246), "S2_train246_test316TA": (s246, tr316),
            "S1b_train316TA+TN_test246": (tr316b, s246)}[name]


def one(task, split, rep, model):
    T = CTX[task]["T"]
    trm, tem = split_masks(T, split)
    pipe, best, p, s = fit_eval(task, rep, model, trm.values, tem.values, split)
    y = T.y.values[tem.values]
    m = metrics(y, p, s)
    lib = T.libsize.values[tem.values]
    m["rho_pred_libsize"] = float(stats.spearmanr(s, lib).statistic) if np.std(s) > 0 else np.nan
    feat = pipe.named_steps["feat"]
    sig = dict(n_up=len(getattr(feat, "up_", [])), n_dn=len(getattr(feat, "dn_", [])))
    preds = pd.DataFrame({"split": split, "rep": rep, "model": model, "key": T.key.values[tem.values], "y": y, "p": p, "score": s})
    return dict(task=task, split=split, rep=rep, model=model, best_params=json.dumps(best), **m, **sig), preds


def perm_one(i, split):
    T = CTX["mash"]["T"]
    trm, tem = split_masks(T, split)
    rng = np.random.default_rng(1000 + i)
    yperm = rng.permutation(T.y.values[trm.values])
    if len(set(yperm)) < 2:
        return np.nan
    _, _, p, s = fit_eval("mash", "R2", "logistic_default", trm.values, tem.values, split, y_override=yperm)
    y = T.y.values[tem.values]
    return roc_auc_score(y, s)


def main():
    build_ctx("mash")
    T = CTX["mash"]["T"]
    T.to_csv(TAB / "sample_manifest.csv", index=False)
    pd.DataFrame([dict(feature=f"R1_{s}", source="Phase 3 gene_sets.csv (MSigDB 2026.1.Mm / custom core)", n_genes_in_universe=len(CTX["mash"]["set_idx"][s]),
                       data_dependent=False) for s in SETS] +
                 [dict(feature="R0_signed_program_score", source="Welch t within TRAINING mice of each split (top 50 up/down, |diff|>=1)",
                       n_genes_in_universe=100, data_dependent=True),
                  dict(feature="PCA_PC1-3", source="PCA on per-sample z-scored log2CPM, fitted on training mice only", n_genes_in_universe=len(CTX["mash"]["genes"]),
                       data_dependent=True)]).to_csv(TAB / "feature_definitions.csv", index=False)
    splits = ["S1_train316TA_test246", "S2_train246_test316TA", "S1b_train316TA+TN_test246"]
    rows = []
    for sp in splits:
        trm, tem = split_masks(T, sp)
        rows.append(dict(split=sp, train_study=",".join(sorted(set(T.study[trm]))), n_train=int(trm.sum()), n_train_pos=int(T.y[trm].sum()),
                         test_study=",".join(sorted(set(T.study[tem]))), n_test=int(tem.sum()), n_test_pos=int(T.y[tem].sum()),
                         overlap_gsm=len(set(T.gsm[trm]) & set(T.gsm[tem]))))
    pd.DataFrame(rows).to_csv(TAB / "train_test_study_splits.csv", index=False)
    jobs = [("mash", sp, rep, mo) for sp in splits for rep in ("R0", "R1", "R2", "PCA") for mo in MODELS
            if not (mo == "R0_rule" and rep != "R0")]
    jobs.sort(key=lambda j: {"random_forest": 0, "tree": 1, "gaussian_nb": 2, "l2_logistic": 3}.get(j[3], 9))
    print(f"{len(jobs)} fits on {N_JOBS} workers; genes={len(CTX['mash']['genes'])}", flush=True)
    out = Parallel(n_jobs=N_JOBS, backend="loky", verbose=5)(delayed(one)(*j) for j in jobs)
    M = pd.DataFrame([o[0] for o in out])
    M.to_csv(RES / "external_metrics.csv", index=False)
    pd.concat([o[1] for o in out]).to_csv(RES / "external_predictions.csv", index=False)
    pj = [(i, sp) for sp in splits[:2] for i in range(200)]
    pr = Parallel(n_jobs=N_JOBS, verbose=2, batch_size=8)(delayed(perm_one)(i, sp) for i, sp in pj)
    perms = {sp: [r for (i, s_), r in zip(pj, pr) if s_ == sp] for sp in splits[:2]}
    pd.DataFrame(perms).to_csv(RES / "permutation_null_auc_R2_logistic.csv", index=False)
    # validity
    prim = M[M.split.isin(splits[:2])]
    v = prim.groupby(["rep", "model"]).apply(lambda d: bool(((d.auc >= .8) & (d.auc_ci_low > .5) & (d.brier < .25)).all() and len(d) == 2),
                                             include_groups=False).rename("valid").reset_index()
    v.to_csv(RES / "validity.csv", index=False)
    print(prim.pivot_table(index=["rep", "model"], columns="split", values=["auc", "brier"]).round(3).to_string())
    print(v[v.valid].to_string(index=False))
    json.dump({sp: dict(mean=float(np.nanmean(x)), q95=float(np.nanquantile(x, .95))) for sp, x in perms.items()},
              open(RES / "permutation_summary.json", "w"), indent=2)
    return M, v


if __name__ == "__main__":
    main()
    sys.exit(0)
