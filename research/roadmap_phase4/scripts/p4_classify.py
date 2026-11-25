"""Task B (DEMONSTRATION ONLY): GAN/MASH-diet vs control mouse classifier.
Train: GSE316023 wk23 TA+TN (16 GAN, 16 chow). Internal: repeated stratified 4-fold (5 repeats) with inner 3-fold tuning.
External: GSE246328 Vehicle 24 w vs Chow 24 w (7 v 8). Descriptive: GSE294180 vehicle mice (CD 4, MD 5).
Features: per-sample log2 CPM over common genes, then per-sample z-score (sample-wise, no leakage);
top-k variance filter + scaler fitted inside the pipeline. Outputs: results/taskB/."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from p4_data import DROOT, P2, REPO, RES, _csv, universe
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold, StratifiedKFold
from sklearn.naive_bayes import GaussianNB
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeClassifier

OUT = RES / "taskB"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261122


class TopVar(BaseEstimator, TransformerMixin):
    def __init__(self, k=200):
        self.k = k

    def fit(self, X, y=None):
        self.idx_ = np.argsort(X.var(axis=0))[::-1][: self.k]
        return self

    def transform(self, X):
        return X[:, self.idx_]


def models():
    tv = {"sel__k": [200, 1000]}
    P = lambda m: Pipeline([("sel", TopVar()), ("sc", StandardScaler()), ("model", m)])  # noqa: E731
    return {
        "dummy": (Pipeline([("model", DummyClassifier(strategy="most_frequent"))]), {}),
        "logistic_l2": (P(LogisticRegression(max_iter=5000)), {**tv, "model__C": [0.01, 0.1, 1]}),
        "gaussian_nb": (P(GaussianNB()), {**tv, "model__var_smoothing": [1e-9, 1e-6, 1e-3]}),
        "tree_gini": (P(DecisionTreeClassifier(criterion="gini", random_state=SEED)), {**tv, "model__max_depth": [1, 2, 3]}),
        "tree_entropy": (P(DecisionTreeClassifier(criterion="entropy", random_state=SEED)), {**tv, "model__max_depth": [1, 2, 3]}),
        "rf": (P(RandomForestClassifier(n_estimators=300, random_state=SEED, n_jobs=1)),
               {**tv, "model__max_features": ["sqrt", 0.1], "model__min_samples_leaf": [1, 3]}),
        # Deviation D2 (post-hoc technical fix): default min_samples_leaf=20 cannot split with 24 training mice,
        # giving constant 0.5 predictions; set to 3 (no other change). Original degenerate run kept in results/taskB/original_hgb_min_leaf20/.
        "hgb": (P(HistGradientBoostingClassifier(max_iter=100, max_leaf_nodes=4, min_samples_leaf=3, random_state=SEED)),
                {**tv, "model__learning_rate": [0.05, 0.2]}),
    }


def load():
    el = pd.read_csv(P2 / "metadata" / "phase3_eligible_samples.csv", dtype={"sample_id": str})
    U = universe("main")

    def mat(acc, ids):
        x = _csv(DROOT / acc / "prepared" / "expression_ready.csv.gz", index_col=0)
        x.columns = x.columns.astype(str)
        return x[ids]

    tr = el[el.contrast_id.isin(["MASH_primary", "MASH_sens_TN"])]
    ex = el[el.contrast_id == "MASH_secondary"]
    A = mat("GSE316023", list(tr.sample_id))
    B = mat("GSE246328", list(ex.sample_id))
    raw = _csv(REPO / "data" / "raw" / "GSE294180_MP_raw_counts_gene_level.csv.gz", index_col=0)
    sm = pd.read_csv(REPO / "metadata" / "sample_metadata.csv")
    veh = sm[sm.microplastic == "Veh"]
    C = raw[list(veh.counts_column)]
    genes = sorted(set(A.index) & set(B.index) & set(C.index) & U)

    def prep(M):
        M = M.loc[genes].astype(float)
        L = np.log2(M / M.sum() * 1e6 + 1)
        return ((L - L.mean()) / L.std()).T.values  # samples x genes, per-sample z

    ya = (tr.arm == "case").astype(int).values
    yb = (ex.arm == "case").astype(int).values
    return prep(A), ya, tr.sample_id.tolist(), prep(B), yb, ex.sample_id.tolist(), prep(C), list(veh.diet), len(genes)


def metrics(y, p):
    yhat = (p >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, yhat, labels=[0, 1]).ravel()
    both = len(set(y)) == 2
    return dict(balanced_accuracy=balanced_accuracy_score(y, yhat), roc_auc=roc_auc_score(y, p) if both else np.nan,
                pr_auc=average_precision_score(y, p) if both else np.nan, sensitivity=tp / (tp + fn) if tp + fn else np.nan,
                specificity=tn / (tn + fp) if tn + fp else np.nan, brier=brier_score_loss(y, p), tn=tn, fp=fp, fn=fn, tp=tp)


def fit_tuned(name, X, y):
    est, grid = models()[name]
    if not grid:
        return est.fit(X, y), {}
    gs = GridSearchCV(est, grid, cv=StratifiedKFold(3, shuffle=True, random_state=SEED), scoring="balanced_accuracy", n_jobs=1)
    gs.fit(X, y)
    return gs.best_estimator_, gs.best_params_


def outer_task(name, rep, k, tr, te, X, y):
    est, bp = fit_tuned(name, X[tr], y[tr])
    p = est.predict_proba(X[te])[:, 1]
    return dict(model=name, repeat=rep, fold=k, **metrics(y[te], p), best_params=json.dumps(bp, default=str))


def main(only=None):
    Xa, ya, ida, Xb, yb, idb, Xc, dietc, ng = load()
    names = only or list(models())
    sfx = "_" + "_".join(only) + "_fix" if only else ""
    cv = RepeatedStratifiedKFold(n_splits=4, n_repeats=5, random_state=SEED)
    splits = [(i // 4, i % 4, tr, te) for i, (tr, te) in enumerate(cv.split(Xa, ya))]
    rows = Parallel(n_jobs=2)(delayed(outer_task)(n, r, k, tr, te, Xa, ya) for n in names for r, k, tr, te in splits)
    pd.DataFrame(rows).to_csv(OUT / f"internal_cv_folds{sfx}.csv", index=False)
    ext, desc = [], []
    for n in names:
        est, bp = fit_tuned(n, Xa, ya)
        pb = est.predict_proba(Xb)[:, 1]
        ext.append(dict(model=n, **metrics(yb, pb), best_params=json.dumps(bp, default=str)))
        pc = est.predict_proba(Xc)[:, 1]
        desc += [dict(model=n, diet=d, p_GAN_state=float(v)) for d, v in zip(dietc, pc)]
    pd.DataFrame(ext).to_csv(OUT / f"external_GSE246328{sfx}.csv", index=False)
    pd.DataFrame(desc).to_csv(OUT / f"descriptive_GSE294180_vehicle{sfx}.csv", index=False)
    summ = pd.DataFrame(rows).groupby("model")[["balanced_accuracy", "roc_auc", "brier"]].agg(["mean", "std"]).round(3)
    summ.to_csv(OUT / f"internal_cv_summary{sfx}.csv")
    json.dump(dict(n_train=len(ya), n_external=len(yb), n_descriptive=len(dietc), n_genes=ng,
                   train_ids=ida, external_ids=idb), open(OUT / f"run_info{sfx}.json", "w"), indent=2)
    print(summ.to_string())
    print(pd.DataFrame(ext)[["model", "balanced_accuracy", "roc_auc", "sensitivity", "specificity", "brier"]].round(3).to_string())
    print(pd.DataFrame(desc).groupby(["model", "diet"]).p_GAN_state.mean().unstack().round(2).to_string())


if __name__ == "__main__":
    import sys
    main(sys.argv[1:] or None)
