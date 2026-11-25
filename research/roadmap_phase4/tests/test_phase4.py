"""Phase 4 leakage and integrity tests. Run: python -m pytest research/roadmap_phase4/tests -q"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.model_selection import GroupKFold

P4 = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P4 / "scripts"))
from p4_data import (  # noqa: E402
    FEATURES_F8,
    FORBIDDEN,
    VARIANTS,
    dl_tau2,
    feature_names,
    feature_table,
    outer_splits,
    target,
)
from p4_models import BASELINES, TUNED, make  # noqa: E402


def test_protocol_lock_unchanged():
    lock = json.loads((P4 / "results" / "protocol_lock.json").read_text())
    for f, h in lock["sha256"].items():
        assert hashlib.sha256((P4 / f).read_bytes()).hexdigest() == h, f


@pytest.mark.parametrize("grouping", ["module", "chromosome"])
def test_outer_splits_disjoint_grouped_and_complete(grouping):
    X = feature_table("primary")
    seen = {}
    for r, k, tr, te, grp in outer_splits(X.index, grouping):
        assert not (tr & te).any() and (tr | te).all()
        assert not set(grp[tr]) & set(grp[te]), "group straddles train/test"
        seen.setdefault(r, np.zeros(len(X), int))
        seen[r] += te
        # inner folds: subsets of outer training and group-disjoint
        idx = np.where(tr)[0]
        for itr, iva in GroupKFold(4).split(idx, groups=grp[tr]):
            assert not set(grp[tr][itr]) & set(grp[tr][iva])
            assert set(idx[iva]) <= set(idx)
    for r, cnt in seen.items():
        assert (cnt == 1).all(), "each gene must be tested exactly once per repeat"


def test_splits_identical_across_calls():
    X = feature_table("primary")
    a = [te.copy() for *_, te, _ in outer_splits(X.index)]
    b = [te.copy() for *_, te, _ in outer_splits(X.index)]
    assert all((x == y).all() for x, y in zip(a, b))


@pytest.mark.parametrize("variant", list(VARIANTS))
def test_no_forbidden_or_target_features(variant):
    feats = feature_names(variant)
    assert not set(feats) & FORBIDDEN
    assert all(f.split("_")[0] in ("MASH", "FAST", "PPARA", "CIRC") for f in feats)
    X = feature_table(variant)
    t = target("PE_in_CD")
    j = X.join(t, how="inner")
    for f in feats:  # no feature is (numerically) the target
        assert abs(np.corrcoef(j[f], j["y"].fillna(0))[0, 1]) < 0.9


def test_tau2_uses_training_rows_only():
    rng = np.random.default_rng(0)
    y, s = rng.normal(size=500), np.full(500, 0.1)
    t_tr = dl_tau2(y[:400], s[:400])
    y2 = y.copy()
    y2[400:] += 100  # changing test rows must not change training tau2
    assert dl_tau2(y2[:400], s[:400]) == t_tr


def test_scaler_fitted_on_training_only():
    X = feature_table("primary")
    est, _ = make("ols_all", FEATURES_F8)
    tr = np.arange(len(X)) < 5000
    Xv = X[FEATURES_F8].values.copy()
    est.fit(Xv[tr], np.zeros(tr.sum()) + Xv[tr, 0])
    mu = est.named_steps["scale"].mean_.copy()
    Xv[~tr] += 1000
    est.fit(Xv[tr], np.zeros(tr.sum()) + Xv[tr, 0])
    assert np.allclose(mu, est.named_steps["scale"].mean_)


def test_single_axis_models_use_one_effect_column():
    for m in [b for b in BASELINES if b.startswith("ols_") and b != "ols_all"]:
        est, _ = make(m, FEATURES_F8)
        X = np.random.default_rng(1).normal(size=(50, 8))
        est.fit(X, X[:, FEATURES_F8.index(m.split("_")[1] + "_eff")])
        assert est.named_steps["model"].coef_.shape == (1,)


def test_model_suite_matches_protocol():
    grid = pd.read_csv(P4 / "model_search_spaces.csv")
    a = set(grid[grid.task == "A"].model)
    for m in TUNED:
        assert m in a


def test_results_complete_if_present():
    f = P4 / "results" / "taskA" / "metrics_primary_positive.csv"
    if not f.exists():
        pytest.skip("not run yet")
    m = pd.read_csv(f)
    m["model"] = m["model"].fillna("null")
    # every model scored on exactly the same held-out genes per split
    n = m.groupby(["run", "repeat", "fold"]).n_test.nunique()
    assert (n == 1).all()
    assert m.groupby("run").model.nunique().min() == len(BASELINES + TUNED)
    assert (m.groupby(["run", "model"]).size() == 15).all()
