"""Phase 4B integrity/leakage tests. Run: python -m pytest research/roadmap_phase4b/tests -q"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

P4B = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P4B / "scripts"))
from p4b_features import A_COLS, B_COLS, C_COLS, F8, RAW_COLS, REPS, RepBuilder, raw_table  # noqa: E402

sys.path.insert(0, str(P4B.parent / "roadmap_phase4" / "scripts"))
from p4_data import feature_table  # noqa: E402


def test_protocol_lock_unchanged():
    lock = json.loads((P4B / "results" / "protocol_lock.json").read_text())
    for f, h in lock["sha256"].items():
        assert hashlib.sha256((P4B / f).read_bytes()).hexdigest() == h, f


def test_rows_identical_to_phase4():
    assert raw_table().index.equals(feature_table("primary").index)


def test_no_pe_or_identifier_columns():
    allowed = ("MASH_", "FAST_", "PPARA_", "CIRC_", "B_HALLMARK_", "B_LIPOGENESIS", "B_CHOL", "B_SANDERSON",
               "FASTclean_", "PPARA5_", "PPARA10_")
    assert all(c.startswith(allowed) for c in RAW_COLS)
    assert not any(t in c.lower() for c in RAW_COLS for t in ("padj", "pvalue", "basemean", "pe_in"))


def test_rep_definitions_and_ablations():
    assert len(REPS) == 10
    full = set(REPS["FULL"])
    assert set(REPS["FULL-A"]) == full - set(A_COLS)
    assert set(REPS["FULL-B"]) == full - set(B_COLS)
    assert set(REPS["FULL-C"]) == full - set(C_COLS)
    assert REPS["F8"] == F8


def test_repbuilder_uses_training_rows_only():
    X = raw_table().values
    tr = np.arange(len(X)) < 6000
    a = RepBuilder("FULL").fit(X[tr])
    X2 = X.copy()
    X2[~tr] *= 50  # perturb test rows
    b = RepBuilder("FULL").fit(X2[tr])
    assert a.tau2_ == b.tau2_ and a.mu_ == b.mu_ and a.sd_ == b.sd_


def test_shrinkage_bounded_and_sign_preserving():
    X = raw_table().values
    rb = RepBuilder("A4").fit(X)
    A = rb.transform(X)
    E = raw_table()[[c[2:].replace("_shr", "_eff") for c in A_COLS]].values
    assert (np.abs(A) <= np.abs(E) + 1e-12).all()
    assert (np.sign(A[np.abs(A) > 0]) == np.sign(E[np.abs(A) > 0])).all()


def test_shuffled_b_keeps_set_sizes_and_breaks_link():
    a, b = raw_table()[B_COLS], raw_table(shuffle_b=True)[B_COLS]
    assert (a.sum() == b.sum()).all()
    assert (a.values != b.values).any()


def test_results_reuse_identical_splits_and_reproduce_phase4():
    f = P4B / "results" / "metrics_main.csv"
    if not f.exists():
        pytest.skip("not run yet")
    m = pd.read_csv(f)
    assert (m.groupby(["target", "repeat", "fold"]).n_test.nunique() == 1).all()
    assert (m.groupby(["rep", "model", "target"]).size() == 15).all()
    p4 = pd.read_csv(P4B.parent / "roadmap_phase4" / "results" / "taskA" / "metrics_primary_positive.csv")
    p4["model"] = p4["model"].fillna("null")
    p4 = p4[p4.run.str.startswith("primary__")].assign(target=lambda d: d.run.str.split("__").str[1])
    p4["model"] = p4.model.replace({"ols_all": "ols"})
    j = m[m.rep == "F8"].merge(p4, on=["target", "model", "repeat", "fold"], suffixes=("", "_p4"))
    assert len(j) == 150
    exact = j[j.model != "catboost"]
    assert np.allclose(exact.rmse, exact.rmse_p4, atol=1e-9)  # ols/ridge/elasticnet/hgb reproduce exactly
    cb = j[j.model == "catboost"]
    # documented: the added scaler shifts CatBoost float quantisation borders (9/30 splits differ, max 0.0016;
    # one near-tied inner hyperparameter choice flipped)
    d = (cb.rmse - cb.rmse_p4).abs()
    assert d.max() < 0.002 and d.mean() < 1e-3
    assert (cb.best_params != cb.best_params_p4).sum() <= 1
