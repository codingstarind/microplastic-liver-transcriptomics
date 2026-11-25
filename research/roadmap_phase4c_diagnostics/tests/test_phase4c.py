"""Phase 4C audit tests. Run: python -m pytest research/roadmap_phase4c_diagnostics/tests -q"""
import hashlib
import json
from pathlib import Path

import pandas as pd

P4C = Path(__file__).resolve().parents[1]
RS = P4C.parent
T = P4C / "tables"


def test_spec_lock_unchanged():
    lock = json.loads((T / "spec_lock.json").read_text())
    assert hashlib.sha256((P4C / lock["file"]).read_bytes()).hexdigest() == lock["sha256"]


def test_fix_experiment_reproduces_phase4b_exactly():
    r = json.loads((T / "fix_experiment_reproducibility.json").read_text())
    assert r["n_overlap"] == 120 and r["max_abs_rmse_diff_vs_phase4b"] < 1e-9


def test_fix_experiment_scope_is_the_predeclared_one():
    m = pd.read_csv(T / "fix_experiment_metrics_by_fold.csv")
    assert set(zip(m.rep, m.model)) == {("F8", "ols"), ("F8+C", "ols"), ("F8+C", "ridge"), ("FULL-A", "ridge")}
    assert (m.groupby(["target", "rep", "model"]).size() == 15).all()
    assert set(m.alpha.dropna().unique()) <= {0.01, 0.1, 1, 10, 100, 1000}


def test_oof_each_gene_once_per_repeat():
    o = pd.read_csv(T / "fix_experiment_oof_predictions.csv.gz")
    c = o.groupby(["target", "rep", "model", "repeat"]).gene.agg(["size", "nunique"])
    assert (c["size"] == c["nunique"]).all() and (c["nunique"] == 11758).all()


def test_collinearity_values():
    s = pd.read_csv(T / "collinearity_summary.csv").set_index("matrix")
    assert s.loc["F8", "max_vif"] < 2 and abs(s.loc["F8", "max_abs_r"] - 0.663) < 0.01
    assert s.loc["FULL", "max_vif"] > 10 and s.loc["F8+C", "max_vif"] < 5
    p = pd.read_csv(T / "correlated_pairs_absr_ge_0.5.csv")
    top = p[p.matrix == "FULL"].iloc[0]
    assert {top.a, top.b} == {"FAST_eff", "A_FAST_shr"} and top.r > 0.99


def test_earlier_phase_files_unchanged():
    for line in (RS / "roadmap_phase4b" / "results" / "input_checksums.sha256").read_text().splitlines():
        h, f = line.split(" *")
        p = (RS / f).resolve()
        assert hashlib.sha256(p.read_bytes()).hexdigest() == h, f


def test_diagnostic_only_bins_labelled():
    s = pd.read_csv(T / "subgroup_errors.csv")
    assert s[s.group.str.startswith("absy_")].grouping.str.contains("DIAGNOSTIC-ONLY").all()
