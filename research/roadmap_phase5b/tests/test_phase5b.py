"""Phase 5B tests. Run: python -m pytest research/roadmap_phase5b/tests -q"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PB = Path(__file__).resolve().parents[1]
RS, REPO = PB.parent, PB.parents[1]
T, R = PB / "tables", PB / "results"
sys.path.insert(0, str(RS / "roadmap_phase4c_diagnostics" / "scripts"))


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def test_protocol_locked_before_scoring():
    lock = json.loads((R / "protocol_lock.json").read_text())
    paths = {"03_resampling_and_nested_cv_protocol.md": PB, "04_targeted_interventions.md": PB, "b04_intervention_targets.R": PB / "scripts",
             "b05_interventions.py": PB / "scripts", "b00_common.py": PB / "scripts", "plan_lock.json": RS / "roadmap_phase5/results",
             "frozen_manifest.json": RS / "roadmap_phase5/frozen"}
    for f, h in lock["sha256"].items():
        assert sha(paths[f] / f) == h, f


def test_count_unit_mismatch_detected_and_harmonised():
    u = pd.read_csv(T / "count_unit_per_sample.csv").groupby("study").counts_per_read_pair.median()
    assert 1.4 < u["HONG"] / u["GSE294180"] < 2.6  # Hong ~2x fragment-equivalent
    assert 0.6 < u["HONG"] * 0.5 < 1.0  # a-priori factor brings it into the fragment range
    h = pd.read_csv(T / "i1_unit_harmonisation_effect.csv").set_index(["unit_factor", "contrast"])
    assert h.loc[(0.5, "PE_vs_CON"), "n_padj05"] < h.loc[(1.0, "PE_vs_CON"), "n_padj05"]


def test_normalisation_equivalence_sign_orientation():
    c = pd.read_csv(T / "equivalence_checks.csv")
    sign = c[c.check.str.contains("exposed minus control|PE minus vehicle")]
    assert len(sign) == 3 and sign.passed.all()


def test_sample_id_overlap_and_shared_control_grouping():
    m = pd.read_csv(RS / "roadmap_phase5/results/external_sample_manifest.csv")
    gsm = m.all_gsm_aliases.str.split(";").explode()
    assert gsm.is_unique
    g294 = set(pd.read_csv(REPO / "metadata/sample_metadata.csv").geo_accession)
    assert not (set(gsm) & g294)
    hong = m[m.study == "HONG"]
    assert (hong.group == "CON").sum() == 3 and len(hong) == 15
    src = (PB / "scripts/b05_interventions.py").read_text()
    assert 'for g in ("CON", "PE", "PS", "PP", "PB")' in src  # bootstrap resamples the 3 shared controls ONCE per replicate


def test_feature_column_order_and_reprediction():
    c = pd.read_csv(T / "equivalence_checks.csv")
    assert c[c.check.str.contains("feature columns|gene order|re-predicts|order sensitivity")].passed.all()


def test_no_train_test_leakage_in_interventions():
    from c00_common import outer_splits, raw_table
    X0 = raw_table()
    for r, k, tr, te, grp in outer_splits(X0.index):
        assert not set(grp[tr]) & set(grp[te])  # whole modules held out
    src = (PB / "scripts/b05_interventions.py").read_text()
    assert src.index("external_de") > src.index("final models (modal inner choice)")  # external data read only after all fitting
    assert "GridSearchCV(p, grid, cv=GroupKFold(4)" in src  # inner tuning on training genes, grouped


def test_cv_reuses_phase4_folds_and_b0_matches_phase4c():
    cv = pd.read_csv(T / "cv_fold_metrics.csv")
    b0 = cv[cv.candidate == "B0_F8_OLS"].sort_values(["repeat", "fold"])
    fx = pd.read_csv(RS / "roadmap_phase4c_diagnostics/tables/fix_experiment_metrics_by_fold.csv")
    f = fx[(fx.target == "PE_in_CD") & (fx.rep == "F8") & (fx.model == "ols")].sort_values(["repeat", "fold"])
    assert len(b0) == 15
    assert abs(b0.r2_vs_train_mean_PE_in_CD.mean() - f.r2_oos.mean()) < 0.01


def test_mash_classifier_hard_rule_not_reported_as_brier():
    m = pd.read_csv(T / "mash_classifier_bootstrap.csv")
    assert (m[m.model == "R0_rule"].metric_prob.str.contains("misclassification")).all()
    assert (m[m.model != "R0_rule"].metric_prob == "Brier").all()


def test_rerun_determinism():
    f = T / "determinism_check.csv"
    if not f.exists():
        pytest.skip("run scripts/b06_determinism.py first")
    d = pd.read_csv(f)
    assert d.identical_shape.all() and (d.max_abs_numeric_diff.fillna(0) < 1e-9).all() and d.non_numeric_equal.all()


def test_earlier_phases_unchanged():
    for sub, lockf in (("roadmap_phase4d", "results/plan_lock.json"), ("roadmap_phase4e", "results/plan_lock.json")):
        p = RS / sub
        if not p.exists():
            pytest.skip(f"{sub} not present")
        lock = json.loads((p / lockf).read_text())
        for f, h in lock["sha256"].items():
            assert sha(p / f) == h
    lk = json.loads((RS / "roadmap_phase5/results/plan_lock.json").read_text())
    assert sha(RS / "roadmap_phase5/locked_plan.md") == lk["sha256"]["roadmap_phase5\\locked_plan.md"]
