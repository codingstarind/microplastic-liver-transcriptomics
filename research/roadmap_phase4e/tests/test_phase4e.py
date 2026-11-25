"""Phase 4E tests. Run: python -m pytest research/roadmap_phase4e/tests -q"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

P4E = Path(__file__).resolve().parents[1]
R, T = P4E / "results", P4E / "tables"


def test_plan_locked():
    lock = json.loads((R / "plan_lock.json").read_text())
    for f, h in lock["sha256"].items():
        assert hashlib.sha256((P4E / f).read_bytes()).hexdigest() == h, f


def test_train_test_studies_disjoint_and_sized():
    s = pd.read_csv(T / "train_test_study_splits.csv").set_index("split")
    assert (s.overlap_gsm == 0).all()
    assert s.loc["S1_train316TA_test246", ["n_train", "n_train_pos", "n_test", "n_test_pos"]].tolist() == [16, 8, 15, 7]
    assert s.loc["S2_train246_test316TA", ["n_train", "n_train_pos", "n_test", "n_test_pos"]].tolist() == [15, 7, 16, 8]
    assert (s.train_study != s.test_study).all()


def test_sample_manifest_unique_and_no_pe_labels_used():
    m = pd.read_csv(T / "sample_manifest.csv")
    assert m.key.is_unique and not m.gsm.duplicated().any()
    assert (m[m.study == "GSE294180"].y == -1).all()  # target mice carry no label in training
    assert len(m[m.study == "GSE294180"]) == 19


def test_feature_definitions_data_dependence_declared():
    f = pd.read_csv(T / "feature_definitions.csv")
    assert f[f.feature.str.startswith("R1_")].data_dependent.eq(False).all() and len(f[f.feature.str.startswith("R1_")]) == 12
    assert f[f.feature.isin(["R0_signed_program_score", "PCA_PC1-3"])].data_dependent.all()


def test_signature_built_on_training_only_code():
    src = (P4E / "scripts" / "e01_mash_state_models.py").read_text()
    assert "welch_signature(Xg[idx]" in src  # signature uses only the rows passed to fit (training fold)
    assert "def fit(self, idx, y)" in src


def test_external_metrics_complete_and_valid_range():
    m = pd.read_csv(R / "external_metrics.csv")
    assert set(m.split) == {"S1_train316TA_test246", "S2_train246_test316TA", "S1b_train316TA+TN_test246"}
    assert m.auc.between(0, 1).all() and m.brier.between(0, 1).all()
    assert (m[m.model == "majority"].auc == 0.5).all()


def test_application_only_for_valid_models_and_contrast_directions():
    a = json.loads((R / "application_decision.json").read_text())
    v = pd.read_csv(R / "validity.csv")
    for rep, mo in a["valid_and_applied"]:
        assert bool(v[(v.rep == rep) & (v.model == mo)].valid.iloc[0])
    s = pd.read_csv(R / "gse294180_frozen_model_scores.csv", index_col=0)
    c = pd.read_csv(R / "gse294180_frozen_model_contrasts.csv")
    g = s.groupby("group")["R2_l2_logistic"].mean()
    est = c[c.score == "R2_l2_logistic"].set_index("contrast").estimate
    assert np.isclose(est["PE_in_MD"], g["MD+PE"] - g["MD+Veh"])
    assert np.isclose(est["Diet_in_Veh"], g["MD+Veh"] - g["CD+Veh"])


def test_earlier_phase4d_results_unchanged():
    lock = json.loads((P4E.parent / "roadmap_phase4d" / "results" / "plan_lock.json").read_text())
    for f, h in lock["sha256"].items():
        assert hashlib.sha256((P4E.parent / "roadmap_phase4d" / f).read_bytes()).hexdigest() == h
