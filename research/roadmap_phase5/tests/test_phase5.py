"""Phase 5 tests. Run: python -m pytest research/roadmap_phase5/tests -q
Integrity of earlier phases is checked only where those directories exist; an archive without them SKIPS (not fails)."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

P5 = Path(__file__).resolve().parents[1]
RS = P5.parent
R, F = P5 / "results", P5 / "frozen"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def test_plan_locked_and_inputs_unchanged():
    lock = json.loads((R / "plan_lock.json").read_text())
    assert sha(P5 / "locked_plan.md") == lock["sha256"]["roadmap_phase5\\locked_plan.md"] or sha(P5 / "locked_plan.md") == lock["sha256"].get("roadmap_phase5/locked_plan.md")
    for k, h in lock["sha256"].items():
        p = Path(k) if Path(k).is_absolute() else RS / k
        if not p.exists():
            pytest.skip(f"{p} not present in this copy")
        assert sha(p) == h, k


def test_biological_sample_uniqueness_and_shared_controls():
    m = pd.read_csv(R / "external_sample_manifest.csv")
    assert m.mouse_id.is_unique
    hong = m[m.study == "HONG"]
    assert len(hong) == 15 and (hong.group == "CON").sum() == 3  # 24 GSMs collapse to 15 mice
    g = pd.read_csv(R / "external_gsm_table.csv")
    assert len(g[g.study == "HONG"]) == 24
    assert all(len(a.split(";")) == 4 for a in hong[hong.group == "CON"].all_gsm_aliases)  # each control = 4 GSM aliases
    assert len(m[m.study == "GSE245069"]) == 16 and (m[m.study == "GSE245069"].group == "CON").sum() == 4
    assert len(m[m.study == "GSE277772"]) == 6


def test_no_overlap_between_independent_study_units():
    m = pd.read_csv(R / "external_sample_manifest.csv")
    gsm = m.all_gsm_aliases.str.split(";").explode()
    assert gsm.is_unique  # no GSM is used by two mice / two study units
    e = pd.read_csv(RS / "roadmap_phase4e/tables/sample_manifest.csv") if (RS / "roadmap_phase4e").exists() else None
    if e is not None:
        assert not set(gsm) & set(e.gsm.astype(str))  # external mice never entered model training or GSE294180


def test_output_rows_match_metadata():
    m = pd.read_csv(R / "external_sample_manifest.csv").set_index("mouse_id")
    s = pd.read_csv(P5 / "per_sample_scores.csv", index_col=0)
    assert set(s.index) == set(m.index)
    assert (s.group == m.loc[s.index, "group"]).all() and (s.study == m.loc[s.index, "study"]).all()


def test_contrast_direction_exposed_minus_control():
    s = pd.read_csv(P5 / "per_sample_scores.csv", index_col=0)
    r = pd.read_csv(P5 / "sample_level_external_results.csv")
    h = s[s.study == "HONG"]
    est = r[(r.study == "HONG") & (r.contrast == "PE_vs_CON") & (r.score == "FAST_singscore")].estimate.iloc[0]
    assert np.isclose(est, h[h.group == "PE"].FAST_singscore.mean() - h[h.group == "CON"].FAST_singscore.mean())
    de = pd.read_csv(R / "external_de/HONG__PE_vs_CON.csv.gz").set_index("gene")
    x = pd.read_csv("D:/mpliver_work/phase5/HONG_expression.csv.gz", index_col=0) if Path("D:/mpliver_work/phase5/HONG_expression.csv.gz").exists() else None
    if x is None:
        pytest.skip("prepared matrix not present")
    g = de[(de.baseMean > 1000)].index[:200]
    cpm = x / x.sum() * 1e6
    naive = np.log2(cpm.loc[g, [c for c in x if "_PE_" in c]].mean(1) + 1) - np.log2(cpm.loc[g, [c for c in x if "_CON_" in c]].mean(1) + 1)
    assert np.corrcoef(naive, de.loc[g, "log2FC"])[0, 1] > 0.9  # DESeq2 log2FC is PE minus CON


def test_gene_set_coverage():
    c = pd.read_csv(R / "score_gene_coverage.csv")
    e = c[c.score == "Phase4E_R0R2"]
    assert (e.R0_signature == 100).all()
    assert (e.background_genes / e.background_frozen >= 0.99).all()
    d = c[c.score.str.startswith("Phase4D")]
    assert (d.n_missing <= 10).all() and (d.background_genes / d.background_frozen >= 0.85).all()
    cc = pd.read_csv(R / "coverage_control_gse294180_rescored.csv")
    assert cc.rel_change.abs().max() < 0.10  # locked rule: reduced backgrounds change GSE294180 estimates by < 10%


def test_frozen_models_serialised_and_reproduce_phase4e():
    man = json.loads((F / "frozen_manifest.json").read_text())
    for fn, h in man["files"].items():
        assert sha(F / fn) == h, fn
    for rep, v in man["phase4e_reproduction"].items():
        assert v["pure_function_max_abs_diff_vs_phase4e"] < 1e-6
    for rep in ("R0", "R2"):
        p = json.loads((F / f"phase4e_{rep}_l2_logistic_params.json").read_text())
        assert p["C"] == 10 and len(p["up"]) == 50 and len(p["dn"]) == 50
        assert not any(k.startswith(("HONG", "G245", "G277")) for k in p["train_mice"])  # no external mouse in training


def test_no_external_label_leakage_in_code():
    for f in ("s02_freeze_models.py",):
        src = (P5 / "scripts" / f).read_text()
        assert "load_prepared" not in src and "external_de" not in src  # freezing never touches external data
    for f in ("s04_gene_level.py", "s05_sample_level.py"):
        src = (P5 / "scripts" / f).read_text()
        for bad in (".fit(", "GridSearchCV", "LogisticRegression(", "Ridge(", "LinearRegression("):
            assert bad not in src, (f, bad)


def test_reproducibility_seed_and_versions():
    lock = json.loads((R / "plan_lock.json").read_text())
    assert lock["seed"] == 20261122 and "sklearn" in lock["versions"]
    g = pd.read_csv(P5 / "gene_level_external_results.csv")
    p = g[(g.study == "HONG") & (g.contrast == "PE_vs_CON") & (g.model == "F8_ols")]
    assert (p.perm_n_labellings == 20).all()  # all C(6,3) relabellings enumerated


def test_earlier_phase4d_unchanged():
    p4d = RS / "roadmap_phase4d"
    if not p4d.exists():
        pytest.skip("Phase 4D folder not present in this copy")
    lock = json.loads((p4d / "results" / "plan_lock.json").read_text())
    for f, h in lock["sha256"].items():
        assert sha(p4d / f) == h
