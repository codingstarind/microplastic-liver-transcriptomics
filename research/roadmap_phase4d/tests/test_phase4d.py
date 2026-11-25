"""Phase 4D tests. Run: python -m pytest research/roadmap_phase4d/tests -q"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

P4D = Path(__file__).resolve().parents[1]
REPO = P4D.parents[1]
R = P4D / "results"


def test_plan_and_manifest_locked():
    lock = json.loads((R / "plan_lock.json").read_text())
    for f, h in lock["sha256"].items():
        assert hashlib.sha256((P4D / f).read_bytes()).hexdigest() == h, f


def test_signature_files_match_manifest_hashes():
    m = pd.read_csv(P4D / "frozen_signature_manifest.csv")
    for _, r in m.iterrows():
        assert hashlib.sha256((P4D / "signatures" / r.file).read_bytes()).hexdigest() == r.sha256


def test_signatures_built_without_pe_data():
    src = (P4D / "scripts" / "s01_build_signatures.py").read_text()
    for bad in ("GSE294180_", "PE_in_CD", "PE_in_MD", "results/deseq2", "provisional_deposited", "Diet_in_Veh", "raw_counts"):
        assert bad not in src, bad
    for f in (P4D / "signatures").glob("*_n*.csv"):
        s = pd.read_csv(f)
        assert set(s.direction) <= {"up", "down"} and s.gene.is_unique
        assert ((s.direction == "up") == (s.eff > 0)).all()  # direction comes from the reference effect sign


def test_signature_genes_in_fixed_universe():
    u = pd.read_csv(REPO / "research/roadmap_phase2/metadata/phase3_gene_universe.csv", index_col=0)
    main = set(u.index[u.universe_main])
    for f in (P4D / "signatures").glob("*_main_n*.csv"):
        assert set(pd.read_csv(f).gene) <= main


def _load_s02_funcs():
    # import singscore/contrast helpers without executing the full pipeline twice
    src = (P4D / "scripts" / "s02_score_and_test.py").read_text()
    ns = {}
    start = src.index("def singscore")
    end = src.index("def ols_table")
    exec("import numpy as np\n" + src[start:end].replace("MD ==", "ns_MD ==").replace("PE ==", "ns_PE =="), ns)
    return ns


def test_singscore_formula_extremes_and_determinism():
    ns = _load_s02_funcs()
    N = 100
    ranks = np.tile(np.arange(1, N + 1)[:, None], (1, 19)).astype(float)
    up, dn = np.arange(90, 100), np.arange(0, 10)  # up genes highest, down genes lowest -> max score 1
    s = ns["singscore"](ranks, None, up, dn)
    assert np.allclose(s, 1.0)
    s2 = ns["singscore"](ranks, None, dn, up)  # reversed -> minimum -1
    assert np.allclose(s2, -1.0)
    assert np.allclose(ns["singscore"](ranks, None, up, dn), s)


def test_contrast_direction_and_sample_alignment():
    sc = pd.read_csv(R / "per_mouse_scores.csv", index_col=0)
    sm = pd.read_csv(REPO / "metadata/sample_metadata.csv").set_index("counts_column")
    assert set(sc.index) == set(sm.index) and len(sc) == 19
    assert (sc.group == sm.loc[sc.index, "group"]).all()
    t = pd.read_csv(R / "contrast_estimates_primary.csv")
    for sig in ("MASH", "FAST"):
        g = sc.groupby("group")[sig].mean()
        est = t[(t.signature == sig)].set_index("contrast").estimate
        assert np.isclose(est["PE_in_CD"], g["CD+PE"] - g["CD+Veh"])
        assert np.isclose(est["PE_in_MD"], g["MD+PE"] - g["MD+Veh"])
        assert np.isclose(est["Diet_in_Veh"], g["MD+Veh"] - g["CD+Veh"])


def test_outputs_complete():
    for f in ("contrast_estimates_primary.csv", "contrast_estimates_sensitivity.csv", "loo_contrasts.csv", "permutation_exact.csv",
              "null_gene_sets.csv", "qc_associations.csv", "per_mouse_scores.csv", "run_summary.json"):
        assert (R / f).exists(), f
    p = pd.read_csv(R / "permutation_exact.csv")
    assert set(p.n_perm) == {126, 252, 31752}
    assert len(pd.read_csv(R / "loo_contrasts.csv")) == 19 * 3 * 4
