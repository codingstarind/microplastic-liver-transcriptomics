"""Lightweight integrity tests for the Phase 2 prepared reference package.

Run:  python -m pytest research/roadmap_phase2/scripts/test_phase2_data.py -q
Tests that need the external matrices are skipped if DATA_ROOT is absent.
"""
from __future__ import annotations

import json

import pandas as pd
import pytest
from p2_common import DATA_ROOT, REPO_P2, file_hashes

DATASETS = ["GSE316023", "GSE246328", "GSE107787", "GSE117134", "GSE282998", "GSE225695"]
EXPECTED_N = {  # contrast role -> expected biological replicate count (verified from GEO)
    ("GSE316023", "MASH_primary:case"): 8, ("GSE316023", "MASH_primary:control"): 8,
    ("GSE316023", "MASH_sens_TN:case"): 8, ("GSE316023", "MASH_sens_TN:control"): 8,
    ("GSE246328", "MASH_secondary:case"): 7, ("GSE246328", "MASH_secondary:control"): 8,
    ("GSE107787", "FASTING_primary:case"): 18, ("GSE107787", "FASTING_primary:control"): 18,
    ("GSE117134", "CIRCADIAN_primary:timecourse"): 18,
    ("GSE282998", "PPARA_primary:case"): 12, ("GSE282998", "PPARA_primary:control"): 4,
    ("GSE282998", "PPARA_specificity_KO:case"): 12, ("GSE282998", "PPARA_specificity_KO:control"): 4,
    ("GSE225695", "FASTING_secondary:case"): 3, ("GSE225695", "FASTING_secondary:control"): 3,
}
VALID_ROLE_PREFIX = ("MASH_primary", "MASH_sens_TN", "MASH_secondary", "FASTING_primary", "FASTING_secondary",
                     "CIRCADIAN_primary", "CIRCADIAN_sens_male", "PPARA_primary", "PPARA_specificity_KO",
                     "PPARA_sens_outbred", "not_in_contrast")
HAVE_DATA = DATA_ROOT.exists()


def meta(acc):
    return pd.read_csv(REPO_P2 / "datasets" / acc / "sample_metadata.csv", dtype={"sample_id": str})


@pytest.mark.parametrize("acc", DATASETS)
def test_no_duplicate_sample_ids(acc):
    m = meta(acc)
    assert not m["sample_id"].duplicated().any()
    assert not m["gsm"].duplicated().any()


@pytest.mark.parametrize("acc", DATASETS)
def test_metadata_matches_geo_soft(acc):
    """Every GSM in the GEO SOFT appears once in our metadata (no unexpected missing samples)."""
    soft = pd.read_csv(REPO_P2 / "metadata" / "geo_raw" / f"{acc}_samples_raw.csv")
    m = meta(acc)
    if acc == "GSE282998":  # only the GW7647 arm was prepared
        soft = soft[soft["ch_treatment"].str.contains("GW7647")]
    assert set(m["gsm"]) == set(soft["gsm"])


@pytest.mark.parametrize("acc", DATASETS)
def test_valid_contrast_labels(acc):
    for r in meta(acc)["contrast_role"].astype(str):
        for part in r.split(";"):
            assert part.startswith(VALID_ROLE_PREFIX), part
            if not part.startswith("not_in_contrast"):
                assert part.split(":")[1] in ("case", "control", "timecourse"), part


@pytest.mark.parametrize("key,n", EXPECTED_N.items())
def test_biological_replicate_counts_and_nonempty_contrasts(key, n):
    acc, role = key
    m = meta(acc)
    got = m["contrast_role"].astype(str).str.split(";").apply(lambda xs: role in xs).sum()
    assert got == n, f"{acc} {role}: {got} != {n}"
    assert got > 0


def test_master_metadata_one_row_per_biological_sample():
    ms = pd.read_csv(REPO_P2 / "metadata" / "sample_metadata_master.csv")
    assert not ms["gsm"].duplicated().any()
    assert set(ms["dataset"]) == set(DATASETS)
    assert len(ms) == sum(len(meta(a)) for a in DATASETS)


def test_no_sample_overlap_between_selected_studies():
    """Selected references must not share GSMs, SRA runs or identical read signatures."""
    runs = pd.read_csv(REPO_P2 / "metadata" / "run_metadata.csv")
    sel = runs[runs["dataset"].isin(DATASETS)].dropna(subset=["read_count"])
    assert not sel["run_accession"].duplicated().any()
    sig = sel["read_count"].astype(int).astype(str) + "_" + sel["base_count"].astype(int).astype(str)
    by = sel.assign(sig=sig).groupby("sig")["dataset"].nunique()
    assert (by == 1).all()
    ms = pd.read_csv(REPO_P2 / "metadata" / "sample_metadata_master.csv")
    assert not ms["gsm"].duplicated().any()


def test_known_overlaps_are_recorded():
    ov = pd.read_csv(REPO_P2 / "qc" / "overlap_checks.csv").set_index("check")
    assert ov.loc["GSE246328 library IDs also in GSE294629", "result"] == 54
    assert ov.loc["Hong-lab control count columns identical to GSE261893", "result"] == 3


@pytest.mark.parametrize("acc", DATASETS)
def test_gene_mapping_one_row_per_original_id(acc):
    g = pd.read_csv(REPO_P2 / "datasets" / acc / "gene_mapping.csv.gz")
    assert not g["original_id"].duplicated().any()
    kept = g[g["kept_in_ready"]]
    assert not kept["vm35_symbol"].duplicated().any() or acc in ("GSE246328", "GSE225695")  # summed many-to-one


@pytest.mark.skipif(not HAVE_DATA, reason="external data root not present")
@pytest.mark.parametrize("acc", DATASETS)
def test_matrix_columns_match_metadata_and_checksums(acc):
    man = json.loads((REPO_P2 / "datasets" / acc / "prepared_manifest.json").read_text())
    for f, info in man["files"].items():
        md5, _ = file_hashes(DATA_ROOT / acc / "prepared" / f)
        assert md5 == info["md5"], f"{acc}/{f} checksum changed"
    m = meta(acc)
    ready = pd.read_csv(DATA_ROOT / acc / "prepared" / "expression_ready.csv.gz", index_col=0, nrows=5)
    assert list(ready.columns.astype(str)) == list(m.loc[m["in_ready_matrix"], "sample_id"])
    orig = pd.read_csv(DATA_ROOT / acc / "prepared" / "expression_original.csv.gz", index_col=0, nrows=5)
    assert set(orig.columns.astype(str)) == set(m["sample_id"])


@pytest.mark.skipif(not HAVE_DATA, reason="external data root not present")
@pytest.mark.parametrize("acc", DATASETS)
def test_ready_matrix_gene_ids_unique_no_missing(acc):
    ready = pd.read_csv(DATA_ROOT / acc / "prepared" / "expression_ready.csv.gz", index_col=0)
    assert ready.index.is_unique
    assert not ready.isna().any().any()
    assert (ready.values >= 0).all()
    man = json.loads((REPO_P2 / "datasets" / acc / "prepared_manifest.json").read_text())
    if "counts" in man["value_type"]:
        assert (ready.values % 1 == 0).all(), "count-type ready matrix must be integer"


@pytest.mark.skipif(not HAVE_DATA, reason="external data root not present")
def test_download_manifest_checksums():
    man = pd.read_csv(REPO_P2 / "metadata" / "download_manifest.csv")
    for _, r in man.iterrows():
        p = DATA_ROOT / r["accession"] / r["file"]
        if p.exists():
            assert file_hashes(p)[0] == r["md5"], r["file"]


def test_inventory_parses_with_fixed_columns():
    inv = pd.read_csv(REPO_P2 / "metadata" / "dataset_master_inventory.csv")
    assert inv.shape[1] == 21 and len(inv) >= 10


@pytest.mark.skipif(not HAVE_DATA, reason="external data root not present")
def test_contrast_plan_ids_exist_in_ready_matrices():
    plan = pd.read_csv(REPO_P2 / "metadata" / "reference_contrast_plan.csv", dtype=str)
    for _, p in plan.iterrows():
        cols = pd.read_csv(DATA_ROOT / p["dataset"] / "prepared" / "expression_ready.csv.gz",
                           nrows=1, index_col=0).columns.astype(str)
        ids = [i for k in ("case_sample_ids", "control_sample_ids", "timecourse_sample_ids")
               for i in str(p[k]).split(";") if i and i != "nan"]
        assert ids, p["contrast_id"]
        assert all(i in cols for i in ids), p["contrast_id"]


def test_phase3_universe_excludes_muscle_and_is_nested():
    u = pd.read_csv(REPO_P2 / "metadata" / "phase3_gene_universe.csv", index_col=0)
    mus = set(pd.read_csv(REPO_P2 / "metadata" / "muscle_exclusion_genes.csv")["vm35_symbol"])
    assert {"Acta1", "Ckm", "Myh1", "Myh4", "Mylpf", "Tnnt3", "Myl1", "Actn3"} <= mus
    assert not (u["universe_main"] & u.index.isin(mus)).any()
    assert not (u["universe_sensitivity"] & ~u["universe_main"]).any()
    assert u.index.is_unique and u["universe_main"].sum() > 10000


def test_phase3_eligible_samples_match_contrast_plan():
    el = pd.read_csv(REPO_P2 / "metadata" / "phase3_eligible_samples.csv", dtype={"sample_id": str})
    plan = pd.read_csv(REPO_P2 / "metadata" / "reference_contrast_plan.csv", dtype=str)
    assert not el.duplicated(["contrast_id", "gsm"]).any()
    for _, p in plan.iterrows():
        e = el[el["contrast_id"] == p["contrast_id"]]
        for arm, col in (("case", "case_gsm"), ("control", "control_gsm"), ("timecourse", "timecourse_gsm")):
            want = set(str(p[col]).split(";")) - {"nan", ""}
            assert set(e.loc[e["arm"] == arm, "gsm"]) == want, (p["contrast_id"], arm)


@pytest.mark.skipif(not HAVE_DATA, reason="external data root not present")
def test_implied_length_checksum():
    man = json.loads((REPO_P2 / "metadata" / "phase3_inputs_manifest.json").read_text())["implied_length"]
    assert file_hashes(DATA_ROOT / "GSE316023" / "prepared" / "implied_length.csv.gz")[0] == man["md5"]
