"""Phase 3 output tests: contrast directions, reference levels, alignment with fixed eligible samples,
coefficient extraction, gene universe, output completeness, plan lock.
Run: python -m pytest research/roadmap_phase3/tests -q
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

P3 = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P3 / "scripts"))
from p3lib import P2, RE, universe  # noqa: E402

REFS = sorted(p.name.replace(".csv.gz", "") for p in RE.glob("*.csv.gz"))
U = universe("main")


def eff(label):
    return pd.read_csv(RE / f"{label}.csv.gz").set_index("gene")


def test_plan_lock_unchanged():
    lock = json.loads((P3 / "results" / "plan_lock.json").read_text())
    assert hashlib.sha256((P3 / "analysis_plan.md").read_bytes()).hexdigest() == lock["sha256"]


@pytest.mark.parametrize("label", REFS)
def test_reference_table_schema_and_universe(label):
    t = pd.read_csv(RE / f"{label}.csv.gz")
    for c in ("gene", "log2FC", "ci95_low", "ci95_high", "stat", "pvalue", "padj", "detected", "dataset", "contrast"):
        assert c in t.columns, (label, c)
    assert t["gene"].is_unique
    assert set(t["gene"]) <= U
    ok = t["log2FC"].notna()
    assert ok.mean() > 0.98
    assert (t.loc[ok, "ci95_low"] <= t.loc[ok, "log2FC"]).all() and (t.loc[ok, "log2FC"] <= t.loc[ok, "ci95_high"]).all()
    p = t["pvalue"].dropna()
    assert ((p >= 0) & (p <= 1)).all()
    # BH within the reported family
    from scipy.stats import false_discovery_control
    pp = t.dropna(subset=["pvalue"])
    assert np.allclose(np.sort(false_discovery_control(pp["pvalue"].values)), np.sort(pp["padj"].values), atol=1e-8)


def test_mash_direction_gan_over_chow():
    for lab in ("MASH1_GSE316023_GANvsChow_TA", "MASH2_GSE246328_GANvsChow_24w"):
        t = eff(lab)
        assert (t["coefficient"] == "diet_GAN_vs_chow").all()
        for g in ("Gpnmb", "Lgals3", "Col3a1", "Cd68", "Trem2"):  # disease-induced markers (in universe) UP in GAN
            assert t.loc[g, "log2FC"] > 0, (lab, g)


def test_fasting_direction():
    t = eff("FAST_GSE107787_FASTvsFED_ZTadj")
    assert t.loc["Hmgcs2", "log2FC"] > 0 and t.loc["Fasn", "log2FC"] < 0  # ketogenesis up, lipogenesis down


def test_ppara_direction_and_ko_attenuation():
    w, k = eff("PPARA_GSE282998_WT_20v0"), eff("PPARA_GSE282998_KO_20v0")
    for g in ("Cyp4a14", "Acot1", "Ehhadh"):
        assert w.loc[g, "log2FC"] > 1 and w.loc[g, "log2FC"] > k.loc[g, "log2FC"]


def test_circadian_contrast_signs():
    t2 = eff("CIRC_GSE117134_WTfemale_T2_ZT12_minus_ZT0")
    assert t2.loc["Dbp", "log2FC"] > 0 and t2.loc["Bmal1", "log2FC"] < 0  # Dbp peaks ~ZT10-12, Bmal1 ~ZT22-0


def test_mash_samples_match_fixed_eligible_list():
    el = pd.read_csv(P2 / "metadata" / "phase3_eligible_samples.csv", dtype={"sample_id": str})
    info = json.loads((P3 / "results" / "mash" / "run_info.json").read_text())
    for key, cid in (("GSE316023_primary", "MASH_primary"), ("GSE246328_primary", "MASH_secondary")):
        e = el[el.contrast_id == cid]
        assert set(info[key]["samples_case"]) == set(e.loc[e.arm == "case", "sample_id"])
        assert set(info[key]["samples_control"]) == set(e.loc[e.arm == "control", "sample_id"])


def test_fasting_clean_subset_is_the_fixed_one():
    el = pd.read_csv(P2 / "metadata" / "phase3_eligible_samples.csv", dtype={"sample_id": str})
    info = json.loads((P3 / "results" / "other_refs_run_info.json").read_text())
    fixed = el[(el.contrast_id == "FASTING_primary") & el.in_clean_sensitivity]
    assert set(info["fasting_clean"]["samples"]) == set(fixed.sample_id)
    assert info["fasting_clean"]["n_fast"] == 4 and info["fasting_clean"]["n_fed"] == 6


def test_positive_control_decision_recomputes():
    pc = json.loads((P3 / "results" / "positive_control" / "positive_control.json").read_text())
    m = pc["M1_vs_target"]
    c1 = m["spearman"] >= .30 and pc["M5_relabel_null"]["observed_rho_meandiff"] > pc["M5_relabel_null"]["null_q99"]
    assert c1 == pc["criteria_results"]["c1_rho_and_null"]
    c = pc["criteria_results"]
    expect = "PASS" if all(c.values()) else ("CONDITIONAL" if c["c1_rho_and_null"] and c["c2_sign_agreement"] else "FAIL")
    assert pc["decision"] == expect


def test_deliverables_exist():
    for f in ["analysis_plan.md", "README.md", "reports/count_offset_validation.md", "reports/mash_reference_report.md",
              "reports/fasting_reference_report.md", "reports/circadian_reference_report.md", "reports/ppara_reference_report.md",
              "reports/reference_comparison_report.md", "reports/baseline_comparison_report.md",
              "reports/phase3_decision_report.md", "figures/mash_positive_control.png",
              "figures/reference_similarity_heatmap.png", "figures/baseline_comparisons.png",
              "figures/fasting_reference.png", "figures/ppara_reference.png", "figures/circadian_reference.png"]:
        assert (P3 / f).exists(), f


def test_posthoc_audit_consistent_and_locked_decision_retained():
    ex = json.loads((P3 / "results" / "posthoc_audit" / "exact_permutation.json").read_text())
    assert ex["n_allocations"] == 12870 and ex["observed_matches_original_montecarlo"]
    assert abs(ex["exact_one_sided_p"] - ex["n_allocations_ge_observed"] / 12870) < 1e-12
    pc = json.loads((P3 / "results" / "positive_control" / "positive_control.json").read_text())
    assert pc["decision"] == "FAIL"  # locked outcome must not be rewritten
    lock = json.loads((P3 / "results" / "posthoc_audit" / "muscle_prereg_lock.json").read_text())
    assert hashlib.sha256((P3 / "posthoc_muscle_rule_prereg.md").read_bytes()).hexdigest() == lock["sha256"]
