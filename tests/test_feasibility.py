import numpy as np
import pandas as pd
import pytest

from mpliver.feasibility import (
    CELLS,
    CONTRASTS,
    design_matrix,
    jackknife_interval,
    leave_one_mouse_out,
    minimum_detectable_lfc,
    moment_dispersion,
    r2_on_design,
)

GROUP_SIZES = [4, 5, 5, 5]


def labels():
    groups = np.repeat(CELLS, GROUP_SIZES)
    return pd.Series(groups, index=[f"s{i}" for i in range(len(groups))])


def test_design_matrix_is_full_rank_for_our_layout():
    g = labels()
    X = design_matrix(g.str.startswith("MD"), g.str.endswith("PE"), g.index)
    assert np.linalg.matrix_rank(X.to_numpy()) == 4
    assert X.loc[g.index[g == "MD+PE"], "diet_MD:PE"].eq(1).all()
    assert X.loc[g.index[g != "MD+PE"], "diet_MD:PE"].eq(0).all()


def test_r2_on_design_flags_a_covariate_that_mirrors_pe():
    g = labels()
    X = design_matrix(g.str.startswith("MD"), g.str.endswith("PE"), g.index)
    mirrors_pe = g.str.endswith("PE").astype(float) + np.linspace(0, 0.01, len(g))
    unrelated = np.random.default_rng(0).normal(size=len(g))
    assert r2_on_design(X, mirrors_pe)["r2_explained_by_design"] > 0.99
    assert r2_on_design(X, unrelated)["vif"] < 3


def test_moment_dispersion_recovers_known_value_on_average():
    rng = np.random.default_rng(1)
    g = labels()
    phi, n_genes = 0.05, 4000
    mu = rng.uniform(200, 2000, size=n_genes)
    sf = pd.Series(rng.uniform(0.8, 1.2, size=len(g)), index=g.index)
    counts = rng.negative_binomial(1 / phi, 1 / (1 + np.outer(mu, sf) * phi))
    norm_counts = pd.DataFrame(counts / sf.to_numpy(), columns=g.index)
    est = moment_dispersion(norm_counts, g, sf)
    assert est.median() == pytest.approx(phi, rel=0.15)


def test_minimum_detectable_lfc_matches_hand_calculation():
    # Two cells of 5, mu = 100, phi = 0.01: var per cell = (0.01 + 0.01) / 5 = 0.004
    se = np.sqrt(2 * 0.004)
    expected = (1.959964 + 0.841621) * se / np.log(2)
    got = minimum_detectable_lfc(100, 0.01, [5, 5, 5, 5], [0, 0, -1, 1], alpha=0.05)
    assert got == pytest.approx(expected, rel=1e-4)


def test_interaction_needs_larger_effect_than_simple_contrast():
    args = dict(mean_count=200, dispersion=0.01, n_cells=GROUP_SIZES, alpha=0.001)
    simple = minimum_detectable_lfc(contrast=CONTRASTS["PE within MD (MD+PE vs MD+Veh)"],
                                    **args)
    interaction = minimum_detectable_lfc(
        contrast=CONTRASTS["Diet x PE interaction (difference of PE effects)"], **args)
    main = minimum_detectable_lfc(contrast=CONTRASTS["PE main effect (average over diets)"],
                                  **args)
    assert main < simple < interaction


def test_leave_one_mouse_out_drops_each_mouse_once():
    g = labels()
    splits = list(leave_one_mouse_out(g))
    assert [m for m, _ in splits] == list(g.index)
    assert all(len(rest) == len(g) - 1 and m not in rest.index for m, rest in splits)


def test_jackknife_interval_matches_textbook_mean_case():
    x = np.array([1.0, 2.0, 4.0, 7.0, 11.0])
    loo_means = np.array([(x.sum() - v) / 4 for v in x])
    low, high = jackknife_interval(loo_means, x.mean())
    se_mean = x.std(ddof=1) / np.sqrt(len(x))   # jackknife SE of a mean equals this
    assert (high - low) / 2 == pytest.approx(1.96 * se_mean)
