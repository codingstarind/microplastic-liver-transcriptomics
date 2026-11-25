"""Shared paths/helpers for the Phase 4C diagnostic audit (read-only on Phases 2-4B)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

P4C = Path(__file__).resolve().parents[1]
RS = P4C.parent
REPO = RS.parent
T, F = P4C / "tables", P4C / "figures"
sys.path.insert(0, str(RS / "roadmap_phase4b" / "scripts"))
sys.path.insert(0, str(RS / "roadmap_phase4" / "scripts"))
from p4_data import _csv, dl_tau2, gene_groups, loo_names, outer_splits, target  # noqa: E402,F401
from p4b_features import B_COLS, C_COLS, F8, RAW_COLS, REPS, RepBuilder, raw_table  # noqa: E402,F401

SEED = 20261122
MATS = {"F8": F8, "A4": REPS["A4"], "B12": B_COLS, "C4": C_COLS, "F8+C": F8 + C_COLS,
        "F8+B+C": F8 + B_COLS + C_COLS, "FULL": REPS["FULL"]}


def nb_interval(d, n_test, n_train):
    J = len(d)
    m, s2 = d.mean(), d.var(ddof=1)
    se = np.sqrt((1 / J + np.mean(n_test) / np.mean(n_train)) * s2)
    t = stats.t.ppf(0.975, J - 1)
    return m, m - t * se, m + t * se


def block_boot(df, stat, groups="module", n=1000, seed=SEED):
    """Module-block bootstrap: resample whole modules with replacement; returns (est, lo, hi)."""
    rng = np.random.default_rng(seed)
    mods = df[groups].unique()
    by = {g: d for g, d in df.groupby(groups)}
    est = stat(df)
    bs = []
    for _ in range(n):
        pick = rng.choice(mods, len(mods), replace=True)
        bs.append(stat(pd.concat([by[g] for g in pick])))
    return est, float(np.nanquantile(bs, .025)), float(np.nanquantile(bs, .975))
