"""Phase 4B representations (feature_engineering_plan.md). Raw gene table = Phase 4 F8 + block-B memberships +
reference-only inputs for block C. RepBuilder turns raw columns into a representation; every data-dependent
quantity (block-A tau2, block-C z-scaling) is learnt in fit() on training rows only."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin

P4B = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P4B.parent / "roadmap_phase4" / "scripts"))
from p4_data import P3, RE3, _csv, feature_table  # noqa: E402

AXES = ["MASH", "FAST", "PPARA", "CIRC"]
F8 = [f"{a}_{k}" for a in AXES for k in ("eff", "se")]
F4 = [f"{a}_eff" for a in AXES]
B_SETS = ["LIPOGENESIS_CORE", "CHOL_SYNTHESIS_CORE", "SANDERSON_PPARA_TARGETS", "HALLMARK_FATTY_ACID_METABOLISM",
          "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
          "HALLMARK_XENOBIOTIC_METABOLISM", "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_PEROXISOME",
          "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_MTORC1_SIGNALING"]
B_COLS = [f"B_{s}" for s in B_SETS]
C_RAW = ["FASTclean_eff", "PPARA5_eff", "PPARA10_eff"]
RAW_COLS = F8 + B_COLS + C_RAW
A_COLS = [f"A_{a}_shr" for a in AXES]
C_COLS = ["C_FASTxPPARA", "C_MASHxFAST", "C_FAST_morning_dev", "C_PPARA_consistent"]
REPS = {
    "F4": F4, "F8": F8, "A4": A_COLS,
    "F8+A": F8 + A_COLS, "F8+B": F8 + B_COLS, "F8+C": F8 + C_COLS,
    "FULL": F8 + A_COLS + B_COLS + C_COLS,
    "FULL-A": F8 + B_COLS + C_COLS, "FULL-B": F8 + A_COLS + C_COLS, "FULL-C": F8 + A_COLS + B_COLS,
}


def raw_table(shuffle_b: bool = False) -> pd.DataFrame:
    X = feature_table("primary")  # identical rows to Phase 4
    gs = pd.read_csv(P3 / "results" / "gene_sets.csv")
    for s, c in zip(B_SETS, B_COLS):
        X[c] = X.index.isin(gs.gene[gs.set == s]).astype(float)
    if shuffle_b:  # control: break gene-membership link, keep set sizes and co-membership pattern
        rng = np.random.default_rng(13)
        X[B_COLS] = X[B_COLS].values[rng.permutation(len(X))]
    for col, f in [("FASTclean_eff", "FASTclean_GSE107787_ZT0ZT4_muscleLow"), ("PPARA5_eff", "PPARA_GSE282998_WT_5v0"),
                   ("PPARA10_eff", "PPARA_GSE282998_WT_10v0")]:
        X[col] = _csv(RE3 / f"{f}.csv.gz").set_index("gene")["log2FC"].reindex(X.index)
    assert not X[RAW_COLS].isna().any().any()
    return X[RAW_COLS]


class RepBuilder(BaseEstimator, TransformerMixin):
    def __init__(self, rep="F8"):
        self.rep = rep

    def fit(self, X, y=None):
        X = np.asarray(X, float)
        ix = {c: i for i, c in enumerate(RAW_COLS)}
        self.tau2_ = {}
        for a in AXES:
            e, s = X[:, ix[f"{a}_eff"]], X[:, ix[f"{a}_se"]]
            self.tau2_[a] = max(0.0, e.var() - (s**2).mean())
        self.mu_ = {a: X[:, ix[f"{a}_eff"]].mean() for a in ("MASH", "FAST", "PPARA")}
        self.sd_ = {a: X[:, ix[f"{a}_eff"]].std() for a in ("MASH", "FAST", "PPARA")}
        return self

    def _all(self, X):
        X = np.asarray(X, float)
        ix = {c: i for i, c in enumerate(RAW_COLS)}
        cols = {c: X[:, ix[c]] for c in F8 + B_COLS}
        for a in AXES:
            e, s = X[:, ix[f"{a}_eff"]], X[:, ix[f"{a}_se"]]
            t2 = self.tau2_[a]
            cols[f"A_{a}_shr"] = e * (t2 / (t2 + s**2)) if t2 > 0 else np.zeros_like(e)
        z = {a: (X[:, ix[f"{a}_eff"]] - self.mu_[a]) / self.sd_[a] for a in ("MASH", "FAST", "PPARA")}
        cols["C_FASTxPPARA"] = z["FAST"] * z["PPARA"]
        cols["C_MASHxFAST"] = z["MASH"] * z["FAST"]
        cols["C_FAST_morning_dev"] = X[:, ix["FASTclean_eff"]] - X[:, ix["FAST_eff"]]
        p20 = X[:, ix["PPARA_eff"]]
        cons = (np.sign(X[:, ix["PPARA5_eff"]]) == np.sign(p20)) & (np.sign(X[:, ix["PPARA10_eff"]]) == np.sign(p20))
        cols["C_PPARA_consistent"] = p20 * cons
        return cols

    def transform(self, X):
        cols = self._all(X)
        return np.column_stack([cols[c] for c in REPS[self.rep]])

    def get_feature_names_out(self, input_features=None):
        return np.array(REPS[self.rep])
