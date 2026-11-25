"""Phase 5C shared paths and frozen-artefact loaders. Reuses Phase 5 code (p5_common) read-only.
Scoring functions score_4d / score_4e are verbatim copies of roadmap_phase5/scripts/s05_sample_level.py
(that file executes on import, so it cannot be imported)."""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

P5C = Path(__file__).resolve().parents[1]
RS = P5C.parent
REPO = RS.parent
P5 = RS / "roadmap_phase5"
P5B = RS / "roadmap_phase5b"
sys.path.insert(0, str(P5 / "scripts"))
from p5_common import SETS12, holm, singscore_signed, singscore_up, welch_ci  # noqa: E402,F401

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

FRZ = P5 / "frozen"
RES, TAB, FIG = P5C / "results", P5C / "tables", P5C / "figures"
WORK = Path("D:/mpliver_work/phase5c")
SEED = 20261122
STUDY = "Z23"  # PRJNA916285, Zhang et al. ACS Nano 2023 (PMID 37486121)

length = pd.read_csv(RS / "roadmap_phase4d/results/cache/gene_lengths.csv", index_col=0)["length"]
BG4D = pd.read_csv(FRZ / "phase4d_background_genes.csv").gene.tolist()
SIGS = {k: pd.read_csv(RS / f"roadmap_phase4d/signatures/{k}_main_n50.csv") for k in ("FAST", "MASH", "PPARA")}
PAR = {r: json.loads((FRZ / f"phase4e_{r}_l2_logistic_params.json").read_text()) for r in ("R0", "R2")}


def sha256(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def verify_frozen() -> dict:
    """Every Phase 5 frozen artefact must match frozen_manifest.json byte-for-byte."""
    man = json.loads((FRZ / "frozen_manifest.json").read_text())["files"]
    res = {f: sha256(FRZ / f) == h for f, h in man.items()}
    assert all(res.values()), res
    return res


def score_4d(X: pd.DataFrame, kind: str, background=None):
    bg = [g for g in (background or BG4D) if g in X.index]
    M = X.loc[bg].astype(float)
    if kind == "counts":
        M = M.div(length[bg], axis=0)
        M = M / M.sum() * 1e6
    R = np.apply_along_axis(stats.rankdata, 1, np.log2(M.T.values + 1))
    gi = {g: i for i, g in enumerate(bg)}
    out, cov = {}, {}
    for k, s in SIGS.items():
        up = np.array([gi[g] for g in s[s.direction == "up"].gene if g in gi])
        dn = np.array([gi[g] for g in s[s.direction == "down"].gene if g in gi])
        out[k] = singscore_signed(R, up, dn)
        cov[k] = dict(n_up=len(up), n_down=len(dn), n_missing=100 - len(up) - len(dn))
    return pd.DataFrame(out, index=X.columns), cov, len(bg)


def score_4e(X: pd.DataFrame, rep: str, background=None):
    par = PAR[rep]
    bg = [g for g in (background or par["background_genes"]) if g in X.index]
    R = np.apply_along_axis(stats.rankdata, 1, X.loc[bg].T.values.astype(float))
    gi = {g: i for i, g in enumerate(bg)}
    F = [singscore_signed(R, np.array([gi[g] for g in par["up"] if g in gi]), np.array([gi[g] for g in par["dn"] if g in gi]))]
    cov = {"R0_signature": sum(g in gi for g in par["up"] + par["dn"])}
    for s, gl in par["sets"].items():
        idx = np.array([gi[g] for g in gl if g in gi])
        F.append(singscore_up(R, idx))
        cov[f"R1_{s}"] = len(idx)
    F = np.column_stack(F)
    Z = (F - np.array(par["scaler_mean"])) / np.array(par["scaler_scale"])
    lo = Z @ np.array(par["coef"]) + par["intercept"]
    return pd.Series(lo, index=X.columns), pd.DataFrame(F, index=X.columns, columns=par["feature_names"]), cov, len(bg)


def load_counts():
    X = pd.read_csv(WORK / "Z23_gene_counts.csv.gz", index_col=0)
    S = pd.read_csv(RES / "z23_sample_manifest.csv").set_index("mouse_id")
    S = S.loc[[c for c in S.index if c in X.columns]]
    return X[list(S.index)], S


def log2cpm(X):
    return np.log2(X / X.sum() * 1e6 + 1)
