"""Phase 5B shared paths. Read-only on Phases 1-5; writes only under research/roadmap_phase5b/."""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PB = Path(__file__).resolve().parents[1]
RS = PB.parent
REPO = RS.parent
P5 = RS / "roadmap_phase5"
TAB, FIG, RES = PB / "tables", PB / "figures", PB / "results"
WORK5 = Path("D:/mpliver_work/phase5")
SEED = 20261122
N_JOBS = int(os.environ.get("P5B_JOBS", "6"))
RSCRIPT = "/root/mpliver-tools/mamba/envs/mpliver-r/bin/Rscript"
sys.path.insert(0, str(RS / "roadmap_phase4c_diagnostics" / "scripts"))
sys.path.insert(0, str(P5 / "scripts"))


def sha256(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def json_dump(o, p):
    Path(p).write_text(json.dumps(o, indent=2, default=str))


def log2cpm(X: pd.DataFrame, prior=1.0) -> pd.DataFrame:
    return np.log2(X / X.sum() * 1e6 + prior)


def dl_tau2(y, s):
    w = 1 / s**2
    yb = (w * y).sum() / w.sum()
    q = (w * (y - yb) ** 2).sum()
    c = w.sum() - (w**2).sum() / w.sum()
    return max(0.0, (q - (len(y) - 1)) / c)


def wr2(y, p, w, ref=None):
    """Weighted R2 vs a constant reference (default: weighted mean of y)."""
    ref = np.average(y, weights=w) if ref is None else ref
    return 1 - np.sum(w * (y - p) ** 2) / np.sum(w * (y - ref) ** 2)


def mod_diff(L: np.ndarray, a: np.ndarray, b: np.ndarray, d0=4.0):
    """Fast mouse-level effect proxy: mean log2CPM difference (a minus b) with a limma-style moderated SE
    (gene variances shrunk toward their median with d0 prior df). L genes x mice; a,b boolean masks."""
    xa, xb = L[:, a], L[:, b]
    na, nb = a.sum(), b.sum()
    d = xa.mean(1) - xb.mean(1)
    ss = ((xa - xa.mean(1, keepdims=True)) ** 2).sum(1) + ((xb - xb.mean(1, keepdims=True)) ** 2).sum(1)
    df = na + nb - 2
    s2 = ss / max(df, 1)
    s0 = np.median(s2[s2 > 0]) if np.any(s2 > 0) else 1e-4
    s2m = (d0 * s0 + df * s2) / (d0 + df)
    return d, np.sqrt(s2m * (1 / na + 1 / nb))
