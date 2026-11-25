"""Shared helpers for Phase 3 comparisons (Python)."""
from __future__ import annotations

import gzip
import io
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, spearmanr

P3 = Path(__file__).resolve().parents[1]
REPO = P3.parents[1]
P2 = REPO / "research" / "roadmap_phase2"
RE = P3 / "reference_effects"
SEED = 20261122
DEG_P, DEG_LFC = 0.05, 0.38


def universe(kind="main"):
    u = pd.read_csv(P2 / "metadata" / "phase3_gene_universe.csv", index_col=0)
    return set(u.index[u["universe_" + kind]])


def target(contrast: str, root: Path | None = None) -> pd.DataFrame:
    root = root or REPO / "results" / "deseq2" / "provisional_deposited"
    t = pd.read_csv(root / f"results_{contrast}.csv.gz").set_index("gene")
    return t[["log2FC", "lfcSE", "padj", "baseMean"]]


def ref(label: str, effect="log2FC", p="padj") -> pd.DataFrame:
    t = pd.read_csv(RE / f"{label}.csv.gz").set_index("gene")
    return t.rename(columns={effect: "log2FC", p: "padj"})[["log2FC", "padj"]]


def compare(a: pd.DataFrame, b: pd.DataFrame, genes: set, name_a="A", name_b="B") -> dict:
    """Effect-level and DEG-level agreement of two signed gene effect tables on a fixed gene set."""
    g = sorted(set(a.index) & set(b.index) & genes)
    x, y = a.loc[g], b.loc[g]
    ok = x["log2FC"].notna() & y["log2FC"].notna()
    x, y = x[ok], y[ok]
    rho = spearmanr(x["log2FC"], y["log2FC"]).statistic
    both_sig = (x["padj"] < DEG_P) & (y["padj"] < DEG_P)
    sign_both = float((np.sign(x.loc[both_sig, "log2FC"]) == np.sign(y.loc[both_sig, "log2FC"])).mean()) if both_sig.any() else np.nan
    bsig = y["padj"] < DEG_P
    sign_ref_sig = float((np.sign(x.loc[bsig, "log2FC"]) == np.sign(y.loc[bsig, "log2FC"])).mean()) if bsig.any() else np.nan
    da = (x["padj"] < DEG_P) & (x["log2FC"].abs() >= DEG_LFC)
    db = (y["padj"] < DEG_P) & (y["log2FC"].abs() >= DEG_LFC)
    n, ka, kb, k = len(x), int(da.sum()), int(db.sum()), int((da & db).sum())
    conc = int((da & db & (np.sign(x["log2FC"]) == np.sign(y["log2FC"]))).sum())
    p_hyp = float(hypergeom.sf(k - 1, n, ka, kb)) if ka and kb else np.nan
    up_a, dn_a = da & (x["log2FC"] > 0), da & (x["log2FC"] < 0)
    up_b, dn_b = db & (y["log2FC"] > 0), db & (y["log2FC"] < 0)
    k_conc_dir = int((up_a & up_b).sum() + (dn_a & dn_b).sum())
    p_dir = float(hypergeom.sf(k_conc_dir - 1, n, int(up_a.sum()), int(up_b.sum())) *
                  0 + hypergeom.sf(int((up_a & up_b).sum()) - 1, n, int(up_a.sum()), int(up_b.sum()))) if up_a.any() and up_b.any() else np.nan
    return dict(a=name_a, b=name_b, n_genes=n, spearman=float(rho),
                n_sig_both=int(both_sig.sum()), sign_agree_sig_both=sign_both,
                n_sig_b=int(bsig.sum()), sign_agree_on_b_sig=sign_ref_sig,
                n_deg_a=ka, n_deg_b=kb, n_deg_overlap=k, jaccard=k / (ka + kb - k) if (ka + kb - k) else np.nan,
                expected_overlap=ka * kb / n if n else np.nan, hypergeom_p=p_hyp,
                n_overlap_concordant=conc, n_overlap_discordant=k - conc,
                up_up_hypergeom_p=p_dir)


def read_gz_csv(path, **kw):
    with gzip.open(path, "rt") as fh:
        return pd.read_csv(io.StringIO(fh.read()), **kw)


def relabel_null(log2norm: pd.DataFrame, n_case: int, target_lfc: pd.Series, n_perm=1000, seed=SEED,
                 strata: list | None = None) -> np.ndarray:
    """Spearman(target, mean-difference) for random relabellings of reference samples.
    strata: optional list of arrays of column positions; relabelling is done within each stratum,
    keeping the per-stratum case count of the observed design."""
    rng = np.random.default_rng(seed)
    g = log2norm.index.intersection(target_lfc.dropna().index)
    X = log2norm.loc[g].values
    t = target_lfc.loc[g].rank().values
    out = np.empty(n_perm)
    m = X.shape[1]
    for i in range(n_perm):
        if strata is None:
            idx = rng.permutation(m)[:n_case]
        else:
            idx = np.concatenate([rng.permutation(s)[:k] for s, k in strata])
        mask = np.zeros(m, bool)
        mask[idx] = True
        d = X[:, mask].mean(1) - X[:, ~mask].mean(1)
        out[i] = np.corrcoef(pd.Series(d).rank().values, t)[0, 1]
    return out
