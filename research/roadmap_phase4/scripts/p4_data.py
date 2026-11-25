"""Phase 4 data layer: gene feature tables, targets, gene groups and frozen outer splits.
Implements evaluation_protocol.md sections 1-3 and 7-8. Read-only on all earlier phases."""
from __future__ import annotations

import gzip
import io
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

P4 = Path(__file__).resolve().parents[1]
RS = P4.parents[0]  # research/
REPO = RS.parent
P2, P3 = RS / "roadmap_phase2", RS / "roadmap_phase3"
RE3 = P3 / "reference_effects"
RES = P4 / "results"
DROOT = Path("D:/mpliver_work/phase2_refs")
SEED = 20261122
OUTER_SEEDS = (101, 102, 103)
N_OUTER = 5
AXES = ("MASH", "FAST", "PPARA", "CIRC")
FEATURES_F8 = [f"{a}_{k}" for a in AXES for k in ("eff", "se")]
FEATURES_F4 = [f"{a}_eff" for a in AXES]
FORBIDDEN = {"gene", "padj", "pvalue", "baseMean", "stat", "target", "y"}


def _csv(path, **kw):
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            return pd.read_csv(io.StringIO(fh.read()), **kw)
    return pd.read_csv(path, **kw)


def _ref(path):
    t = _csv(path).set_index("gene")
    se = "lfcSE" if "lfcSE" in t.columns else "SE"
    return t[["log2FC", se]].rename(columns={"log2FC": "eff", se: "se"})


REF_FILES = {
    "primary": {"MASH": RE3 / "MASH1_GSE316023_GANvsChow_TA.csv.gz",
                "FAST": RE3 / "FAST_GSE107787_FASTvsFED_ZTadj.csv.gz",
                "PPARA": RE3 / "PPARA_GSE282998_WT_20v0.csv.gz",
                "CIRC": RE3 / "CIRC_GSE117134_WTfemale_T1_ZT6_minus_ZT2.csv.gz"},
}
_alt = RES / "alt_universe_refs"
VARIANTS = {
    "primary": dict(),
    "V1_MASH2": dict(MASH=RE3 / "MASH2_GSE246328_GANvsChow_24w.csv.gz"),
    "V2_FASTclean": dict(FAST=RE3 / "FASTclean_GSE107787_ZT0ZT4_muscleLow.csv.gz"),
    "V3_PPARA10v0": dict(PPARA=RE3 / "PPARA_GSE282998_WT_10v0.csv.gz"),
    "V4_CIRC_T2": dict(CIRC=RE3 / "CIRC_GSE117134_WTfemale_T2_ZT12_minus_ZT0.csv.gz"),
    "V5_strict_universe": dict(universe="sensitivity"),
    "V6_alt_muscle": dict(MASH=P3 / "results" / "posthoc_audit" / "MASH1_altuniverse.csv.gz",
                          FAST=_alt / "FAST.csv.gz", PPARA=_alt / "PPARA_20v0.csv.gz", CIRC=_alt / "CIRC_T1.csv.gz",
                          universe="alt"),
    "V7_chromosome_groups": dict(grouping="chromosome"),
    "V8_effects_only": dict(features="F4"),
}


@lru_cache(None)
def universe(kind: str) -> frozenset:
    if kind == "alt":
        return frozenset(_csv(P3 / "results" / "posthoc_audit" / "alt_universe.csv")["vm35_symbol"])
    u = pd.read_csv(P2 / "metadata" / "phase3_gene_universe.csv", index_col=0)
    return frozenset(u.index[u["universe_" + kind]])


def target(name: str) -> pd.DataFrame:
    """Return DataFrame[y, s] indexed by gene. Names: PE_in_CD, PE_in_MD, Diet_in_Veh,
    LOO:<dropped>:<contrast>, DEPTH_CD, DEPTH_MD, PERM:<contrast>."""
    if name.startswith("LOO:"):
        _, dropped, c = name.split(":")
        p = REPO / "results" / "deseq2" / "loo" / f"provisional_deposited__without_{dropped}" / f"results_{c}.csv.gz"
    elif name.startswith("DEPTH_"):
        p = RES / "depth_thinning" / f"{name}.csv.gz"
    elif name.startswith("PERM:"):
        t = target(name.split(":")[1]).copy()
        rng = np.random.default_rng(7)
        t[["y", "s"]] = t[["y", "s"]].values[rng.permutation(len(t))]
        return t
    else:
        p = REPO / "results" / "deseq2" / "provisional_deposited" / f"results_{name}.csv.gz"
    t = _csv(p).set_index("gene")
    out = pd.DataFrame({"y": t["log2FC"], "s": t["lfcSE"]})
    if "padj" in t:
        out["deg"] = (t["padj"] < 0.05) & (t["log2FC"].abs() >= 0.38)
    return out


def loo_names() -> list[str]:
    d = REPO / "results" / "deseq2" / "loo"
    return sorted(p.name.split("without_")[1] for p in d.glob("provisional_deposited__without_*"))


def feature_table(variant: str = "primary") -> pd.DataFrame:
    v = VARIANTS[variant]
    files = {**REF_FILES["primary"], **{k: v[k] for k in AXES if k in v}}
    U = universe(v.get("universe", "main"))
    parts = []
    for a in AXES:
        r = _ref(files[a])
        parts.append(r.rename(columns={"eff": f"{a}_eff", "se": f"{a}_se"}))
    X = pd.concat(parts, axis=1, join="inner")
    X = X.loc[X.index.intersection(list(U))].dropna()
    return X.sort_index()


def feature_names(variant: str) -> list[str]:
    return FEATURES_F4 if VARIANTS[variant].get("features") == "F4" else FEATURES_F8


@lru_cache(None)
def gene_groups() -> pd.DataFrame:
    f = RES / "gene_groups.csv"
    if f.exists():
        return pd.read_csv(f, index_col=0)
    from sklearn.cluster import AgglomerativeClustering
    genes = sorted(universe("main") | universe("alt"))
    x = _csv(DROOT / "GSE117134" / "prepared" / "expression_ready.csv.gz", index_col=0)
    x = np.log2(x.reindex(genes) + 1)
    sd = x.std(axis=1)
    ok = sd > 0
    z = x[ok].sub(x[ok].mean(axis=1), axis=0).div(sd[ok], axis=0)
    lab = AgglomerativeClustering(n_clusters=20, linkage="ward").fit_predict(z.values)
    mod = pd.Series(-1, index=genes, name="module")
    mod[z.index] = lab
    ann = _csv(REPO / "metadata" / "gene_annotation_gencode_vM35.csv.gz").set_index("gene")
    chrom = ann["chromosome"].reindex(genes).fillna("unknown").str.split(";").str[0]
    g = pd.DataFrame({"module": mod, "chromosome": chrom})
    g.to_csv(f)
    return g


def outer_splits(genes: pd.Index, grouping: str = "module"):
    """Yield (repeat, fold, train_mask, test_mask, groups) with whole groups assigned to folds."""
    grp = gene_groups().loc[genes, grouping].astype(str).values
    uniq = np.array(sorted(set(grp)))
    for r, seed in enumerate(OUTER_SEEDS):
        rng = np.random.default_rng(seed)
        perm = rng.permutation(uniq)
        fold_of = {g: i % N_OUTER for i, g in enumerate(perm)}
        f = np.array([fold_of[g] for g in grp])
        for k in range(N_OUTER):
            yield r, k, f != k, f == k, grp


def dl_tau2(y, s):
    """DerSimonian-Laird between-gene variance (training rows only)."""
    w = 1 / s**2
    ybar = (w * y).sum() / w.sum()
    q = (w * (y - ybar) ** 2).sum()
    k = len(y)
    c = w.sum() - (w**2).sum() / w.sum()
    return max(0.0, (q - (k - 1)) / c)
