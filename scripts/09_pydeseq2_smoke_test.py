"""Step 9 (optional): check that pydeseq2 can fit our factorial design.

Requires `pip install pydeseq2` (not in requirements.txt yet; the team has not
chosen between R DESeq2 and pydeseq2).

Two checks, neither of which tests real genes:
    1. Synthetic data with the real group sizes (4/5/5/5) and known diet, PE and
       interaction effects: does `~ diet + microplastic + diet:microplastic`
       fit, and are the planted effects recovered?
    2. Real counts: do pydeseq2's size factors match our Phase 1 median-of-ratios
       values? (Normalisation only; no dispersion fit, no Wald tests.)
"""

import sys

import numpy as np
import pandas as pd

from mpliver.config import PROCESSED_DIR, REPO_ROOT, SAMPLE_METADATA_FILE, TABLES_DIR

OUT = REPO_ROOT / "results" / "methodology" / "pydeseq2_smoke_test.md"


def synthetic_counts(meta: pd.DataFrame, n_genes: int = 3000, seed: int = 0):
    """10% of genes carry a planted effect (median-of-ratios assumes most genes do not change)."""
    rng = np.random.default_rng(seed)
    md = (meta["diet"] == "MD").to_numpy().astype(float)
    pe = (meta["microplastic"] == "PE").to_numpy().astype(float)
    base = rng.lognormal(mean=6, sigma=1.2, size=n_genes)
    effects = np.zeros((n_genes, 3))            # log2 diet, PE, interaction
    effects[:100, 0] = 1.5                      # diet genes
    effects[100:200, 1] = 1.0                   # PE genes
    effects[200:300, 2] = 1.5                   # interaction genes
    log2_mu = (np.log2(base)[:, None] + effects[:, [0]] * md + effects[:, [1]] * pe
               + effects[:, [2]] * md * pe)
    depth = rng.uniform(0.8, 1.25, size=len(meta))
    mu = 2 ** log2_mu * depth
    dispersion = 0.015
    counts = rng.negative_binomial(1 / dispersion, 1 / (1 + mu * dispersion))
    genes = [f"g{i}" for i in range(n_genes)]
    return pd.DataFrame(counts.T, index=meta["counts_column"], columns=genes), effects


def main() -> int:
    try:
        import pydeseq2
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.default_inference import DefaultInference
    except ImportError:
        print("pydeseq2 is not installed: pip install pydeseq2")
        return 1

    meta = pd.read_csv(SAMPLE_METADATA_FILE)
    # Reference levels must be set explicitly: the formula engine otherwise orders
    # levels alphabetically, which would make PE (not vehicle) the reference.
    design_meta = meta.set_index("counts_column")[["diet", "microplastic"]].copy()
    design_meta["diet"] = pd.Categorical(design_meta["diet"], categories=["CD", "MD"])
    design_meta["microplastic"] = pd.Categorical(design_meta["microplastic"],
                                                 categories=["Veh", "PE"])
    lines = [f"# pydeseq2 smoke test (pydeseq2 {pydeseq2.__version__})", ""]

    # 1. Synthetic factorial data -------------------------------------------------
    counts, effects = synthetic_counts(meta)
    dds = DeseqDataSet(counts=counts, metadata=design_meta,
                       design="~ diet + microplastic + diet:microplastic",
                       refit_cooks=True, inference=DefaultInference(n_cpus=1), quiet=True)
    dds.deseq2()
    lfc = dds.varm["LFC"]                       # natural-log coefficients
    lfc.columns = [str(c) for c in lfc.columns]
    lines += ["## 1. Synthetic data, real group sizes (4/5/5/5)", "",
              f"Design columns: `{list(lfc.columns)}`", "",
              "| Planted effect (log2) | Estimated coefficient | Mean estimate on planted genes "
              "(log2) | Mean on null genes |", "|---|---|---|---|"]
    targets = {0: "diet", 1: "microplastic", 2: ":"}
    for idx, key in targets.items():
        col = next(c for c in lfc.columns if key in c and (idx == 2 or ":" not in c))
        est = lfc[col].to_numpy() / np.log(2)
        planted = effects[:, idx] != 0
        lines.append(f"| {effects[planted, idx][0]} | `{col}` | {est[planted].mean():.2f} | "
                     f"{est[effects.sum(axis=1) == 0].mean():.2f} |")

    # 2. Size factors on the real data -----------------------------------------------
    order = meta["counts_column"]
    real = pd.read_csv(PROCESSED_DIR / "counts_filtered.csv.gz", index_col=0)[order].T
    dds_real = DeseqDataSet(counts=real.astype(int), metadata=design_meta.loc[order],
                            design="~ diet + microplastic + diet:microplastic", quiet=True)
    dds_real.fit_size_factors()
    ours = pd.read_csv(TABLES_DIR / "size_factors.csv", index_col=0)["size_factor"]
    sf = (dds_real.obs["size_factors"] if "size_factors" in dds_real.obs
          else dds_real.obsm["size_factors"])
    theirs = pd.Series(np.asarray(sf, dtype=float), index=order)
    diff = (theirs / ours.loc[order] - 1).abs().max()
    lines += ["", "## 2. Real data: size factors only (no model fit, no tests)", "",
              f"Maximum relative difference between pydeseq2 and Phase 1 size factors: "
              f"**{diff:.2e}** across 19 samples (13,777 filtered genes).", ""]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
