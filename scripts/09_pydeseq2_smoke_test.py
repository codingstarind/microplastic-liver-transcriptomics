"""Step 9 (optional): check that pydeseq2 can fit our factorial design.

Requires `pip install -r requirements-methodology.txt` (pydeseq2 is not yet adopted
for analysis; the team has not chosen between R DESeq2 and pydeseq2).

What this checks:
    1. Synthetic data with the real group sizes (4/5/5/5) and planted diet, PE
       and interaction effects: does `~ diet + microplastic + diet:microplastic`
       fit with the intended reference levels, and are the planted coefficients
       recovered on average?
    2. Real counts (skipped with --synthetic-only): do pydeseq2's size factors
       match our Phase 1 median-of-ratios values? Normalisation only.

What this does NOT check: whether p-values are calibrated, whether the FDR is
controlled at the nominal level, or how the model behaves with outliers,
contamination or the real dispersion structure. Those need a dedicated
simulation study (many replicate datasets, null and non-null genes) and are
listed as open validation work in reports/methodology_proposal.md.
"""

import argparse
import sys

import numpy as np
import pandas as pd

from mpliver.config import PROCESSED_DIR, REPO_ROOT, SAMPLE_METADATA_FILE, TABLES_DIR

OUT = REPO_ROOT / "results" / "methodology" / "pydeseq2_smoke_test.md"
DESIGN = "~ diet + microplastic + diet:microplastic"
EXPECTED_COLUMNS = ["Intercept", "diet[T.MD]", "microplastic[T.PE]",
                    "diet[T.MD]:microplastic[T.PE]"]
TOLERANCE = 0.15  # log2 units, for the mean estimate over 100 planted genes


def synthetic_counts(meta: pd.DataFrame, n_genes: int = 3000, seed: int = 0):
    """10% of genes carry a planted effect (median-of-ratios assumes most genes do not change)."""
    rng = np.random.default_rng(seed)
    md = (meta["diet"] == "MD").to_numpy().astype(float)
    pe = (meta["microplastic"] == "PE").to_numpy().astype(float)
    base = rng.lognormal(mean=6, sigma=1.2, size=n_genes)
    effects = np.zeros((n_genes, 3))            # log2 diet, PE, interaction
    effects[:100, 0] = 1.5
    effects[100:200, 1] = 1.0
    effects[200:300, 2] = 1.5
    log2_mu = (np.log2(base)[:, None] + effects[:, [0]] * md + effects[:, [1]] * pe
               + effects[:, [2]] * md * pe)
    mu = 2 ** log2_mu * rng.uniform(0.8, 1.25, size=len(meta))
    dispersion = 0.015
    counts = rng.negative_binomial(1 / dispersion, 1 / (1 + mu * dispersion))
    genes = [f"g{i}" for i in range(n_genes)]
    return pd.DataFrame(counts.T, index=meta["counts_column"], columns=genes), effects


def design_metadata(meta: pd.DataFrame) -> pd.DataFrame:
    # Reference levels must be explicit: the formula engine otherwise orders
    # levels alphabetically, which would make PE (not vehicle) the reference.
    df = meta.set_index("counts_column")[["diet", "microplastic"]].copy()
    df["diet"] = pd.Categorical(df["diet"], categories=["CD", "MD"])
    df["microplastic"] = pd.Categorical(df["microplastic"], categories=["Veh", "PE"])
    return df


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--synthetic-only", action="store_true",
                        help="skip the real-data size-factor check (for CI)")
    parser.add_argument("--strict", action="store_true",
                        help="exit non-zero if reference levels or recovery checks fail")
    args = parser.parse_args()
    try:
        import pydeseq2
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.default_inference import DefaultInference
    except ImportError:
        print("pydeseq2 is not installed: pip install -r requirements-methodology.txt")
        return 1

    meta = pd.read_csv(SAMPLE_METADATA_FILE)
    design_meta = design_metadata(meta)
    failures = []
    lines = [f"# pydeseq2 smoke test (pydeseq2 {pydeseq2.__version__})", "",
             "**Scope:** checks that the factorial design fits with the intended reference "
             "levels and that planted coefficients are recovered on average. It does **not** "
             "validate p-value calibration, FDR control, or behaviour under outliers and "
             "contamination.", ""]

    # 1. Synthetic factorial data -------------------------------------------------
    counts, effects = synthetic_counts(meta)
    dds = DeseqDataSet(counts=counts, metadata=design_meta, design=DESIGN,
                       inference=DefaultInference(n_cpus=1), quiet=True)
    dds.deseq2()
    lfc = dds.varm["LFC"]
    lfc.columns = [str(c) for c in lfc.columns]
    if list(lfc.columns) != EXPECTED_COLUMNS:
        failures.append(f"unexpected design columns {list(lfc.columns)}")
    lines += ["## 1. Synthetic data, real group sizes (4/5/5/5): coefficient recovery", "",
              f"Design columns: `{list(lfc.columns)}`", "",
              "| Planted effect (log2) | Coefficient | Mean estimate, planted genes (log2) "
              "| Mean estimate, null genes (log2) |", "|---|---|---|---|"]
    null = effects.sum(axis=1) == 0
    for idx, col in enumerate(EXPECTED_COLUMNS[1:]):
        if col not in lfc.columns:
            continue
        est = lfc[col].to_numpy() / np.log(2)   # pydeseq2 stores natural-log coefficients
        planted = effects[:, idx] != 0
        target = effects[planted, idx][0]
        if abs(est[planted].mean() - target) > TOLERANCE or abs(est[null].mean()) > TOLERANCE:
            failures.append(f"{col}: planted {target}, recovered {est[planted].mean():.2f}")
        lines.append(f"| {target} | `{col}` | {est[planted].mean():.2f} | "
                     f"{est[null].mean():.2f} |")

    # 2. Size factors on the real data -----------------------------------------------
    if not args.synthetic_only:
        order = meta["counts_column"]
        real = pd.read_csv(PROCESSED_DIR / "counts_filtered.csv.gz", index_col=0)[order].T
        dds_real = DeseqDataSet(counts=real.astype(int), metadata=design_meta.loc[order],
                                design=DESIGN, quiet=True)
        dds_real.fit_size_factors()
        sf = (dds_real.obs["size_factors"] if "size_factors" in dds_real.obs
              else dds_real.obsm["size_factors"])
        theirs = pd.Series(np.asarray(sf, dtype=float), index=order)
        ours = pd.read_csv(TABLES_DIR / "size_factors.csv", index_col=0)["size_factor"]
        diff = (theirs / ours.loc[order] - 1).abs().max()
        lines += ["", "## 2. Real data: size factors only (no dispersion fit, no tests)", "",
                  f"Maximum relative difference between pydeseq2 and Phase 1 size factors: "
                  f"**{diff:.2e}** across 19 samples (13,777 filtered genes).", ""]
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))
    for failure in failures:
        print(f"[failed] {failure}")
    return 1 if (args.strict and failures) else 0


if __name__ == "__main__":
    sys.exit(main())
