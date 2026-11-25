"""Step 1: build and freeze signed signatures from INDEPENDENT reference effects only (no GSE294180 PE data).

Deterministic rules (also in analysis_plan.md):
  universe U   = Phase 2 main universe (11,888; locked 309 muscle exclusion) [sensitivity: alt 88-gene exclusion]
  MASH         = padj<0.05 AND |log2FC|>=1 in BOTH GSE316023 (MASH1) and GSE246328 (MASH2), same sign;
                 rank by min(|Wald stat| MASH1, |Wald stat| MASH2), descending
  FAST         = GSE107787 ZT-adjusted padj<0.05 AND |log2FC|>=1, AND clean ZT0/ZT4 subset padj<0.05 same sign;
                 rank by |t| of the ZT-adjusted model
  PPARA (secondary) = Phase 3 class 'PPARa_dependent' AND |WT 20v0 log2FC|>=1 AND sign(5v0)=sign(10v0)=sign(20v0);
                 rank by |Wald stat| 20v0
  size         = up to 50 up + 50 down per signature (sensitivity: 25); ties broken by gene symbol
Outputs: signatures/*.csv (gene, direction, rank, reference stats), frozen_signature_manifest.csv with SHA-256.
"""
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

P4D = Path(__file__).resolve().parents[1]
RS = P4D.parent
RE = RS / "roadmap_phase3" / "reference_effects"
SIG = P4D / "signatures"


def rd(name, base=RE):
    return pd.read_csv(base / f"{name}.csv.gz").set_index("gene")


def universe(kind):
    if kind == "alt":
        return set(pd.read_csv(RS / "roadmap_phase3/results/posthoc_audit/alt_universe.csv").vm35_symbol)
    u = pd.read_csv(RS / "roadmap_phase2/metadata/phase3_gene_universe.csv", index_col=0)
    return set(u.index[u.universe_main])


def refs(kind):
    if kind == "main":
        return dict(m1=rd("MASH1_GSE316023_GANvsChow_TA"), m2=rd("MASH2_GSE246328_GANvsChow_24w"),
                    fa=rd("FAST_GSE107787_FASTvsFED_ZTadj"), fc=rd("FASTclean_GSE107787_ZT0ZT4_muscleLow"),
                    p20=rd("PPARA_GSE282998_WT_20v0"), p5=rd("PPARA_GSE282998_WT_5v0"), p10=rd("PPARA_GSE282998_WT_10v0"))
    # alt-universe refits (Phase 3 post-hoc + Phase 4 V6); clean-subset & 5v0/10v0 exist on main universe only
    pa = RS / "roadmap_phase3/results/posthoc_audit"
    alt4 = RS / "roadmap_phase4/results/alt_universe_refs"
    f = pd.read_csv(alt4 / "FAST.csv.gz").set_index("gene").rename(columns={"SE": "lfcSE"})
    f["stat"] = f.log2FC / f.lfcSE
    p = pd.read_csv(alt4 / "PPARA_20v0.csv.gz").set_index("gene").rename(columns={"SE": "lfcSE"})
    p["stat"] = p.log2FC / p.lfcSE
    return dict(m1=rd("MASH1_altuniverse", pa), m2=rd("MASH2_altuniverse", pa), fa=f, fc=rd("FASTclean_GSE107787_ZT0ZT4_muscleLow"),
                p20=p, p5=rd("PPARA_GSE282998_WT_5v0"), p10=rd("PPARA_GSE282998_WT_10v0"))


def pick(df, n, key):
    out = []
    for d, sgn in (("up", 1), ("down", -1)):
        c = df[np.sign(df.eff) == sgn].copy()
        c["_k"] = -c[key]
        c = c.reset_index().sort_values(["_k", "gene"]).head(n)
        c["direction"], c["rank"] = d, np.arange(1, len(c) + 1)
        out.append(c.drop(columns="_k"))
    return pd.concat(out)


def build(kind="main", n=50):
    U = universe(kind)
    R = refs(kind)
    cl = pd.read_csv(RS / "roadmap_phase3/results/ppara/genotype_classification.csv.gz").set_index("gene")["class"]
    sigs, cand = {}, {}
    g = sorted(U & set(R["m1"].index) & set(R["m2"].index))
    a, b = R["m1"].loc[g], R["m2"].loc[g]
    ok = (a.padj < .05) & (b.padj < .05) & (np.sign(a.log2FC) == np.sign(b.log2FC)) & (a.log2FC.abs() >= 1) & (b.log2FC.abs() >= 1)
    df = pd.DataFrame({"eff": a.log2FC, "eff2": b.log2FC, "key": np.minimum(a.stat.abs(), b.stat.abs())})[ok]
    sigs["MASH"], cand["MASH"] = pick(df, n, "key"), df
    g = sorted(U & set(R["fa"].index) & set(R["fc"].index))
    a, b = R["fa"].loc[g], R["fc"].loc[g]
    ok = (a.padj < .05) & (a.log2FC.abs() >= 1) & (b.padj < .05) & (np.sign(a.log2FC) == np.sign(b.log2FC))
    df = pd.DataFrame({"eff": a.log2FC, "eff2": b.log2FC, "key": a.stat.abs()})[ok]
    sigs["FAST"], cand["FAST"] = pick(df, n, "key"), df
    g = sorted(U & set(R["p20"].index))
    a = R["p20"].loc[g]
    s5, s10 = np.sign(R["p5"].log2FC.reindex(g)), np.sign(R["p10"].log2FC.reindex(g))
    ok = (cl.reindex(g) == "PPARa_dependent") & (a.log2FC.abs() >= 1) & (s5 == np.sign(a.log2FC)) & (s10 == np.sign(a.log2FC))
    df = pd.DataFrame({"eff": a.log2FC, "eff2": R["p5"].log2FC.reindex(g), "key": a.stat.abs()})[ok]
    sigs["PPARA"], cand["PPARA"] = pick(df, n, "key"), df
    return sigs, cand


if __name__ == "__main__":
    SIG.mkdir(exist_ok=True)
    man = []
    for kind, n in [("main", 50), ("main", 25), ("alt", 50)]:
        sigs, cand = build(kind, n)
        for name, s in sigs.items():
            tag = f"{name}_{kind}_n{n}"
            f = SIG / f"{tag}.csv"
            s.to_csv(f, index=False)
            cand[name].to_csv(SIG / f"candidates_{name}_{kind}.csv")
            man.append(dict(signature=tag, role={"MASH": "primary", "FAST": "primary", "PPARA": "secondary"}[name],
                            universe=kind, n_per_direction_cap=n, n_up=int((s.direction == "up").sum()),
                            n_down=int((s.direction == "down").sum()), candidates_up=int((cand[name].eff > 0).sum()),
                            candidates_down=int((cand[name].eff < 0).sum()), file=f.name,
                            sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
    M = pd.DataFrame(man)
    M.to_csv(P4D / "frozen_signature_manifest.csv", index=False)
    mash, fast, pp = (set(pd.read_csv(SIG / f"{k}_main_n50.csv").gene) for k in ("MASH", "FAST", "PPARA"))
    print(M.to_string(index=False))
    print("overlap MASH-FAST", len(mash & fast), "MASH-PPARA", len(mash & pp), "FAST-PPARA", len(fast & pp))
    sys.exit(0)
