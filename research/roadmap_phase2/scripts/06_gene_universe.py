"""Step 06: candidate gene universes and key-gene coverage (reporting only).

Universe U0 = GSE294180 filtered genes (13,777; >=10 counts in >=4 samples).
For each reference: genes present in its ready matrix, and genes passing the
pre-specified detection rule in its primary-contrast samples.
Writes metadata/gene_universe_membership.csv.gz (one row per U0 gene) and
metadata/key_gene_coverage.csv. Phase 3 must fix the final universe before
looking at PE-vs-reference results.
"""
import pandas as pd
from gene_mapping import detected_genes, own_universe
from p2_common import DATA_ROOT, REPO_P2

REFS = {
    "GSE316023": ("MASH_primary", "counts"), "GSE246328": ("MASH_secondary", "counts"),
    "GSE107787": ("FASTING_primary", "rpkm"), "GSE225695": ("FASTING_secondary", "counts"),
    "GSE117134": ("CIRCADIAN_primary", "rpkm"), "GSE282998": ("PPARA_primary", "counts"),
}
KEY = {
    "lipogenesis/SREBP1": ["Fasn", "Acly", "Acaca", "Me1", "Elovl6", "Scd1", "Srebf1", "Pnpla3", "Mlxipl"],
    "cholesterol/SREBP2": ["Hmgcr", "Hmgcs1", "Sqle", "Ldlr", "Cyp51", "Fdps", "Idi1", "Pcsk9", "Insig1"],
    "PPARa targets": ["Cyp4a14", "Cyp4a10", "Acot1", "Ehhadh", "Cpt1a", "Hmgcs2", "Pdk4", "Angptl4", "Fgf21",
                      "Plin2", "Plin4", "Cd36", "Acox1", "Ppara"],
    "clock": ["Bmal1", "Clock", "Per1", "Per2", "Per3", "Cry1", "Nr1d1", "Nr1d2", "Dbp", "Tef", "Nfil3"],
    "MASH/fibrosis/inflammation": ["Col1a1", "Col3a1", "Timp1", "Lgals3", "Cd68", "Gpnmb", "Trem2", "Ccl2", "Mmp12"],
    "skeletal muscle (technical)": ["Acta1", "Ckm", "Myh1", "Myh4", "Mylpf", "Tnnt3"],
}

if __name__ == "__main__":
    own_all, own_filt = own_universe()
    u0 = sorted(own_filt)
    mem = pd.DataFrame(index=pd.Index(u0, name="vm35_symbol"))
    cov = []
    for acc, (prefix, kind) in REFS.items():
        m = pd.read_csv(REPO_P2 / "datasets" / acc / "sample_metadata.csv", dtype={"sample_id": str})
        x = pd.read_csv(DATA_ROOT / acc / "prepared" / "expression_ready.csv.gz", index_col=0)
        x.columns = x.columns.astype(str)
        s = list(m.loc[m["contrast_role"].astype(str).str.contains(prefix), "sample_id"])
        det = set(detected_genes(x, s, kind))
        mem[f"{acc}_present"] = mem.index.isin(x.index)
        mem[f"{acc}_detected"] = mem.index.isin(det)
        for cls, genes in KEY.items():
            for g in genes:
                cov.append(dict(gene=g, program=cls, dataset=acc, present=g in x.index, detected=g in det,
                                in_GSE294180_filtered=g in own_filt))
    core = ["GSE316023", "GSE107787", "GSE117134", "GSE282998"]
    mem["in_all_primary_present"] = mem[[f"{a}_present" for a in core]].all(axis=1)
    mem["in_all_primary_detected"] = mem[[f"{a}_detected" for a in core]].all(axis=1)
    mem.to_csv(REPO_P2 / "metadata" / "gene_universe_membership.csv.gz")
    cv = pd.DataFrame(cov)
    piv = cv.pivot_table(index=["program", "gene", "in_GSE294180_filtered"], columns="dataset",
                         values="detected", aggfunc="first").reset_index()
    piv.to_csv(REPO_P2 / "metadata" / "key_gene_coverage.csv", index=False)
    print("U0 (GSE294180 filtered):", len(u0))
    for c in mem.columns:
        print(f"  {c}: {int(mem[c].sum())}")
    print(piv.to_string())
