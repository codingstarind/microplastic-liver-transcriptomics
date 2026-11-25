"""Direction C: compare every sensitivity variant with the primary fit.

Gene level (per contrast): number of FDR < 0.05 genes, Spearman correlation of Wald
statistics with the primary fit, sign agreement and retention of the primary FDR genes.
Pathway level: fry direction and p-value for the pre-specified headline pathways.
Gene classes: share of primary FDR genes that are muscle-tracking, mitochondrial,
circadian/clock-output or among the 50 most abundant transcripts.

Usage: python analysis/robustness/summarise_robustness.py <primary_label>
"""

import json
import sys

import numpy as np
import pandas as pd

from mpliver.config import REPO_ROOT, TABLES_DIR

CONTRASTS = ["PE_in_CD", "PE_in_MD", "Diet_in_Veh", "Interaction"]
HEADLINE_PATHWAYS = ["SANDERSON_PPARA_TARGETS", "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                     "HALLMARK_MTORC1_SIGNALING", "HALLMARK_FATTY_ACID_METABOLISM",
                     "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_XENOBIOTIC_METABOLISM",
                     "HALLMARK_INFLAMMATORY_RESPONSE", "WP_PPAR_SIGNALING_PATHWAY"]
# Core clock and clock-output genes (fixed list, chosen before looking at results).
CIRCADIAN = ["Arntl", "Bmal1", "Clock", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2",
             "Nr1d1", "Nr1d2", "Dbp", "Tef", "Hlf", "Nfil3", "Bhlhe40", "Bhlhe41", "Ciart",
             "Rorc"]
DESEQ = REPO_ROOT / "results" / "deseq2"
PATHWAYS = REPO_ROOT / "results" / "pathways"


def load(label, contrast):
    path = DESEQ / label / f"results_{contrast}.csv.gz"
    return pd.read_csv(path).set_index("gene") if path.exists() else None


def compare_genes(primary, variant):
    common = primary.index.intersection(variant.index)
    p, v = primary.loc[common], variant.loc[common]
    sig = p.index[p["padj"] < 0.05]
    stat_col = "stat" if "stat" in v.columns else ("t" if "t" in v.columns else "F")
    lfc_col = "log2FC" if "log2FC" in v.columns else "logFC"
    retained = (v.loc[sig, "padj"] < 0.05).mean() if len(sig) else np.nan
    return {"n_fdr_genes": int((v["padj"] < 0.05).sum()),
            "spearman_stat_vs_primary": float(p["stat"].corr(v[stat_col], method="spearman"))
            if stat_col != "F" else np.nan,
            "spearman_lfc_vs_primary": float(p["log2FC"].corr(v[lfc_col], method="spearman")),
            "sign_agreement_primary_fdr_genes":
                float((np.sign(p.loc[sig, "log2FC"]) == np.sign(v.loc[sig, lfc_col])).mean())
                if len(sig) else np.nan,
            "retention_primary_fdr_genes": float(retained)}


def pathway_calls(label):
    path = PATHWAYS / label / "limma_fry_camera.csv"
    if not path.exists():
        return None
    t = pd.read_csv(path)
    t = t[t["pathway"].isin(HEADLINE_PATHWAYS)]
    return t.set_index(["pathway", "contrast"])[["fry_direction", "fry_pvalue", "fry_fdr"]]


def main() -> int:
    primary_label = sys.argv[1]
    out = REPO_ROOT / "results" / "robustness" / primary_label
    out.mkdir(parents=True, exist_ok=True)
    variants = sorted({p.parent.name if p.parent.parent.name == "deseq2" else
                       f"loo/{p.parent.name}"
                       for p in DESEQ.glob(f"**/{primary_label}__*/contrast_summary.csv")})

    gene_rows = []
    for contrast in CONTRASTS:
        primary = load(primary_label, contrast)
        for v in variants:
            res = load(v, contrast)
            if res is not None:
                gene_rows.append({"variant": v, "contrast": contrast,
                                  **compare_genes(primary, res)})
        for method in ["edgeR_QL", "limma_voom"]:
            path = DESEQ / primary_label / f"{method}_{contrast}.csv.gz"
            if path.exists():
                gene_rows.append({"variant": method, "contrast": contrast,
                                  **compare_genes(primary, pd.read_csv(path).set_index("gene"))})
        gene_rows.append({"variant": "PRIMARY", "contrast": contrast,
                          "n_fdr_genes": int((primary["padj"] < 0.05).sum())})
    genes = pd.DataFrame(gene_rows)
    genes.round(4).to_csv(out / "gene_level_robustness.csv", index=False)

    primary_calls = pathway_calls(primary_label)
    path_rows = []
    for v in ["PRIMARY", *variants]:
        calls = primary_calls if v == "PRIMARY" else pathway_calls(v)
        if calls is None:
            continue
        for (pw, contrast), row in calls.iterrows():
            ref = primary_calls.loc[(pw, contrast)]
            path_rows.append({"variant": v, "pathway": pw, "contrast": contrast,
                              "fry_direction": row["fry_direction"],
                              "fry_pvalue": row["fry_pvalue"], "fry_fdr": row["fry_fdr"],
                              "same_direction_as_primary":
                                  row["fry_direction"] == ref["fry_direction"]})
    pathways = pd.DataFrame(path_rows)
    pathways.round(5).to_csv(out / "pathway_robustness.csv", index=False)

    # Gene classes among the primary FDR genes.
    muscle = set(pd.read_csv(REPO_ROOT / "results" / "methodology" / ".." / "tables" /
                             "muscle_tracking_genes.csv").iloc[:, 0])
    top50 = set(pd.read_csv(TABLES_DIR / "top_expressed_genes.csv").iloc[:30, 0])
    class_rows = []
    for contrast in CONTRASTS:
        res = load(primary_label, contrast).dropna(subset=["pvalue"])
        sig = res.index[res["padj"] < 0.05]
        for name, members in {"muscle_tracking": muscle,
                              "mitochondrial": {g for g in res.index if g.startswith("mt-")},
                              "circadian_clock": set(CIRCADIAN),
                              "top30_abundant": top50}.items():
            class_rows.append({"contrast": contrast, "gene_class": name,
                               "n_in_fdr_genes": len(set(sig) & members),
                               "frac_of_fdr_genes": len(set(sig) & members) / max(len(sig), 1),
                               "frac_of_tested_genes": len(set(res.index) & members) / len(res)})
    pd.DataFrame(class_rows).round(4).to_csv(out / "gene_class_checks.csv", index=False)

    clock = pd.concat({c: load(primary_label, c).reindex(CIRCADIAN)[["log2FC", "lfcSE", "padj"]]
                       for c in CONTRASTS}, axis=1).dropna(how="all")
    clock.round(4).to_csv(out / "circadian_genes_by_contrast.csv")

    summary = {"primary": primary_label, "n_variants": len(variants),
               "variants": variants}
    (out / "robustness_index.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(genes[genes["variant"].isin(["PRIMARY", "edgeR_QL", "limma_voom"]) |
                ~genes["variant"].str.startswith("loo/")].round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
