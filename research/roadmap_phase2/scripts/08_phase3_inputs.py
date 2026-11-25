"""Step 08: fix the Phase 3 inputs BEFORE any reference signature or PE comparison.

Nothing here looks at GSE294180 PE effects. Outputs:
  metadata/muscle_exclusion_genes.csv   fixed muscle rule (independent tissue data)
  metadata/phase3_gene_universe.csv     main + sensitivity universes
  metadata/phase3_eligible_samples.csv  one row per reference sample x contrast, with flags
  DATA_ROOT/GSE316023/prepared/implied_length.csv.gz  per-sample effective-length matrix
      (count / (FPKM * N / 1e9)), for an optional tximport-style length offset

Muscle rule (pre-specified): using GSE107787's own gastrocnemius samples (36) and the
livers whose canonical muscle score is < 2 RPKM, a gene is muscle-enriched if
median muscle RPKM >= 10 AND >= 20 x (median clean-liver RPKM + 1). The 8 canonical
sarcomere markers are always excluded. Thresholds were set before inspecting which
genes they select and are not tuned against any PE result.
"""
import json

import numpy as np
import pandas as pd

from gene_mapping import collapse, map_ids, own_universe
from p2_common import DATA_ROOT, REPO_P2, file_hashes

CANON = ["Acta1", "Ckm", "Myh1", "Myh4", "Mylpf", "Tnnt3", "Myl1", "Actn3"]
PRIMARY = {"GSE316023": "MASH_primary", "GSE107787": "FASTING_primary",
           "GSE117134": "CIRCADIAN_primary", "GSE282998": "PPARA_primary"}


def muscle_genes() -> pd.DataFrame:
    x = pd.read_excel(DATA_ROOT / "GSE107787" / "GSE107787_RNAseq_expression_results.xlsx").set_index("Gene")
    mp = map_ids(x.index, "symbol")
    key = x.index.astype(str) + np.where(x.index.duplicated(keep=False),
                                          "#row" + pd.Series(range(len(x))).astype(str).values, "")
    x.index, mp["original_id"] = key, key
    X, _ = collapse(x, mp, "max_mean")
    mus = X[[c for c in X if c.startswith("Muscle")]]
    liv = X[[c for c in X if c.startswith("Liver")]]
    clean = liv.columns[liv.loc[CANON].sum() < 2]
    mm, lm = mus.median(axis=1), liv[clean].median(axis=1)
    rule = (mm >= 10) & (mm >= 20 * (lm + 1))
    genes = sorted(set(rule[rule].index) | set(CANON))
    return pd.DataFrame({"vm35_symbol": genes,
                         "muscle_median_rpkm": mm.reindex(genes).round(2).values,
                         "clean_liver_median_rpkm": lm.reindex(genes).round(3).values,
                         "reason": ["canonical+rule" if g in CANON and rule.get(g, False) else
                                    "canonical" if g in CANON else "rule" for g in genes]}), len(clean)


def implied_length() -> dict:
    c = pd.read_csv(DATA_ROOT / "GSE316023" / "GSE316023_COUNT_NCOMMS-25-40823-T.csv.gz", index_col=0)
    f = pd.read_csv(DATA_ROOT / "GSE316023" / "GSE316023_FPKM_NCOMMS-25-40823-T.csv.gz", index_col=0)
    L = (c / (f * c.sum() / 1e9)).replace([np.inf, -np.inf], np.nan)
    L = L[~c.index.duplicated(keep=False)]           # unique original IDs only
    L = L.apply(lambda r: r.fillna(r.mean()), axis=1)  # zero-count cells: gene's mean length
    ready = pd.read_csv(DATA_ROOT / "GSE316023" / "prepared" / "expression_ready.csv.gz", index_col=0)
    mp = pd.read_csv(REPO_P2 / "datasets" / "GSE316023" / "gene_mapping.csv.gz")
    kept = mp[mp["kept_in_ready"]].set_index("original_id")["vm35_symbol"]
    kept = kept[kept.index.isin(L.index)]
    out = L.loc[kept.index]
    out.index = kept.values
    out = out.reindex(ready.index)[ready.columns]
    p = DATA_ROOT / "GSE316023" / "prepared" / "implied_length.csv.gz"
    out.round(2).to_csv(p)
    md5, sha = file_hashes(p)
    return dict(path=str(p), md5=md5, sha256=sha, genes_with_length=int(out.notna().all(axis=1).sum()),
                genes_missing_length=int(out.isna().any(axis=1).sum()))


if __name__ == "__main__":
    mus, n_clean = muscle_genes()
    mus.to_csv(REPO_P2 / "metadata" / "muscle_exclusion_genes.csv", index=False)
    mset = set(mus["vm35_symbol"])
    _, own_filt = own_universe()
    mem = pd.read_csv(REPO_P2 / "metadata" / "gene_universe_membership.csv.gz", index_col=0)
    u = pd.DataFrame(index=mem.index)
    u["muscle_excluded"] = u.index.isin(mset)
    u["universe_main"] = mem["in_all_primary_present"] & ~u["muscle_excluded"]
    u["universe_sensitivity"] = mem["in_all_primary_detected"] & ~u["muscle_excluded"]
    u.to_csv(REPO_P2 / "metadata" / "phase3_gene_universe.csv")

    ms = pd.read_csv(REPO_P2 / "metadata" / "sample_metadata_master.csv", dtype={"sample_id": str})
    qc = pd.read_csv(REPO_P2 / "qc" / "qc_sample_metrics.csv", dtype={"sample_id": str})
    qc = qc.set_index(["dataset", "sample_id"])
    rows = []
    for _, r in ms.iterrows():
        for role in str(r["contrast_role"]).split(";"):
            if role.startswith("not_in_contrast"):
                continue
            cid, arm = role.split(":")
            q = qc.loc[(r["dataset"], r["sample_id"])] if (r["dataset"], r["sample_id"]) in qc.index else None
            ms_score = float(q["muscle_score"]) if q is not None else np.nan
            muscle_high = bool(ms_score > (50 if r["dataset"] == "GSE282998" else 10))
            clean_sens = True
            if cid == "FASTING_primary":      # key sensitivity: ZT0 + ZT4, muscle-low livers only
                clean_sens = r["zt"] in ("ZT0", "ZT4") and not muscle_high
            elif cid.startswith("PPARA") or cid.startswith("CIRCADIAN_sens"):
                clean_sens = not muscle_high
            rows.append(dict(contrast_id=cid, arm=arm, dataset=r["dataset"], gsm=r["gsm"],
                             sample_id=r["sample_id"], group=r["group"], zt=r["zt"], dose=r["dose"],
                             muscle_score=round(ms_score, 2), muscle_high=muscle_high,
                             corr_outlier_flag=bool(q["corr_outlier_flag"]) if q is not None else False,
                             eligible_main=True, in_clean_sensitivity=clean_sens))
    el = pd.DataFrame(rows)
    el.to_csv(REPO_P2 / "metadata" / "phase3_eligible_samples.csv", index=False)
    il = implied_length()
    summary = dict(clean_livers_for_rule=n_clean, muscle_genes=len(mset),
                   muscle_genes_in_GSE294180_filtered=int(len(mset & own_filt)),
                   universe_main=int(u["universe_main"].sum()), universe_sensitivity=int(u["universe_sensitivity"].sum()),
                   eligible_rows=len(el), implied_length=il)
    (REPO_P2 / "metadata" / "phase3_inputs_manifest.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(el.groupby(["contrast_id", "arm"]).agg(n=("gsm", "size"), muscle_high=("muscle_high", "sum"),
                                                 clean_sens=("in_clean_sensitivity", "sum")).to_string())
