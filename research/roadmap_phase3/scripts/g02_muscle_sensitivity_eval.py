"""POST-HOC EXPLORATORY (step 3/3): positive-control metrics on the pre-registered alternative universe.
Also re-states the locked primary numbers for side-by-side reporting. Output: results/posthoc_audit/muscle_sensitivity.json"""
import json

import pandas as pd
from p3lib import P3, compare, ref, target, universe

OUT = P3 / "results" / "posthoc_audit"
alt = set(pd.read_csv(OUT / "alt_universe.csv").vm35_symbol)
tgt = target("Diet_in_Veh")


def alt_ref(lab):
    return pd.read_csv(OUT / f"{lab}_altuniverse.csv.gz").set_index("gene")[["log2FC", "padj"]]


m1a, m2a = alt_ref("MASH1"), alt_ref("MASH2")
m1, m2 = ref("MASH1_GSE316023_GANvsChow_TA"), ref("MASH2_GSE246328_GANvsChow_24w")
keys = ["n_genes", "spearman", "n_sig_both", "sign_agree_sig_both", "n_deg_overlap", "jaccard"]
res = {
    "status": "POST-HOC EXPLORATORY; does not replace locked results",
    "primary_locked_309_list": {
        "MASH1_vs_Diet": {k: compare(tgt, m1, universe("main"), "t", "r")[k] for k in keys},
        "MASH2_vs_Diet": {k: compare(tgt, m2, universe("main"), "t", "r")[k] for k in keys},
        "MASH1_vs_MASH2": {k: compare(m1, m2, universe("main"), "a", "b")[k] for k in keys}},
    "alt_88_muscle_restricted": {
        "MASH1_vs_Diet": {k: compare(tgt, m1a, alt, "t", "r")[k] for k in keys},
        "MASH2_vs_Diet": {k: compare(tgt, m2a, alt, "t", "r")[k] for k in keys},
        "MASH1_vs_MASH2": {k: compare(m1a, m2a, alt, "a", "b")[k] for k in keys}},
}
watch = ["Pdk4", "Plin4", "Col1a1", "Col1a2", "Lpl", "Cpt1b", "Cav1", "Vldlr"]
res["readmitted_watch_genes"] = {g: dict(MASH1=float(m1a.loc[g, "log2FC"]) if g in m1a.index else None,
                                         MASH2=float(m2a.loc[g, "log2FC"]) if g in m2a.index else None,
                                         Diet_in_Veh=float(tgt.loc[g, "log2FC"]) if g in tgt.index else None) for g in watch}
(OUT / "muscle_sensitivity.json").write_text(json.dumps(res, indent=2))
print(json.dumps(res, indent=1))
