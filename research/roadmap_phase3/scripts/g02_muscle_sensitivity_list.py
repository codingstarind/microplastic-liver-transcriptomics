"""POST-HOC EXPLORATORY: build the muscle-restricted list and alternative universe exactly as
pre-registered in posthoc_muscle_rule_prereg.md (hash in results/posthoc_audit/muscle_prereg_lock.json).
Step 1 of 3 (2: g02_muscle_sensitivity_fit.R, 3: g02_muscle_sensitivity_eval.py)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

P3 = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P3.parent / "roadmap_phase2" / "scripts"))
from gene_mapping import collapse, map_ids  # noqa: E402

P2 = P3.parent / "roadmap_phase2"
OUT = P3 / "results" / "posthoc_audit"
CANON = ["Acta1", "Ckm", "Myh1", "Myh4", "Mylpf", "Tnnt3", "Myl1", "Actn3"]

x = pd.read_excel("D:/mpliver_work/phase2_refs/GSE107787/GSE107787_RNAseq_expression_results.xlsx").set_index("Gene")
mp = map_ids(x.index, "symbol")
key = x.index.astype(str) + np.where(x.index.duplicated(keep=False), "#row" + pd.Series(range(len(x))).astype(str).values, "")
x.index, mp["original_id"] = key, key
X, _ = collapse(x, mp, "max_mean")
mus = X[[c for c in X if c.startswith("Muscle")]]
liv = X[[c for c in X if c.startswith("Liver")]]
clean = liv.columns[liv.loc[CANON].sum() < 2]
assert len(clean) == 21
mm, lm = mus.median(axis=1), liv[clean].median(axis=1)
rule = (mm >= 10) & (mm >= 100 * (lm + 1)) & (lm < 1)
genes = sorted(set(rule[rule].index) | set(CANON))
pd.DataFrame({"vm35_symbol": genes, "muscle_median_rpkm": mm.reindex(genes).round(2).values,
              "clean_liver_median_rpkm": lm.reindex(genes).round(3).values}).to_csv(OUT / "muscle_restricted_genes.csv", index=False)
mem = pd.read_csv(P2 / "metadata" / "gene_universe_membership.csv.gz", index_col=0)
orig = pd.read_csv(P2 / "metadata" / "phase3_gene_universe.csv", index_col=0)
alt = mem.index[mem["in_all_primary_present"] & ~mem.index.isin(genes)]
readmitted = sorted(set(alt) - set(orig.index[orig.universe_main]))
pd.Series(sorted(alt), name="vm35_symbol").to_csv(OUT / "alt_universe.csv", index=False)
pd.Series(readmitted, name="vm35_symbol").to_csv(OUT / "readmitted_genes.csv", index=False)
print("muscle-restricted:", len(genes), "| alt universe:", len(alt), "| re-admitted:", len(readmitted))
print("watch genes re-admitted:", {g: g in readmitted for g in ["Pdk4", "Plin4", "Col1a1", "Col1a2", "Lpl", "Cpt1b", "Cav1", "Vldlr", "Mgp"]})
