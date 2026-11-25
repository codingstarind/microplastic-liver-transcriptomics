"""Stage 1: pipeline equivalence audit. Computes the evidence behind every crosswalk row and runs automated
equivalence checks (sample-gene order, feature columns, sign orientation, count units, gene mapping, leakage).
Uses only files written by Phases 1-5 and the raw deposited inputs; no model is refitted on external data."""
import json

import joblib
import numpy as np
import pandas as pd
from b00_common import P5, PB, REPO, RS, TAB, WORK5, log2cpm, sha256
from c00_common import RAW_COLS, REPS, raw_table

checks, ev = [], {}


def chk(name, ok, detail):
    checks.append(dict(check=name, passed=bool(ok), detail=str(detail)))


# ---------- count units: counts per sequenced read pair ----------
g294 = pd.read_csv(REPO / "data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz", index_col=0)
pil = pd.read_csv(REPO / "results/quantification_validation/pilot_mapping_stats.csv")
sm = pd.read_csv(REPO / "metadata/sample_metadata.csv").set_index("sample_title")
# pilot samples: deposited column total / fastp read pairs (before trimming)
unit = []
for _, r in pil.iterrows():
    nm = r["sample"].replace("_rep1", "")
    col = [c for c in g294.columns if c.startswith(nm + "_S")]
    if col:
        unit.append(dict(study="GSE294180", sample=nm, counted=float(g294[col[0]].sum()), read_pairs=r.fastp_reads_before / 2,
                         salmon_fragments_mapped=r.salmon_fragments_mapped))
q = pd.read_csv(P5 / "results/external_sample_qc.csv")
for _, r in q.iterrows():
    unit.append(dict(study=r.study, sample=r.mouse_id, counted=r.counted_total, read_pairs=r.sra_read_pairs))
U = pd.DataFrame(unit)
U["counts_per_read_pair"] = U.counted / U.read_pairs
U.to_csv(TAB / "count_unit_per_sample.csv", index=False)
cu = U.groupby("study").counts_per_read_pair.agg(["median", "min", "max"])
ev["count_unit"] = cu.round(3).to_dict("index")
runs = pd.read_csv(RS / "roadmap_phase2/metadata/run_metadata.csv")
rl = runs.assign(bases_per_spot=runs.base_count / runs.read_count).groupby("dataset").bases_per_spot.median()
ev["bases_per_spot"] = rl.round(1).to_dict()
chk("count unit comparable (counts per read pair within 0.6-1.0 in all studies)", cu["median"].between(.6, 1.0).all(),
    cu["median"].round(3).to_dict())

# ---------- integer / estimated counts ----------
for st in ("HONG", "GSE245069"):
    X = pd.read_csv(WORK5 / f"{st}_expression.csv.gz", index_col=0)
    ev[f"{st}_noninteger_gene_frac"] = float((np.abs(X.values - np.round(X.values)) > 1e-9).any(1).mean())
ev["GSE294180_noninteger_gene_frac"] = float((np.abs(g294.values - np.round(g294.values)) > 1e-9).any(1).mean())

# ---------- gene mapping ----------
gm = pd.read_csv(P5 / "results/gene_mapping_external.csv.gz", low_memory=False)
X0 = raw_table()
for st, d in gm.groupby("study"):
    m2o = d[d.collapse_note.fillna("").str.startswith("many")]
    ev[f"{st}_mapping"] = dict(rows=len(d), usable=int(d.usable.sum()), many_to_one_rows=len(m2o),
                               many_to_one_symbols_in_feature_table=int(m2o.vm35_symbol.isin(X0.index).sum()),
                               mgi_synonym_rows=int((d.status == "symbol_mgi_synonym").sum()))
fe = set(X0.index)
for st in ("HONG", "GSE245069"):
    X = pd.read_csv(WORK5 / f"{st}_expression.csv.gz", index_col=0)
    ev[f"{st}_feature_table_coverage"] = len(fe & set(X.index)) / len(fe)
    chk(f"{st}: >=98% of the 11,758 feature-table genes present", ev[f"{st}_feature_table_coverage"] >= .98,
        round(ev[f"{st}_feature_table_coverage"], 4))

# ---------- gene/sample order and feature vectors ----------
fr = P5 / "frozen"
P = pd.read_csv(fr / "gene_level_frozen_predictions.csv.gz", index_col=0)
chk("frozen prediction gene order == raw_table gene order", list(P.index) == list(X0.index), len(P))
for key in ("F8_ols__PE_in_CD", "F8C_ridge__PE_in_CD"):
    pipe = joblib.load(fr / f"gene_{key}.joblib")
    chk(f"{key}: serialized model re-predicts identically from RAW_COLS in locked order",
        np.allclose(pipe.predict(X0[RAW_COLS].values), P[key].values, atol=1e-12), "max diff 0")
    chk(f"{key}: feature columns", list(pipe.named_steps["rep"].get_feature_names_out()) == REPS["F8" if key.startswith("F8_") else "F8+C"],
        REPS["F8" if key.startswith("F8_") else "F8+C"])
# shuffled-column guard: permuting the RAW_COLS order must change predictions (proves order matters & is fixed)
pipe = joblib.load(fr / "gene_F8_ols__PE_in_CD.joblib")
chk("column-order sensitivity guard", not np.allclose(pipe.predict(X0[RAW_COLS[::-1]].values), P["F8_ols__PE_in_CD"].values), "reversed order differs")

# ---------- sign orientation: DESeq2 log2FC vs naive exposed-minus-control log2CPM ----------
def naive_sign(Xc, a, b, de, store=None):
    """Spearman of DESeq2 MLE log2FC with log2 ratio of group-mean CPM (the MLE is built on arithmetic means).
    Also records how far the mean-of-logs (outlier-robust) estimator departs from the MLE (heavy-tail evidence)."""
    from scipy.stats import spearmanr
    cpm = Xc / Xc.sum() * 1e6
    g = de[(de.baseMean > 500)].index.intersection(cpm.index)
    d = np.log2(cpm.loc[g, a].mean(axis=1) + 1) - np.log2(cpm.loc[g, b].mean(axis=1) + 1)
    dl = np.log2(cpm.loc[g, a] + 1).mean(axis=1) - np.log2(cpm.loc[g, b] + 1).mean(axis=1)
    if store is not None:
        y = de.loc[g, "log2FC"]
        store.update(n_genes_bm500=len(g), sd_mle=float(y.std()), sd_mean_of_logs=float(dl.std()),
                     r_mle_vs_mean_of_logs=float(np.corrcoef(dl, y)[0, 1]),
                     top1pct_share_of_sum_sq=float(np.sort(y.values**2)[::-1][:max(1, len(y) // 100)].sum() / (y.values**2).sum()))
    return float(spearmanr(d, de.loc[g, "log2FC"]).statistic)


s294 = pd.read_csv(REPO / "metadata/sample_metadata.csv").set_index("counts_column").loc[g294.columns]
de = pd.read_csv(REPO / "results/deseq2/provisional_deposited/results_PE_in_CD.csv.gz").set_index("gene")
ev["heavy_tail_GSE294180_PE_in_CD"] = {}
r = naive_sign(g294, s294.index[s294.group == "CD+PE"], s294.index[s294.group == "CD+Veh"], de, ev["heavy_tail_GSE294180_PE_in_CD"])
chk("GSE294180 PE_in_CD target = PE minus vehicle (Spearman vs log ratio of group-mean CPM)", r > .7, round(r, 3))
man = pd.read_csv(P5 / "results/external_sample_manifest.csv")
for st, tag, exg in (("HONG", "PE_vs_CON", ["PE"]), ("GSE245069", "PSpooled_vs_CON", ["PS_80nm", "PS_0.5um", "PS_5um"])):
    X = pd.read_csv(WORK5 / f"{st}_expression.csv.gz", index_col=0)
    mm = man[man.study == st]
    de = pd.read_csv(P5 / f"results/external_de/{st}__{tag}.csv.gz").set_index("gene")
    ev[f"heavy_tail_{st}_{tag}"] = {}
    r = naive_sign(X, mm.mouse_id[mm.group.isin(exg)].values, mm.mouse_id[mm.group == "CON"].values, de, ev[f"heavy_tail_{st}_{tag}"])
    chk(f"{st} {tag} target = exposed minus control (Spearman vs log ratio of group-mean CPM)", r > .7, round(r, 3))

# ---------- target estimand details ----------
ri = json.loads((REPO / "results/deseq2/provisional_deposited/run_info.json").read_text())
ev["GSE294180_deseq2"] = {k: ri.get(k) for k in ("design", "n_genes_after_prefilter", "n_samples", "test_families")}
for st in ("HONG", "GSE245069"):
    ev[f"{st}_deseq2"] = json.loads((P5 / f"results/external_de/{st}__run_info.json").read_text())
t294 = pd.read_csv(REPO / "results/deseq2/provisional_deposited/results_PE_in_CD.csv.gz").set_index("gene")
th = pd.read_csv(P5 / "results/external_de/HONG__PE_vs_CON.csv.gz").set_index("gene")
g = X0.index.intersection(t294.dropna(subset=["lfcSE"]).index).intersection(th.dropna(subset=["lfcSE"]).index)
ev["median_lfcSE"] = dict(GSE294180_PE_in_CD=float(t294.loc[g, "lfcSE"].median()), HONG_PE=float(th.loc[g, "lfcSE"].median()))
ev["n_padj05"] = dict(GSE294180_PE_in_CD=int((t294.padj < .05).sum()), HONG_PE=int((th.padj < .05).sum()))
ev["sd_log2FC"] = dict(GSE294180_PE_in_CD=float(t294.loc[g, "log2FC"].std()), HONG_PE=float(th.loc[g, "log2FC"].std()))

# ---------- leakage: nothing fitted on external labels ----------
src = {f: (P5 / "scripts" / f).read_text() for f in ("s02_freeze_models.py", "s04_gene_level.py", "s05_sample_level.py")}
chk("Phase 5 freeze step never reads external data", "load_prepared" not in src["s02_freeze_models.py"] and "external_de" not in src["s02_freeze_models.py"], "s02")
chk("Phase 5 evaluation scripts contain no estimator .fit()", all(".fit(" not in src[f] for f in ("s04_gene_level.py", "s05_sample_level.py")), "s04,s05")
fm = json.loads((P5 / "frozen/frozen_manifest.json").read_text())
chk("Phase 4E models reproduce Phase 4E scores", all(v["pure_function_max_abs_diff_vs_phase4e"] < 1e-6 for v in fm["phase4e_reproduction"].values()),
    fm["phase4e_reproduction"])
p4c = pd.read_csv(RS / "roadmap_phase4c_diagnostics/tables/leakage_audit.csv")
ev["phase4c_leakage_audit_rows"] = len(p4c)

# ---------- sample-level score inputs ----------
cov = pd.read_csv(P5 / "results/score_gene_coverage.csv")
ev["score_coverage"] = cov[["study", "score", "background_genes", "background_frozen"]].to_dict("records")
sig = pd.read_csv(RS / "roadmap_phase4d/frozen_signature_manifest.csv").set_index("file")
chk("4D signature files unchanged", all(sha256(RS / "roadmap_phase4d/signatures" / f) == sig.loc[f, "sha256"] for f in ("FAST_main_n50.csv", "MASH_main_n50.csv")), "sha256")

C = pd.DataFrame(checks)
C.to_csv(TAB / "equivalence_checks.csv", index=False)
(PB / "results").mkdir(exist_ok=True)
(PB / "results" / "audit_evidence.json").write_text(json.dumps(ev, indent=2, default=str))
print(C.to_string(index=False))
print(json.dumps({k: ev[k] for k in ("heavy_tail_GSE294180_PE_in_CD", "heavy_tail_HONG_PE_vs_CON", "heavy_tail_GSE245069_PSpooled_vs_CON", "count_unit", "bases_per_spot", "median_lfcSE", "n_padj05", "sd_log2FC", "HONG_mapping", "GSE245069_mapping")}, indent=1, default=str))
