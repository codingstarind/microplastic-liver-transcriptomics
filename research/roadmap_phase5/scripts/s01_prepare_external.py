"""Step 1: audit + prepare independent microplastic-liver datasets (no effect estimation here).

* Hong lab: 4 GEO series (PE/PS/PP/PB) sharing ONE control group. Verifies byte-identical control columns and
  identical control SRA read counts, then builds ONE 15-mouse matrix (3 CON + 4 x 3 exposed).
* GSE245069: 16 mice, deposited read_count (estimated counts; kept as float here, rounded only inside DESeq2).
* GSE277772: 6 mice, normalised values only (TMM/FPKM-like). Kept for rank-based exploratory scoring only.
Identifiers are mapped to GENCODE vM35 symbols with the Phase 2 mapper (MGI synonym rescue).
Outputs: D:/mpliver_work/phase5/<study>_expression.csv.gz, results/external_sample_manifest.csv,
results/external_sample_qc.csv, results/gene_mapping_external.csv, ../dataset_provenance.csv."""
import gzip
import re
import sys

import numpy as np
import pandas as pd
from p5_common import DROOT, P5, RES, RS, STUDIES, WORK, sha256

sys.path.insert(0, str(RS / "roadmap_phase2" / "scripts"))
from gene_mapping import collapse, map_ids  # noqa: E402

WORK.mkdir(parents=True, exist_ok=True)
RES.mkdir(exist_ok=True)
runs = pd.read_csv(RS / "roadmap_phase2/metadata/run_metadata.csv")


def soft_titles(acc):
    out, cur = {}, None
    with gzip.open(DROOT / acc / f"{acc}_family.soft.gz", "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("^SAMPLE"):
                cur = line.split("=")[1].strip()
            elif line.startswith("!Sample_title") and cur:
                out[cur] = line.split("=", 1)[1].strip()
    return out


manifest, mapping_rows, prov = [], [], []

# ---------------- Hong lab ----------------
HS = {"GSE261893": "PE", "GSE261725": "PS", "GSE286116": "PP", "GSE262399": "PB"}
tabs = {acc: pd.read_csv(DROOT / acc / f"{acc}_readcount_Con_{pol}.txt.gz", sep="\t", index_col=0) for acc, pol in HS.items()}
con = {acc: t[["CON1", "CON2", "CON3"]] for acc, t in tabs.items()}
ref = con["GSE261893"]
identical = {acc: bool(c.index.equals(ref.index) and np.array_equal(c.values, ref.values)) for acc, c in con.items()}
assert all(identical.values()), identical  # shared controls: same 3 libraries deposited 4 times
ctrl_reads = {}
for acc in HS:
    r = runs[runs.dataset == acc]
    t = soft_titles(acc)
    cg = sorted([g for g, ti in t.items() if "normal food" in ti], key=lambda g: t[g])
    ctrl_reads[acc] = [float(r[r.gsm == g].read_count.iloc[0]) for g in cg]
assert len({tuple(v) for v in ctrl_reads.values()}) == 1, ctrl_reads  # identical read counts per rep across series
H = pd.concat([ref.rename(columns=lambda c: f"HONG_CON_{c[-1]}")] +
              [t[[c for c in t.columns if not c.startswith("CON")]].rename(columns=lambda c, p=pol: f"HONG_{p}_{c.split('-')[1]}")
               for (acc, pol), t in zip(HS.items(), tabs.values())], axis=1)
assert np.allclose(H.values, np.round(H.values))
m = map_ids(H.index, "symbol")
Hc, mm = collapse(H, m, "sum")
mapping_rows.append(mm.assign(study="HONG"))
Hc.to_csv(WORK / "HONG_expression.csv.gz")
for acc, pol in HS.items():
    t = soft_titles(acc)
    r = runs[runs.dataset == acc]
    for g, ti in sorted(t.items(), key=lambda kv: kv[1]):
        rep = re.search(r"rep(\d)", ti).group(1)
        is_ctrl = "normal food" in ti
        manifest.append(dict(study="HONG", mouse_id=f"HONG_CON_{rep}" if is_ctrl else f"HONG_{pol}_{rep}", gsm=g, series=acc,
                             title=ti, group="CON" if is_ctrl else pol, polymer="none" if is_ctrl else pol,
                             size="NOT STATED", dose="NOT STATED", route="in food", duration="71 wk", strain="C57BL/6",
                             sex="male", control_diet="normal chow (grain-based assumed; NOT STATED)",
                             sra_read_pairs=float(r[r.gsm == g].read_count.iloc[0]) if (r.gsm == g).any() else np.nan,
                             shared_control=is_ctrl, role=STUDIES["HONG"]["role"]))
prov.append(dict(study_unit="HONG", geo_series="GSE261893;GSE261725;GSE286116;GSE262399", lab="Hong lab (no linked paper)",
                 polymer="PE (primary); PS, PP, PB (same lab, same controls)", species_strain="Mus musculus C57BL/6 male",
                 exposure="MP in food, 71 wk; particle size and dose NOT STATED", control="normal chow (3 mice shared by 4 series)",
                 n_biological_mice=15, n_gsm=24, n_control_mice=3,
                 input_type=STUDIES["HONG"]["input"], quantification="HISAT2/StringTie/Ballgown, GRCm39 RefSeq",
                 integer_counts=True, genes_raw=len(H), genes_vm35=len(Hc),
                 shared_control_evidence="CON1-3 count columns byte-identical across all 4 series; control SRA read counts identical (48.6/50.1/67.0 M)",
                 comparable_to_GSE294180="partly: PE vs control on a control diet, but chow (not purified), 71 wk in food (not 8 wk gavage), unknown size/dose, older mice",
                 role="PRIMARY independent PE test (one study)",
                 files=";".join(f"{a}_readcount_Con_{p}.txt.gz" for a, p in HS.items()),
                 sha256=";".join(sha256(DROOT / a / f"{a}_readcount_Con_{p}.txt.gz")[:16] for a, p in HS.items())))

# ---------------- GSE245069 ----------------
f = DROOT / "GSE245069" / "GSE245069_gene_expression.txt.gz"
g = pd.read_csv(f, sep="\t")
rc = [c for c in g.columns if c.startswith("read_count_Liver_")]
tp = [c for c in g.columns if c.startswith("tpm_Liver_")]
t = soft_titles("GSE245069")
title2gsm = {v.replace("0.5um", "05um"): k for k, v in t.items()}
ren = {c: "G245_" + c.replace("read_count_Liver_", "") for c in rc}
G = g.set_index("gene_symbol")[rc].rename(columns=ren)
GT = g.set_index("gene_symbol")[tp].rename(columns=lambda c: "G245_" + c.replace("tpm_Liver_", ""))
nonint = float((np.abs(G.values - np.round(G.values)) > 1e-9).any(axis=1).mean())
m = map_ids(G.index, "symbol")
Gc, mm = collapse(G, m, "sum")
GTc, _ = collapse(GT, m, "sum")
mapping_rows.append(mm.assign(study="GSE245069"))
Gc.to_csv(WORK / "GSE245069_expression.csv.gz")
GTc.to_csv(WORK / "GSE245069_tpm.csv.gz")
r = runs[runs.dataset == "GSE245069"]
for c in Gc.columns:
    lab = c.replace("G245_", "")
    gsm = title2gsm["Liver_" + lab]
    grp = lab.rsplit("_", 1)[0]
    manifest.append(dict(study="GSE245069", mouse_id=c, gsm=gsm, series="GSE245069", title=t[gsm],
                         group={"H2O": "CON", "80nm": "PS_80nm", "05um": "PS_0.5um", "5um": "PS_5um"}[grp],
                         polymer="none" if grp == "H2O" else "PS", size={"H2O": "-", "80nm": "80 nm", "05um": "0.5 um", "5um": "5 um"}[grp],
                         dose="1 mg/L drinking water", route="drinking water", duration="12 wk", strain="ICR (outbred)", sex="male",
                         control_diet="standard chow (NOT STATED in GEO)",
                         sra_read_pairs=float(r[r.gsm == gsm].read_count.iloc[0]) if (r.gsm == gsm).any() else np.nan,
                         shared_control=grp == "H2O", role=STUDIES["GSE245069"]["role"]))
prov.append(dict(study_unit="GSE245069", geo_series="GSE245069", lab="Li et al., Sci Rep 2025 (PMID 39833454)",
                 polymer="PS 80 nm / 0.5 um / 5 um", species_strain="Mus musculus ICR (outbred) male, 6 wk at start",
                 exposure="1 mg/L in drinking water, 12 wk", control="water (4 mice shared by 3 size contrasts)",
                 n_biological_mice=16, n_gsm=16, n_control_mice=4, input_type=STUDIES["GSE245069"]["input"],
                 quantification="deposited read_count / FPKM / TPM (GRCm38.p6, NCBI)", integer_counts=False,
                 genes_raw=len(G), genes_vm35=len(Gc), shared_control_evidence="design: one water group for 3 sizes",
                 comparable_to_GSE294180="no: different polymer (PS), outbred strain, water route, 12 wk; cross-polymer transfer only",
                 role="SECONDARY cross-polymer (one study, sizes correlated through shared control)",
                 files=f.name, sha256=sha256(f)[:16], noninteger_gene_fraction=round(nonint, 4),
                 aibd_artifacts="NOT FOUND on this machine (repo, Desktop, Downloads, Documents, D: searched); analysis done fresh here"))

# ---------------- GSE277772 ----------------
t = soft_titles("GSE277772")
parts = {}
for fp in sorted((DROOT / "GSE277772" / "raw_extract").glob("GSM*.txt.gz")):
    gsm = fp.name.split("_")[0]
    parts[gsm] = pd.read_csv(fp, sep="\t", index_col=0).iloc[:, 0]
E = pd.DataFrame(parts)
nonint77 = float((np.abs(E.values - np.round(E.values)) > 1e-9).mean())
m = map_ids(E.index, "ensembl")
Ec, mm = collapse(E, m, "sum")
mapping_rows.append(mm.assign(study="GSE277772"))
ids = {gsm: ("G277_CON_" if "normal" in t[gsm] else "G277_PS_") + re.search(r"rep(\d)", t[gsm]).group(1) for gsm in Ec.columns}
Ec = Ec.rename(columns=ids)
Ec.to_csv(WORK / "GSE277772_expression.csv.gz")
r = runs[runs.dataset == "GSE277772"]
for gsm, mid in ids.items():
    manifest.append(dict(study="GSE277772", mouse_id=mid, gsm=gsm, series="GSE277772", title=t[gsm],
                         group="CON" if "CON" in mid else "PS", polymer="none" if "CON" in mid else "PS", size="NOT STATED",
                         dose="0.4 mg twice weekly, oral", route="oral gavage", duration="8 wk", strain="C57BL/6", sex="male",
                         control_diet="NOT STATED", sra_read_pairs=float(r[r.gsm == gsm].read_count.iloc[0]),
                         shared_control=False, role=STUDIES["GSE277772"]["role"]))
prov.append(dict(study_unit="GSE277772", geo_series="GSE277772", lab="Environ Pollut 2025 (PMID 39956511)", polymer="PS (size NOT STATED)",
                 species_strain="Mus musculus C57BL/6 male, 6 wk at start", exposure="0.4 mg PS-MP oral, twice weekly, 8 wk",
                 control="water gavage (3 mice)", n_biological_mice=6, n_gsm=6, n_control_mice=3,
                 input_type=STUDIES["GSE277772"]["input"], quantification="HISAT2/featureCounts; deposited values normalised (no raw count file on GEO)",
                 integer_counts=False, genes_raw=len(E), genes_vm35=len(Ec), noninteger_gene_fraction=round(nonint77, 4),
                 shared_control_evidence="n/a", comparable_to_GSE294180="no (PS); counts unavailable",
                 role="EXPLORATORY rank-only (sample-level singscores; gene-level limma on log values descriptive only). NOT a validation",
                 files="GSE277772_RAW.tar (6 per-sample txt)", sha256=sha256(DROOT / "GSE277772" / "GSE277772_RAW.tar")[:16]))

M = pd.DataFrame(manifest)
M = M[~M.duplicated("mouse_id")].reset_index(drop=True)  # Hong controls: keep ONE row per mouse
gsm_all = pd.DataFrame(manifest)
ctrl_alias = gsm_all[(gsm_all.study == "HONG") & gsm_all.shared_control].groupby("mouse_id").gsm.apply(lambda s: ";".join(sorted(s)))
M["all_gsm_aliases"] = M.mouse_id.map(ctrl_alias).fillna(M.gsm)
M.to_csv(RES / "external_sample_manifest.csv", index=False)
gsm_all.to_csv(RES / "external_gsm_table.csv", index=False)
pd.concat(mapping_rows).to_csv(RES / "gene_mapping_external.csv.gz", index=False)

# ---------------- sample QC (technical covariates; no labels used) ----------------
mus = set(pd.read_csv(RS / "roadmap_phase2/metadata/muscle_exclusion_genes.csv").vm35_symbol)
qc = []
for st, X in (("HONG", Hc), ("GSE245069", Gc), ("GSE277772", Ec)):
    cpm = X / X.sum() * 1e6
    mg = [x for x in mus if x in X.index]
    for s in X.columns:
        qc.append(dict(study=st, mouse_id=s, counted_total=float(X[s].sum()), detected_genes=int((X[s] > 0).sum()),
                       muscle_cpm_sum=float(cpm.loc[mg, s].sum()), Alb_cpm=float(cpm.loc["Alb", s]) if "Alb" in X.index else np.nan,
                       Xist_cpm=float(cpm.loc["Xist", s]) if "Xist" in X.index else np.nan,
                       Ddx3y_cpm=float(cpm.loc["Ddx3y", s]) if "Ddx3y" in X.index else np.nan,
                       top100_share=float(np.sort(X[s].values)[::-1][:100].sum() / X[s].sum())))
Q = pd.DataFrame(qc).merge(M[["mouse_id", "group", "sra_read_pairs"]], on="mouse_id")
Q.to_csv(RES / "external_sample_qc.csv", index=False)
P = pd.DataFrame(prov)
P.to_csv(P5 / "dataset_provenance.csv", index=False)
with pd.option_context("display.width", 250, "display.max_columns", 30):
    print(P[["study_unit", "n_biological_mice", "n_gsm", "genes_raw", "genes_vm35", "integer_counts"]].to_string(index=False))
    print(Q.groupby(["study", "group"])[["counted_total", "sra_read_pairs", "detected_genes", "muscle_cpm_sum", "Xist_cpm", "Ddx3y_cpm"]].mean().round(1).to_string())
print("shared-control identity:", identical, "control read counts:", ctrl_reads["GSE261893"])
