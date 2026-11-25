"""Step 03: build per-dataset sample metadata, original and analysis-ready
expression matrices, and gene-mapping tables for the selected references.

Nothing is pooled across studies. Each dataset is written to its own folder:
  external : DATA_ROOT/<GSE>/prepared/expression_original.csv.gz
             DATA_ROOT/<GSE>/prepared/expression_ready.csv.gz
  repo     : datasets/<GSE>/sample_metadata.csv
             datasets/<GSE>/gene_mapping.csv.gz
             datasets/<GSE>/prepared_manifest.json
  repo     : metadata/sample_metadata_master.csv, metadata/gene_mapping_summary.csv

No differential expression or reference-signature estimation is done here.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from gene_mapping import collapse, detected_genes, map_ids, mapping_summary, own_universe
from p2_common import DATA_ROOT, REPO_P2, file_hashes

RAW = REPO_P2 / "metadata" / "geo_raw"
OWN_ALL, OWN_FILT = own_universe()

BASE_COLS = ["dataset", "gsm", "sample_id", "title", "tissue", "species", "strain", "sex", "age",
             "genotype", "treatment", "dose", "route", "exposure_duration", "housing_temp_C", "zt",
             "feeding_state", "fasting_before_sacrifice", "group", "biological_unit",
             "contrast_role", "notes"]


def soft(acc: str) -> pd.DataFrame:
    return pd.read_csv(RAW / f"{acc}_samples_raw.csv")


def finish_meta(acc: str, m: pd.DataFrame) -> pd.DataFrame:
    m.insert(0, "dataset", acc)
    for c in BASE_COLS:
        if c not in m:
            m[c] = ""
    return m[BASE_COLS + [c for c in m.columns if c not in BASE_COLS]]


def write(acc, original, ready, meta, mapping, kind, contrast_samples, note):
    pdir = DATA_ROOT / acc / "prepared"
    pdir.mkdir(parents=True, exist_ok=True)
    rdir = REPO_P2 / "datasets" / acc
    rdir.mkdir(parents=True, exist_ok=True)
    # matrix columns must equal metadata sample_id order
    assert list(ready.columns) == list(meta.loc[meta["in_ready_matrix"], "sample_id"]), acc
    original.to_csv(pdir / "expression_original.csv.gz")
    ready.to_csv(pdir / "expression_ready.csv.gz")
    meta.to_csv(rdir / "sample_metadata.csv", index=False)
    mapping.to_csv(rdir / "gene_mapping.csv.gz", index=False)
    det = detected_genes(ready, contrast_samples, "counts" if kind.endswith("counts") else "rpkm")
    summ = mapping_summary(acc, mapping, ready, OWN_ALL, OWN_FILT, det)
    files = {}
    for f in ("expression_original.csv.gz", "expression_ready.csv.gz"):
        md5, sha = file_hashes(pdir / f)
        files[f] = dict(path=str(pdir / f), md5=md5, sha256=sha)
    man = dict(dataset=acc, value_type=kind, note=note, original_shape=list(original.shape),
               ready_shape=list(ready.shape), ready_index="GENCODE vM35 gene symbol (GSE294180 annotation)",
               files=files, mapping_summary=summ)
    (rdir / "prepared_manifest.json").write_text(json.dumps(man, indent=2))
    print(f"{acc}: original {original.shape}, ready {ready.shape}, detected {len(det)}")
    return summ


# ---------------------------------------------------------------- GSE316023
def gse316023():
    acc = "GSE316023"
    d = soft(acc)
    m = pd.DataFrame({
        "gsm": d["gsm"], "sample_id": d["title"], "title": d["title"], "tissue": "liver",
        "species": "Mus musculus", "strain": "C57BL/6J (JAX)", "sex": "male",
        "age": "7-8 wk at diet start", "genotype": "WT",
        "treatment": np.where(d["ch_treatment"].str.contains("GAN"), "GAN_D09100310", "chow_LabDiet5001"),
        "route": "ad libitum diet",
        "exposure_duration": d["ch_timepoint_duration_on_diet"].str.replace("wk", "") + " wk on diet",
        "housing_temp_C": np.where(d["ch_housing_temperature"] == "TA", "21", "29"),
        "zt": "NOT STATED (lights on 06:00)", "feeding_state": "ad libitum",
        "fasting_before_sacrifice": "5 h",
        "biological_unit": "individual mouse",
    })
    m["group"] = (np.where(m["treatment"].str.startswith("GAN"), "GAN", "Chow") + "_" +
                  d["ch_timepoint_duration_on_diet"] + "_" + d["ch_housing_temperature"])
    m["contrast_role"] = m["group"].map({
        "GAN_wk23_TA": "MASH_primary:case", "Chow_wk23_TA": "MASH_primary:control",
        "GAN_wk23_TN": "MASH_sens_TN:case", "Chow_wk23_TN": "MASH_sens_TN:control"}).fillna("not_in_contrast")
    m["notes"] = "Chow controls exist only at wk23; A=TA(21C), B=TN(29C) per GSM characteristics"
    m = finish_meta(acc, m)
    c = pd.read_csv(DATA_ROOT / acc / "GSE316023_COUNT_NCOMMS-25-40823-T.csv.gz", index_col=0)
    assert set(c.columns) == set(m["sample_id"]), "GSE316023 column/sample mismatch"
    c = c[m["sample_id"]]
    mp = map_ids(c.index, "symbol")
    # duplicated original IDs: give each row a unique key so mapping stays one row per matrix row
    key = c.index.astype(str) + np.where(c.index.duplicated(keep=False),
                                          "#row" + pd.Series(range(len(c))).astype(str).values, "")
    c.index, mp["original_id"] = key, key
    ready_float, mapping = collapse(c, mp, "max_mean")
    ready = np.rint(ready_float).astype(int)  # RSEM-style expected counts -> integer (see report)
    m["in_ready_matrix"] = True
    m["library_size"] = m["sample_id"].map(c.sum().round().astype(int))
    note = ("Deposited non-integer estimated counts (6% of cells non-integer; FPKM file implies a gene length "
            "that is ~constant across samples, median CV 1%). Ready matrix: tab-suffix IDs recovered, "
            "duplicate symbols resolved by highest mean, values rounded to integers for count models.")
    contrast = list(m.loc[m.contrast_role.str.startswith("MASH_primary"), "sample_id"])
    return write(acc, c, ready, m, mapping, "estimated_counts", contrast, note), m


# ---------------------------------------------------------------- GSE246328
def gse246328():
    acc = "GSE246328"
    d = soft(acc)
    lib = d["title"].str.extract(r"\[(\d+)\]")[0]
    trt = d["ch_treatment"]
    total = {"Vehicle 8 w": "42 wk GAN (34 pre + 8)", "Vehicle 16 w": "50 wk GAN (34 + 16)",
             "Vehicle 24 w": "58 wk GAN (34 + 24)", "Chow-reversal 8 w": "34 wk GAN then 8 wk chow",
             "Chow-reversal 16 w": "34 wk GAN then 16 wk chow", "Chow-reversal 24 w": "34 wk GAN then 24 wk chow",
             "Chow Vehicle 24 w": "58 wk chow"}
    m = pd.DataFrame({
        "gsm": d["gsm"], "sample_id": lib, "title": d["title"], "tissue": "liver", "species": "Mus musculus",
        "strain": "C57BL/6J (Janvier; C57BL/6JRj per GSE294629)", "sex": "male (per GSE294629; not stated in GSE246328)",
        "age": "5-6 wk at diet start", "genotype": "WT", "treatment": trt,
        "route": "ad libitum diet + s.c. vehicle dosing", "exposure_duration": trt.map(total),
        "housing_temp_C": "21", "zt": "NOT STATED (lights on 03:00)", "feeding_state": "ad libitum",
        "fasting_before_sacrifice": "4 h (terminal blood)", "biological_unit": "individual mouse",
        "group": trt.str.replace(" ", "_"),
    })
    m["contrast_role"] = trt.map({"Vehicle 24 w": "MASH_secondary:case",
                                  "Chow Vehicle 24 w": "MASH_secondary:control"}).fillna("not_in_contrast")
    m["notes"] = ("GAN groups had liver biopsy before intervention; chow group biopsy status NOT STATED; "
                  "same libraries re-deposited in GSE294629")
    m = finish_meta(acc, m)
    c = pd.read_csv(DATA_ROOT / acc / "GSE246328_GAN-DIO-NASH-Senescence_counts.txt.gz", sep="\t", index_col=0)
    c.columns = c.columns.astype(str)
    assert set(c.columns) == set(m["sample_id"]), "GSE246328 column/sample mismatch"
    c = c[m["sample_id"]]
    mp = map_ids(c.index, "ensembl")
    ready, mapping = collapse(c, mp, "sum")
    m["in_ready_matrix"] = True
    m["library_size"] = m["sample_id"].map(c.sum())
    note = "Deposited integer STAR gene counts (Ensembl GRCm38 v96); Ensembl->vM35 symbol, many-to-one summed."
    contrast = list(m.loc[m.contrast_role.str.startswith("MASH_secondary"), "sample_id"])
    return write(acc, c, ready.astype(int), m, mapping, "raw_counts", contrast, note), m


# ---------------------------------------------------------------- GSE107787
def gse107787():
    acc = "GSE107787"
    d = soft(acc)
    t = d["title"].str.extract(r"(Liver|Muscle) (FAST|FED) ZT(\d+) Replicate (\d+)")
    sid = t[0] + "_" + t[1] + "_ZT" + t[2] + "_Rep" + t[3]
    m = pd.DataFrame({
        "gsm": d["gsm"], "sample_id": sid, "title": d["title"],
        "tissue": t[0].str.lower().replace({"muscle": "gastrocnemius"}), "species": "Mus musculus",
        "strain": "C57BL/6J (JAX)", "sex": "male", "age": "8 wk", "genotype": "WT",
        "treatment": t[1].map({"FAST": "24h_fast", "FED": "ad_libitum_chow"}),
        "route": "", "exposure_duration": t[1].map({"FAST": "24 h fast ending at harvest", "FED": ""}),
        "housing_temp_C": "NOT STATED", "zt": "ZT" + t[2].astype(int).astype(str),
        "feeding_state": t[1].str.lower(), "fasting_before_sacrifice": t[1].map({"FAST": "24 h", "FED": "none"}),
        "biological_unit": "individual mouse (single-housed 1 wk before harvest)",
        "group": t[1] + "_ZT" + t[2],
    })
    liver = m["tissue"] == "liver"
    m["contrast_role"] = np.where(liver, "FASTING_primary:" + np.where(t[1] == "FAST", "case", "control"),
                                  "not_in_contrast(muscle)")
    m.loc[liver & (t[1] == "FED"), "contrast_role"] += ";CIRCADIAN_sens_male:timecourse"
    m["notes"] = "Fasted mice moved to new cages without food 24 h before sacrifice; FAST and FED sampled at the same ZTs"
    m = finish_meta(acc, m)
    x = pd.read_excel(DATA_ROOT / acc / "GSE107787_RNAseq_expression_results.xlsx").set_index("Gene")
    assert set(x.columns) == set(m["sample_id"]), "GSE107787 column/sample mismatch"
    x = x[m["sample_id"]]
    mp = map_ids(x.index, "symbol")
    key = x.index.astype(str) + np.where(x.index.duplicated(keep=False),
                                          "#row" + pd.Series(range(len(x))).astype(str).values, "")
    x.index, mp["original_id"] = key, key
    m["in_ready_matrix"] = m["tissue"] == "liver"
    ready, mapping = collapse(x[m.loc[m.in_ready_matrix, "sample_id"]], mp, "max_mean")
    m["library_size"] = "not available (RPKM only)"
    note = "Deposited RPKM only (ELAND/RefSeq in-house pipeline). Ready = liver samples, vM35 symbols, duplicates by highest mean. No counts."
    contrast = list(m.loc[m.in_ready_matrix, "sample_id"])
    return write(acc, x, ready, m, mapping, "RPKM", contrast, note), m


# ---------------------------------------------------------------- GSE117134
def gse117134():
    acc = "GSE117134"
    d = soft(acc)
    geno = d["title"].str.extract(r"^(\S+) ZT(\d+) Replicate (\d+)")
    sid = geno[0] + "_ZT" + geno[1] + "_BR" + geno[2]
    m = pd.DataFrame({
        "gsm": d["gsm"], "sample_id": sid, "title": d["title"], "tissue": "liver", "species": "Mus musculus",
        "strain": "C57BL/6J background", "sex": "female", "age": "8-12 wk",
        "genotype": d["ch_genotype_variation"], "treatment": "ad_libitum_chow", "route": "",
        "exposure_duration": "", "housing_temp_C": "NOT STATED", "zt": "ZT" + geno[1].astype(int).astype(str),
        "feeding_state": "ad libitum", "fasting_before_sacrifice": "none", "biological_unit": "individual mouse",
        "group": geno[0] + "_ZT" + geno[1],
    })
    m["contrast_role"] = np.where(geno[0] == "WT", "CIRCADIAN_primary:timecourse", "not_in_contrast(clock mutant)")
    m["notes"] = "Female only; WT used; KO and Liver-RE excluded from reference"
    m = finish_meta(acc, m)
    x = pd.read_csv(DATA_ROOT / acc / "GSE117134_casava_gene_expression.tsv.gz", sep="\t", index_col=0)
    assert set(x.columns) == set(m["sample_id"]), "GSE117134 column/sample mismatch"
    x = x[m["sample_id"]]
    mp = map_ids(x.index, "symbol")
    key = x.index.astype(str) + np.where(x.index.duplicated(keep=False),
                                          "#row" + pd.Series(range(len(x))).astype(str).values, "")
    x.index, mp["original_id"] = key, key
    m["in_ready_matrix"] = m["genotype"] == "WT"
    ready, mapping = collapse(x[m.loc[m.in_ready_matrix, "sample_id"]], mp, "max_mean")
    m["library_size"] = "not available (RPKM only)"
    note = "Deposited RPKM only (same in-house pipeline as GSE107787). Ready = WT samples only."
    contrast = list(m.loc[m.in_ready_matrix, "sample_id"])
    return write(acc, x, ready, m, mapping, "RPKM", contrast, note), m


# ---------------------------------------------------------------- GSE282998
def gse282998():
    acc = "GSE282998"
    d = soft(acc)
    d = d[d["ch_treatment"].str.contains("GW7647")].copy()
    dose = d["ch_treatment"].str.extract(r"^(\d+)")[0].astype(int)
    strain = d["ch_strain"]
    geno = np.where(strain == "PPARa KO", "Ppara-null", "WT")
    m = pd.DataFrame({
        "gsm": d["gsm"], "sample_id": d["title"], "title": d["title"], "tissue": "liver (FFPE)",
        "species": "Mus musculus",
        "strain": strain.map({"PPARa KO": "B6;129S4-Ppara<tm1Gonz>/J", "B6129SF2/J": "B6129SF2/J",
                              "J:ARC(S)": "J:ARC(S) Swiss outbred"}),
        "sex": "male", "age": "9-12 wk", "genotype": geno, "treatment": "GW7647",
        "dose": dose.astype(str) + " mg/kg/day", "route": "oral gavage (vehicle 1% CMC, per paper)",
        "exposure_duration": "5 days; necropsy 24 h after last dose", "housing_temp_C": "20-26",
        "zt": "NOT STATED", "feeding_state": "ad libitum (not stated otherwise)", "fasting_before_sacrifice": "NOT STATED",
        "biological_unit": "individual mouse",
        "group": strain.str.replace(r"[^A-Za-z0-9]", "", regex=True) + "_GW" + dose.astype(str),
    })
    m["contrast_role"] = "not_in_contrast"
    wt = m["strain"] == "B6129SF2/J"
    m.loc[wt & (dose > 0), "contrast_role"] = "PPARA_primary:case"
    m.loc[wt & (dose == 0), "contrast_role"] = "PPARA_primary:control"
    ko = m["genotype"] == "Ppara-null"
    m.loc[ko & (dose > 0), "contrast_role"] = "PPARA_specificity_KO:case"
    m.loc[ko & (dose == 0), "contrast_role"] = "PPARA_specificity_KO:control"
    arc = m["strain"].str.startswith("J:ARC")
    m.loc[arc & (dose > 0), "contrast_role"] = "PPARA_sens_outbred:case"
    m.loc[arc & (dose == 0), "contrast_role"] = "PPARA_sens_outbred:control"
    m["notes"] = np.where(ko, "GEO genotype field says 'WT' for this KO-strain sample; strain field 'PPARa KO' used", "")
    m = finish_meta(acc, m)
    x = pd.read_csv(DATA_ROOT / acc / "GSE282998_GW7647_gene_counts.csv.gz", index_col=0)
    assert set(x.columns) == set(m["sample_id"]), "GSE282998 column/sample mismatch"
    x = x[m["sample_id"]]
    mp = map_ids(x.index, "tempo_probe")
    ready, mapping = collapse(x, mp, "max_mean")
    m["in_ready_matrix"] = True
    m["library_size"] = m["sample_id"].map(x.sum())
    note = ("TempO-Seq whole-transcriptome probe counts (30,146 probes, ~21.4k genes), FFPE. Ready = one probe per gene "
            "(highest mean across the arm, as in Heintz et al. 2025). Integer counts.")
    contrast = list(m.loc[m.contrast_role.str.startswith("PPARA_primary"), "sample_id"])
    return write(acc, x, ready.astype(int), m, mapping, "targeted_probe_counts", contrast, note), m


# ---------------------------------------------------------------- GSE225695
def gse225695():
    acc = "GSE225695"
    f = DATA_ROOT / acc / "GSE225695_mRNA_Expression_Profiling.xlsx"
    if not f.exists():
        print("GSE225695 xlsx not downloaded; skipped")
        return None, None
    d = soft(acc)
    grp = d["title"].str.extract(r"Liver, (Fed|Fast|Refed-2h)-(\d)")
    m = pd.DataFrame({
        "gsm": d["gsm"], "title": d["title"], "tissue": "liver", "species": "Mus musculus",
        "strain": "C57BL/6J (SLAC)", "sex": "male", "age": "10 wk (8 wk + 2 wk acclimation)", "genotype": "WT",
        "treatment": grp[0], "zt": "NOT STATED", "housing_temp_C": "21",
        "feeding_state": grp[0].str.lower(), "fasting_before_sacrifice": grp[0].map(
            {"Fed": "none", "Fast": "24 h", "Refed-2h": "24 h then 2 h refeed"}),
        "biological_unit": "individual mouse", "group": grp[0],
    })
    m["sample_id"] = grp[0] + "-" + grp[1]
    m["contrast_role"] = grp[0].map({"Fast": "FASTING_secondary:case", "Fed": "FASTING_secondary:control"}).fillna(
        "not_in_contrast(refed)")
    m["notes"] = "Sacrifice clock time per group NOT STATED; n=3"
    m = finish_meta(acc, m)
    raw = pd.read_excel(f, header=25)
    cnt_cols = [c for c in raw.columns if str(c).endswith("_count")]
    x = raw.set_index("gene_id")[cnt_cols]
    x.columns = [c[: -len("_count")] for c in cnt_cols]
    assert set(x.columns) == set(m["sample_id"]), "GSE225695 column/sample mismatch"
    x = x[m["sample_id"]]
    x.index = x.index.astype(str)
    mp = map_ids(x.index, "ensembl")
    key = x.index + np.where(x.index.duplicated(keep=False),
                             "#row" + pd.Series(range(len(x))).astype(str).values, "")
    x.index, mp["original_id"] = key, key
    nonint = int((x.values % 1 != 0).sum())
    ready, mapping = collapse(x, mp, "sum")
    m["in_ready_matrix"] = True
    m["library_size"] = m["sample_id"].map(x.sum().round().astype(int))
    note = (f"Deposited HISAT2/StringTie 'count' columns from xlsx (non-integer cells: {nonint}); "
            "StringTie-derived counts are estimates. Ensembl->vM35, many-to-one summed. n=3/group.")
    contrast = list(m.loc[m.contrast_role.str.startswith("FASTING_secondary"), "sample_id"])
    if nonint:
        ready = np.rint(ready)
    return write(acc, x, ready.astype(int), m, mapping, "estimated_counts", contrast, note), m


if __name__ == "__main__":
    summaries, metas = [], []
    for fn in (gse316023, gse246328, gse107787, gse117134, gse282998, gse225695):
        s, mm = fn()
        summaries.append(s)
        metas.append(mm)
    pd.DataFrame(summaries).to_csv(REPO_P2 / "metadata" / "gene_mapping_summary.csv", index=False)
    master = pd.concat(metas, ignore_index=True)
    master.to_csv(REPO_P2 / "metadata" / "sample_metadata_master.csv", index=False)
    print("master metadata rows:", len(master))
