"""Gene-identifier mapping to the GSE294180 annotation (GENCODE vM35).

Strategy (every decision is recorded per original identifier):
  1. Original identifiers are preserved verbatim (column `original_id`).
  2. Malformed identifiers (embedded tab / quote / whitespace) are flagged
     `malformed`; the text before the first tab is used as a candidate symbol
     but the row is NOT used in analysis-ready matrices (its values may be
     misaligned).
  3. Ensembl gene IDs: version suffix stripped, matched to vM35 gene_ids
     (status `ensembl_match`). Unmatched Ensembl IDs -> `ensembl_not_in_vM35`
     (retired / pre-M35 genes).
  4. Symbols: exact match to vM35 gene_name (`symbol_exact`); otherwise a
     case-insensitive unique match (`symbol_case`); otherwise
     `symbol_not_in_vM35` (typically renamed genes or old Gm/Rik models).
     No synonym tables are used, so renamed genes are reported, not rescued.
  5. TempO-Seq probe IDs `Symbol_probeNumber`: symbol part mapped as in 4.
  6. Many-to-one: several original rows mapping to one vM35 gene are
     resolved by a stated rule per data type (see `collapse`).
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
ANNOT = REPO / "metadata" / "gene_annotation_gencode_vM35.csv.gz"


def load_vm35() -> tuple[pd.DataFrame, dict, dict, dict]:
    a = pd.read_csv(ANNOT)
    ens2sym = {}
    for sym, ids in zip(a["gene"], a["gene_ids"]):
        for gid in str(ids).split(";"):
            ens2sym[gid.split(".")[0].strip()] = sym
    sym_exact = {s: s for s in a["gene"]}
    low = a["gene"].str.lower()
    uniq_low = low[~low.duplicated(keep=False)]
    sym_case = dict(zip(uniq_low, a.loc[uniq_low.index, "gene"]))
    return a, ens2sym, sym_exact, sym_case


VM35, ENS2SYM, SYM_EXACT, SYM_CASE = load_vm35()


def load_mgi_synonyms() -> tuple[dict, str]:
    """Old/alternative symbol -> current vM35 symbol, from MGI MRK_List1.rpt.

    Used only for identifiers that do not match a vM35 symbol exactly or
    case-insensitively. A synonym is used only if it points to exactly one
    vM35 symbol. Withdrawn MGI symbols ('withdrawn, = X') map to X.
    """
    import os
    root = Path(os.environ.get("MPLIVER_P2_DATA", "D:/mpliver_work/phase2_refs"))
    f = root / "MGI" / "MRK_List1.rpt"
    if not f.exists():
        return {}, "MGI file absent: no synonym rescue"
    t = pd.read_csv(f, sep="\t", dtype=str, usecols=lambda c: c in (
        "Marker Symbol", "Status", "Marker Name", "Marker Synonyms (pipe-separated)"))
    cand: dict[str, set] = {}
    off = t[t["Status"] == "O"]
    for sym, syns in zip(off["Marker Symbol"], off["Marker Synonyms (pipe-separated)"].fillna("")):
        if sym not in SYM_EXACT:
            continue
        for a in syns.split("|"):
            a = a.strip()
            if a and a not in SYM_EXACT:
                cand.setdefault(a, set()).add(sym)
    wd = t[t["Status"] == "W"]
    for sym, name in zip(wd["Marker Symbol"], wd["Marker Name"].fillna("")):
        mm = re.match(r"withdrawn, = (\S+)$", name)
        if mm and mm.group(1) in SYM_EXACT and sym not in SYM_EXACT:
            cand.setdefault(sym, set()).add(mm.group(1))
    syn = {k: next(iter(v)) for k, v in cand.items() if len(v) == 1}
    return syn, f"MGI MRK_List1.rpt ({pd.Timestamp(f.stat().st_mtime, unit='s'):%Y-%m-%d}); {len(syn)} unambiguous synonyms"


MGI_SYN, MGI_NOTE = load_mgi_synonyms()
ENS_RE = re.compile(r"^ENSMUSG\d{11}(\.\d+)?$")


def map_ids(ids, id_type: str) -> pd.DataFrame:
    """id_type: 'ensembl' | 'symbol' | 'tempo_probe'."""
    rows = []
    for orig in ids:
        s = str(orig)
        status, cand, sym = "", s, None
        if re.fullmatch(r"[^\t\s\"]+\t1", s):
            # GSE316023 export artefact '<symbol>\t1'; the row's values are aligned
            status = "tab_suffix"
            cand = s.split("\t")[0]
        elif "\t" in s or '"' in s or s != s.strip():
            status = "malformed"
            cand = s.split("\t")[0].strip().strip('"')
        if id_type == "tempo_probe":
            cand = cand.rsplit("_", 1)[0]
        if id_type == "ensembl":
            base = cand.split(".")[0]
            if not ENS_RE.match(cand):
                status = status or "invalid_id"
            elif base in ENS2SYM:
                sym, status = ENS2SYM[base], status or "ensembl_match"
            else:
                status = status or "ensembl_not_in_vM35"
        else:
            if cand in SYM_EXACT:
                sym = cand
                status = "tab_suffix_recovered" if status == "tab_suffix" else (status or "symbol_exact")
            elif cand.lower() in SYM_CASE:
                sym = SYM_CASE[cand.lower()]
                status = "tab_suffix_recovered" if status == "tab_suffix" else (status or "symbol_case")
            elif cand in MGI_SYN:
                sym = MGI_SYN[cand]
                status = "tab_suffix_recovered" if status == "tab_suffix" else (status or "symbol_mgi_synonym")
            elif cand in ("", "nan", "None"):
                status = status or "invalid_id"
            else:
                status = "tab_suffix_unmapped" if status == "tab_suffix" else (status or "symbol_not_in_vM35")
        rows.append(dict(original_id=s, candidate=cand, vm35_symbol=sym, status=status))
    m = pd.DataFrame(rows)
    m["usable"] = m["vm35_symbol"].notna() & (m["status"] != "malformed")
    # a recovered tab-suffix row whose symbol also exists as a clean row is not used
    clean = set(m.loc[m["usable"] & (m["status"] != "tab_suffix_recovered"), "vm35_symbol"])
    conflict = (m["status"] == "tab_suffix_recovered") & m["vm35_symbol"].isin(clean)
    m.loc[conflict, "usable"] = False
    m.loc[conflict, "status"] = "tab_suffix_conflict_with_clean_row"
    # exact duplicate original identifiers: flag (resolved by collapse rule)
    m["original_id_duplicated"] = m["original_id"].duplicated(keep=False)
    return m


def collapse(mat: pd.DataFrame, mapping: pd.DataFrame, rule: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collapse usable rows onto vM35 symbols.

    rule 'sum'          : sum rows (counts from Ensembl IDs merged into one symbol)
    rule 'max_mean'     : keep the row with the highest mean (probes, RPKM)
    Returns (collapsed matrix, mapping with `kept_in_ready` and `collapse_note`).
    """
    m = mapping.set_index("original_id").copy()
    m["kept_in_ready"] = False
    m["collapse_note"] = ""
    use = m.index[m["usable"]]
    sub = mat.loc[use]
    sym = m.loc[use, "vm35_symbol"]
    dup = sym[sym.duplicated(keep=False)]
    m.loc[dup.index, "collapse_note"] = "many_to_one"
    if rule == "sum":
        out = sub.groupby(sym.values).sum()
        m.loc[use, "kept_in_ready"] = True
        m.loc[dup.index, "collapse_note"] = "many_to_one_summed"
    elif rule == "max_mean":
        means = sub.mean(axis=1)
        order = means.sort_values(ascending=False).index
        first = ~sym.loc[order].duplicated(keep="first")
        keep = order[first.values]
        out = sub.loc[keep]
        out.index = sym.loc[keep].values
        m.loc[keep, "kept_in_ready"] = True
        dropped = dup.index.difference(keep)
        m.loc[dropped, "collapse_note"] = "many_to_one_dropped_lower_mean"
        m.loc[dup.index.intersection(keep), "collapse_note"] = "many_to_one_kept_highest_mean"
    else:
        raise ValueError(rule)
    out = out.sort_index()
    out.index.name = "vm35_symbol"
    return out, m.reset_index()


def mapping_summary(name: str, mapping: pd.DataFrame, ready: pd.DataFrame, own_all: set, own_filt: set,
                    detected: pd.Index | None = None) -> dict:
    s = dict(
        dataset=name,
        total_rows=len(mapping),
        malformed=int(mapping["status"].isin(["malformed"]).sum()),
        tab_suffix_ids=int(mapping["status"].str.startswith("tab_suffix").sum()),
        duplicated_original_ids=int(mapping["original_id_duplicated"].sum()),
        invalid_id=int((mapping["status"] == "invalid_id").sum()),
        valid_ids=int((~mapping["status"].isin(["malformed", "invalid_id"])).sum()),
        mapped_to_vM35=int(mapping["usable"].sum()),
        unmapped=int((~mapping["usable"] & ~mapping["status"].isin(["malformed", "invalid_id"])).sum()),
        unique_vM35_genes=int(ready.shape[0]),
        rows_lost_to_collapse=int(mapping["usable"].sum() - ready.shape[0]),
        shared_with_GSE294180_all=int(len(set(ready.index) & own_all)),
        shared_with_GSE294180_filtered=int(len(set(ready.index) & own_filt)),
    )
    if detected is not None:
        s["detected_in_contrast_samples"] = int(len(detected))
        s["detected_and_in_GSE294180_filtered"] = int(len(set(detected) & own_filt))
    for k, v in mapping["status"].value_counts().items():
        s[f"status_{k}"] = int(v)
    return s


def own_universe() -> tuple[set, set]:
    raw = pd.read_csv(REPO / "data" / "raw" / "GSE294180_MP_raw_counts_gene_level.csv.gz", usecols=[0])
    allg = set(raw.iloc[:, 0].astype(str))
    filt = pd.read_csv(REPO / "data" / "processed" / "gene_filter.csv")
    fset = set(filt.loc[filt["passes_filter"].astype(bool), "gene"].astype(str))
    return allg, fset


def detected_genes(ready: pd.DataFrame, samples: list[str], kind: str) -> pd.Index:
    """Pre-specified detection rule used only for reporting coverage:
    counts: >= 10 counts in >= min_group_size samples (here: half the samples, min 3);
    RPKM/FPKM: >= 1 in >= half the samples. Not a final analysis filter."""
    x = ready[samples]
    k = max(3, len(samples) // 2)
    thr = 10 if kind == "counts" else 1.0
    return ready.index[(x >= thr).sum(axis=1) >= k]


__all__ = ["map_ids", "collapse", "mapping_summary", "own_universe", "detected_genes", "np"]
