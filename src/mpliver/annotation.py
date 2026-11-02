"""Gene-level annotation from a GENCODE GTF (gene_id -> symbol, biotype, chromosome)."""

import gzip
import re
from pathlib import Path

import pandas as pd

_ATTR_RE = re.compile(r'(\S+) "([^"]*)"')


def parse_gtf_attributes(field: str) -> dict:
    return dict(_ATTR_RE.findall(field))


def read_gene_annotation(gtf_path: Path) -> pd.DataFrame:
    """Return one row per ``gene`` feature in the GTF."""
    rows = []
    opener = gzip.open if gtf_path.suffix == ".gz" else open
    with opener(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parse_gtf_attributes(parts[8])
            rows.append({
                "gene_id": attrs["gene_id"],
                "gene_name": attrs.get("gene_name"),
                "gene_type": attrs.get("gene_type"),
                "chromosome": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
            })
    return pd.DataFrame(rows).set_index("gene_id")


def annotation_by_symbol(genes: pd.DataFrame) -> pd.DataFrame:
    """Collapse a gene-level annotation to one row per gene symbol.

    The deposited counts are indexed by symbol. A few symbols belong to more than
    one Ensembl gene, so their IDs, biotypes and chromosomes are joined with ';'
    and ``n_gene_ids`` records how many genes share the symbol.
    """
    def join_unique(values):
        return ";".join(sorted(set(map(str, values))))

    grouped = genes.reset_index().groupby("gene_name")
    out = pd.DataFrame({
        "gene_ids": grouped["gene_id"].agg(join_unique),
        "n_gene_ids": grouped["gene_id"].size(),
        "gene_type": grouped["gene_type"].agg(join_unique),
        "chromosome": grouped["chromosome"].agg(join_unique),
        "start": grouped["start"].min(),
        "end": grouped["end"].max(),
    })
    out.index.name = "gene_symbol"
    return out
