"""Build transcript-to-gene tables from the GENCODE M35 transcript FASTA headers.

Header format: >ENSMUST...|ENSMUSG...|OTTMUSG|OTTMUST|transcript_name|gene_name|length|biotype|
Salmon's --gencode option keeps only the first field, so the transcript ID is the key.

Two tables are written:
    tx2gene_symbol.tsv   transcript -> gene symbol (matches the deposited matrix)
    tx2gene_ensembl.tsv  transcript -> Ensembl gene ID (one row per real gene)

Usage: python pipelines/quantification/make_tx2gene.py <transcripts.fa.gz> <out_dir>
"""

import gzip
import sys
from pathlib import Path


def main() -> int:
    fasta, out_dir = Path(sys.argv[1]), Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    with gzip.open(fasta, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                fields = line[1:].rstrip("\n").split("|")
                rows.append((fields[0], fields[1], fields[5], fields[7]))
    with open(out_dir / "tx2gene_symbol.tsv", "w") as sym, \
            open(out_dir / "tx2gene_ensembl.tsv", "w") as ens:
        sym.write("TXNAME\tGENEID\n")
        ens.write("TXNAME\tGENEID\tGENE_SYMBOL\tTX_BIOTYPE\n")
        for tx, gene_id, symbol, biotype in rows:
            sym.write(f"{tx}\t{symbol}\n")
            ens.write(f"{tx}\t{gene_id}\t{symbol}\t{biotype}\n")
    print(f"{len(rows)} transcripts, {len({r[2] for r in rows})} symbols, "
          f"{len({r[1] for r in rows})} Ensembl genes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
