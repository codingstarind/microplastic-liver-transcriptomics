"""Step 3: audit the deposited count matrix without changing it.

Outputs (all in results/tables/ unless noted):
    audit_summary.json           integrity checks and value-type evidence
    counts_preview.csv           first rows of the file exactly as deposited
    per_sample_audit.csv         per-sample totals, detection, mito/biotype fractions
    gene_biotype_summary.csv     genes and counts by GENCODE biotype
    top_expressed_genes.csv      genes taking the largest share of reads
    metadata/gene_annotation_gencode_vM35.csv.gz   annotation of genes in the matrix
"""

import json
import sys

import pandas as pd

from mpliver.annotation import annotation_by_symbol, read_gene_annotation
from mpliver.config import (
    COUNTS_FILE,
    GENE_ANNOTATION_FILE,
    GTF_FILE,
    SAMPLE_METADATA_FILE,
    TABLES_DIR,
)
from mpliver.counts import (
    audit_counts,
    check_sample_ids,
    per_sample_summary,
    read_counts_file,
    split_annotation_and_counts,
    value_type_evidence,
)

SEX_MARKERS = ["Xist", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty"]


def main() -> int:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    raw = read_counts_file(COUNTS_FILE)
    metadata = pd.read_csv(SAMPLE_METADATA_FILE)

    raw.head(10).to_csv(TABLES_DIR / "counts_preview.csv", index=False)
    annotation_cols, counts = split_annotation_and_counts(raw)

    summary = {
        "file": COUNTS_FILE.name,
        "raw_shape_rows_x_columns": list(raw.shape),
        "raw_columns": [str(c) for c in raw.columns],
        "raw_column_dtypes": {str(c): str(t) for c, t in raw.dtypes.items()},
        "annotation_columns": [annotation_cols.index.name, *annotation_cols.columns],
        "first_gene_ids": counts.index[:5].tolist(),
        "integrity": audit_counts(counts),
        "value_type_evidence": value_type_evidence(counts),
        "sample_id_check": check_sample_ids(counts.columns, metadata["counts_column"]),
    }

    # The deposited file is indexed by gene symbol, so match on GENCODE gene_name.
    gencode_genes = read_gene_annotation(GTF_FILE)
    by_symbol = annotation_by_symbol(gencode_genes)
    matched = counts.index.isin(by_symbol.index)
    absent = by_symbol.index.difference(counts.index)
    summary["gencode_vM35_match"] = {
        "n_gencode_genes": int(len(gencode_genes)),
        "n_gencode_unique_symbols": int(len(by_symbol)),
        "n_counts_ids_matching_gencode_symbol": int(matched.sum()),
        "n_counts_ids_unmatched": int((~matched).sum()),
        "unmatched_examples": counts.index[~matched][:10].tolist(),
        "n_counts_symbols_shared_by_multiple_gene_ids": int(
            (by_symbol.reindex(counts.index)["n_gene_ids"] > 1).sum()),
        "n_gencode_symbols_absent_from_counts": int(len(absent)),
        "absent_symbols_by_type": by_symbol.loc[absent, "gene_type"].value_counts()
                                           .head(10).to_dict(),
        "absent_symbols_by_chromosome": by_symbol.loc[absent, "chromosome"].value_counts()
                                                 .head(10).to_dict(),
    }
    ann = by_symbol.reindex(counts.index)
    ann.insert(0, "gene_name", counts.index)
    ann.index.name = "gene"
    ann.to_csv(GENE_ANNOTATION_FILE)

    # Per-sample audit on the original values.
    per_sample = per_sample_summary(counts)
    is_mito = ann["chromosome"].eq("chrM")
    per_sample["fraction_counts_chrM"] = counts[is_mito.to_numpy()].sum() / counts.sum()
    per_sample["top_gene_symbol"] = ann["gene_name"].reindex(per_sample["top_gene"]).to_numpy()
    per_sample["fraction_counts_top_gene"] = per_sample["max_gene_count"] / counts.sum()
    per_sample = per_sample.join(
        metadata.set_index("counts_column")[["geo_accession", "group", "library_name"]])
    per_sample.index.name = "counts_column"
    per_sample.to_csv(TABLES_DIR / "per_sample_audit.csv")

    biotype = pd.DataFrame({
        "n_genes": ann["gene_type"].fillna("unannotated").value_counts(),
        "fraction_of_all_counts": counts.sum(axis=1).groupby(
            ann["gene_type"].fillna("unannotated")).sum() / counts.to_numpy().sum(),
    }).sort_values("n_genes", ascending=False)
    biotype.index.name = "gene_type"
    biotype.to_csv(TABLES_DIR / "gene_biotype_summary.csv")

    share = counts.div(counts.sum(), axis=1)
    top = pd.DataFrame({
        "gene_name": ann["gene_name"],
        "gene_type": ann["gene_type"],
        "mean_fraction_of_library": share.mean(axis=1),
        "min_fraction": share.min(axis=1),
        "max_fraction": share.max(axis=1),
    }).sort_values("mean_fraction_of_library", ascending=False).head(30)
    top.to_csv(TABLES_DIR / "top_expressed_genes.csv")

    # Sex check: the paper reports male mice, GEO does not record sex.
    sex_genes = [g for g in SEX_MARKERS if g in counts.index]
    sex_table = counts.loc[sex_genes].T
    sex_table["group"] = per_sample["group"]
    sex_table.index.name = "counts_column"
    sex_table.to_csv(TABLES_DIR / "sex_marker_counts.csv")
    y_genes = [g for g in sex_genes if g != "Xist"]
    summary["sex_check"] = {
        "markers_present": sex_genes,
        "xist_count_range": [float(counts.loc["Xist"].min()), float(counts.loc["Xist"].max())]
        if "Xist" in counts.index else None,
        "y_marker_total_range": [float(counts.loc[y_genes].sum().min()),
                                 float(counts.loc[y_genes].sum().max())],
        "chrY_symbols_in_counts": int(ann["chromosome"].eq("chrY").sum()),
        "chrY_total_count_range": [float(counts[ann["chromosome"].eq("chrY").to_numpy()]
                                         .sum().min()),
                                   float(counts[ann["chromosome"].eq("chrY").to_numpy()]
                                         .sum().max())],
    }

    # Identical count profiles: mostly genes with a handful of reads, or genuinely
    # ambiguous loci (e.g. identical gene copies) - check which.
    expressed = counts[counts.sum(axis=1) > 0]
    dup_mask = expressed.duplicated(keep=False)
    mean_counts = expressed.mean(axis=1)
    summary["duplicate_profiles"] = {
        "n_rows_in_duplicate_groups": int(dup_mask.sum()),
        "n_with_mean_count_ge_10": int((dup_mask & (mean_counts >= 10)).sum()),
        "examples_mean_ge_10": expressed[dup_mask & (mean_counts >= 10)].index[:20].tolist(),
    }

    (TABLES_DIR / "audit_summary.json").write_text(
        json.dumps(summary, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, default=str))
    return 0 if summary["sample_id_check"]["all_match"] else 1


if __name__ == "__main__":
    sys.exit(main())
