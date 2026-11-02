"""Step 2: build metadata/sample_metadata.csv from the GEO SOFT family file.

Each row is one GSM record. Group labels come from the GEO ``treatment``
characteristic; the counts-column name is used only as a cross-check.
Study-level facts that GEO does not record per sample (strain, sex, age, dose,
diet catalogue numbers) live in metadata/study_design.csv with their source.
"""

import sys

from mpliver.config import GROUP_ORDER, SAMPLE_METADATA_FILE, SOFT_FILE
from mpliver.geo_metadata import build_sample_metadata

COLUMN_ORDER = [
    "geo_accession", "sample_title", "counts_column", "group", "diet", "microplastic",
    "replicate", "geo_treatment", "geo_batch", "library_name", "tissue", "genotype",
    "organism", "molecule", "library_strategy", "instrument", "biosample", "sra_experiment",
    "submission_date", "group_from_column_name",
]


def main() -> int:
    metadata = build_sample_metadata(SOFT_FILE)

    mismatched = metadata[metadata["group"] != metadata["group_from_column_name"]]
    if not mismatched.empty:
        print("[error] GEO treatment and counts-column name disagree:")
        print(mismatched[["geo_accession", "geo_treatment", "counts_column"]])
        return 1

    metadata["group"] = metadata["group"].astype("category").cat.set_categories(GROUP_ORDER)
    metadata = metadata.sort_values(["group", "replicate"])[COLUMN_ORDER]
    metadata.to_csv(SAMPLE_METADATA_FILE, index=False)

    print(f"Wrote {len(metadata)} samples to {SAMPLE_METADATA_FILE}")
    print(metadata.groupby("group", observed=True).size().rename("n").to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
