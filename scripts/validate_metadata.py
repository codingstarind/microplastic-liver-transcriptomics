"""Consistency checks on the committed metadata files (runs in CI, no download needed).

Checks that the sample table still describes the published design (19 samples,
4 / 5 / 5 / 5 per group, unique IDs, labels agreeing with column names) and that
every raw file has a provenance record with checksums.
"""

import json
import sys

import pandas as pd

from mpliver.config import GROUP_ORDER, PROVENANCE_FILE, SAMPLE_METADATA_FILE, SOURCE_FILES
from mpliver.geo_metadata import group_from_column_name

EXPECTED_GROUP_SIZES = {"CD+Veh": 4, "CD+PE": 5, "MD+Veh": 5, "MD+PE": 5}


def main() -> int:
    errors = []
    meta = pd.read_csv(SAMPLE_METADATA_FILE)

    if len(meta) != 19:
        errors.append(f"expected 19 samples, found {len(meta)}")
    for col in ["geo_accession", "counts_column", "library_name", "biosample"]:
        if meta[col].duplicated().any():
            errors.append(f"duplicate values in {col}")
    sizes = meta["group"].value_counts().to_dict()
    if sizes != EXPECTED_GROUP_SIZES:
        errors.append(f"group sizes {sizes} differ from {EXPECTED_GROUP_SIZES}")
    if set(meta["group"]) - set(GROUP_ORDER):
        errors.append("unknown group labels")
    inferred = meta["counts_column"].map(group_from_column_name)
    if (inferred != meta["group"]).any():
        errors.append("group label disagrees with counts column name")

    provenance = json.loads(PROVENANCE_FILE.read_text(encoding="utf-8"))
    for key in SOURCE_FILES:
        record = provenance.get("files", {}).get(key)
        if not record or not record.get("sha256") or not record.get("downloaded_utc"):
            errors.append(f"missing provenance/checksum for {key}")

    for error in errors:
        print(f"[error] {error}")
    if not errors:
        print("Metadata and provenance checks passed.")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
