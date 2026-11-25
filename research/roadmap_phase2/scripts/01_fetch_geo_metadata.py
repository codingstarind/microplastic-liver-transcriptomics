"""Step 01: fetch GEO family SOFT files and supplementary-file listings
for every Phase 2 candidate series. Small files only; no expression
matrices and no FASTQ. Rerunnable: existing files are not re-downloaded.

Outputs
  DATA_ROOT/<GSE>/<GSE>_family.soft.gz            (external, checksummed)
  metadata/geo_raw/<GSE>_samples_raw.csv          (all SOFT sample fields)
  metadata/geo_raw/<GSE>_series.txt               (series-level fields)
  metadata/geo_raw/suppl_listing.csv              (file, size per series)
"""
from __future__ import annotations

import csv
import sys

import pandas as pd
from p2_common import REPO_P2, fetch, list_suppl, parse_soft_samples, parse_soft_series

CANDIDATES = [
    # MASH references
    "GSE316023", "GSE246328", "GSE294629",
    # fasting / circadian
    "GSE107787", "GSE117134", "GSE225697",
    # PPARa
    "GSE282998", "GSE8295",
    # external microplastic (Track B status only)
    "GSE245069", "GSE277772", "GSE261893", "GSE261725", "GSE286116", "GSE262399",
]


def main(accs: list[str]) -> None:
    out = REPO_P2 / "metadata" / "geo_raw"
    out.mkdir(parents=True, exist_ok=True)
    listing_rows = []
    for acc in accs:
        soft = fetch(acc, "soft/", f"{acc}_family.soft.gz")
        s = pd.DataFrame(parse_soft_samples(soft))
        s.to_csv(out / f"{acc}_samples_raw.csv", index=False)
        ser = parse_soft_series(soft)
        with open(out / f"{acc}_series.txt", "w", encoding="utf-8") as fh:
            for k, v in ser.items():
                fh.write(f"{k}: {v}\n")
        try:
            for name, size in list_suppl(acc):
                listing_rows.append(dict(accession=acc, file=name, size=size))
        except Exception as exc:  # noqa: BLE001
            listing_rows.append(dict(accession=acc, file=f"LISTING_FAILED: {exc}", size=""))
        print(f"  {acc}: {len(s)} samples")
    with open(out / "suppl_listing.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "file", "size"])
        w.writeheader()
        w.writerows(listing_rows)


if __name__ == "__main__":
    main(sys.argv[1:] or CANDIDATES)
