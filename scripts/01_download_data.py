"""Step 1: download the GEO files for GSE294180 and record provenance.

Usage:
    python scripts/01_download_data.py            # download missing files, verify existing ones
    python scripts/01_download_data.py --refresh  # re-download and overwrite the provenance entry

Raw files are never modified after download. If a file already exists, its
SHA-256 is checked against metadata/data_provenance.json; a mismatch stops the
script instead of silently replacing the data.

If NCBI refuses scripted downloads, download the files in a browser from the
URLs in src/mpliver/config.py, put them in data/raw/ under the same names and
run this script: it checks they are complete and records their checksums.
"""

import argparse
import json
import sys

from mpliver.config import GEO_ACCESSION, PROVENANCE_FILE, RAW_DIR, REPO_ROOT, SOURCE_FILES
from mpliver.download import (
    download_file,
    is_valid_gzip,
    load_provenance,
    provenance_record,
    verify_against_provenance,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--refresh", action="store_true", help="re-download all files")
    args = parser.parse_args()

    provenance = load_provenance(PROVENANCE_FILE)
    files = provenance.setdefault("files", {})
    provenance["accession"] = GEO_ACCESSION
    provenance["geo_page"] = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={GEO_ACCESSION}"

    failures = []
    for key, spec in SOURCE_FILES.items():
        dest = RAW_DIR / spec["filename"]
        if dest.exists() and key in files and not args.refresh:
            if verify_against_provenance(dest, files[key]):
                print(f"[ok] {dest.name}: checksum matches provenance record")
                continue
            print(f"[error] {dest.name} differs from the recorded checksum. "
                  "Investigate before re-downloading with --refresh.")
            return 1
        manual = dest.exists() and key not in files and not args.refresh
        if manual:
            # File placed by hand (e.g. downloaded in a browser when NCBI refuses scripts).
            if not is_valid_gzip(dest):
                print(f"[error] {dest.name} is not a complete gzip file; download it again.")
                return 1
            print(f"[register] {dest.name}: manually downloaded file, recording checksums")
        else:
            print(f"[download] {spec['url']}")
            try:
                download_file(spec["url"], dest)
            except RuntimeError as error:
                print(f"[failed] {error}")
                failures.append(key)
                continue
        files[key] = provenance_record(key, spec["url"], dest, spec["description"])
        files[key]["local_path"] = dest.relative_to(REPO_ROOT).as_posix()
        files[key]["acquisition"] = ("manual browser download from url, placed in data/raw"
                                     if manual else "scripts/01_download_data.py")
        record = files[key]
        print(f"[saved] {dest.name} ({record['size_bytes']} bytes, "
              f"sha256 {record['sha256'][:12]}...)")

    PROVENANCE_FILE.parent.mkdir(parents=True, exist_ok=True)
    PROVENANCE_FILE.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(f"Provenance written to {PROVENANCE_FILE}")
    if failures:
        print(f"[error] not downloaded: {', '.join(failures)} (re-run the script later)")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
