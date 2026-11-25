"""Shared helpers for Phase 2 data preparation.

Large downloaded files live outside the repository in DATA_ROOT
(default D:/mpliver_work/phase2_refs, override with env MPLIVER_P2_DATA).
Every download is skipped if the file already exists, and every file's
size, MD5 and SHA-256 are appended to a per-dataset manifest in the repo.
"""
from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import os
import re
import time
import urllib.request
from pathlib import Path

REPO_P2 = Path(__file__).resolve().parents[1]
DATA_ROOT = Path(os.environ.get("MPLIVER_P2_DATA", "D:/mpliver_work/phase2_refs"))
USER_AGENT = "mpliver-phase2-prep/1.0 (academic research script)"


def geo_stub(acc: str) -> str:
    """GSE316023 -> GSE316nnn ; GSE8295 -> GSE8nnn."""
    num = acc[3:]
    return f"{acc[:3]}{num[:-3] if len(num) > 3 else ''}nnn"


def geo_ftp(acc: str, sub: str = "") -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_stub(acc)}/{acc}/{sub}"


def ds_dir(acc: str) -> Path:
    d = DATA_ROOT / acc
    d.mkdir(parents=True, exist_ok=True)
    return d


def file_hashes(path: Path) -> tuple[str, str]:
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            md5.update(chunk)
            sha.update(chunk)
    return md5.hexdigest(), sha.hexdigest()


def download(url: str, dest: Path, retries: int = 3) -> bool:
    """Download url to dest unless dest exists. Returns True if downloaded now."""
    if dest.exists() and dest.stat().st_size > 0:
        return False
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as out:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    out.write(b)
            tmp.replace(dest)
            return True
        except Exception as exc:  # noqa: BLE001
            print(f"  retry {attempt + 1} for {url}: {exc}")
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"failed to download {url}")


def record_manifest(acc: str, url: str, path: Path) -> None:
    """Append/refresh one row in metadata/download_manifest.csv."""
    man = REPO_P2 / "metadata" / "download_manifest.csv"
    rows = []
    if man.exists():
        with open(man, newline="", encoding="utf-8") as fh:
            rows = [r for r in csv.DictReader(fh) if r["local_path"] != str(path)]
    md5, sha = file_hashes(path)
    mtime = dt.datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
    rows.append(
        dict(
            accession=acc,
            file=path.name,
            url=url,
            local_path=str(path),
            bytes=path.stat().st_size,
            md5=md5,
            sha256=sha,
            downloaded=mtime,
        )
    )
    rows.sort(key=lambda r: (r["accession"], r["file"]))
    with open(man, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def fetch(acc: str, sub: str, name: str) -> Path:
    url = geo_ftp(acc, sub) + name
    dest = ds_dir(acc) / name
    new = download(url, dest)
    print(f"{'downloaded' if new else 'present   '} {acc}/{name} ({dest.stat().st_size:,} B)")
    record_manifest(acc, url, dest)
    return dest


def list_suppl(acc: str) -> list[tuple[str, str]]:
    """Return [(filename, size_text)] from the GEO suppl directory listing."""
    req = urllib.request.Request(geo_ftp(acc, "suppl/"), headers={"User-Agent": USER_AGENT})
    html = urllib.request.urlopen(req, timeout=60).read().decode("utf-8", "replace")
    out = []
    for m in re.finditer(r'<a href="([^"/?][^"]*)">[^<]*</a>\s+(\S+ \S+)\s+(\S+)', html):
        out.append((m.group(1), m.group(3)))
    return out


def parse_soft_samples(soft_gz: Path) -> list[dict]:
    """Parse ^SAMPLE blocks of a GEO family SOFT file into flat dicts.

    Repeated keys (characteristics, supplementary files, relations) are
    joined with ' | '. Characteristics 'key: value' pairs also become
    columns named ch_<key>.
    """
    samples, cur = [], None
    with gzip.open(soft_gz, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                cur = {"gsm": line.split("=", 1)[1].strip()}
                samples.append(cur)
            elif line.startswith("^"):
                cur = None
            elif cur is not None and line.startswith("!Sample_"):
                k, _, v = line[8:].partition(" = ")
                if k.startswith("characteristics_ch") and ":" in v:
                    ck, _, cv = v.partition(":")
                    col = "ch_" + re.sub(r"\W+", "_", ck.strip().lower()).strip("_")
                    cur[col] = (cur[col] + " | " + cv.strip()) if col in cur else cv.strip()
                cur[k] = (cur[k] + " | " + v) if k in cur else v
    return samples


def parse_soft_series(soft_gz: Path) -> dict:
    out: dict = {}
    with gzip.open(soft_gz, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("^SAMPLE"):
                break
            if line.startswith("!Series_"):
                k, _, v = line[8:].rstrip("\n").partition(" = ")
                out[k] = (out[k] + " | " + v) if k in out else v
    return out
