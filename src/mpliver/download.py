"""Download source files and record their provenance (URL, date, size, checksums)."""

import gzip
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

# NCBI's FTP-over-HTTPS server rejects the default python-requests user agent.
HEADERS = {"User-Agent": "mpliver/0.1 (BT3041 course project; GEO data download)"}


def file_checksums(path: Path) -> dict:
    """Return MD5 and SHA-256 hex digests of a file."""
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            md5.update(chunk)
            sha256.update(chunk)
    return {"md5": md5.hexdigest(), "sha256": sha256.hexdigest()}


def download_file(url: str, dest: Path, retries: int = 5, timeout: int = 60) -> Path:
    """Download ``url`` to ``dest``, retrying incomplete transfers.

    The file is written to a temporary name and only renamed once the number of
    received bytes matches the server's Content-Length, so a partial download
    never masquerades as the real file.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_error = None
    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, headers=HEADERS, stream=True, timeout=timeout) as response:
                response.raise_for_status()
                expected = response.headers.get("Content-Length")
                received = 0
                with open(tmp, "wb") as handle:
                    for chunk in response.iter_content(chunk_size=1 << 16):
                        handle.write(chunk)
                        received += len(chunk)
            if expected is not None and received != int(expected):
                raise IOError(f"incomplete download: {received} of {expected} bytes")
            tmp.replace(dest)
            return dest
        except (requests.RequestException, IOError) as error:
            last_error = error
            time.sleep(2 * attempt)
    tmp.unlink(missing_ok=True)
    raise RuntimeError(f"failed to download {url} after {retries} attempts: {last_error}")


def is_valid_gzip(path: Path) -> bool:
    """True if the whole file decompresses without error (catches truncated downloads)."""
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(1 << 20):
                pass
        return True
    except (OSError, EOFError):
        return False


def provenance_record(key: str, url: str, path: Path, description: str) -> dict:
    """Build a provenance entry for one downloaded file."""
    return {
        "key": key,
        "url": url,
        "local_path": path.as_posix(),
        "description": description,
        "downloaded_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "size_bytes": path.stat().st_size,
        **file_checksums(path),
    }


def load_provenance(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return {}


def verify_against_provenance(path: Path, record: dict) -> bool:
    """True if ``path`` has the SHA-256 recorded at download time."""
    return path.exists() and file_checksums(path)["sha256"] == record["sha256"]
