"""Download the remaining GSE294180 FASTQs from ENA and quantify each sample as soon as
both of its files are verified (download and quantification overlap).

- Each file is fetched as parallel HTTP byte-range segments (ENA supports ranges).
  Segments are retried individually; finished segments are recorded so a restart resumes.
- Every file is checked against the ENA MD5 in the manifest before use.
- When both mates of a run pass, fastp -> Salmon runs in WSL via quantify_run.sh,
  at most QUANT_WORKERS samples at a time. FASTQs are kept.

Usage (Windows, repository venv):
    python pipelines/quantification/download_and_quantify.py [connections] [quant_workers]
"""

import hashlib
import json
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests

from mpliver.config import REPO_ROOT

WORK = Path("D:/mpliver_work")
OUT = WORK / "fastq_full"
LOG = WORK / "logs" / "full_requant.log"
SEGMENT = 32 * 1024 * 1024
HEADERS = {"User-Agent": "mpliver/0.1 (BT3041 course project; ENA FASTQ download)"}
_local = threading.local()
_lock = threading.Lock()


def log(msg: str) -> None:
    line = f"{time.strftime('%H:%M:%S')} {msg}"
    with _lock:
        print(line, flush=True)
        with open(LOG, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


def session() -> requests.Session:
    if not hasattr(_local, "s"):
        _local.s = requests.Session()
        _local.s.headers.update(HEADERS)
    return _local.s


class RemoteFile:
    def __init__(self, url, size, md5, path, on_verified):
        self.url, self.size, self.md5, self.path = url, size, md5, path
        self.done_marker = path.with_suffix(path.suffix + ".segments")
        self.on_verified = on_verified
        self.segments = [(o, min(o + SEGMENT, size) - 1) for o in range(0, size, SEGMENT)]
        finished = (set(json.loads(self.done_marker.read_text()))
                    if self.done_marker.exists() else set())
        self.remaining = {i for i in range(len(self.segments)) if i not in finished}
        self.finished = finished
        if not path.exists() or path.stat().st_size != size:
            with open(path, "wb") as fh:
                fh.truncate(size)
            self.finished, self.remaining = set(), set(range(len(self.segments)))

    def fetch(self, i: int) -> None:
        start, end = self.segments[i]
        for attempt in range(12):
            try:
                r = session().get(self.url, headers={"Range": f"bytes={start}-{end}"},
                                  timeout=(20, 120))
                if r.status_code != 206 or len(r.content) != end - start + 1:
                    raise IOError(f"status {r.status_code}, {len(r.content)} bytes")
                with open(self.path, "r+b") as fh:
                    fh.seek(start)
                    fh.write(r.content)
                break
            except (requests.RequestException, IOError):
                time.sleep(min(2 ** attempt, 30))
        else:
            raise RuntimeError(f"segment {i} of {self.path.name} failed")
        with _lock:
            self.finished.add(i)
            self.remaining.discard(i)
            self.done_marker.write_text(json.dumps(sorted(self.finished)))
            complete = not self.remaining
        if complete:
            self.verify()

    def verify(self) -> None:
        h = hashlib.md5()
        with open(self.path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 24), b""):
                h.update(chunk)
        if h.hexdigest() == self.md5:
            self.done_marker.unlink(missing_ok=True)
            log(f"verified {self.path.name}")
            self.on_verified(self)
        else:
            log(f"MD5 MISMATCH {self.path.name}; re-downloading")
            self.done_marker.unlink(missing_ok=True)
            with _lock:
                self.finished, self.remaining = set(), set(range(len(self.segments)))
            for i in range(len(self.segments)):
                self.fetch(i)


def to_wsl(p: Path) -> str:
    return "/mnt/" + p.drive[0].lower() + p.as_posix()[2:]


def quantify(sample: str, run: str, r1: Path, r2: Path) -> None:
    repo = to_wsl(REPO_ROOT)
    cmd = (f"cd {repo} && FASTQ_R1={to_wsl(r1)} FASTQ_R2={to_wsl(r2)} THREADS=6 "
           f"bash pipelines/quantification/quantify_run.sh {run} /mnt/d/mpliver_work {sample}")
    log(f"quantifying {sample} ({run})")
    res = subprocess.run(["wsl", "-d", "Ubuntu-24.04", "-e", "bash", "-c", cmd],
                         capture_output=True, text=True)
    log(f"{sample}: exit {res.returncode} {res.stdout.strip().splitlines()[-1:]}")


def main() -> int:
    connections = int(sys.argv[1]) if len(sys.argv) > 1 else 96
    quant_workers = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    OUT.mkdir(parents=True, exist_ok=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(REPO_ROOT / "pipelines" / "quantification" /
                           "full_requant_manifest.tsv", sep="\t")
    todo = manifest[~manifest["already_quantified_in_pilot"]]
    quant_pool = ThreadPoolExecutor(max_workers=quant_workers)
    verified: dict[str, set] = {}
    quant_jobs = []

    def make_callback(sample, run, r1, r2):
        def cb(remote):
            with _lock:
                verified.setdefault(run, set()).add(remote.path)
                ready = len(verified[run]) == 2
            if ready:
                quant_jobs.append(quant_pool.submit(quantify, sample, run, r1, r2))
        return cb

    files, tasks = [], []
    for _, row in todo.iterrows():
        r1 = OUT / f"{row['run']}_1.fastq.gz"
        r2 = OUT / f"{row['run']}_2.fastq.gz"
        cb = make_callback(row["sample"], row["run"], r1, r2)
        for url, size, md5, path in [(row.r1_url, row.r1_bytes, row.r1_md5, r1),
                                     (row.r2_url, row.r2_bytes, row.r2_md5, r2)]:
            f = RemoteFile(url, int(size), md5, path, cb)
            files.append(f)
            if not f.remaining:
                f.verify()
            tasks.extend((f, i) for i in sorted(f.remaining))
    total = sum(f.size for f in files)
    log(f"{len(todo)} runs, {len(files)} files, {total / 1e9:.1f} GB, "
        f"{len(tasks)} segments, {connections} connections")

    start = time.time()

    def progress():
        while any(f.remaining for f in files):
            got = sum(len(f.finished) * SEGMENT for f in files)
            el = time.time() - start
            log(f"progress {min(got, total) / 1e9:.1f}/{total / 1e9:.1f} GB, "
                f"{got / 1e6 / max(el, 1):.1f} MB/s")
            time.sleep(60)
    threading.Thread(target=progress, daemon=True).start()

    with ThreadPoolExecutor(max_workers=connections) as pool:
        for fut in [pool.submit(f.fetch, i) for f, i in tasks]:
            fut.result()
    log(f"all downloads finished in {(time.time() - start) / 60:.1f} min; "
        "waiting for quantification")
    quant_pool.shutdown(wait=True)
    log("ALL DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
