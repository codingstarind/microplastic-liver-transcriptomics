"""Phase 5C step 1: download PRJNA916285 liver RNA-seq runs (full depth) and quantify them with the SAME
pipeline as the GSE294180 training data: fastq-dump -> fastp -> Salmon 1.10 (GENCODE vM35, k31, -l A,
--seqBias, --useVBOpt), i.e. pipelines/quantification/quantify_run.sh, unchanged (SRA_LOCAL mode).

Source: NCBI SRA open-data mirror on AWS (.sra files). ENA FASTQ was rate-limited to ~0-25 KB/s on
2026-11-23; the link itself caps at ~1.1 MB/s, so each file is fetched as 32 parallel HTTP byte-range
segments (resumable). Integrity: byte size = S3 Content-Length; spots dumped must equal the ENA read_count.
No labels are used here. Order: PBS1-4 + MPs1-4 (primary) first, then EGCG1-4 + EGCG_MPs1-4.

Usage: python c01_fetch_quant.py [primary|all]   (env: CONN, PAR_RUNS, QWORKERS, QTHREADS)
Outputs: D:/mpliver_work/phase5c/sra/<run>.sra ; D:/mpliver_work/phase5c/requant/<label>/salmon/quant.sf ;
         manifests/download_manifest.csv ; logs/c01_fetch_quant.log
"""
import csv
import hashlib
import json
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

P5C = Path(__file__).resolve().parents[1]
WORK = Path("D:/mpliver_work/phase5c")
SRA = WORK / "sra"
WHICH = sys.argv[1] if len(sys.argv) > 1 else "primary"
SEG = 16 * 1024 * 1024
import os
CONN = int(os.environ.get("CONN", 24))          # connections per run
PAR_RUNS = int(os.environ.get("PAR_RUNS", 3))   # runs downloaded at the same time
QWORKERS = int(os.environ.get("QWORKERS", 2))   # quantifications at the same time
QTHREADS = int(os.environ.get("QTHREADS", 10))  # Salmon/fastp threads per quantification
WSL = ["wsl.exe", "-d", "Ubuntu-24.04", "-u", "root", "-e", "bash", "-lc"]
QSCRIPT = "/mnt/c/Users/Soham/Desktop/microplastic/pipelines/quantification/quantify_run.sh"
HDR = {"User-Agent": "mpliver/0.1 (BT3041 course project)"}
LOG = P5C / "logs" / "c01_fetch_quant.log"
MAN = P5C / "manifests" / "download_manifest.csv"
_lock = threading.Lock()


def log(m):
    line = f"{time.strftime('%H:%M:%S')} {m}"
    with _lock:
        print(line, flush=True)
        with open(LOG, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


def label(title):
    g, i = title[:-1], title[-1]
    return {"PBS": "PBS", "MPs": "MPs", "EGCG": "EGCG", "EGCG_MPs": "EGCGMPs"}[g] + "_" + i


def fetch(url, out):
    size = int(requests.head(url, headers=HDR, timeout=60).headers["Content-Length"])
    done_f = out.with_suffix(".segments.json")
    segs = [(o, min(o + SEG, size) - 1) for o in range(0, size, SEG)]
    finished = set(json.loads(done_f.read_text())) if done_f.exists() and out.exists() and out.stat().st_size == size else set()
    if not out.exists() or out.stat().st_size != size:
        with open(out, "wb") as fh:
            fh.truncate(size)
    todo = [i for i in range(len(segs)) if i not in finished]
    s_local = threading.local()

    def get(i):
        s, e = segs[i]
        if not hasattr(s_local, "s"):
            s_local.s = requests.Session()
        for a in range(60):
            try:
                r = s_local.s.get(url, headers={**HDR, "Range": f"bytes={s}-{e}"}, timeout=(20, 600))
                if r.status_code == 206 and len(r.content) == e - s + 1:
                    with _lock:
                        with open(out, "r+b") as fh:
                            fh.seek(s)
                            fh.write(r.content)
                        finished.add(i)
                        done_f.write_text(json.dumps(sorted(finished)))
                    return
                raise IOError(r.status_code)
            except Exception as ex:  # noqa: BLE001
                time.sleep(min(2 ** a, 30))
                if a > 3:
                    log(f"retry {out.name} seg {i} ({ex})")
        raise RuntimeError(f"failed {url} seg {i}")

    with ThreadPoolExecutor(CONN) as ex:
        list(ex.map(get, todo))
    assert len(finished) == len(segs) and out.stat().st_size == size
    return size


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def quant(r, lab, sra, size):
    w = "/mnt/d/mpliver_work/phase5c"
    cmd = f"SRA_LOCAL={w}/sra/{sra.name} THREADS={QTHREADS} bash {QSCRIPT} {r['run_accession']} {w} {lab}"
    t0 = time.time()
    p = subprocess.run(WSL + [cmd], capture_output=True, text=True)
    log(f"{lab} quant rc={p.returncode} {time.time()-t0:.0f}s {p.stdout.strip()[-200:]} {p.stderr.strip()[-300:]}")
    if p.returncode != 0:
        return
    fj = json.loads((WORK / "requant" / lab / "fastp.json").read_text())
    pairs_in = fj["summary"]["before_filtering"]["total_reads"] // 2
    ok = pairs_in == int(r["read_count"])
    log(f"{lab} spots dumped {pairs_in} vs ENA read_count {r['read_count']} -> {'OK' if ok else 'MISMATCH'}")
    new = not MAN.exists()
    with _lock, open(MAN, "a", newline="") as fh:
        wr = csv.writer(fh)
        if new:
            wr.writerow(["run", "label", "sample_accession", "sample_title", "ena_read_pairs", "pairs_dumped", "read_count_match",
                         "sra_bytes", "sra_sha256", "source_url", "fetched"])
        wr.writerow([r["run_accession"], lab, r["sample_accession"], r["sample_title"], r["read_count"], pairs_in, ok, size, sha256(sra),
                     f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{r['run_accession']}/{r['run_accession']}", time.strftime("%Y-%m-%d %H:%M")])


def main():
    runs = list(csv.DictReader(open(P5C / "manifests/ena_filereport_PRJNA916285.tsv"), delimiter="\t"))
    order = {"PBS": 0, "MPs": 1, "EGCG": 2, "EGCG_MPs": 3}
    runs.sort(key=lambda r: (order[r["sample_title"][:-1]], r["sample_title"][-1]))
    # interleave PBS/MPs so a usable 2v2, 3v3 ... exists early
    prim = [r for r in runs if r["sample_title"][:-1] in ("PBS", "MPs")]
    prim = [x for pair in zip(prim[:4], prim[4:]) for x in pair]
    sec = [r for r in runs if r["sample_title"][:-1] in ("EGCG", "EGCG_MPs")]
    sec = [x for pair in zip(sec[:4], sec[4:]) for x in pair]
    runs = prim + (sec if WHICH == "all" else [])
    SRA.mkdir(parents=True, exist_ok=True)
    done = {x["run"] for x in csv.DictReader(open(MAN))} if MAN.exists() else set()
    qpool = ThreadPoolExecutor(QWORKERS)
    futs = []

    def one(r):
        lab = label(r["sample_title"])
        out = SRA / f"{r['run_accession']}.sra"
        t0 = time.time()
        try:
            size = fetch(f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{r['run_accession']}/{r['run_accession']}", out)
        except Exception as ex:  # noqa: BLE001  (one failed run must not stop the others; rerun resumes it)
            log(f"{lab} {r['run_accession']} DOWNLOAD FAILED ({ex}); rerun the script to resume")
            return
        log(f"{lab} {r['run_accession']} downloaded {size/1e9:.2f} GB in {(time.time()-t0)/60:.1f} min ({size/1e6/max(1,time.time()-t0):.2f} MB/s)")
        with _lock:
            futs.append(qpool.submit(quant, r, lab, out, size))

    with ThreadPoolExecutor(PAR_RUNS) as dl:
        list(dl.map(one, [r for r in runs if r["run_accession"] not in done]))
    for f in list(futs):
        f.result()
    log("ALL DONE")


if __name__ == "__main__":
    main()
