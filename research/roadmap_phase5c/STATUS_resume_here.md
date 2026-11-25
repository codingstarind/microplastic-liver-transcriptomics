# Phase 5C: status (paused 2026-11-23, resume next session)

**Nothing has been scored on the new data yet.** The protocol is not yet written or locked, and no external result exists.

## Accession audit (done)
- **PRJNA1513170, PRJNA1536617, PRJNA1305997, PRJNA1476926** (the four discovery-report leads):
  - NCBI BioProject: "is not public";
  - E-utilities: no summary;
  - ENA: 0 runs.
  - **No data can be downloaded.** They are embargoed or withdrawn at the source.
- **PRJNA916285** (Zhang et al. 2023, *ACS Nano* 17:15125, PMID 37486121) is the **selected study**:
  - polystyrene 5 µm, oral; C57BL/6 male, 10 weeks old; liver bulk RNA-seq;
  - 16 distinct BioSamples / 16 runs (biological replicates 1–16): PBS 4, MPs 4, EGCG 4, EGCG+MPs 4;
  - never used in any earlier phase, so it is **untouched**;
  - dose and duration are NOT VERIFIED (the paper is paywalled).
- **GSE291616 / PRJNA1234520:** snRNA-seq, PBS vs 50 nm vs 500 nm PS, 2 mice per arm, about 20 GB per sample. Excluded (2 v 2, nuclei, about 120 GB).
- **PRJNA913771:** airborne nanoplastics, 3 runs, no control arm deposited. Excluded.
- **PRJNA1248512** is GSE294180 itself.
- **Already seen (comparators only):** GSE245069, GSE277772 and the Hong series.

## Download
- **Script:** `scripts/c01_fetch_quant.py all`. It downloads `.sra` files from the AWS SRA mirror with 32 range connections and is resumable. It then runs fastp → Salmon vM35, the same pipeline as GSE294180.
- **Speed:** the link is about 1.1 MB/s, so each run takes about 35 min; all 16 runs (about 40 GB) take about 9–10 h. ENA was rate-limited.
- **Files:**
  - progress: `logs/c01_fetch_quant.log`
  - finished runs: `manifests/download_manifest.csv`
  - data: `D:/mpliver_work/phase5c/`
- **To resume after an interruption,** rerun the same command; finished segments and runs are skipped.

## Written but NOT yet locked or run
- `c00_common.py`
- `c02_prepare.py`
- `c02_de.R`
- `c03_score.py` (locked rules: primary F8-OLS PE_in_CD and FAST singscore; MPs vs PBS 4 v 4, 70 labellings)
- `c04_biology.py`

## Next steps
1. Write `03_locked_protocol.md` and a `c00_lock.py` hash lock.
2. Run `c03_score.py --dry-run` (label-blind).
3. Once the 8 primary mice are quantified: run `c02_prepare` → `c02_de.R` (WSL) → `c03` → `c04`.
4. Then figures, tests, the report and the Phase 6 memo.
