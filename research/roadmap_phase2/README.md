# Phase 2: reference dataset discovery, verification and preparation

**Branch:** `research/phase2-datasets` (not merged; for team review). **Prepared:** 2026-11-22.
**Scope:** dataset audit and data preparation only. No differential expression, signatures, regression, classifiers or biological interpretation were run in this phase.

## What is where

| Path | Content |
|---|---|
| `reports/dataset_selection_report.md` | Which MASH, fasting, circadian and PPARα references were selected or rejected, and why |
| `reports/gene_mapping_report.md` | Annotation versions, mapping rules, gene counts at each step, overlap with GSE294180 |
| `reports/qc_report.md` | QC, technical-vs-treatment checks, sample-overlap checks, plots in `qc/` |
| `reports/external_microplastic_status.md` | GSE245069, GSE277772, Hong-lab series (Track B status only) |
| `reports/phase3_readiness_report.md` | What is ready, what is not, valid and misleading comparisons, Phase 3 starting point |
| `metadata/dataset_master_inventory.csv` | One row per dataset (selected, secondary, rejected, Track B) |
| `metadata/sample_metadata_master.csv` | One row per biological sample of the 6 prepared references (metadata index, **not** an expression matrix) |
| `metadata/reference_contrast_plan.csv` | Planned contrasts with eligible GSM and matrix column IDs (generated from the metadata) |
| `metadata/gene_mapping_summary.csv`, `gene_universe_membership.csv.gz`, `key_gene_coverage.csv` | Gene-mapping counts and candidate universes |
| `metadata/download_manifest.csv` | URL, local path, size, MD5, SHA-256 and date of every downloaded file |
| `metadata/run_metadata.csv` | ENA run table (reads per run) for all audited series |
| `metadata/geo_raw/` | Parsed GEO SOFT sample and series fields for every audited series |
| `datasets/<GSE>/` | Per-dataset `sample_metadata.csv`, `gene_mapping.csv.gz` (one row per original identifier), `prepared_manifest.json` (format, notes, checksums of prepared matrices, mapping summary) |
| `qc/` | `qc_sample_metrics.csv`, `qc_technical_vs_treatment.csv`, `overlap_checks.csv`, one `<GSE>_qc.png` per dataset |
| `scripts/` | Rerunnable pipeline and tests |

**Large files are not in git.** Downloads and prepared matrices live in `D:/mpliver_work/phase2_refs/<GSE>/` (override with env `MPLIVER_P2_DATA`):
- `<GSE>/<deposited files>`: as downloaded from GEO, checksummed in `metadata/download_manifest.csv`.
- `<GSE>/prepared/expression_original.csv.gz`: the deposited matrix, all samples, original identifiers, values unchanged. Two exceptions: GSE316023 duplicated identifiers get a `#row<n>` suffix so that every row is addressable, and GSE107787/GSE117134 are read from xlsx/tsv.
- `<GSE>/prepared/expression_ready.csv.gz`: analysis-ready matrix. Rows are GENCODE vM35 symbols (the GSE294180 annotation); columns are exactly the `in_ready_matrix` rows of `datasets/<GSE>/sample_metadata.csv`, in order. Counts stay integer counts; RPKM stays RPKM. **Nothing is normalised, transformed or pooled.**

## How to regenerate

From `research/roadmap_phase2/scripts/` (Python ≥ 3.10 with pandas, numpy, scipy, matplotlib, openpyxl, pytest):

```powershell
python 01_fetch_geo_metadata.py
python 02_download_processed.py
python 03_prepare_references.py
python 04_fetch_run_metadata.py
python 05_qc.py
python 06_gene_universe.py
python 07_contrast_plan.py
python 08_phase3_inputs.py      # fixed muscle rule, gene universe, eligible samples, GSE316023 length offset
python -m pytest test_phase2_data.py -q
```

Every download step skips files that already exist and refreshes their checksums. Total download is about 150 MB (no FASTQ). NCBI throughput from this network was about 20–40 KB/s, so a first run takes about 15–20 minutes.

## Tests (`scripts/test_phase2_data.py`)

The tests check:
- duplicate sample IDs and GSMs;
- every GEO GSM is present (no unexpected missing samples);
- valid contrast labels;
- biological replicate counts per contrast, and that no contrast is empty;
- one metadata row per biological sample;
- no shared GSM, SRA run or identical read signature between the selected studies;
- known overlaps (GSE246328/GSE294629, Hong-lab controls) are recorded;
- one mapping row per original identifier;
- matrix columns equal the metadata, in order;
- unique, non-missing gene IDs, with integer values for count-type data;
- prepared-file and download checksums.
