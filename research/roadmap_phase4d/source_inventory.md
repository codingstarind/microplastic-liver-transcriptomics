# Source-of-truth inventory (Phase 4D, 2026-11-23)

Local copies were read; nothing earlier was modified. "Verified" means recomputed or read from the file in this phase. "Reported" means taken from a report.

| Path | Used for | Status |
|---|---|---|
| `data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz` | 19-mouse expression | **Verified:** 55,836 × 19 integers, no NA, no duplicates; columns = `metadata/sample_metadata.csv` (4/5/5/5) |
| `data/processed/gene_filter.csv` | scoring background | **Verified:** 13,777 genes pass |
| `results/deseq2/provisional_deposited/vst_matrix.csv.gz` | input sensitivity | **Verified:** 13,777 × 19 |
| `results/deseq2/covariates.csv` | QC | **Verified:** depth, quality_proxy (yield), muscle, clock_phase |
| `data/raw/gencode.vM35.annotation.gtf.gz` | gene length (TPM-proxy) | parsed, cached |
| `research/roadmap_phase2/` | universes (11,888 / alt 12,010), muscle list | Verified files |
| `research/roadmap_phase3/reference_effects/`, `results/ppara/` | signature construction | **Verified** (checksums match Phase 4B record) |
| `research/roadmap_phase3/results/positive_control/`, `posthoc_audit/` | locked MASH positive control: **FAIL on c1** (exact p = 0.036) | Verified (tests) |
| `research/roadmap_phase4/`, `roadmap_phase4b/`, `roadmap_phase4c_diagnostics/` | gene-level results | Verified (tests, metrics files) |
| `research/literature/`, `research/cross_study_replication/` | context | Reported |
| `research/*.zip` snapshots | **not used** (local folders are the source of truth) | — |
| Independent GPT audit memo | **NOT FOUND locally**; the user's summary was used | — |
| R `singscore` / `GSVA` | **not installed**; singscore implemented in Python from the primary formula (Foroutan 2018, Methods); **GSVA not run** | — |
