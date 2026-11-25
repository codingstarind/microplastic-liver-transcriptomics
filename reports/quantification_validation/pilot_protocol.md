# Re-quantification pilot: protocol and decision rules (written before the pilot was run)

Date fixed: 2026-11-02, before any pilot output existed. Results are reported in
`reports/quantification_validation/quantification_validation_report.md`. They are
judged against these rules, which will not be changed after seeing the data.

## Question

Which tximport aggregation produced the deposited integer matrix
`GSE294180_MP_raw_counts_gene_level.csv.gz`? Can it be inferred reliably, and does
our planned DESeq2 analysis require full re-quantification?

## Samples (one per group, chosen to include atypical libraries)

| Sample | Run | Group | Why chosen |
|---|---|---|---|
| CD_1 | SRR33050963 | CD+Veh | deepest library (19.1 M pairs); **quantified twice** to measure run-to-run noise |
| CD_MP_5 | SRR33050950 | CD+PE | shallowest library (13.5 M pairs), lowest counts-per-read-pair proxy in its group |
| MD_1 | SRR33050959 | MD+Veh | heaviest skeletal-muscle contamination |
| MD_MP_4 | SRR33050946 | MD+PE | highest cholangiocyte/stellate marker reads; low proxy |

## Processing (reproducing the authors' reported choices)

- **FASTQ source:** ENA. Streamed, not stored.
- **fastp 0.23.4:** `--detect_adapter_for_pe --trim_poly_g --trim_poly_x`; otherwise defaults.
- **Salmon 1.10.3:** `quant -l A --validateMappings --seqBias --useVBOpt`. Index from the GENCODE M35 transcript FASTA, `-k 31 --gencode`, no decoys, default duplicate removal.
- **Assumptions the authors did not state:** decoys, k-mer size, duplicate handling, fastp options beyond those named, tx2gene key. Each is listed in the report.

## Comparisons

1. For each pilot sample, run tximport three times: `countsFromAbundance` = `"no"`, `"scaledTPM"` and `"lengthScaledTPM"`, with tx2gene by **gene symbol** to match the deposited matrix. Then round to integers.
   - `"lengthScaledTPM"` uses the average transcript length *across the samples imported together*, so a 4-sample import cannot reproduce a 19-sample import exactly. This is noted as a structural limit of the pilot.
2. Compare each version with the deposited column for the same sample, on the genes present in both.

## Pre-specified metrics (genes with deposited count ≥ 20 unless stated)

| Metric | Purpose |
|---|---|
| **M1** column-total ratio (ours / deposited, all genes) | Gross scale agreement |
| **M2** median \|log2((deposited+1)/(ours+1))\| (MALR) | Typical gene-level disagreement |
| **M3** fraction of genes within ±5% | Closeness |
| **M4** fraction of genes with exactly equal integers (all genes with count > 0 in either) | Exact reproduction |
| **M5** slope of log2 ratio on log2 gene effective length (from Salmon, abundance-weighted) | **Main discriminator:** `"no"` counts scale with effective length while `"scaledTPM"` counts do not, so comparing against the wrong setting leaves a length trend |
| **M6** overlap of the gene-symbol set (which of the 1,061 symbols missing from the deposited matrix are also missing from ours) | Structural fingerprint of the index and tx2gene construction |

**Noise floor:** the same metrics computed between the two independent quantifications of CD_1 (`"no"` setting). This captures run-to-run variation from multithreading, `--seqBias` sampling and VB optimisation.

## Decision rules

- A setting is **consistent** with the deposited matrix if, in **all four** samples:
  - M2 ≤ 2 × the noise-floor M2 (or ≤ 0.02 if the floor is smaller);
  - |M5| ≤ 0.02 log2 per log2-length unit;
  - M1 within 1 ± 0.02.
- The original setting is **identified** only if exactly one setting is consistent *and* every other setting has M2 at least 2× larger or |M5| > 0.05.
- If none is consistent (for example because the authors used a different index, decoys or other trimming), the setting is **not identifiable**, whatever the correlations. **High correlation alone is never taken as recovery.**
- **Consequence for the analysis:**
  - If `"no"` is identified, or not identifiable: the deposited matrix lacks its required offset. Full re-quantification (19 samples) is needed for the confirmatory DESeq2 analysis. Deposited-matrix results stay provisional.
  - If `"scaledTPM"` or `"lengthScaledTPM"` is identified: the deposited matrix may be used directly as counts (the vignette's "bias-corrected counts without an offset" route). Full re-quantification is optional, used as a sensitivity check.
