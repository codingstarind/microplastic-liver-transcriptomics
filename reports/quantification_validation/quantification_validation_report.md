# RNA-seq quantification validation: can the deposited count matrix be used, and how?

**Status: IN PROGRESS.** Sections 1–3 are complete. Section 4 (pilot results) is pending: the four pilot runs are being downloaded manually because scripted downloads are slow. Section 5 lists what the outcome will decide. Decision rules were fixed in advance in [pilot_protocol.md](pilot_protocol.md).

## 1. The problem

- The deposited matrix (`GSE294180_MP_raw_counts_gene_level.csv.gz`) has integer gene-level counts from Salmon v1.10.3 → tximport, summarised to **gene symbols**.
- The tximport vignette (Bioconductor, read 2026-11-02) advises against passing *original* gene-level estimated counts (`countsFromAbundance = "no"`) to count-based methods **without** the average-transcript-length offset. The exception is data without length bias (3′-tagged or unfragmented), which this full-length poly(A) library is not.
- Two routes are acceptable:
  - original counts plus offset (`DESeqDataSetFromTximport`);
  - bias-corrected counts (`"scaledTPM"` or `"lengthScaledTPM"`) used directly.
- The length matrix was not deposited, and the setting was not reported.

## 2. Can the setting be recovered from published sources?

| Source | Checked | What it gives |
|---|---|---|
| GEO series/sample records | Yes | Tool versions, Salmon flags (`--validateMapping -l A --seqBias --useVBOpt`), GENCODE GRCm39 M35, "Raw counts". **No tximport setting.** |
| Paper methods (Europe PMC full text) | Yes | Same as GEO, plus "After estimateSizeFactors were calculated, the DESeq function was used". **No tximport setting**, no `DESeqDataSetFromTximport` mention. |
| Zenodo record 19603085 (authors') | Yes (API metadata) | Spatial data only (h5, cell metadata, UMAP). No bulk processing code or files. |
| Supplementary tables S1–S18 | **Not obtained.** PMC serves a proof-of-work bot challenge to scripts; manual download requested. | Unknown. Their contents are not described or used here. |
| Code repository | None referenced in the paper | – |

**Conclusion:** the setting **cannot be recovered from published information**. The only route is empirical re-quantification.

## 3. What the deposited matrix itself can and cannot tell us

**Established:**

- integers, so rounding happened after tximport;
- all 55,836 identifiers are GENCODE M35 gene symbols;
- 228 symbols map to more than one Ensembl gene;
- 1,061 annotated symbols are absent, 407 of them on chrY;
- column totals are 10.4–15.5 M.

**Not discriminating:**

- Column totals are preserved under all three tximport settings.
- The settings differ only through per-sample effective lengths, which were not deposited.

**Testable prediction (from M6 of the protocol):** Salmon's default index removes transcripts with identical sequences. Genes made only of such duplicates (e.g. multi-copy chrY genes, rRNA and snoRNA copies) then have no rows after tximport. If our default index loses the same symbols as the deposited matrix, that structurally links the two pipelines.

## 4. Pilot (pending)

**Infrastructure built:**

- micromamba environments in WSL: `mpliver-quant` (Salmon 1.10.3, fastp 0.23.4, sra-tools 3.2.1) and `mpliver-r` (R 4.4.3, tximport 1.34.0, DESeq2 1.46.0, …). Specs in `envs/`.
- GENCODE M35 transcript FASTA downloaded to `D:\mpliver_work\reference`. The MD5 matches GENCODE's published checksum; 149,138 transcripts.
- Salmon index: k = 31, `--gencode`, no decoys, default duplicate handling. Built in 3.7 min; 645 MB; 0.7 GB peak RAM.
- tx2gene tables by symbol and by Ensembl ID: 56,897 symbols, 57,132 genes.

**Samples:** CD_1 (quantified twice for the noise floor), CD_MP_5, MD_1, MD_MP_4. One per group, chosen to include atypical libraries.

**Infrastructure issues encountered (documented):**

1. **C: drive full.** The Windows system drive had ~10 GB free, and WSL's virtual disk filled it while the R environment and a temporary index were being built. Recovered by removing caches and failed builds; the user also freed space. Large files now live on D:.
2. **Slow downloads.** ENA (~0.5 MB/s per stream) and the NCBI SRA AWS mirror (~2.5 MB/s) are both slow here; NCBI FTP refused scripted requests. The pilot `.sra` files (1.3–1.8 GB each) are being downloaded by the user in a browser.
   - `pipelines/quantification/run_pilot_from_downloads.sh` quantifies each file when complete and then deletes it.
   - The pipeline streams `fastq-dump` → `fastp` → FIFOs → Salmon, so no FASTQ is ever written.

## 5. What the pilot will decide (per the pre-specified rules)

| Outcome | Consequence |
|---|---|
| `"no"` identified, **or** no setting identifiable | The deposited matrix lacks its required offset. Full re-quantification of 19 runs (~31 GB of `.sra`) is needed for confirmatory DESeq2 (`DESeqDataSetFromTximport`). All current results stay provisional. |
| `"scaledTPM"` or `"lengthScaledTPM"` identified | The deposited matrix is a valid bias-corrected count matrix. Current DESeq2 results can be promoted from provisional; re-quantification is an optional sensitivity check. |

Either way, the pilot also provides the **real Salmon mapping rates and fastp statistics** for 4 of the 19 samples. These are needed to interpret the PE-associated library-quality proxy.
