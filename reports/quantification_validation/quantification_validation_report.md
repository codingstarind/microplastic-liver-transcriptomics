# RNA-seq quantification validation: can the deposited count matrix be used, and how?

**Status: PILOT COMPLETE (2026-11-02).** Decision rules were fixed in advance in [pilot_protocol.md](pilot_protocol.md).

**Verdict:**

- By the pre-registered rules, the original tximport setting is **not identifiable**. No setting meets all the closeness thresholds in all samples.
- `scaledTPM` is clearly **excluded**. The evidence strongly favours `countsFromAbundance = "no"` (original estimated counts, rounded) from a Salmon run configured slightly differently from ours, and that matrix needs a length offset the deposited data lacks.
- Per the pre-specified consequence, **full re-quantification of all 19 samples is required** for confirmatory DESeq2.
- Results from the deposited matrix remain **provisional**. See §4–§5.

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

**Samples:**

- **Planned:** CD_1 (twice), CD_MP_5, MD_1, MD_MP_4, chosen to include atypical libraries.
- **Actually run (deviation):** CD_1 (SRR33050963, quantified twice), CD_MP_3 (SRR33050952), MD_4 (SRR33050956) and MD_MP_3 (SRR33050947).
- The team downloaded these from ENA because of the bandwidth limits. All 8 FASTQ.gz files passed ENA MD5 checks (`D:\mpliver_work\fastq_pilot\validation.json`).
- **Added afterwards:** CD_MP_5 (SRR33050950), the shallowest library with the lowest counts-per-read-pair proxy, from the original plan. Both files were MD5-verified independently against ENA.
- The pilot now covers one mouse per group, the deepest library (CD_1) and the shallowest (CD_MP_5). It does not include the most contaminated sample (MD_1) or MD_MP_4.

**Infrastructure issues encountered (documented):**

1. **C: drive full.** The Windows system drive had ~10 GB free, and WSL's virtual disk filled it while the R environment and a temporary index were being built. Recovered by removing caches and failed builds; the user also freed space. Large files now live on D:.
2. **Slow downloads.** ENA (~0.5 MB/s per stream) and the NCBI SRA AWS mirror (~2.5 MB/s) are both slow here; NCBI FTP refused scripted requests. The pilot `.sra` files (1.3–1.8 GB each) are being downloaded by the user in a browser.
   - `pipelines/quantification/run_pilot_from_downloads.sh` quantifies each file when complete and then deletes it.
   - The pipeline streams `fastq-dump` → `fastp` → FIFOs → Salmon, so no FASTQ is ever written.

### 4.1 Run statistics (the real values the deposited data lacked)

| Sample | Group | fastp reads in → out | Fragments processed | Mapped | Mapping rate | Library type |
|---|---|---|---|---|---|---|
| CD_1 | CD+Veh | 38.11 M → 37.38 M | 18.69 M | 16.08 M | **86.1%** | ISR |
| CD_MP_3 | CD+PE | 31.27 M → 30.45 M | 15.23 M | 12.64 M | **83.0%** | ISR |
| CD_MP_5 | CD+PE | 27.07 M → 26.33 M | 13.16 M | 11.11 M | **84.4%** | ISR |
| MD_4 | MD+Veh | 33.99 M → 33.26 M | 16.63 M | 14.24 M | **85.6%** | ISR |
| MD_MP_3 | MD+PE | 29.98 M → 29.23 M | 14.61 M | 12.43 M | **85.0%** | ISR |

- Library type ISR means stranded, read 1 antisense, consistent with the stranded kit.
- fastp removed 1.9–2.7% of reads.
- PE samples map 0.6–3.0 points lower than the vehicle sample of the same diet. That is consistent with the lower gene-counts-per-read-pair proxy, but with 2 vehicle samples it is only indicative.
- Output: `results/quantification_validation/pilot_mapping_stats.csv`.

### 4.2 Pre-specified metrics (`results/quantification_validation/pilot_metrics.csv`)

Noise floor: CD_1 quantified twice. M2 = 0.000; 99.5% of genes give identical integers; M1 = 1.000. **Salmon is reproducible here**, so the M2 threshold is the protocol's 0.02.

| Setting | M1 total ratio (ours/deposited) | M2 median \|log2 ratio\| | M3 within ±5% | M4 exact integers | M5 length slope | Spearman |
|---|---|---|---|---|---|---|
| `"no"` | 1.045–1.072 | **0.029–0.037** | 0.68–0.73 | **0.30–0.33** | −0.024 to −0.005 | 0.995–0.996 |
| `"lengthScaledTPM"` | 1.045–1.072 | 0.066–0.072 | 0.49–0.52 | 0.21–0.22 | −0.016 to 0.021 | 0.991–0.993 |
| `"scaledTPM"` | 1.045–1.072 | 1.11–1.32 | 0.02–0.03 | 0.05 | **0.96–0.97** | 0.86–0.88 |

Ranges are over the 5 samples.

**Applying the decision rules:**

- **`scaledTPM`: inconsistent.** It has a strong length trend (M5 ≈ 0.96): the deposited counts scale with transcript length and scaledTPM counts do not.
- **`"no"`: closest, but not consistent.**
  - M2 exceeds 0.02 in every sample (0.029–0.037).
  - |M5| > 0.02 in CD_1 (−0.024).
  - **M1 is outside 1 ± 0.02 in every sample** (+4.5% to +7.2%).
- **`lengthScaledTPM`: not consistent.** M2 is about 2× that of `"no"` (1.9–2.5×) and M3/M4 are lower. Part of its excess deviation is structural, because its average-length scaling depends on the samples imported together (5 here vs 19 originally).
- **Formal outcome: "not identifiable".** High correlation (≥ 0.99) was *not* taken as recovery.

### 4.3 Structural fingerprint (M6)

- **The gene sets are identical.** Our default-index pipeline drops exactly the same 1,061 GENCODE M35 symbols as the deposited matrix, with 0 differences in either direction.
- 1,582 duplicate-sequence transcripts were removed by Salmon's default index.
- This shows the authors used the GENCODE M35 transcript FASTA with Salmon's default duplicate removal and a symbol-level tx2gene, as we did.

### 4.4 Where our counts differ from the deposited ones

Shown for CD_1; the other samples are similar.

- Our total is 4.5% higher (0.69 M fragments). The excess is concentrated:
  - **mitochondrial genes**: 0.29 M (mt-Atp6 +243%, mt-Atp8 +648%, mt-Co2 +144%, mt-Co3 +213%, mt-Nd4l +206%);
  - a few genes with very large excesses (e.g. Gm28439, +32,000 counts where the deposited matrix has ~1);
  - the top 200 genes account for 0.56 M.
- The median gene is +1.8% (log2 ratio +0.026). 81% of genes are higher in ours.
- Processed pseudogenes are *lower* in ours (mean log2 ratio −1.1), lncRNAs higher (+0.17).

This pattern means the authors' Salmon run assigned ambiguous reads differently: mitochondrial overlapping and NUMT-like genes, pseudogene vs parent gene. The most likely, **untested**, explanation is a **decoy-aware index** (genome decoys), which the paper does not mention. Decoys change mapping and assignment without changing the gene set. Testing it would need a gentrome index for the mouse genome (≈ 20 GB RAM), beyond the 16 GB available to WSL here.

### 4.5 Does the missing length offset matter?

Effective transcript length per gene (from our Salmon runs, genes with ≥ 20 counts; `pilot_effective_length_variation.csv`):

- **Variation between samples:**
  - median per-gene SD of log2 effective length = 0.069 (about 5%);
  - 38% of genes vary by > 0.1 log2;
  - 1% vary by > 0.55 log2.
- **PE vs vehicle** mean difference:
  - 12% of genes differ by > 0.2 log2;
  - 29% differ by > 0.1 log2.

  This is comparable to the effect sizes of interest, so the offset is **not negligible** at the gene level.
- **However,** the PE-vs-vehicle length shift is **uncorrelated** with the provisional PE_in_CD log2FC (Spearman 0.015), and PE-significant genes do not have larger shifts (median 0.058 vs 0.068 for all genes).
- So the missing offset adds gene-specific noise and potential bias, but there is no sign that it creates the provisional PE signal. Caveat: 5 samples; effective lengths include isoform-estimation noise.

### 4.6 Attempted approximate offset from pilot lengths: invalid (documented negative result)

- **What we tried:** we attached group-averaged pilot effective lengths (CD+Veh from 1 mouse, CD+PE from 2, MD+Veh from 1, MD+PE from 1) as a tximport-style `avgTxLength` offset to all 19 deposited samples (`analysis/deseq2/run_deseq2_factorial.R --length_offset=`).
- **What happened:** significant genes jumped (PE_in_CD 858 → 2,180; interaction 87 → 1,396). The change in every gene's log2FC equalled minus the injected group length difference (correlation 1.000). Newly significant genes were those whose single-mouse length estimate happened to differ by about 0.3 log2.
- **Conclusion:** with 1–2 pilot mice per group, a group-level offset turns one mouse's isoform-estimation noise into an artificial group effect. **This approach is not a valid approximation**, and its results (`results/deseq2/provisional_deposited__pilot_length_offset/`) must not be interpreted. A valid offset needs per-sample lengths for all 19 mice.

### 4.7 Are effective-length differences systematic by group, or per-mouse noise?

From the two CD+PE pilot mice and the single vehicle mouse of each diet (10,531 genes with ≥ 20 counts):

| Check | Observed | Expected if lengths vary only randomly per mouse |
|---|---|---|
| Correlation of (CD_MP_3 − CD_1) with (CD_MP_5 − CD_1) | 0.501 | 0.50 (from the shared CD_1); higher if PE shifts lengths consistently |
| Correlation of the PE − Veh length shift in CD with the same shift in MD | 0.035 | 0 |
| Robust SD within the CD+PE pair (per-mouse SD) | 0.069 (0.049) | – |
| Robust SD of (PE mean − Veh) in CD | 0.091 | 0.060 |

**Reading:**

- There is no sign of a PE-associated, consistent shift in transcript effective lengths. The two PE mice do not share a common deviation, and the PE-vs-vehicle shifts do not replicate across diets.
- The between-group spread is somewhat larger than pure per-mouse noise predicts. The likeliest source is the single vehicle mouse (CD_1), which cannot be resolved with n = 1.
- Implication: the missing offset mainly adds gene-specific noise to the 19-sample analysis, not a systematic PE bias. This is consistent with the length shifts being uncorrelated with the provisional PE effects (§4.5).

## 5. Consequences (per the pre-specified rules)

- **Outcome "not identifiable" → full re-quantification (19 samples) is required** for confirmatory DESeq2 using `DESeqDataSetFromTximport` with the length offset.
- The re-quantified data will be **kept separate** from the deposited matrix (`D:\mpliver_work\requant\`, results labelled `requant_*`).
- The deposited-matrix results stay labelled PROVISIONAL and become a comparison/sensitivity layer.
- **Download manifest** for the 14 remaining runs (ENA URLs, sizes, MD5): `pipelines/quantification/full_requant_manifest.tsv` (33.2 GB).
### 5.1 Team decision (2026-11-02): proceed with the deposited matrix under stated assumptions

The team decided that the pilot serves as the verification step and that the analysis proceeds on the **deposited matrix**. Full re-quantification of the remaining 14 runs is **not done**, for lack of time. This departs from the pre-specified consequence above and is recorded here and in `reports/phase2/analysis_plan.md`.

The working assumptions, each with the pilot evidence behind it:

| # | Assumption | Evidence | Remaining risk |
|---|---|---|---|
| A1 | The deposited values are Salmon/tximport estimated counts (`"no"`-type, rounded), not TPM-scaled | `scaledTPM` excluded (length slope ≈ 0.96); `"no"` fits best (§4.2); identical gene set (§4.3) | Exact setting formally not identifiable |
| A2 | The missing length offset adds gene-level noise but no systematic PE or group bias | Length shifts uncorrelated with PE effects (ρ = 0.015); no shared PE length shift (§4.7) | Based on 5 mice; within-group length variation of individual genes still inflates their variance |
| A3 | Differences between the authors' and our Salmon runs (totals +4.5–7.2%; mitochondrial and pseudogene assignment) do not affect protein-coding conclusions | The median protein-coding gene differs by about +2%, similar across groups | **Mitochondrially encoded genes and processed pseudogenes are flagged as less reliable**, including the PE_in_CD enrichment of mt- genes |

Results from the deposited matrix are therefore labelled **"deposited counts, pilot-verified under assumptions A1–A3"**, not "provisional". Any claim that hinges on mitochondrial or pseudogene expression, or on single genes with large within-group length variability, should carry this caveat. Full re-quantification remains the way to remove these assumptions if time allows. The manifest and pipeline are ready.

### Earlier pre-specified outcome table

| Outcome | Consequence |
|---|---|
| `"no"` identified, **or** no setting identifiable | The deposited matrix lacks its required offset. Full re-quantification of 19 runs (~31 GB of `.sra`) is needed for confirmatory DESeq2 (`DESeqDataSetFromTximport`). All current results stay provisional. |
| `"scaledTPM"` or `"lengthScaledTPM"` identified | The deposited matrix is a valid bias-corrected count matrix. Current DESeq2 results can be promoted from provisional; re-quantification is an optional sensitivity check. |

Either way, the pilot also provides the **real Salmon mapping rates and fastp statistics** for 4 of the 19 samples. These are needed to interpret the PE-associated library-quality proxy.
