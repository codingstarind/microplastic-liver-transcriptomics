# Data dictionary

Every data object in the project, what each field means, and whether it is
**original** (as deposited by the authors) or **derived** (computed by us).
Scripts that create each object are given in brackets.

## 1. Original data (never modified)

### `data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz` (GEO supplementary file)

| Field | Type | Meaning |
|---|---|---|
| first column (unnamed in the file) | text | Gene identifier. **This is a gene symbol** from GENCODE M35 (e.g. `Alb`, `0610005C13Rik`). 477 rows use an unversioned Ensembl gene ID instead (`ENSMUSG…`) because GENCODE gives those genes no other name. Loaded as the index `gene`. |
| `CD_1_S39` … `MD_MP_5_S57` (19 columns) | integer | Gene-level read counts for one liver sample. The name encodes diet (`CD`/`MD`), microplastic (`_MP` = PE gavage), replicate number and sequencing library (`S39`–`S57`). Values are whole numbers ≥ 0. They were produced by Salmon and tximport and are almost certainly *rounded* estimated counts (see the audit report). |

Shape: 55,836 rows (genes) × 20 columns (1 identifier + 19 samples). The order of
the columns in the file is CD+PE, CD+Veh, MD+PE, MD+Veh. Every analysis reorders
them to match `sample_metadata.csv`.

### `data/raw/GSE294180_family.soft.gz` (GEO SOFT)
Text records for the series, the platform and the 19 GSM samples. It is the source of `sample_metadata.csv`.

### `data/raw/gencode.vM35.annotation.gtf.gz` (GENCODE)
Reference gene annotation, the same release the authors used for Salmon. It is used only to attach biotype and chromosome to each symbol.

## 2. Metadata (curated, committed)

### `metadata/sample_metadata.csv` [02]
One row per sample. Group labels come from the GEO `treatment` characteristic, not from the sample names.

| Column | Meaning | Source |
|---|---|---|
| `geo_accession` | GSM accession | GEO |
| `sample_title` | GEO title, e.g. `Liver_MD_MP_3` | GEO |
| `counts_column` | Matching column name in the count file | GEO `description` field |
| `group` | `CD+Veh`, `CD+PE`, `MD+Veh` or `MD+PE` | derived from `geo_treatment` |
| `diet` | `CD` = control diet (D09100304), `MD` = MASH diet (D09100310) | derived from `geo_treatment` |
| `microplastic` | `Veh` = vehicle gavage, `PE` = 2 mg/day polyethylene gavage | derived from `geo_treatment` |
| `replicate` | Replicate number within group (1–5) | derived (order of GSM records) |
| `geo_treatment` | Treatment text as submitted (e.g. `MASH Diet + Microplastics`) | GEO |
| `geo_batch` | Batch label as submitted (`B1` for all samples) | GEO |
| `library_name` | Sequencing library / sample-sheet index (`S39`–`S57`) | GEO |
| `tissue`, `genotype`, `organism`, `molecule`, `library_strategy`, `instrument` | As submitted | GEO |
| `biosample`, `sra_experiment` | NCBI BioSample and SRA experiment IDs (FASTQ files) | GEO |
| `submission_date` | GEO submission date | GEO |
| `group_from_column_name` | Group inferred from `counts_column`; only a cross-check, and identical to `group` for all 19 | derived |

### `metadata/study_design.csv`
Study-level facts (strain, sex, age, diet composition, dose, protocols). Each row cites its source (GEO or a section of the paper) and carries notes on discrepancies.

### `metadata/data_provenance.json` [01]
For each raw file: URL, local path, how it was obtained, UTC download time, size in bytes, MD5 and SHA-256.

### `metadata/gene_annotation_gencode_vM35.csv.gz` [03]
One row per gene in the count file (index `gene`).

| Column | Meaning |
|---|---|
| `gene_name` | Same as the index (the symbol used in the count file) |
| `gene_ids` | GENCODE Ensembl gene ID(s) carrying this symbol, `;`-separated |
| `n_gene_ids` | Number of Ensembl genes sharing the symbol (> 1 for 228 symbols) |
| `gene_type` | GENCODE biotype(s), e.g. `protein_coding`, `lncRNA` |
| `chromosome` | Chromosome(s); `chrM` = mitochondrial genome |
| `start`, `end` | Genomic span on GRCm39 (1-based, smallest start / largest end) |

## 3. Derived matrices (`data/processed/`, not committed, regenerate with [04])

All are genes × 19 samples, with columns in `sample_metadata.csv` order.

| File | Scale | Transformation from the original | Used for |
|---|---|---|---|
| `counts_original.csv.gz` | counts | none (columns reordered only) | QC, library sizes, marker scores |
| `gene_filter.csv` | – | per-gene `mean_count`, `n_samples_ge_min_count`, `passes_filter` | records which genes were kept |
| `counts_filtered.csv.gz` | counts | rows kept if ≥ 10 counts in ≥ 4 samples (13,777 genes) | input to normalisation |
| `normalized_counts_filtered.csv.gz` | counts | each column ÷ its median-of-ratios size factor | mean–variance plot |
| `log2_normalized_filtered.csv.gz` | log2 | log2(normalised count + 1) | EDA, PCA, clustering |
| `log2_cpm_filtered.csv.gz` | log2 | log2(count ÷ library size × 10⁶ + 1) | PCA sensitivity check only |

No matrix has been gene-wise standardised (z-scored); PCA centres genes internally but does not scale them.

## 4. Derived tables (`results/tables/`, committed)

| File | Content |
|---|---|
| `audit_summary.json` | Integrity checks, value-type evidence, ID matching, sex check, duplicate profiles |
| `counts_preview.csv` | First 10 rows of the count file exactly as deposited |
| `per_sample_audit.csv` | Library size, genes detected, chrM fraction, top gene per sample |
| `gene_biotype_summary.csv` | Genes and share of counts by biotype |
| `top_expressed_genes.csv` | 30 genes with the largest mean share of the library |
| `sex_marker_counts.csv` | Xist and Y-chromosome marker counts |
| `filter_summary.json` | Low-expression filter rule and outcome |
| `size_factors.csv` | Library sizes and median-of-ratios size factors |
| `library_summary_by_group.csv`, `gene_level_summary.csv` | Descriptive summaries |
| `sample_correlation_summary.csv` | Mean / minimum correlation of each sample to the others |
| `variance_share_by_factor.csv` | Per-gene share of variance (η²) tied to diet / PE / group, vs shuffled labels |
| `tissue_marker_cpm.csv` | Marker reads per million for 7 tissue / cell-type sets |
| `muscle_tracking_genes.csv` | Genes removed in the muscle-sensitivity PCA |
| `pca_variance_explained.csv` | % variance per PC for each feature set |
| `pca_pc_associations.csv` | η² of each PC with diet / PE / group; Spearman ρ with library size, library order, muscle score |
| `pca_primary_scores.csv`, `pca_primary_loadings.csv` | PC scores per sample and gene loadings |
| `pca_outlier_screen.csv` | Distance of each sample to its group centroid (PC1–4) |
| `clustering_label_agreement.csv` | Adjusted Rand index of tree cuts (k = 2–6) with group / diet / PE |

### Key derived quantities

- **Library size**: sum of all counts in a sample (units: reads assigned to genes).
- **Size factor** (median-of-ratios): for sample *j*, the median over genes of count<sub>ij</sub> / geometric mean<sub>i</sub>. Genes with a zero in any sample are skipped. A value of 1.2 means the sample looks 20% "deeper" than a typical sample.
- **CPM**: count / library size × 10⁶.
- **η² (eta squared)**: between-group sum of squares / total sum of squares, from 0 to 1. It is a descriptive effect-size measure; no p-value is attached.
- **Marker CPM**: summed CPM of a predefined marker gene set (see `src/mpliver/qc_markers.py`).
