# Compare pilot re-quantification with the deposited GSE294180 matrix.
#
# Implements the pre-specified metrics M1-M6 of
# reports/quantification_validation/pilot_protocol.md for the three tximport
# countsFromAbundance settings, plus the run-to-run noise floor (CD_1 x2).
#
# Usage (inside WSL, mpliver-r environment):
#   Rscript pipelines/quantification/compare_tximport.R <work_dir> <repo_dir>

suppressPackageStartupMessages({
  library(tximport)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
work <- args[1]
repo <- args[2]
out_dir <- file.path(repo, "results", "quantification_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pilot <- c(CD_1_rep1 = "CD_1_S39", CD_MP_5 = "CD_MP_5_S52",
           MD_1 = "MD_1_S43", MD_MP_4 = "MD_MP_4_S56")
quant_file <- function(label) file.path(work, "requant", label, "salmon", "quant.sf")

tx2gene <- read.delim(file.path(work, "reference", "tx2gene_symbol.tsv"))
deposited <- read.csv(gzfile(file.path(repo, "data", "raw",
                                       "GSE294180_MP_raw_counts_gene_level.csv.gz")),
                      row.names = 1, check.names = FALSE)

import <- function(labels, setting) {
  files <- setNames(vapply(labels, quant_file, ""), labels)
  suppressMessages(tximport(files, type = "salmon", tx2gene = tx2gene,
                            countsFromAbundance = setting))
}

# Metrics comparing one integer vector "ours" with "ref" (named by gene symbol).
compare <- function(ours, ref, eff_len) {
  genes <- intersect(names(ours), names(ref))
  o <- ours[genes]; r <- ref[genes]
  keep <- r >= 20
  lr <- log2((r[keep] + 1) / (o[keep] + 1))
  len <- eff_len[genes][keep]
  ok <- is.finite(log2(len)) & len > 0
  slope <- unname(coef(lm(lr[ok] ~ log2(len[ok])))[2])
  either <- (o > 0) | (r > 0)
  c(M1_total_ratio = sum(o) / sum(r),
    M2_median_abs_log2_ratio = median(abs(lr)),
    M3_frac_within_5pct = mean(abs(o[keep] / r[keep] - 1) <= 0.05),
    M4_frac_exact = mean(o[either] == r[either]),
    M5_length_slope = slope,
    spearman_log = cor(log1p(o[keep]), log1p(r[keep]), method = "spearman"),
    n_genes_compared = sum(keep))
}

settings <- c("no", "scaledTPM", "lengthScaledTPM")
base <- import(names(pilot), "no")          # effective lengths come from this import
rows <- list()
gene_tables <- list()
for (setting in settings) {
  txi <- if (setting == "no") base else import(names(pilot), setting)
  counts <- round(txi$counts)
  for (label in names(pilot)) {
    m <- compare(counts[, label], setNames(deposited[[pilot[label]]], rownames(deposited)),
                 base$length[, label])
    rows[[length(rows) + 1]] <- data.frame(sample = label, setting = setting, t(m))
    if (setting == "no") {
      gene_tables[[label]] <- data.frame(
        gene = rownames(counts), sample = label,
        ours_no = counts[, label],
        effective_length = base$length[, label])
    }
  }
  for (label in names(pilot)) {
    gene_tables[[label]][[paste0("ours_", setting)]] <- counts[gene_tables[[label]]$gene, label]
  }
}

# Noise floor: two independent quantifications of the same run (CD_1).
reps <- import(c("CD_1_rep1", "CD_1_rep2"), "no")
rep_counts <- round(reps$counts)
floor_m <- compare(rep_counts[, "CD_1_rep2"], rep_counts[, "CD_1_rep1"],
                   reps$length[, "CD_1_rep1"])
rows[[length(rows) + 1]] <- data.frame(sample = "CD_1 rep2 vs rep1", setting = "noise floor (no)",
                                       t(floor_m))
metrics <- do.call(rbind, rows)
write.csv(metrics, file.path(out_dir, "pilot_metrics.csv"), row.names = FALSE)

# M6: which gene symbols exist in each matrix?
ours_genes <- rownames(base$counts)
dep_genes <- rownames(deposited)
all_symbols <- unique(tx2gene$GENEID)
missing_dep <- setdiff(all_symbols, dep_genes)
missing_ours <- setdiff(all_symbols, ours_genes)
dups <- file.path(work, "index", "salmon_index_gencode_vM35_k31", "duplicate_clusters.tsv")
n_dup <- if (file.exists(dups)) nrow(read.delim(dups)) else NA
m6 <- list(
  n_symbols_in_annotation = length(all_symbols),
  n_symbols_deposited = length(dep_genes),
  n_symbols_ours = length(ours_genes),
  n_missing_deposited = length(missing_dep),
  n_missing_ours = length(missing_ours),
  n_missing_in_both = length(intersect(missing_dep, missing_ours)),
  n_missing_only_deposited = length(setdiff(missing_dep, missing_ours)),
  n_missing_only_ours = length(setdiff(missing_ours, missing_dep)),
  examples_missing_only_deposited = head(setdiff(missing_dep, missing_ours), 15),
  examples_missing_only_ours = head(setdiff(missing_ours, missing_dep), 15),
  n_duplicate_transcripts_removed_by_salmon_index = n_dup
)
write_json(m6, file.path(out_dir, "pilot_gene_set_overlap.json"), auto_unbox = TRUE,
           pretty = TRUE)

# Salmon / fastp run statistics (the real mapping rates the deposited data lacked).
stats <- do.call(rbind, lapply(c(names(pilot), "CD_1_rep2"), function(label) {
  meta <- fromJSON(file.path(work, "requant", label, "salmon", "aux_info", "meta_info.json"))
  libfmt <- fromJSON(file.path(work, "requant", label, "salmon", "lib_format_counts.json"))
  fp <- fromJSON(file.path(work, "requant", label, "fastp.json"))
  data.frame(sample = label,
             fastp_reads_before = fp$summary$before_filtering$total_reads,
             fastp_reads_after = fp$summary$after_filtering$total_reads,
             salmon_fragments_processed = meta$num_processed,
             salmon_fragments_mapped = meta$num_mapped,
             salmon_percent_mapped = meta$percent_mapped,
             library_type = libfmt$expected_format,
             sum_numreads = sum(read.delim(quant_file(label))$NumReads))
}))
write.csv(stats, file.path(out_dir, "pilot_mapping_stats.csv"), row.names = FALSE)

genes_long <- do.call(rbind, gene_tables)
genes_long$deposited <- deposited[cbind(match(genes_long$gene, rownames(deposited)),
                                        match(pilot[genes_long$sample], colnames(deposited)))]
write.csv(genes_long, gzfile(file.path(work, "requant", "pilot_gene_comparison.csv.gz")),
          row.names = FALSE)

print(metrics, digits = 4)
print(stats)
str(m6)
