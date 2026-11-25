# DESeq2 factorial model for GSE294180: ~ diet + microplastic + diet:microplastic
#
# Usage (inside WSL, mpliver-r environment):
#   Rscript analysis/deseq2/run_deseq2_factorial.R <repo_dir> <input> <out_label> [options]
#
#   <input>  deposited       deposited GEO integer matrix, no length offset (PROVISIONAL:
#                            the tximport setting is unverified, see the validation report)
#            requant:<dir>   Salmon quant.sf files in <dir>/<sample>/salmon, imported with
#                            tximport and DESeqDataSetFromTximport (length offset included)
#   options  --drop=<sample>       drop one sample (sensitivity analyses only)
#            --covariate=<column>  add a numeric covariate from <repo>/results/deseq2/covariates.csv
#
# Reference levels are set explicitly (diet: CD, microplastic: Veh) and checked.
# Each contrast is its own test family; DESeq2 applies BH-FDR within each results()
# call after independent filtering. Cook's-distance outlier *replacement* is not
# possible with 4-5 replicates per group (it needs >= 7); flagged genes get NA p-values.

suppressPackageStartupMessages({
  library(DESeq2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
repo <- args[1]
input <- args[2]
label <- args[3]
opts <- args[-(1:3)]
drop <- sub("^--drop=", "", grep("^--drop=", opts, value = TRUE))
covariate <- sub("^--covariate=", "", grep("^--covariate=", opts, value = TRUE))
# --length_offset=<csv>: group-level average transcript lengths (columns named by group,
# e.g. from the re-quantification pilot), attached as DESeq2's avgTxLength assay so the
# deposited counts get an approximate tximport length offset.
length_offset <- sub("^--length_offset=", "", grep("^--length_offset=", opts, value = TRUE))

provisional <- input == "deposited"
out_dir <- file.path(repo, "results", "deseq2", label)
fig_dir <- file.path(repo, "figures", "deseq2", label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(file.path(repo, "metadata", "sample_metadata.csv"))
rownames(meta) <- meta$counts_column
meta$diet <- factor(meta$diet, levels = c("CD", "MD"))
meta$microplastic <- factor(meta$microplastic, levels = c("Veh", "PE"))
if (length(drop)) meta <- meta[!startsWith(meta$counts_column, paste0(drop, "_S")), ]

design <- ~ diet + microplastic + diet:microplastic
if (length(covariate)) {
  cov <- read.csv(file.path(repo, "results", "deseq2", "covariates.csv"), row.names = 1)
  meta[[covariate]] <- as.numeric(scale(cov[rownames(meta), covariate]))
  design <- as.formula(paste("~", covariate, "+ diet + microplastic + diet:microplastic"))
}

# ---- Counts -------------------------------------------------------------------
if (provisional) {
  counts <- read.csv(gzfile(file.path(repo, "data", "raw",
                                      "GSE294180_MP_raw_counts_gene_level.csv.gz")),
                     row.names = 1, check.names = FALSE)[, rownames(meta)]
  dds <- DESeqDataSetFromMatrix(as.matrix(counts), meta, design)
  if (length(length_offset)) {
    glen <- read.csv(length_offset, row.names = 1, check.names = FALSE)
    stopifnot(all(rownames(dds) %in% rownames(glen)))
    grp_of <- paste(meta$diet, meta$microplastic, sep = "+")
    lenmat <- as.matrix(glen[rownames(dds), grp_of])
    lenmat[!is.finite(lenmat) | lenmat <= 0] <- 1
    dimnames(lenmat) <- dimnames(dds)
    assay(dds, "avgTxLength") <- lenmat   # used by estimateSizeFactors (tximport-style)
  }
} else {
  suppressPackageStartupMessages(library(tximport))
  qdir <- sub("^requant:", "", input)
  tx2gene <- read.delim(file.path(qdir, "..", "reference", "tx2gene_symbol.tsv"))
  files <- setNames(file.path(qdir, sub("_S[0-9]+$", "", rownames(meta)), "salmon", "quant.sf"),
                    rownames(meta))
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)  # countsFromAbundance = "no"
  dds <- DESeqDataSetFromTximport(txi, meta, design)
}

# Same pre-filter as Phase 1: >= 10 counts in >= 4 samples (smallest group size).
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep, ]
dds <- DESeq(dds, quiet = TRUE)

expected <- c("Intercept", "diet_MD_vs_CD", "microplastic_PE_vs_Veh", "dietMD.microplasticPE")
coef_names <- resultsNames(dds)
if (!all(expected %in% coef_names)) stop("unexpected coefficients: ", paste(coef_names, collapse = ", "))

contrasts <- list(
  PE_in_CD = list(contrast = list("microplastic_PE_vs_Veh"),
                  meaning = "CD+PE vs CD+Veh (PE effect under control diet)"),
  PE_in_MD = list(contrast = list(c("microplastic_PE_vs_Veh", "dietMD.microplasticPE")),
                  meaning = "MD+PE vs MD+Veh (PE effect under MASH diet)"),
  Diet_in_Veh = list(contrast = list("diet_MD_vs_CD"),
                     meaning = "MD+Veh vs CD+Veh (diet effect without PE)"),
  Interaction = list(contrast = list("dietMD.microplasticPE"),
                     meaning = "(MD+PE - MD+Veh) - (CD+PE - CD+Veh), log2 scale")
)

norm <- counts(dds, normalized = TRUE)
grp <- paste(meta[colnames(dds), "diet"], meta[colnames(dds), "microplastic"], sep = "+")
lgm <- sapply(split(seq_len(ncol(norm)), grp), function(i) rowMeans(log2(norm[, i, drop = FALSE] + 0.5)))
naive <- list(PE_in_CD = lgm[, "CD+PE"] - lgm[, "CD+Veh"],
              PE_in_MD = lgm[, "MD+PE"] - lgm[, "MD+Veh"],
              Diet_in_Veh = lgm[, "MD+Veh"] - lgm[, "CD+Veh"],
              Interaction = (lgm[, "MD+PE"] - lgm[, "MD+Veh"]) - (lgm[, "CD+PE"] - lgm[, "CD+Veh"]))

summary_rows <- list()
for (name in names(contrasts)) {
  res <- results(dds, contrast = contrasts[[name]]$contrast, alpha = 0.05)
  shr <- lfcShrink(dds, contrast = contrasts[[name]]$contrast, res = res, type = "ashr",
                   quiet = TRUE)
  tab <- data.frame(gene = rownames(res), baseMean = res$baseMean,
                    log2FC = res$log2FoldChange, lfcSE = res$lfcSE,
                    ci95_low = res$log2FoldChange - qnorm(0.975) * res$lfcSE,
                    ci95_high = res$log2FoldChange + qnorm(0.975) * res$lfcSE,
                    stat = res$stat, pvalue = res$pvalue, padj = res$padj,
                    log2FC_ashr = shr$log2FoldChange, lfcSE_ashr = shr$lfcSE,
                    naive_log2_group_mean_diff = naive[[name]][rownames(res)])
  write.csv(tab, gzfile(file.path(out_dir, paste0("results_", name, ".csv.gz"))), row.names = FALSE)

  well <- tab$baseMean > 100 & !is.na(tab$pvalue)
  sign_agree <- mean(sign(tab$log2FC[well]) == sign(tab$naive_log2_group_mean_diff[well]) |
                       abs(tab$naive_log2_group_mean_diff[well]) < 0.05)
  summary_rows[[name]] <- data.frame(
    contrast = name, meaning = contrasts[[name]]$meaning,
    n_genes_tested = sum(!is.na(res$pvalue)),
    n_outlier_flagged_NA = sum(is.na(res$pvalue) & res$baseMean > 0),
    n_independent_filtered = sum(!is.na(res$pvalue) & is.na(res$padj)),
    n_padj_lt_0.05 = sum(res$padj < 0.05, na.rm = TRUE),
    n_padj_lt_0.05_absLFC_ge_0.38 = sum(res$padj < 0.05 & abs(res$log2FoldChange) >= 0.38, na.rm = TRUE),
    pi0_storey_lambda0.5 = min(1, mean(res$pvalue > 0.5, na.rm = TRUE) / 0.5),
    sign_check_vs_group_means = round(sign_agree, 4),
    corr_with_group_means = round(cor(tab$log2FC[well], tab$naive_log2_group_mean_diff[well]), 4))

  png(file.path(fig_dir, paste0("ma_", name, ".png")), width = 1600, height = 1200, res = 220)
  plotMA(res, ylim = c(-4, 4), alpha = 0.05,
         main = paste0(name, if (provisional) " [PROVISIONAL: deposited counts]" else ""),
         colSig = "#d62728")
  mtext(contrasts[[name]]$meaning, side = 3, line = 0.2, cex = 0.7)
  dev.off()
}
summary_tab <- do.call(rbind, summary_rows)
write.csv(summary_tab, file.path(out_dir, "contrast_summary.csv"), row.names = FALSE)

# ---- Diagnostics -----------------------------------------------------------------
png(file.path(fig_dir, "dispersion_trend.png"), width = 1600, height = 1200, res = 220)
plotDispEsts(dds, main = paste("Mean-dispersion trend", if (provisional) "[PROVISIONAL]" else ""))
dev.off()

png(file.path(fig_dir, "pvalue_histograms.png"), width = 2400, height = 1400, res = 220)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (name in names(contrasts)) {
  tab <- read.csv(gzfile(file.path(out_dir, paste0("results_", name, ".csv.gz"))))
  hist(tab$pvalue[tab$baseMean > 10], breaks = 40, col = "#2a78d6", border = "white",
       main = name, xlab = "raw p-value (baseMean > 10)")
}
dev.off()

cooks <- assays(dds)[["cooks"]]
png(file.path(fig_dir, "cooks_distance_by_sample.png"), width = 2400, height = 1200, res = 220)
par(mar = c(8, 4, 3, 1))
boxplot(log10(cooks + 1e-8), las = 2, outline = FALSE, col = "#d9d9d9",
        ylab = "log10 Cook's distance", main = "Cook's distance per sample (all genes)")
dev.off()
cutoff <- qf(0.99, length(coef_names), ncol(dds) - length(coef_names))
influence <- data.frame(sample = colnames(dds),
                        n_genes_cooks_above_cutoff = colSums(cooks > cutoff, na.rm = TRUE),
                        median_cooks = apply(cooks, 2, median, na.rm = TRUE),
                        size_factor = if (is.null(sizeFactors(dds)))
                          colMeans(normalizationFactors(dds)) else sizeFactors(dds))
write.csv(influence, file.path(out_dir, "sample_influence.csv"), row.names = FALSE)

vsd <- vst(dds, blind = FALSE)
pca <- prcomp(t(assay(vsd)[order(-apply(assay(vsd), 1, var))[1:500], ]))
write.csv(data.frame(sample = colnames(vsd), pca$x[, 1:4]),
          file.path(out_dir, "vst_pca_scores.csv"), row.names = FALSE)
saveRDS(dds, file.path(out_dir, "dds.rds"))
write.csv(data.frame(gene = rownames(vsd), assay(vsd), check.names = FALSE),
          gzfile(file.path(out_dir, "vst_matrix.csv.gz")), row.names = FALSE)

info <- list(
  input = input, provisional = provisional, design = deparse(design),
  dropped_samples = drop, covariate = covariate, length_offset = length_offset,
  n_samples = ncol(dds), n_genes_after_prefilter = nrow(dds),
  coefficient_names = coef_names,
  reference_levels = list(diet = levels(dds$diet)[1], microplastic = levels(dds$microplastic)[1]),
  cooks_cutoff_F99 = cutoff, min_replicates_for_replace = 7,
  versions = list(R = R.version.string, DESeq2 = as.character(packageVersion("DESeq2")),
                  ashr = as.character(packageVersion("ashr"))),
  test_families = "one BH-FDR family per contrast (results() call); independent filtering on"
)
write_json(info, file.path(out_dir, "run_info.json"), auto_unbox = TRUE, pretty = TRUE)
print(summary_tab)
