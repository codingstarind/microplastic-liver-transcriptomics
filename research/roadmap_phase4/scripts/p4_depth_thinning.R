# Negative/stress control N2: depth-only effect vectors.
# For each diet, every vehicle library is binomially thinned to (median PE total / median vehicle total) of that diet.
# Paired DESeq2 ~ mouse + thinned (original vs thinned copy of the same mouse) -> log2FC attributable only to depth.
# Uses GSE294180 deposited counts read-only. Seed 11. Writes to research/roadmap_phase4/results/depth_thinning/.
source("/mnt/c/Users/Soham/Desktop/microplastic/research/roadmap_phase3/scripts/common.R")
P4 <- file.path(REPO, "research/roadmap_phase4"); OUT <- file.path(P4, "results/depth_thinning")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(11)
cnt <- rd(file.path(REPO, "data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz"))
M <- as.matrix(cnt[, -1]); rownames(M) <- cnt[[1]]; storage.mode(M) <- "integer"
sm <- fread(file.path(REPO, "metadata/sample_metadata.csv")); sm <- sm[match(colnames(M), sm$counts_column)]
filt <- fread(file.path(REPO, "data/processed/gene_filter.csv"))[passes_filter == TRUE, gene]
info <- list()
for (d in c("CD", "MD")) {
  veh <- sm$counts_column[sm$diet == d & sm$microplastic == "Veh"]
  pe <- sm$counts_column[sm$diet == d & sm$microplastic == "PE"]
  p <- median(colSums(M[, pe])) / median(colSums(M[, veh]))
  thin <- apply(M[, veh], 2, function(x) rbinom(length(x), x, p))
  rownames(thin) <- rownames(M)
  X <- cbind(M[filt, veh], thin[filt, ]); colnames(X) <- c(paste0(veh, "_orig"), paste0(veh, "_thin"))
  cd <- data.frame(mouse = factor(rep(veh, 2)), thinned = factor(rep(c("orig", "thin"), each = length(veh)), levels = c("orig", "thin")))
  dds <- DESeq(DESeqDataSetFromMatrix(X, cd, ~ mouse + thinned), quiet = TRUE)
  r <- results(dds, name = "thinned_thin_vs_orig", independentFiltering = FALSE)
  fwrite(data.table(gene = rownames(r), log2FC = r$log2FoldChange, lfcSE = r$lfcSE, pvalue = r$pvalue, padj = r$padj),
         file.path(OUT, paste0("DEPTH_", d, ".csv.gz")))
  info[[d]] <- list(thinning_p = p, n_vehicle = length(veh), n_padj05 = sum(r$padj < .05, na.rm = TRUE),
                    sd_log2FC = sd(r$log2FoldChange, na.rm = TRUE))
}
write_json(info, file.path(OUT, "run_info.json"), auto_unbox = TRUE, pretty = TRUE)
print(info)
