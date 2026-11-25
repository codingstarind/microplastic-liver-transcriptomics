# Stage 2 (diagnostic, pre-lock): uncertainty of the Hong PE-vs-CON external target under mouse resampling.
# Same estimator as Phase 5 (DESeq2 Wald MLE, 5-group model ~group, prefilter >=10 counts in >= smallest group).
# Fits: full; leave-one-mouse-out for each of the 6 PE/CON mice (other polymers kept); all 9 two-vs-two subsets.
suppressPackageStartupMessages({library(DESeq2); library(data.table); library(BiocParallel)})
P5 <- "/mnt/c/Users/Soham/Desktop/microplastic/research/roadmap_phase5"
PB <- "/mnt/c/Users/Soham/Desktop/microplastic/research/roadmap_phase5b"
out <- file.path(PB, "results", "hong_resampling"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
x <- fread(cmd = "zcat /mnt/d/mpliver_work/phase5/HONG_expression.csv.gz", header = TRUE)
C <- as.matrix(x[, -1]); rownames(C) <- x[[1]]
man <- fread(file.path(P5, "results", "external_sample_manifest.csv"))[study == "HONG"]
C <- round(C[, man$mouse_id]); storage.mode(C) <- "integer"
grp <- relevel(factor(man$group), ref = "CON"); names(grp) <- man$mouse_id
fit_one <- function(keep_mice, tag) {
  g <- droplevels(grp[keep_mice]); Ck <- C[, keep_mice]
  keep <- rowSums(Ck >= 10) >= min(table(g))
  dds <- DESeq(DESeqDataSetFromMatrix(Ck[keep, ], data.frame(group = g, row.names = keep_mice), ~ group), quiet = TRUE)
  r <- results(dds, contrast = c("group", "PE", "CON"), independentFiltering = FALSE)
  data.table(gene = rownames(r), fit = tag, log2FC = r$log2FoldChange, lfcSE = r$lfcSE)
}
pe <- man$mouse_id[man$group == "PE"]; co <- man$mouse_id[man$group == "CON"]; all <- man$mouse_id
jobs <- list(list(all, "full"))
for (m in c(pe, co)) jobs[[length(jobs) + 1]] <- list(setdiff(all, m), paste0("LOO_", m))
for (dp in pe) for (dc in co) jobs[[length(jobs) + 1]] <- list(setdiff(all, c(dp, dc)), paste0("2v2_drop_", dp, "_", dc))
res <- bplapply(jobs, function(j) fit_one(j[[1]], j[[2]]), BPPARAM = MulticoreParam(8))
fwrite(rbindlist(res), file.path(out, "hong_pe_target_resampling.csv.gz"))
cat("fits:", length(jobs), "\n")
