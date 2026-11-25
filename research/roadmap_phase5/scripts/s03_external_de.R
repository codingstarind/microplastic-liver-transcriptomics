# Step 3: external plastic-vs-control effect estimates + camera pathway tests (one study per call).
# Usage: Rscript s03_external_de.R <HONG|GSE245069|GSE277772>
# Counts: DESeq2 Wald, MLE log2FC + lfcSE (same quantity as the GSE294180 training target), independentFiltering off,
#         prefilter >= 10 counts in >= n_min samples. GSE245069 estimated counts rounded ONLY here (DESeq2 needs integers).
# GSE277772: normalised values only -> limma-trend on log2(x+1); descriptive.
# Pathways: limma camera (inter-gene correlation adjusted) on voom/log values, fixed Phase 3 gene sets.
suppressPackageStartupMessages({library(DESeq2); library(limma); library(edgeR); library(data.table)})
args <- commandArgs(TRUE); st <- args[1]
P5 <- "/mnt/c/Users/Soham/Desktop/microplastic/research/roadmap_phase5"
WORK <- "/mnt/d/mpliver_work/phase5"
out <- file.path(P5, "results", "external_de"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) { x <- fread(cmd = paste("zcat", shQuote(f)), header = TRUE); m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]; m }
X <- rd(file.path(WORK, paste0(st, "_expression.csv.gz")))
man <- fread(file.path(P5, "results", "external_sample_manifest.csv"))[study == st]
X <- X[, man$mouse_id]
grp <- gsub("[^A-Za-z0-9]", "", man$group)          # CON, PE, PS, PP, PB, PS80nm, PS05um, PS5um
grp <- relevel(factor(grp), ref = "CON")
gs <-fread(file.path(dirname(P5), "roadmap_phase3", "results", "gene_sets.csv"))
sets <- split(gs$gene, gs$set)
set.seed(20261122)
save_res <- function(res, name, extra = list()) {
  d <- data.table(gene = rownames(res), baseMean = res$baseMean, log2FC = res$log2FoldChange, lfcSE = res$lfcSE,
                  stat = res$stat, pvalue = res$pvalue, padj = res$padj)
  fwrite(d, file.path(out, paste0(st, "__", name, ".csv.gz")))
}
cam_rows <- list()
run_camera <- function(y, design, contrast, name, weights = NULL) {
  idx <- ids2indices(sets, rownames(y))
  cm <- camera(y, idx, design, contrast = contrast, weights = weights, inter.gene.cor = 0.01)
  cam_rows[[name]] <<- data.table(study = st, contrast = name, set = rownames(cm), n_genes = cm$NGenes,
                                  camera_direction = cm$Direction, camera_p = cm$PValue, camera_fdr = cm$FDR)
}
info <- list(study = st, n = ncol(X), groups = as.list(table(grp)))
if (st %in% c("HONG", "GSE245069")) {
  C <- round(X); storage.mode(C) <- "integer"
  nmin <- min(table(grp))
  keep <- rowSums(C >= 10) >= nmin
  info$genes_after_prefilter <- sum(keep); info$noninteger_cells_rounded <- sum(abs(X - round(X)) > 1e-9)
  dds <- DESeqDataSetFromMatrix(C[keep, ], data.frame(group = grp, row.names = colnames(C)), ~ group)
  dds <- DESeq(dds, quiet = TRUE)
  rn <- resultsNames(dds); info$resultsNames <- rn
  others <- setdiff(levels(grp), "CON")
  for (g in others) save_res(results(dds, contrast = c("group", g, "CON"), independentFiltering = FALSE), paste0(g, "_vs_CON"))
  # voom log-CPM for camera, same design
  dge <- calcNormFactors(DGEList(C[keep, ]))
  des <- model.matrix(~ grp); colnames(des) <- c("Intercept", paste0(levels(grp)[-1]))
  v <- voom(dge, des)
  for (g in others) { cv <- setNames(rep(0, ncol(des)), colnames(des)); cv[g] <- 1; run_camera(v, des, cv, paste0(g, "_vs_CON")) }
  if (st == "GSE245069") {
    num <- setNames(rep(0, length(rn)), rn); num[grep("^group_PS", rn)] <- 1 / 3
    save_res(results(dds, contrast = num, independentFiltering = FALSE), "PSpooled_vs_CON")
    cv <- setNames(rep(0, ncol(des)), colnames(des)); cv[others] <- 1 / 3
    run_camera(v, des, cv, "PSpooled_vs_CON")
    # sensitivity: limma-voom on UNROUNDED estimated counts
    dge2 <- calcNormFactors(DGEList(X[keep, ])); v2 <- voom(dge2, des); fit <- eBayes(contrasts.fit(lmFit(v2, des), cbind(pooled = cv)))
    tt <- topTable(fit, number = Inf, sort.by = "none")
    fwrite(data.table(gene = rownames(tt), log2FC = tt$logFC, lfcSE = tt$logFC / tt$t, pvalue = tt$P.Value, padj = tt$adj.P.Val),
           file.path(out, paste0(st, "__PSpooled_vs_CON__limmavoom_unrounded.csv.gz")))
  }
  if (st == "HONG") {
    # sensitivity: PE vs CON fitted on its own 6 mice
    s6 <- grp %in% c("CON", "PE"); g6 <- droplevels(grp[s6])
    k6 <- rowSums(C[, s6] >= 10) >= 3
    d6 <- DESeq(DESeqDataSetFromMatrix(C[k6, s6], data.frame(group = g6, row.names = colnames(C)[s6]), ~ group), quiet = TRUE)
    save_res(results(d6, contrast = c("group", "PE", "CON"), independentFiltering = FALSE), "PE_vs_CON__6mice_only")
  }
} else {
  L <- log2(X + 1)
  keep <- rowSums(X >= 1) >= 3
  info$genes_after_prefilter <- sum(keep)
  des <- model.matrix(~ grp); colnames(des) <- c("Intercept", "PS")
  fit <- eBayes(lmFit(L[keep, ], des), trend = TRUE)
  tt <- topTable(fit, coef = "PS", number = Inf, sort.by = "none")
  fwrite(data.table(gene = rownames(tt), baseMean = rowMeans(X[keep, ][rownames(tt), ]), log2FC = tt$logFC, lfcSE = tt$logFC / tt$t,
                    stat = tt$t, pvalue = tt$P.Value, padj = tt$adj.P.Val), file.path(out, paste0(st, "__PS_vs_CON.csv.gz")))
  run_camera(L[keep, ], des, c(0, 1), "PS_vs_CON")
}
fwrite(rbindlist(cam_rows), file.path(out, paste0(st, "__camera.csv")))
writeLines(jsonlite::toJSON(info, auto_unbox = TRUE, pretty = TRUE), file.path(out, paste0(st, "__run_info.json")))
cat(st, "done\n")
