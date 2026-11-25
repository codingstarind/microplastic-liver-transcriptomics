# Phase 5C step 2b: tximport (countsFromAbundance = "no", the GSE294180 setting) + DESeq2 contrasts + camera.
# Same estimand as the GSE294180 training target: DESeq2 Wald MLE log2FC + lfcSE, independentFiltering off,
# prefilter >= 10 counts in >= 4 samples (smallest group). Each 2-group contrast is fitted on its own 8 mice.
#   C1 MPs_vs_PBS         PRIMARY  (PBS 1-4 vs MPs 1-4)
#   C2 EGCGMPs_vs_EGCG    secondary (MP effect on an EGCG background; 8 different mice)
#   C3 MPmain_additive    secondary (16 mice, ~ EGCG + MP; MP coefficient)
#   C4 EGCG_vs_PBS        specificity comparator: a NON-plastic perturbation
# Usage (WSL, mpliver-r): Rscript c02_de.R
suppressPackageStartupMessages({library(tximport); library(DESeq2); library(limma); library(edgeR); library(data.table); library(jsonlite)})
P5C <- "/mnt/c/Users/Soham/Desktop/microplastic/research/roadmap_phase5c"
W <- "/mnt/d/mpliver_work/phase5c"
out <- file.path(P5C, "results", "de"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
man <- fread(file.path(P5C, "results", "z23_sample_manifest.csv"))
tx2 <- read.delim("/mnt/d/mpliver_work/reference/tx2gene_symbol.tsv")
files <- setNames(file.path(W, "requant", man$label, "salmon", "quant.sf"), man$mouse_id)
txi <- suppressMessages(tximport(files, type = "salmon", tx2gene = tx2, countsFromAbundance = "no"))
C <- round(txi$counts); storage.mode(C) <- "integer"
fwrite(data.table(gene = rownames(C), C), file.path(W, "Z23_gene_counts.csv.gz"))
fwrite(data.table(gene = rownames(txi$length), txi$length), file.path(W, "Z23_gene_efflength.csv.gz"))
gs <- fread("/mnt/c/Users/Soham/Desktop/microplastic/research/roadmap_phase3/results/gene_sets.csv"); sets <- split(gs$gene, gs$set)
info <- list(n_mice = ncol(C), genes = nrow(C), tximport = "countsFromAbundance=no, rounded")
cam <- list()
save_res <- function(res, name) fwrite(data.table(gene = rownames(res), baseMean = res$baseMean, log2FC = res$log2FoldChange, lfcSE = res$lfcSE,
                                                  stat = res$stat, pvalue = res$pvalue, padj = res$padj), file.path(out, paste0("Z23__", name, ".csv.gz")))
two_group <- function(ctrl, expo, name) {
  m <- man[group %in% c(ctrl, expo)]
  if (sum(m$group == ctrl) < 3 || sum(m$group == expo) < 3) { info[[name]] <<- "SKIPPED: mice not available"; return(invisible()) }
  g <- relevel(factor(m$group), ref = ctrl); Cm <- C[, m$mouse_id]
  keep <- rowSums(Cm >= 10) >= min(table(g))
  dds <- DESeq(DESeqDataSetFromMatrix(Cm[keep, ], data.frame(g = g, row.names = m$mouse_id), ~ g), quiet = TRUE)
  save_res(results(dds, contrast = c("g", expo, ctrl), independentFiltering = FALSE), name)
  des <- model.matrix(~ g); v <- voom(calcNormFactors(DGEList(Cm[keep, ])), des)
  cm <- camera(v, ids2indices(sets, rownames(v)), des, contrast = 2, inter.gene.cor = 0.01)
  cam[[name]] <<- data.table(contrast = name, set = rownames(cm), n_genes = cm$NGenes, direction = cm$Direction, p = cm$PValue, fdr = cm$FDR)
  info[[name]] <<- list(n_ctrl = sum(g == ctrl), n_exp = sum(g == expo), genes_after_prefilter = sum(keep), n_padj05 = sum(results(dds, contrast = c("g", expo, ctrl), independentFiltering = FALSE)$padj < 0.05, na.rm = TRUE))
}
set.seed(20261122)
two_group("PBS", "MPs", "MPs_vs_PBS")
two_group("EGCG", "EGCGMPs", "EGCGMPs_vs_EGCG")
two_group("PBS", "EGCG", "EGCG_vs_PBS")
if (nrow(man) == 16) {
  d <- data.frame(EGCG = factor(man$egcg), MP = factor(man$mp_exposed), row.names = man$mouse_id)
  keep <- rowSums(C >= 10) >= 4
  dds <- DESeq(DESeqDataSetFromMatrix(C[keep, man$mouse_id], d, ~ EGCG + MP), quiet = TRUE)
  save_res(results(dds, name = "MP_TRUE_vs_FALSE", independentFiltering = FALSE), "MPmain_additive")
  des <- model.matrix(~ EGCG + MP, d); v <- voom(calcNormFactors(DGEList(C[keep, man$mouse_id])), des)
  cm <- camera(v, ids2indices(sets, rownames(v)), des, contrast = 3, inter.gene.cor = 0.01)
  cam[["MPmain_additive"]] <- data.table(contrast = "MPmain_additive", set = rownames(cm), n_genes = cm$NGenes, direction = cm$Direction, p = cm$PValue, fdr = cm$FDR)
  info$MPmain_additive <- list(genes_after_prefilter = sum(keep))
} else info$MPmain_additive <- "SKIPPED: needs all 16 mice"
fwrite(rbindlist(cam), file.path(out, "Z23__camera.csv"))
writeLines(toJSON(info, auto_unbox = TRUE, pretty = TRUE), file.path(out, "Z23__run_info.json"))
cat("c02_de done\n")
