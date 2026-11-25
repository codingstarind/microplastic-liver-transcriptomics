# Locked interventions, R part.
# I1: Hong count-unit harmonisation (factor 0.5 primary, 0.472 sensitivity), identical Phase 5 DESeq2 pipeline,
#     all 4 polymer contrasts vs the shared CON.
# I2b: GSE294180 PE main effect averaged over diets from the saved 19-mouse dds (numeric contrast 0,0,1,0.5).
suppressPackageStartupMessages({library(DESeq2); library(data.table); library(BiocParallel)})
REPO <- "/mnt/c/Users/Soham/Desktop/microplastic"
P5 <- file.path(REPO, "research/roadmap_phase5"); PB <- file.path(REPO, "research/roadmap_phase5b")
out <- file.path(PB, "results", "intervention_targets"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
x <- fread(cmd = "zcat /mnt/d/mpliver_work/phase5/HONG_expression.csv.gz", header = TRUE)
C0 <- as.matrix(x[, -1]); rownames(C0) <- x[[1]]
man <- fread(file.path(P5, "results", "external_sample_manifest.csv"))[study == "HONG"]
C0 <- C0[, man$mouse_id]
grp <- relevel(factor(man$group), ref = "CON")
run_unit <- function(f) {
  C <- round(C0 * f); storage.mode(C) <- "integer"
  keep <- rowSums(C >= 10) >= 3
  dds <- DESeq(DESeqDataSetFromMatrix(C[keep, ], data.frame(group = grp, row.names = colnames(C)), ~ group), quiet = TRUE)
  rbindlist(lapply(c("PE", "PS", "PP", "PB"), function(g) {
    r <- results(dds, contrast = c("group", g, "CON"), independentFiltering = FALSE)
    data.table(gene = rownames(r), contrast = paste0(g, "_vs_CON"), unit_factor = f, baseMean = r$baseMean, log2FC = r$log2FoldChange,
               lfcSE = r$lfcSE, padj = r$padj)
  }))
}
res <- bplapply(c(1, 0.5, 0.472), run_unit, BPPARAM = MulticoreParam(3))
fwrite(rbindlist(res), file.path(out, "hong_unit_harmonised_targets.csv.gz"))
dds <- readRDS(file.path(REPO, "results/deseq2/provisional_deposited/dds.rds"))
stopifnot(identical(resultsNames(dds), c("Intercept", "diet_MD_vs_CD", "microplastic_PE_vs_Veh", "dietMD.microplasticPE")))
r <- results(dds, contrast = c(0, 0, 1, 0.5), independentFiltering = FALSE)
fwrite(data.table(gene = rownames(r), baseMean = r$baseMean, log2FC = r$log2FoldChange, lfcSE = r$lfcSE),
       file.path(out, "gse294180_PE_avg_over_diets.csv.gz"))
r1 <- results(dds, contrast = c(0, 0, 1, 0), independentFiltering = FALSE)  # reproduces PE_in_CD (sanity)
fwrite(data.table(gene = rownames(r1), log2FC = r1$log2FoldChange), file.path(out, "gse294180_PE_in_CD_recomputed.csv.gz"))
cat("done\n")
