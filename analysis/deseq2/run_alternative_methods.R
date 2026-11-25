# Robustness: same factorial design and contrasts with edgeR quasi-likelihood and limma-voom.
# Agreement with DESeq2 shows software dependence only; it is not independent validation.
#
# Usage: Rscript analysis/deseq2/run_alternative_methods.R <repo_dir> <label>   (deposited counts)

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
repo <- args[1]
label <- args[2]
out_dir <- file.path(repo, "results", "deseq2", label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(file.path(repo, "metadata", "sample_metadata.csv"))
rownames(meta) <- meta$counts_column
meta$diet <- factor(meta$diet, levels = c("CD", "MD"))
meta$microplastic <- factor(meta$microplastic, levels = c("Veh", "PE"))
counts <- read.csv(gzfile(file.path(repo, "data", "raw",
                                    "GSE294180_MP_raw_counts_gene_level.csv.gz")),
                   row.names = 1, check.names = FALSE)[, rownames(meta)]
counts <- as.matrix(counts[rowSums(counts >= 10) >= 4, ])  # same pre-filter as DESeq2 run

design <- model.matrix(~ diet * microplastic, data = meta)
stopifnot(identical(colnames(design),
                    c("(Intercept)", "dietMD", "microplasticPE", "dietMD:microplasticPE")))
L <- cbind(PE_in_CD = c(0, 0, 1, 0), PE_in_MD = c(0, 0, 1, 1),
           Diet_in_Veh = c(0, 1, 0, 0), Interaction = c(0, 0, 0, 1))

y <- DGEList(counts)
y <- normLibSizes(y)                       # TMM
y <- estimateDisp(y, design)
qfit <- glmQLFit(y, design, robust = TRUE)
v <- voom(y, design)
lfit <- eBayes(contrasts.fit(lmFit(v, design), L))

for (k in colnames(L)) {
  qt <- topTags(glmQLFTest(qfit, contrast = L[, k]), n = Inf, sort.by = "none")$table
  lt <- topTable(lfit, coef = k, number = Inf, sort.by = "none")
  write.csv(data.frame(gene = rownames(qt), logFC = qt$logFC, F = qt$F, pvalue = qt$PValue,
                       padj = qt$FDR),
            gzfile(file.path(out_dir, paste0("edgeR_QL_", k, ".csv.gz"))), row.names = FALSE)
  write.csv(data.frame(gene = rownames(lt), logFC = lt$logFC, t = lt$t, pvalue = lt$P.Value,
                       padj = lt$adj.P.Val),
            gzfile(file.path(out_dir, paste0("limma_voom_", k, ".csv.gz"))), row.names = FALSE)
}
cat("edgeR", as.character(packageVersion("edgeR")), "limma", as.character(packageVersion("limma")), "\n")
