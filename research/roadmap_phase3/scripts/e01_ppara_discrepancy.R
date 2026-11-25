# Why does Jung 2026 report Ppara log2FC +0.53 (CD+PE vs CD+Veh, padj 0.011) when our joint model gives +0.16?
# Same deposited counts. We vary one analysis choice at a time.
source("/mnt/c/Users/Soham/Desktop/microplastic/research/roadmap_phase3/scripts/common.R")
OUT <- file.path(P3, "results/ppara_discrepancy"); dir.create(OUT, showWarnings = FALSE)
cnt <- rd(file.path(REPO, "data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz"))
M <- as.matrix(cnt[, -1]); rownames(M) <- cnt[[1]]; storage.mode(M) <- "integer"
sm <- fread(file.path(REPO, "metadata/sample_metadata.csv"))
sm <- sm[match(colnames(M), sm$counts_column)]
grp <- factor(sm$group, levels = c("CD+Veh", "CD+PE", "MD+Veh", "MD+PE"))
gene <- "Ppara"; rows <- list()
add <- function(label, dds, name = NULL, contrast = NULL, note = "") {
  r <- if (is.null(contrast)) results(dds, name = name) else results(dds, contrast = contrast)
  rows[[length(rows) + 1]] <<- data.table(approach = label, log2FC = r[gene, "log2FoldChange"], lfcSE = r[gene, "lfcSE"],
                                          pvalue = r[gene, "pvalue"], padj = r[gene, "padj"], note = note)
}
cd <- grp %in% c("CD+Veh", "CD+PE")
# 1 joint 4-group model, all genes (our factorial is equivalent for this contrast), prefilter as ours
keep <- rowSums(M >= 10) >= 4
d1 <- DESeq(DESeqDataSetFromMatrix(M[keep, ], data.frame(g = grp), ~ g), quiet = TRUE)
add("A joint 4-group fit, 19 samples, our prefilter (>=10 in >=4)", d1, contrast = c("g", "CD+PE", "CD+Veh"))
# 2 pairwise CD-only fit, no prefilter (DESeq2 default)
d2 <- DESeq(DESeqDataSetFromMatrix(M[, cd], data.frame(g = droplevels(grp[cd])), ~ g), quiet = TRUE)
add("B pairwise CD-only fit (9 samples), no prefilter", d2, name = "g_CD.PE_vs_CD.Veh")
# 3 pairwise with our prefilter
d3 <- DESeq(DESeqDataSetFromMatrix(M[keep, cd], data.frame(g = droplevels(grp[cd])), ~ g), quiet = TRUE)
add("C pairwise CD-only, our prefilter", d3, name = "g_CD.PE_vs_CD.Veh")
# 4 pairwise, size factors taken from the joint 19-sample fit
d4 <- DESeqDataSetFromMatrix(M[keep, cd], data.frame(g = droplevels(grp[cd])), ~ g)
sizeFactors(d4) <- sizeFactors(d1)[cd]; d4 <- DESeq(d4, quiet = TRUE)
add("D pairwise CD-only, joint-fit size factors", d4, name = "g_CD.PE_vs_CD.Veh")
# 5 joint fit, dispersions from joint but no Cook's replacement
d5 <- DESeq(DESeqDataSetFromMatrix(M[keep, ], data.frame(g = grp), ~ g), quiet = TRUE, minReplicatesForReplace = Inf)
add("E joint fit, no outlier replacement", d5, contrast = c("g", "CD+PE", "CD+Veh"))
# 6 pairwise apeglm / ashr shrinkage
s2 <- lfcShrink(d2, coef = "g_CD.PE_vs_CD.Veh", type = "ashr", quiet = TRUE)
rows[[length(rows) + 1]] <- data.table(approach = "F pairwise CD-only, ashr-shrunk", log2FC = s2[gene, "log2FoldChange"],
                                       lfcSE = s2[gene, "lfcSE"], pvalue = NA, padj = NA, note = "")
# 7 existing factorial result (published Phase 2)
ex <- rd(file.path(REPO, "results/deseq2/provisional_deposited/results_PE_in_CD.csv.gz"))[gene == "Ppara"]
rows[[length(rows) + 1]] <- data.table(approach = "G existing Phase 2 factorial (~diet*PE), MLE", log2FC = ex$log2FC, lfcSE = ex$lfcSE,
                                       pvalue = ex$pvalue, padj = ex$padj, note = "results_PE_in_CD.csv.gz")
# 8 MD contrast pairwise (Jung: +0.34, padj 0.23)
md <- grp %in% c("MD+Veh", "MD+PE")
d8 <- DESeq(DESeqDataSetFromMatrix(M[, md], data.frame(g = droplevels(grp[md])), ~ g), quiet = TRUE)
add("H pairwise MD-only (Jung reports +0.34)", d8, name = "g_MD.PE_vs_MD.Veh")
res <- rbindlist(rows)
fwrite(res, file.path(OUT, "ppara_estimates_by_approach.csv"))
# per-sample normalised counts under the two size-factor sets
pn <- data.table(sample = colnames(M), group = as.character(grp), raw = M[gene, ],
                 sf_joint = sizeFactors(d1), norm_joint = counts(d1, normalized = TRUE)[gene, ])
pn[cd, sf_pairwise := sizeFactors(d2)]; pn[cd, norm_pairwise := counts(d2, normalized = TRUE)[gene, ]]
fwrite(pn, file.path(OUT, "ppara_per_sample.csv"))
print(res); print(pn)
