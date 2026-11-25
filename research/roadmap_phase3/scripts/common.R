# Shared helpers for Phase 3 R scripts (run inside WSL env mpliver-r).
suppressPackageStartupMessages({
  library(DESeq2); library(limma); library(edgeR); library(data.table); library(jsonlite)
})
REPO <- "/mnt/c/Users/Soham/Desktop/microplastic"
P2 <- file.path(REPO, "research/roadmap_phase2")
P3 <- file.path(REPO, "research/roadmap_phase3")
DROOT <- "/mnt/d/mpliver_work/phase2_refs"
SEED <- 20261122
# fread cannot read .gz without R.utils (not installed); stream through zcat instead
rd <- function(f, ...) if (endsWith(f, ".gz")) fread(cmd = paste("zcat", shQuote(f)), header = TRUE, ...) else fread(f, header = TRUE, ...)
set.seed(SEED)

universe <- function() {
  u <- fread(file.path(P2, "metadata/phase3_gene_universe.csv"))
  list(main = u$vm35_symbol[u$universe_main], sens = u$vm35_symbol[u$universe_sensitivity])
}
eligible <- function() fread(file.path(P2, "metadata/phase3_eligible_samples.csv"),
                             colClasses = c(sample_id = "character"))
ready <- function(acc) {
  x <- rd(file.path(DROOT, acc, "prepared/expression_ready.csv.gz"))
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]; m
}

# DESeq2 Wald results for every gene in `genes` (independent filtering OFF => padj for all)
deseq_table <- function(dds, coef_name = NULL, contrast = NULL, genes = rownames(dds), n_detect = 4) {
  r <- if (!is.null(contrast)) results(dds, contrast = contrast, independentFiltering = FALSE)
       else results(dds, name = coef_name, independentFiltering = FALSE)
  r <- r[genes, ]
  cts <- counts(dds)[genes, ]
  out <- data.table(gene = genes, baseMean = r$baseMean, log2FC = r$log2FoldChange, lfcSE = r$lfcSE,
                    ci95_low = r$log2FoldChange - 1.959964 * r$lfcSE,
                    ci95_high = r$log2FoldChange + 1.959964 * r$lfcSE,
                    stat = r$stat, pvalue = r$pvalue)
  out[, padj := p.adjust(pvalue, "BH")]
  out[, detected := rowSums(cts >= 10) >= n_detect]
  cooks <- assays(dds)[["cooks"]]
  out[, max_cooks := if (!is.null(cooks)) apply(cooks[genes, , drop = FALSE], 1, max) else NA_real_]
  out[]
}

add_ashr <- function(tab, dds, coef_name = NULL, contrast = NULL) {
  s <- if (!is.null(contrast)) lfcShrink(dds, contrast = contrast, type = "ashr", quiet = TRUE)
       else lfcShrink(dds, coef = coef_name, type = "ashr", quiet = TRUE)
  tab[, log2FC_ashr := s[tab$gene, "log2FoldChange"]]
  tab[]
}

diag_plots <- function(dds, tab, file, title) {
  png(file, width = 1100, height = 420)
  par(mfrow = c(1, 3))
  plotDispEsts(dds, main = paste(title, "dispersion"))
  with(tab, plot(log10(baseMean + 1), log2FC, pch = 16, cex = .3,
                 col = ifelse(!is.na(padj) & padj < 0.05, "red", "grey50"),
                 main = paste(title, "MA (MLE, red padj<0.05)"), xlab = "log10 baseMean"))
  abline(h = 0)
  ck <- assays(dds)[["cooks"]]
  if (!is.null(ck)) boxplot(log10(ck + 1e-6), las = 2, main = "Cook's distance per sample (log10)", cex.axis = .6)
  dev.off()
}

gene_sets <- function() {
  f <- file.path(P3, "results/gene_sets.csv")
  if (file.exists(f)) {
    g <- fread(f)
    return(split(g$gene, g$set))
  }
  suppressPackageStartupMessages(library(msigdbr))
  h <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
  c2 <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "M2")
  ppar <- c2[c2$gs_name %in% c("REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA",
                               "WP_PPAR_SIGNALING_PATHWAY", "SANDERSON_PPARA_TARGETS"), ]
  sets <- c(split(h$gene_symbol, h$gs_name), split(ppar$gene_symbol, ppar$gs_name))
  # custom sets, defined in analysis_plan.md before any fit
  sets$LIPOGENESIS_CORE <- c("Fasn", "Acly", "Acaca", "Me1", "Elovl6", "Scd1", "Srebf1", "Mlxipl",
                             "Acss2", "Thrsp", "Pnpla3", "Gpam")
  sets$CHOL_SYNTHESIS_CORE <- c("Hmgcr", "Hmgcs1", "Mvk", "Pmvk", "Mvd", "Idi1", "Fdps", "Fdft1", "Sqle",
                                "Lss", "Cyp51", "Msmo1", "Nsdhl", "Hsd17b7", "Tm7sf2", "Dhcr7", "Dhcr24",
                                "Sc5d", "Ldlr", "Pcsk9", "Insig1")
  sets <- lapply(sets, unique)
  fwrite(data.table(set = rep(names(sets), lengths(sets)), gene = unlist(sets)), f)
  sets
}
CORE_SETS <- c("LIPOGENESIS_CORE", "CHOL_SYNTHESIS_CORE", "SANDERSON_PPARA_TARGETS",
               "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_INFLAMMATORY_RESPONSE",
               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

# fry (self-contained) + camera (competitive) on a log-expression matrix with its own design
pathway_tests <- function(expr, design, contrast_vec, label, min_size = 5) {
  sets <- gene_sets()
  idx <- ids2indices(sets, rownames(expr))
  idx <- idx[lengths(idx) >= min_size]
  fr <- fry(expr, index = idx, design = design, contrast = contrast_vec, sort = "none")
  cm <- camera(expr, index = idx, design = design, contrast = contrast_vec, sort = FALSE)
  out <- data.table(contrast = label, set = rownames(fr), n_genes = fr$NGenes,
                    fry_direction = fr$Direction, fry_p = fr$PValue,
                    camera_direction = cm[rownames(fr), "Direction"], camera_p = cm[rownames(fr), "PValue"])
  out[, family := ifelse(grepl("^HALLMARK", set), "Hallmark", "PPAR_custom")]
  out[, fry_fdr := p.adjust(fry_p, "BH"), by = family]
  out[, camera_fdr := p.adjust(camera_p, "BH"), by = family]
  out[, core := set %in% CORE_SETS]
  out[]
}

versions <- function() list(R = R.version.string, DESeq2 = as.character(packageVersion("DESeq2")),
                            limma = as.character(packageVersion("limma")),
                            edgeR = as.character(packageVersion("edgeR")),
                            ashr = as.character(packageVersion("ashr")))
