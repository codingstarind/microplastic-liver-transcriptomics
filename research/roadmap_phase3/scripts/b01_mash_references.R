# Stage B/D: independent MASH reference signatures (GSE316023 primary, GSE246328 secondary)
# plus pathway tests for the GSE294180 Diet_in_Veh target (existing VST, existing design; new custom sets).
# Primary GSE316023 model = NO length offset (Stage A decision); offset model = approximate sensitivity.
source("/mnt/c/Users/Soham/Desktop/microplastic/research/roadmap_phase3/scripts/common.R")
RE <- file.path(P3, "reference_effects"); RS <- file.path(P3, "results/mash"); FG <- file.path(P3, "figures/diagnostics")
dir.create(RS, recursive = TRUE, showWarnings = FALSE); dir.create(FG, recursive = TRUE, showWarnings = FALSE)
U <- universe(); el <- eligible()
info <- list(versions = versions(), seed = SEED)

fit_two_group <- function(cts, cond, ref, case, offset_len = NULL) {
  cd <- data.frame(diet = factor(cond, levels = c(ref, case)), row.names = colnames(cts))
  dds <- DESeqDataSetFromMatrix(cts, cd, ~ diet)
  if (!is.null(offset_len)) assays(dds)[["avgTxLength"]] <- offset_len[rownames(cts), colnames(cts)]
  DESeq(dds, quiet = TRUE)
}
coef_of <- function(dds) { nm <- resultsNames(dds)[2]; stopifnot(grepl("^diet_", nm)); nm }

run_ref <- function(acc, cid, ref_lab, case_lab, label, offset_len = NULL, loo = TRUE) {
  e <- el[contrast_id == cid]
  cts <- ready(acc)
  genes <- intersect(U$main, rownames(cts))
  s <- e$sample_id
  cts <- cts[genes, s]
  storage.mode(cts) <- "integer"
  cond <- ifelse(e$arm == "case", case_lab, ref_lab)
  dds <- fit_two_group(cts, cond, ref_lab, case_lab, offset_len)
  cn <- coef_of(dds)
  stopifnot(cn == paste0("diet_", case_lab, "_vs_", ref_lab))
  tab <- deseq_table(dds, coef_name = cn, genes = genes)
  tab <- add_ashr(tab, dds, coef_name = cn)
  tab[, in_universe_sensitivity := gene %in% U$sens]
  tab[, `:=`(dataset = acc, contrast = label, coefficient = cn)]
  fwrite(tab, file.path(RE, paste0(label, ".csv.gz")))
  diag_plots(dds, tab, file.path(FG, paste0(label, "_diagnostics.png")), label)
  # normalised log counts for the relabelling null (python)
  nl <- log2(counts(dds, normalized = TRUE) + 1)
  fwrite(data.table(gene = rownames(nl), nl), file.path(RS, paste0(label, "_log2norm.csv.gz")))
  # sample influence: Cook's > F(0.99) cutoff count per sample
  p <- length(resultsNames(dds)); m <- ncol(dds)
  cut <- qf(0.99, p, m - p)
  ck <- assays(dds)[["cooks"]][genes, ]
  infl <- data.table(sample = colnames(ck), n_genes_cooks_gt_F99 = colSums(ck > cut, na.rm = TRUE),
                     size_factor = sizeFactors(dds) %||% NA)
  fwrite(infl, file.path(RS, paste0(label, "_sample_influence.csv")))
  lo <- NULL
  if (loo) {
    lo <- sapply(seq_along(s), function(i) {
      d2 <- fit_two_group(cts[, -i], cond[-i], ref_lab, case_lab, offset_len)
      results(d2, name = cn, independentFiltering = FALSE)[genes, "log2FoldChange"]
    })
    colnames(lo) <- paste0("without_", s); rownames(lo) <- genes
    fwrite(data.table(gene = genes, lo), file.path(RS, paste0(label, "_loo_log2FC.csv.gz")))
  }
  # pathway tests on voom of the same samples/design
  v <- voom(DGEList(cts) |> calcNormFactors(), model.matrix(~ factor(cond, levels = c(ref_lab, case_lab))))
  pw <- pathway_tests(v$E, v$design, c(0, 1), label)
  fwrite(pw, file.path(RS, paste0(label, "_pathways.csv")))
  list(label = label, coefficient = cn, n_case = sum(cond == case_lab), n_control = sum(cond == ref_lab),
       samples_case = s[cond == case_lab], samples_control = s[cond == ref_lab], n_genes = length(genes),
       n_padj05 = sum(tab$padj < 0.05, na.rm = TRUE),
       n_deg = sum(tab$padj < 0.05 & abs(tab$log2FC) >= 0.38, na.rm = TRUE),
       n_outlier_replaced = if (!is.null(assays(dds)[["replaceCounts"]])) sum(mcols(dds)$replace, na.rm = TRUE) else 0,
       n_na_p = sum(is.na(tab$pvalue)), offset = !is.null(offset_len),
       design = "~ diet", reference_level = ref_lab)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- GSE316023 ----------------------------------------------------------------------------
lenf <- rd(file.path(P3, "results/offset_audit/gse316023_length_matrix_universe16.csv.gz"))
L <- as.matrix(lenf[, -1]); rownames(L) <- lenf[[1]]
allg <- intersect(U$main, rownames(ready("GSE316023")))
L2 <- matrix(1000, length(allg), ncol(L), dimnames = list(allg, colnames(L)))  # genes w/o length: no offset
L2[intersect(allg, rownames(L)), ] <- L[intersect(allg, rownames(L)), ]
info$GSE316023_primary <- run_ref("GSE316023", "MASH_primary", "chow", "GAN", "MASH1_GSE316023_GANvsChow_TA")
info$GSE316023_offset <- run_ref("GSE316023", "MASH_primary", "chow", "GAN", "MASH1_GSE316023_GANvsChow_TA_lengthOffset",
                                 offset_len = L2, loo = FALSE)
info$GSE316023_TN <- run_ref("GSE316023", "MASH_sens_TN", "chow", "GAN", "MASH1sens_GSE316023_GANvsChow_TN", loo = FALSE)
info$n_genes_no_length_offset_fallback <- length(setdiff(allg, rownames(L)))
# ---- GSE246328 ----------------------------------------------------------------------------
info$GSE246328_primary <- run_ref("GSE246328", "MASH_secondary", "chow", "GAN", "MASH2_GSE246328_GANvsChow_24w")

# ---- GSE294180 Diet_in_Veh pathway tests (existing VST + existing factorial design) -------
vs <- rd(file.path(REPO, "results/deseq2/provisional_deposited/vst_matrix.csv.gz"))
V <- as.matrix(vs[, -1]); rownames(V) <- vs[[1]]
sm <- fread(file.path(REPO, "metadata/sample_metadata.csv"))
sm <- sm[match(colnames(V), sm$counts_column)]
diet <- factor(sm$diet, levels = c("CD", "MD")); mp <- factor(sm$microplastic, levels = c("Veh", "PE"))
X <- model.matrix(~ diet * mp)
V <- V[intersect(rownames(V), U$main), ]
pw <- rbind(pathway_tests(V, X, c(0, 1, 0, 0), "GSE294180_Diet_in_Veh"),
            pathway_tests(V, X, c(0, 0, 1, 0), "GSE294180_PE_in_CD"),
            pathway_tests(V, X, c(0, 0, 1, 1), "GSE294180_PE_in_MD"))
fwrite(pw, file.path(RS, "GSE294180_pathways.csv"))
info$GSE294180_pathways <- list(source = "results/deseq2/provisional_deposited/vst_matrix.csv.gz",
                                design = "~ diet * microplastic", genes = nrow(V))
write_json(info, file.path(RS, "run_info.json"), auto_unbox = TRUE, pretty = TRUE)
cat("done\n")
