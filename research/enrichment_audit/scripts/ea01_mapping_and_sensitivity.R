# Enrichment audit, step 1 (computational only; reads existing inputs, writes only to research/enrichment_audit/).
#
#  (1) Gene-mapping audit: MSigDB 2026.1.Mm Hallmark (50) + 3 PPAR sets vs the GSE294180 deposited symbols
#      (raw matrix and the 13,777-gene DESeq2 universe), with an Ensembl-ID rescue check through GENCODE vM35.
#  (2) Reproduction check: re-computes fry/camera for the primary fit in memory and compares with the stored
#      results/pathways/provisional_deposited/limma_fry_camera.csv (nothing in results/ is rewritten).
#  (3) Sensitivities missing from results/pathways/ (same design, contrasts, gene sets and tests as
#      analysis/pathways/pathway_analysis.R Method 3 + fgsea):
#        S_muscle_excluded : the 309 Phase-2 muscle-exclusion genes removed before testing
#        S_cov_depth_muscle: depth and muscle covariates jointly (existing runs had one covariate at a time)
#        S_depth_thinned   : vehicle libraries binomially thinned to the median PE/vehicle total-count ratio
#                            of their diet (as research/roadmap_phase4/scripts/p4_depth_thinning.R, seed 11),
#                            then the full factorial DESeq2 + VST + fry/camera/fgsea refit.
#
# Run (WSL): /root/mpliver-tools/mamba/envs/mpliver-r/bin/Rscript ea01_mapping_and_sensitivity.R
suppressPackageStartupMessages({
  library(msigdbr); library(fgsea); library(limma); library(DESeq2); library(data.table); library(jsonlite)
})
REPO <- "/mnt/c/Users/Soham/Desktop/microplastic"
OUT <- file.path(REPO, "research/enrichment_audit")
dir.create(file.path(OUT, "results/sensitivity"), recursive = TRUE, showWarnings = FALSE)
rd <- function(f, ...) if (endsWith(f, ".gz")) fread(cmd = paste("zcat", shQuote(f)), header = TRUE, ...) else fread(f, header = TRUE, ...)

PPAR_SETS <- c("REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA", "WP_PPAR_SIGNALING_PATHWAY", "SANDERSON_PPARA_TARGETS")
hall <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
c2 <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "M2")
ppar <- c2[c2$gs_name %in% PPAR_SETS, ]
stopifnot(length(unique(hall$gs_name)) == 50, length(unique(ppar$gs_name)) == 3)
to_list <- function(df) lapply(split(df$gene_symbol, df$gs_name), unique)
FAMS <- list(Hallmark = to_list(hall), PPAR = to_list(ppar))

# ---------------------------------------------------------------- (1) mapping audit
raw <- rd(file.path(REPO, "data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz"))
raw_genes <- raw[[1]]
vst0 <- rd(file.path(REPO, "results/deseq2/provisional_deposited/vst_matrix.csv.gz"))
univ <- vst0$gene
gtf <- fread(cmd = paste("zcat", shQuote(file.path(REPO, "data/raw/gencode.vM35.annotation.gtf.gz")),
                         "| awk -F'\\t' '$3==\"gene\"' | grep -o 'gene_id \"[^\"]*\"; gene_type \"[^\"]*\"; gene_name \"[^\"]*\"'"),
             header = FALSE, sep = ";")
gtf <- data.table(ens = sub("\\..*", "", gsub('.*"(.*)"', "\\1", gtf$V1)), vm35_name = gsub('.*"(.*)"', "\\1", gtf$V3))
map_rows <- list(); miss_rows <- list()
for (fam in names(FAMS)) {
  src <- if (fam == "Hallmark") hall else ppar
  for (s in names(FAMS[[fam]])) {
    sym <- FAMS[[fam]][[s]]
    sub <- unique(src[src$gs_name == s, c("gene_symbol", "ensembl_gene")])
    in_raw <- sym %in% raw_genes; in_univ <- sym %in% univ
    miss <- sym[!in_raw]
    # rescue: MSigDB Ensembl id -> vM35 gene_name present in the deposited matrix under another symbol
    m_ens <- sub[sub$gene_symbol %in% miss & !is.na(sub$ensembl_gene), ]
    resc <- merge(as.data.table(m_ens), gtf, by.x = "ensembl_gene", by.y = "ens")
    resc <- resc[vm35_name %in% raw_genes & !(vm35_name %in% sym)]
    ci <- miss[toupper(miss) %in% toupper(raw_genes)]
    map_rows[[length(map_rows) + 1]] <- data.table(
      family = fam, pathway = s, n_msigdb_symbols = length(sym), n_in_raw_matrix = sum(in_raw),
      n_in_deseq2_universe = sum(in_univ), pct_in_universe = round(100 * mean(in_univ), 1),
      n_not_in_raw = length(miss), n_rescuable_by_ensembl = uniqueN(resc$gene_symbol),
      rescuable_symbols = paste(unique(paste0(resc$gene_symbol, "->", resc$vm35_name)), collapse = ";"),
      n_case_only_mismatch = length(ci))
    if (length(miss)) miss_rows[[length(miss_rows) + 1]] <- data.table(family = fam, pathway = s, symbol = miss,
      in_gencode_vm35_names = miss %in% gtf$vm35_name)
  }
}
map_tab <- rbindlist(map_rows); fwrite(map_tab, file.path(OUT, "results/gene_mapping_audit.csv"))
fwrite(rbindlist(miss_rows), file.path(OUT, "results/gene_mapping_unmatched_symbols.csv"))
dup <- sum(duplicated(raw_genes))
cat("mapping: median % in universe", median(map_tab$pct_in_universe), "; rescuable total", sum(map_tab$n_rescuable_by_ensembl),
    "; duplicated raw symbols", dup, "\n")

# ---------------------------------------------------------------- shared helpers (mirror pathway_analysis.R)
meta0 <- read.csv(file.path(REPO, "metadata/sample_metadata.csv"))
covs <- read.csv(file.path(REPO, "results/deseq2/covariates.csv"), row.names = 1)
L4 <- rbind(PE_in_CD = c(0, 0, 1, 0), PE_in_MD = c(0, 0, 1, 1), Diet_in_Veh = c(0, 1, 0, 0), Interaction = c(0, 0, 0, 1))
prep_meta <- function(samples, covariates = character(0)) {
  m <- meta0[match(samples, meta0$counts_column), ]
  m$diet <- factor(m$diet, levels = c("CD", "MD")); m$microplastic <- factor(m$microplastic, levels = c("Veh", "PE"))
  for (cv in covariates) m[[cv]] <- as.numeric(scale(covs[m$counts_column, cv]))
  m
}
contrast_mat <- function(n_cov) if (n_cov) cbind(L4[, 1], matrix(0, 4, n_cov), L4[, -1]) else L4
run_fry_camera <- function(expr, meta, covariates, lfc_list, variant) {
  f <- as.formula(paste("~", paste(c(covariates, "diet * microplastic"), collapse = " + ")))
  design <- model.matrix(f, data = meta); L <- contrast_mat(length(covariates))
  stopifnot(ncol(L) == ncol(design))
  genes <- rownames(expr); rows <- list()
  for (fam in names(FAMS)) {
    idx <- ids2indices(FAMS[[fam]], genes, remove.empty = TRUE); idx <- idx[lengths(idx) >= 10]
    for (k in rownames(L)) {
      fr <- fry(expr, index = idx, design = design, contrast = L[k, ], sort = "none")
      cm <- camera(expr, index = idx, design = design, contrast = L[k, ], sort = FALSE)
      lfc <- lfc_list[[k]]
      rows[[length(rows) + 1]] <- data.table(variant = variant, family = fam, pathway = rownames(fr), contrast = k,
        n_genes = fr$NGenes,
        mean_gene_log2FC_ashr = if (is.null(lfc)) NA_real_ else vapply(rownames(fr), function(s) mean(lfc[intersect(names(lfc), genes[idx[[s]]])], na.rm = TRUE), 0),
        fry_direction = fr$Direction, fry_pvalue = fr$PValue, fry_pvalue_mixed = fr$PValue.Mixed,
        camera_direction = cm[rownames(fr), "Direction"], camera_pvalue = cm[rownames(fr), "PValue"])
    }
  }
  out <- rbindlist(rows)
  out[, `:=`(fry_fdr = p.adjust(fry_pvalue, "BH"), fry_fdr_mixed = p.adjust(fry_pvalue_mixed, "BH"),
             camera_fdr = p.adjust(camera_pvalue, "BH")), by = .(family, contrast)]
  out[]
}
run_fgsea <- function(stat_list, variant) {
  set.seed(20261102); rows <- list()
  for (k in names(stat_list)) for (fam in names(FAMS)) {
    fg <- suppressWarnings(fgseaMultilevel(FAMS[[fam]], stat_list[[k]], minSize = 15, maxSize = 500, eps = 0))
    fg$leadingEdge <- vapply(fg$leadingEdge, function(x) paste(head(x, 15), collapse = ";"), "")
    rows[[length(rows) + 1]] <- data.table(variant = variant, family = fam, contrast = k, fg)
  }
  rbindlist(rows)
}
load_lfc <- function(label) {
  setNames(lapply(rownames(L4), function(k) {
    r <- rd(file.path(REPO, "results/deseq2", label, paste0("results_", k, ".csv.gz"))); setNames(r$log2FC_ashr, r$gene)
  }), rownames(L4))
}
expr0 <- as.matrix(vst0[, -1]); rownames(expr0) <- vst0$gene
meta_p <- prep_meta(colnames(expr0))
lfc0 <- load_lfc("provisional_deposited")

# ---------------------------------------------------------------- (2) reproduction check
rep <- run_fry_camera(expr0, meta_p, character(0), lfc0, "reproduced_primary")
stored <- fread(file.path(REPO, "results/pathways/provisional_deposited/limma_fry_camera.csv"))
cmp <- merge(rep, stored, by = c("family", "pathway", "contrast"), suffixes = c("_new", "_stored"))
repro <- data.table(n_rows_compared = nrow(cmp),
  max_abs_diff_fry_p = max(abs(cmp$fry_pvalue_new - cmp$fry_pvalue_stored)),
  max_abs_diff_camera_p = max(abs(cmp$camera_pvalue_new - cmp$camera_pvalue_stored)),
  max_abs_diff_fry_fdr = max(abs(cmp$fry_fdr_new - cmp$fry_fdr_stored)),
  direction_mismatches = sum(cmp$fry_direction_new != cmp$fry_direction_stored),
  n_genes_mismatches = sum(cmp$n_genes_new != cmp$n_genes_stored))
fwrite(repro, file.path(OUT, "results/reproduction_check_primary_fry_camera.csv")); print(repro)

# ---------------------------------------------------------------- (3a) muscle genes excluded
musc <- fread(file.path(REPO, "research/roadmap_phase2/metadata/muscle_exclusion_genes.csv"))$vm35_symbol
keep <- !(rownames(expr0) %in% musc)
sm <- run_fry_camera(expr0[keep, ], meta_p, character(0), lfc0, "S_muscle_excluded")
sm_fg <- run_fgsea(lapply(setNames(rownames(L4), rownames(L4)), function(k) {
  r <- rd(file.path(REPO, "results/deseq2/provisional_deposited", paste0("results_", k, ".csv.gz")))
  r <- r[!is.na(stat) & !(gene %in% musc)]; setNames(r$stat, r$gene) }), "S_muscle_excluded")
musc_info <- list(n_muscle_genes_list = length(musc), n_removed_from_universe = sum(!keep),
                  removed_from_sets = sapply(FAMS$Hallmark, function(g) sum(g %in% musc[musc %in% rownames(expr0)])))

# ---------------------------------------------------------------- (3b) depth + muscle covariates jointly
sdm <- run_fry_camera(expr0, prep_meta(colnames(expr0), c("depth", "muscle")), c("depth", "muscle"), NULL, "S_cov_depth_muscle")

# ---------------------------------------------------------------- (3c) depth thinning, full refit
set.seed(11)
M <- as.matrix(raw[, -1]); rownames(M) <- raw_genes; storage.mode(M) <- "integer"
mt <- prep_meta(colnames(M)); thin_info <- list()
for (d in c("CD", "MD")) {
  veh <- mt$counts_column[mt$diet == d & mt$microplastic == "Veh"]; pe <- mt$counts_column[mt$diet == d & mt$microplastic == "PE"]
  p <- median(colSums(M[, pe])) / median(colSums(M[, veh]))
  for (s in veh) M[, s] <- rbinom(nrow(M), M[, s], p)
  thin_info[[d]] <- list(thinning_p = p, vehicle = veh)
}
rownames(mt) <- mt$counts_column
dds <- DESeqDataSetFromMatrix(M, mt, ~ diet + microplastic + diet:microplastic)
dds <- dds[rowSums(counts(dds) >= 10) >= 4, ]
dds <- DESeq(dds, quiet = TRUE, minReplicatesForReplace = 7)
cons <- list(PE_in_CD = list("microplastic_PE_vs_Veh"), PE_in_MD = list(c("microplastic_PE_vs_Veh", "dietMD.microplasticPE")),
             Diet_in_Veh = list("diet_MD_vs_CD"), Interaction = list("dietMD.microplasticPE"))
stat_t <- list(); lfc_t <- list(); de_sum <- list()
for (k in names(cons)) {
  r <- results(dds, contrast = cons[[k]], alpha = 0.05)
  sh <- lfcShrink(dds, contrast = cons[[k]], res = r, type = "ashr", quiet = TRUE)
  fwrite(data.table(gene = rownames(r), log2FC = r$log2FoldChange, stat = r$stat, pvalue = r$pvalue, padj = r$padj,
                    log2FC_ashr = sh$log2FoldChange), file.path(OUT, "results/sensitivity", paste0("depth_thinned_deseq2_", k, ".csv.gz")))
  ok <- !is.na(r$stat); stat_t[[k]] <- setNames(r$stat[ok], rownames(r)[ok]); lfc_t[[k]] <- setNames(sh$log2FoldChange, rownames(r))
  de_sum[[k]] <- sum(r$padj < 0.05, na.rm = TRUE)
}
vt <- assay(vst(dds, blind = FALSE))
st <- run_fry_camera(vt, prep_meta(colnames(vt)), character(0), lfc_t, "S_depth_thinned")
st_fg <- run_fgsea(stat_t, "S_depth_thinned")
lib_after <- colSums(M)

fwrite(rbindlist(list(sm, sdm, st), fill = TRUE), file.path(OUT, "results/sensitivity/new_sensitivity_fry_camera.csv"))
fwrite(rbindlist(list(sm_fg, st_fg)), file.path(OUT, "results/sensitivity/new_sensitivity_fgsea.csv"))
write_json(list(msigdb_version = unique(hall$db_version), msigdbr = as.character(packageVersion("msigdbr")),
                limma = as.character(packageVersion("limma")), fgsea = as.character(packageVersion("fgsea")),
                DESeq2 = as.character(packageVersion("DESeq2")), R = R.version.string,
                duplicated_raw_symbols = dup, muscle = musc_info,
                depth_thinning = c(thin_info, list(seed = 11, n_genes_after_prefilter = nrow(dds), n_deg_padj05 = de_sum,
                  library_sizes_after = as.list(lib_after)))),
           file.path(OUT, "results/sensitivity/run_info.json"), auto_unbox = TRUE, pretty = TRUE)
cat("done\n")
