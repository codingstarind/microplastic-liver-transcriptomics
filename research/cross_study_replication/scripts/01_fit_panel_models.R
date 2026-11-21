# Cross-study replication: model fits that the existing DESeq2 outputs do not already contain.
#
# Usage (inside WSL, mpliver-r environment, from the repository root):
#   Rscript research/cross_study_replication/scripts/01_fit_panel_models.R <repo_dir> [n_boot]
#
# Reads the PRIMARY factorial fit (results/deseq2/provisional_deposited/dds.rds, design
# ~ diet + microplastic + diet:microplastic, reference levels CD and Veh) and writes to
# research/cross_study_replication/intermediate/:
#
#   additive_model_panel.csv   ~ diet + microplastic fitted genome-wide; Wald PE main effect
#                              (BH over all tested genes) for the panel genes
#   interaction_lrt_panel.csv  LRT of full (factorial) vs reduced (additive) model, genome-wide
#                              BH, panel rows only
#   ols_vst_panel.csv          course-style two-way ANOVA (lm on VST values) per panel gene:
#                              additive vs factorial R^2, adjusted R^2, AIC, interaction F test
#   bootstrap_panel.csv        exploratory stratified bootstrap (mice resampled with
#                              replacement within each of the 4 groups)
#
# Nothing in results/ is modified.

suppressPackageStartupMessages({
  library(DESeq2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
repo <- args[1]
n_boot <- if (length(args) >= 2) as.integer(args[2]) else 2000L
set.seed(20261121)

panel <- c("Fasn", "Acly", "Acaca", "Ppara", "Acot3", "Abcc3", "Nr1i3", "Nr4a1")
out_dir <- file.path(repo, "research", "cross_study_replication", "intermediate")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dds <- readRDS(file.path(repo, "results", "deseq2", "provisional_deposited", "dds.rds"))
stopifnot(levels(dds$diet)[1] == "CD", levels(dds$microplastic)[1] == "Veh")
stopifnot(all(panel %in% rownames(dds)))
z <- qnorm(0.975)

# ---- 1. Additive model, genome-wide ----------------------------------------------------
dds_add <- dds
design(dds_add) <- ~ diet + microplastic
dds_add <- DESeq(dds_add, quiet = TRUE)
res_add <- results(dds_add, name = "microplastic_PE_vs_Veh", alpha = 0.05)
add_tab <- data.frame(gene = panel,
                      term = "PE main effect (additive model; averaged over diets)",
                      log2FC = res_add[panel, "log2FoldChange"],
                      lfcSE = res_add[panel, "lfcSE"],
                      pvalue = res_add[panel, "pvalue"],
                      padj_genomewide = res_add[panel, "padj"])
add_tab$ci95_low <- add_tab$log2FC - z * add_tab$lfcSE
add_tab$ci95_high <- add_tab$log2FC + z * add_tab$lfcSE
res_add_diet <- results(dds_add, name = "diet_MD_vs_CD", alpha = 0.05)
add_tab$diet_log2FC_additive <- res_add_diet[panel, "log2FoldChange"]
add_tab$n_genes_tested_genomewide <- sum(!is.na(res_add$pvalue))
write.csv(add_tab, file.path(out_dir, "additive_model_panel.csv"), row.names = FALSE)

# ---- 2. Likelihood-ratio test: factorial vs additive -----------------------------------
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ diet + microplastic, quiet = TRUE)
res_lrt <- results(dds_lrt, alpha = 0.05)
lrt_tab <- data.frame(gene = panel,
                      lrt_stat = res_lrt[panel, "stat"],
                      lrt_df = 1,
                      lrt_pvalue = res_lrt[panel, "pvalue"],
                      lrt_padj_genomewide = res_lrt[panel, "padj"],
                      n_genes_tested_genomewide = sum(!is.na(res_lrt$pvalue)),
                      n_genes_lrt_padj_lt_0.05 = sum(res_lrt$padj < 0.05, na.rm = TRUE))
write.csv(lrt_tab, file.path(out_dir, "interaction_lrt_panel.csv"), row.names = FALSE)

# ---- 3. Course-style OLS two-way ANOVA on VST values ----------------------------------
vsd <- vst(dds, blind = FALSE)
cd <- as.data.frame(colData(dds))
ols_rows <- list()
for (g in panel) {
  y <- assay(vsd)[g, ]
  d <- data.frame(y = y, diet = cd$diet, pe = cd$microplastic)
  m_add <- lm(y ~ diet + pe, data = d)
  m_int <- lm(y ~ diet * pe, data = d)
  a <- anova(m_add, m_int)
  ci <- confint(m_int)["dietMD:pePE", ]
  ols_rows[[g]] <- data.frame(
    gene = g,
    r2_additive = summary(m_add)$r.squared, r2_factorial = summary(m_int)$r.squared,
    adj_r2_additive = summary(m_add)$adj.r.squared, adj_r2_factorial = summary(m_int)$adj.r.squared,
    aic_additive = AIC(m_add), aic_factorial = AIC(m_int),
    interaction_estimate_vst = coef(m_int)[["dietMD:pePE"]],
    interaction_ci95_low_vst = ci[[1]], interaction_ci95_high_vst = ci[[2]],
    interaction_F = a$F[2], interaction_F_pvalue = a$`Pr(>F)`[2],
    residual_df_factorial = m_int$df.residual,
    shapiro_pvalue_residuals = shapiro.test(residuals(m_int))$p.value)
}
write.csv(do.call(rbind, ols_rows), file.path(out_dir, "ols_vst_panel.csv"), row.names = FALSE)

# ---- 4. Exploratory stratified bootstrap ------------------------------------------------
# Mice are resampled with replacement within each diet x PE group, so every replicate keeps
# the 2 x 2 structure and group sizes (4/5/5/5). Each mouse keeps its own size factor.
# Dispersions are FIXED at the primary fit's values: re-estimating them on resamples that
# contain duplicated mice would underestimate them. Only coefficients are re-estimated.
sub <- dds[panel, ]
grp <- paste(sub$diet, sub$microplastic, sep = "+")
idx_by_grp <- split(seq_len(ncol(sub)), grp)
boot_rows <- vector("list", n_boot)
for (b in seq_len(n_boot)) {
  idx <- unlist(lapply(idx_by_grp, function(i) i[sample.int(length(i), length(i), replace = TRUE)]))
  bs <- DESeqDataSetFromMatrix(counts(sub)[, idx], colData(sub)[idx, c("diet", "microplastic")],
                               ~ diet + microplastic + diet:microplastic)
  colnames(bs) <- make.unique(colnames(sub)[idx])
  sizeFactors(bs) <- sizeFactors(sub)[idx]
  dispersions(bs) <- dispersions(sub)
  bs <- suppressMessages(nbinomWaldTest(bs, quiet = TRUE))
  co <- coef(bs)   # DESeq2 coefficients are on the log2 scale
  boot_rows[[b]] <- data.frame(
    replicate = b, gene = panel,
    PE_in_CD = co[, "microplastic_PE_vs_Veh"],
    PE_in_MD = co[, "microplastic_PE_vs_Veh"] + co[, "dietMD.microplasticPE"],
    Diet_in_Veh = co[, "diet_MD_vs_CD"],
    Interaction = co[, "dietMD.microplasticPE"],
    n_distinct_mice = length(unique(idx)))
}
boot <- do.call(rbind, boot_rows)
write.csv(boot, gzfile(file.path(out_dir, "bootstrap_panel.csv.gz")), row.names = FALSE)

# Reference fit with the same fixed-dispersion machinery on the ORIGINAL data, so the
# bootstrap is compared with estimates computed identically.
ref <- DESeqDataSetFromMatrix(counts(sub), colData(sub)[, c("diet", "microplastic")],
                              ~ diet + microplastic + diet:microplastic)
sizeFactors(ref) <- sizeFactors(sub)
dispersions(ref) <- dispersions(sub)
ref <- suppressMessages(nbinomWaldTest(ref, quiet = TRUE))
co <- coef(ref)
write.csv(data.frame(gene = panel,
                     PE_in_CD = co[, "microplastic_PE_vs_Veh"],
                     PE_in_MD = co[, "microplastic_PE_vs_Veh"] + co[, "dietMD.microplasticPE"],
                     Diet_in_Veh = co[, "diet_MD_vs_CD"],
                     Interaction = co[, "dietMD.microplasticPE"]),
          file.path(out_dir, "bootstrap_reference_fit.csv"), row.names = FALSE)

write_json(list(n_boot = n_boot, seed = 20261121, panel = panel,
                primary_fit = "results/deseq2/provisional_deposited/dds.rds",
                versions = list(R = R.version.string,
                                DESeq2 = as.character(packageVersion("DESeq2")))),
           file.path(out_dir, "run_info_01.json"), auto_unbox = TRUE, pretty = TRUE)
cat("done\n")
