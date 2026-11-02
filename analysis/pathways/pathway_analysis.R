# Pathway-level analysis for the factorial model (plan: reports/phase2/analysis_plan.md).
#
# Method 1: per-mouse pathway scores (mean z of VST expression over set genes),
#   fitted with lm(score ~ diet * microplastic); contrasts with SE, 95% CI, t-test, BH per family.
#   Diluted for sets whose genes move in opposite directions (see the pathway report).
# Method 2 (supporting): preranked fgsea on each contrast's DESeq2 Wald statistic
#   (gene = sampling unit; ignores inter-gene correlation, so anti-conservative).
# Method 3 (added after Method 1 failed a positive-control check, see deviations):
#   limma fry (self-contained rotation test, mouse = unit, directional and mixed) and
#   camera (competitive, corrected for inter-gene correlation) on the VST matrix with the
#   same factorial design. Categories of the PE response use fry.
#
# Usage: Rscript analysis/pathways/pathway_analysis.R <repo_dir> <deseq2_label> [fast]
#   fast: skip the exploratory Reactome family (used for the 19 leave-one-out refits)

suppressPackageStartupMessages({
  library(msigdbr)
  library(fgsea)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
repo <- args[1]
label <- args[2]
fast <- length(args) >= 3 && args[3] == "fast"
res_dir <- file.path(repo, "results", "deseq2", label)
out_dir <- file.path(repo, "results", "pathways", label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20261102)

FOCUS <- c(
  "HALLMARK_FATTY_ACID_METABOLISM" = "Fatty-acid metabolism",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS" = "Cholesterol",
  "HALLMARK_PEROXISOME" = "Peroxisome / bile acids",
  "HALLMARK_BILE_ACID_METABOLISM" = "Peroxisome / bile acids",
  "HALLMARK_ADIPOGENESIS" = "Lipid storage",
  "HALLMARK_INFLAMMATORY_RESPONSE" = "Inflammation",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "Inflammation",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING" = "Inflammation",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE" = "Inflammation",
  "HALLMARK_XENOBIOTIC_METABOLISM" = "Xenobiotic",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY" = "Oxidative / cellular stress",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE" = "Oxidative / cellular stress",
  "HALLMARK_P53_PATHWAY" = "Oxidative / cellular stress",
  "HALLMARK_HYPOXIA" = "Oxidative / cellular stress",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" = "Mitochondrial",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "Fibrosis / remodelling",
  "HALLMARK_TGF_BETA_SIGNALING" = "Fibrosis / remodelling",
  "HALLMARK_MYOGENESIS" = "Contamination flag")
# Deviation from the plan: the Reactome PPARA set is named
# REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA in MSigDB 2026.1.Mm, and
# SANDERSON_PPARA_TARGETS (mouse liver PPARalpha targets) was added from the set list.
PPAR_SETS <- c("REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA",
               "WP_PPAR_SIGNALING_PATHWAY", "SANDERSON_PPARA_TARGETS")

# ---- Gene sets (MSigDB mouse-native release) ---------------------------------------
hall <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
c2 <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "M2")
ppar <- c2[c2$gs_name %in% PPAR_SETS, ]
react <- c2[c2$gs_subcollection == "CP:REACTOME", ]
db_version <- unique(hall$db_version)
to_list <- function(df) split(df$gene_symbol, df$gs_name)

vst <- as.data.table(read.csv(gzfile(file.path(res_dir, "vst_matrix.csv.gz")), check.names = FALSE))
genes <- vst$gene
expr <- as.matrix(vst[, -1]); rownames(expr) <- genes
meta <- read.csv(file.path(repo, "metadata", "sample_metadata.csv"))
meta <- meta[match(colnames(expr), meta$counts_column), ]
meta$diet <- factor(meta$diet, levels = c("CD", "MD"))
meta$microplastic <- factor(meta$microplastic, levels = c("Veh", "PE"))

# Use the same covariate as the DESeq2 fit (sensitivity variants), so pathway models match.
run_info <- jsonlite::fromJSON(file.path(res_dir, "run_info.json"))
covariate <- if (length(run_info$covariate)) run_info$covariate else character(0)
if (length(covariate)) {
  cov <- read.csv(file.path(repo, "results", "deseq2", "covariates.csv"), row.names = 1)
  meta$cov <- as.numeric(scale(cov[meta$counts_column, covariate]))
}
design_formula <- if (length(covariate)) ~ cov + diet * microplastic else ~ diet * microplastic

# ---- Method 1: per-mouse scores + factorial linear model ----------------------------
z <- t(scale(t(expr)))
score_sets <- function(sets) {
  sets <- lapply(sets, intersect, genes)
  sets <- sets[lengths(sets) >= 10]
  scores <- sapply(sets, function(g) colMeans(z[g, , drop = FALSE]))
  list(scores = scores, sizes = lengths(sets))
}
# Contrast vectors over coefficients (Intercept, [cov,] dietMD, microplasticPE, dietMD:microplasticPE)
L <- rbind(PE_in_CD = c(0, 0, 1, 0), PE_in_MD = c(0, 0, 1, 1),
           Diet_in_Veh = c(0, 1, 0, 0), Interaction = c(0, 0, 0, 1))
if (length(covariate)) L <- cbind(L[, 1], 0, L[, -1])
stopifnot(ncol(L) == ncol(model.matrix(design_formula, data = meta)))
fit_scores <- function(sc, family) {
  rows <- list()
  for (set in colnames(sc$scores)) {
    meta$score <- sc$scores[, set]
    fit <- lm(update(design_formula, score ~ .), data = meta)
    b <- coef(fit); V <- vcov(fit); df <- fit$df.residual
    for (k in rownames(L)) {
      est <- sum(L[k, ] * b); se <- sqrt(drop(t(L[k, ]) %*% V %*% L[k, ]))
      tq <- qt(0.975, df)
      rows[[length(rows) + 1]] <- data.frame(
        family = family, pathway = set, n_genes = sc$sizes[[set]], contrast = k,
        estimate = est, se = se, ci95_low = est - tq * se, ci95_high = est + tq * se,
        t = est / se, df = df, pvalue = 2 * pt(-abs(est / se), df))
    }
  }
  out <- do.call(rbind, rows)
  out$padj <- ave(out$pvalue, out$contrast, FUN = function(p) p.adjust(p, "BH"))
  out
}
hall_sc <- score_sets(to_list(hall))
ppar_sc <- score_sets(to_list(ppar))
scores_tab <- rbind(fit_scores(hall_sc, "Hallmark"), fit_scores(ppar_sc, "PPAR"))
if (!fast) {
  react_sc <- score_sets(to_list(react[react$gs_name %in% names(which(table(react$gs_name) <= 500)), ]))
  scores_tab <- rbind(scores_tab, fit_scores(react_sc, "Reactome (exploratory)"))
}
scores_tab$theme <- ifelse(scores_tab$pathway %in% names(FOCUS), FOCUS[scores_tab$pathway],
                           ifelse(scores_tab$family == "PPAR", "PPAR", ""))
fwrite(scores_tab, file.path(out_dir, "pathway_score_effects.csv.gz"))
fwrite(data.table(sample = rownames(hall_sc$scores), hall_sc$scores),
       file.path(out_dir, "hallmark_scores_per_mouse.csv"))

# ---- Categorise the PE response across diets (Hallmark + PPAR) ----------------------
excl <- function(lo, hi) lo > 0 | hi < 0
categorise <- function(cd, md, int) {
  cd_on <- excl(cd$ci95_low, cd$ci95_high); md_on <- excl(md$ci95_low, md$ci95_high)
  int_on <- excl(int$ci95_low, int$ci95_high)
  if (!cd_on && !md_on) return(if (int_on) "No clear PE effect in either diet (interaction CI excludes 0)"
                               else "No clear PE effect")
  if (cd_on && md_on && sign(cd$estimate) != sign(md$estimate)) return("Opposing")
  if (cd_on && md_on) {
    if (!int_on) return("Similar in both diets")
    return(if (sign(int$estimate) == sign(cd$estimate)) "Amplified under MD" else "Attenuated under MD")
  }
  only <- if (cd_on) "CD-only" else "MD-only"
  paste0(only, if (int_on) ", interaction supported" else ", interaction not supported")
}
cats <- list()
for (fam in c("Hallmark", "PPAR")) {
  sub <- scores_tab[scores_tab$family == fam, ]
  for (set in unique(sub$pathway)) {
    s <- sub[sub$pathway == set, ]
    g <- function(k) s[s$contrast == k, ]
    cats[[length(cats) + 1]] <- data.frame(
      family = fam, pathway = set, theme = s$theme[1], n_genes = s$n_genes[1],
      PE_in_CD = g("PE_in_CD")$estimate, PE_in_CD_low = g("PE_in_CD")$ci95_low,
      PE_in_CD_high = g("PE_in_CD")$ci95_high,
      PE_in_MD = g("PE_in_MD")$estimate, PE_in_MD_low = g("PE_in_MD")$ci95_low,
      PE_in_MD_high = g("PE_in_MD")$ci95_high,
      Interaction = g("Interaction")$estimate, Interaction_low = g("Interaction")$ci95_low,
      Interaction_high = g("Interaction")$ci95_high, Interaction_padj = g("Interaction")$padj,
      Diet_in_Veh = g("Diet_in_Veh")$estimate, Diet_in_Veh_low = g("Diet_in_Veh")$ci95_low,
      Diet_in_Veh_high = g("Diet_in_Veh")$ci95_high,
      category = categorise(g("PE_in_CD"), g("PE_in_MD"), g("Interaction")))
  }
}
cat_tab <- do.call(rbind, cats)
fwrite(cat_tab, file.path(out_dir, "pe_response_categories.csv"))

# ---- Method 2: preranked fgsea --------------------------------------------------------
fg_rows <- list()
for (k in rownames(L)) {
  res <- as.data.table(read.csv(gzfile(file.path(res_dir, paste0("results_", k, ".csv.gz")))))
  res <- res[!is.na(stat)]
  ranks <- setNames(res$stat, res$gene)
  fams <- list(list("Hallmark", to_list(hall)), list("PPAR", to_list(ppar)))
  if (!fast) fams[[3]] <- list("Reactome (exploratory)", to_list(react))
  for (fam in fams) {
    fg <- suppressWarnings(fgseaMultilevel(fam[[2]], ranks, minSize = 15, maxSize = 500, eps = 0))
    fg$leadingEdge <- vapply(fg$leadingEdge, function(x) paste(head(x, 15), collapse = ";"), "")
    fg$family <- fam[[1]]; fg$contrast <- k
    fg_rows[[length(fg_rows) + 1]] <- fg
  }
}
fg_tab <- rbindlist(fg_rows)
fwrite(fg_tab, file.path(out_dir, "fgsea_results.csv.gz"))

# ---- Method 3: limma fry / camera (mouse = unit, inter-gene correlation accounted for) --
suppressPackageStartupMessages(library(limma))
design <- model.matrix(design_formula, data = meta)
fit_lfc <- lapply(rownames(L), function(k) {
  res <- read.csv(gzfile(file.path(res_dir, paste0("results_", k, ".csv.gz"))))
  setNames(res$log2FC_ashr, res$gene)
})
names(fit_lfc) <- rownames(L)
lim_rows <- list()
for (fam in list(list("Hallmark", to_list(hall)), list("PPAR", to_list(ppar)))) {
  idx <- ids2indices(fam[[2]], genes, remove.empty = TRUE)
  idx <- idx[lengths(idx) >= 10]
  for (k in rownames(L)) {
    fr <- fry(expr, index = idx, design = design, contrast = L[k, ], sort = "none")
    cm <- camera(expr, index = idx, design = design, contrast = L[k, ], sort = FALSE)
    lfc <- fit_lfc[[k]]
    lim_rows[[length(lim_rows) + 1]] <- data.frame(
      family = fam[[1]], pathway = rownames(fr), contrast = k, n_genes = fr$NGenes,
      mean_gene_log2FC_ashr = vapply(rownames(fr), function(s)
        mean(lfc[intersect(names(lfc), genes[idx[[s]]])], na.rm = TRUE), 0),
      frac_genes_up = vapply(rownames(fr), function(s)
        mean(lfc[intersect(names(lfc), genes[idx[[s]]])] > 0, na.rm = TRUE), 0),
      fry_direction = fr$Direction, fry_pvalue = fr$PValue, fry_pvalue_mixed = fr$PValue.Mixed,
      camera_direction = cm[rownames(fr), "Direction"], camera_pvalue = cm[rownames(fr), "PValue"])
  }
}
lim_tab <- do.call(rbind, lim_rows)
grp_key <- paste(lim_tab$family, lim_tab$contrast)
lim_tab$fry_fdr <- ave(lim_tab$fry_pvalue, grp_key, FUN = function(p) p.adjust(p, "BH"))
lim_tab$fry_fdr_mixed <- ave(lim_tab$fry_pvalue_mixed, grp_key, FUN = function(p) p.adjust(p, "BH"))
lim_tab$camera_fdr <- ave(lim_tab$camera_pvalue, grp_key, FUN = function(p) p.adjust(p, "BH"))
lim_tab$theme <- ifelse(lim_tab$pathway %in% names(FOCUS), FOCUS[lim_tab$pathway],
                        ifelse(lim_tab$family == "PPAR", "PPAR", ""))
fwrite(lim_tab, file.path(out_dir, "limma_fry_camera.csv"))

# Categories from fry (directional test, FDR < 0.05 within family and contrast).
fry_call <- function(row) if (row$fry_fdr < 0.05) ifelse(row$fry_direction == "Up", 1, -1) else 0
fry_cats <- list()
for (fam in unique(lim_tab$family)) for (set in unique(lim_tab$pathway[lim_tab$family == fam])) {
  s <- lim_tab[lim_tab$family == fam & lim_tab$pathway == set, ]
  g <- function(k) s[s$contrast == k, ]
  cd <- fry_call(g("PE_in_CD")); md <- fry_call(g("PE_in_MD")); it <- fry_call(g("Interaction"))
  int_mixed <- g("Interaction")$fry_fdr_mixed < 0.05
  category <- if (cd == 0 && md == 0) {
    if (it != 0) "No clear PE effect in either diet (interaction detected)" else "No clear PE effect"
  } else if (cd != 0 && md != 0 && cd != md) "Opposing"
  else if (cd != 0 && md != 0) {
    if (it == 0) "Similar in both diets" else if (it == cd) "Amplified under MD" else "Attenuated under MD"
  } else paste0(if (cd != 0) "CD-only" else "MD-only",
                if (it != 0) ", interaction supported" else ", interaction not supported")
  fry_cats[[length(fry_cats) + 1]] <- data.frame(
    family = fam, pathway = set, theme = s$theme[1], n_genes = s$n_genes[1],
    PE_in_CD_call = cd, PE_in_MD_call = md, Interaction_call = it,
    Interaction_mixed_fdr = g("Interaction")$fry_fdr_mixed,
    Diet_in_Veh_call = fry_call(g("Diet_in_Veh")),
    PE_in_CD_meanLFC = g("PE_in_CD")$mean_gene_log2FC_ashr,
    PE_in_MD_meanLFC = g("PE_in_MD")$mean_gene_log2FC_ashr,
    Diet_in_Veh_meanLFC = g("Diet_in_Veh")$mean_gene_log2FC_ashr,
    category_fry = category)
}
fry_cat_tab <- do.call(rbind, fry_cats)
fry_cat_tab <- merge(fry_cat_tab, cat_tab[, c("family", "pathway", "category")],
                     by = c("family", "pathway"), all.x = TRUE)
names(fry_cat_tab)[names(fry_cat_tab) == "category"] <- "category_meanz_score"
fwrite(fry_cat_tab, file.path(out_dir, "pe_response_categories_fry.csv"))
print(table(fry_cat_tab$category_fry))

writeLines(jsonlite::toJSON(list(msigdb_version = db_version,
                                 msigdbr = as.character(packageVersion("msigdbr")),
                                 fgsea = as.character(packageVersion("fgsea")),
                                 n_hallmark_scored = ncol(hall_sc$scores),
                                 n_reactome_scored = if (fast) 0 else ncol(react_sc$scores),
                                 ppar_sets_found = colnames(ppar_sc$scores)),
                            auto_unbox = TRUE, pretty = TRUE),
           file.path(out_dir, "run_info.json"))
print(table(cat_tab$category))
print(cat_tab[cat_tab$theme != "" | cat_tab$family == "PPAR",
              c("pathway", "PE_in_CD", "PE_in_MD", "Interaction", "Diet_in_Veh", "category")],
      digits = 2)
