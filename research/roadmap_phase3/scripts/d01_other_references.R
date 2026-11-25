# Stages E-G: fasting (GSE107787, GSE225695), circadian (GSE117134; GSE107787 FED male sens.),
# PPARalpha (GSE282998 TempO-Seq). Models exactly as in analysis_plan.md sections 6-8.
source("/mnt/c/Users/Soham/Desktop/microplastic/research/roadmap_phase3/scripts/common.R")
RE <- file.path(P3, "reference_effects"); FG <- file.path(P3, "figures/diagnostics")
U <- universe(); el <- eligible(); info <- list(versions = versions(), seed = SEED)
dir.create(file.path(P3, "results/fasting"), showWarnings = FALSE)
dir.create(file.path(P3, "results/circadian"), showWarnings = FALSE)
dir.create(file.path(P3, "results/ppara"), showWarnings = FALSE)

limma_table <- function(fit, coef, genes, label, extra = list()) {
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none", confint = TRUE)
  tt <- tt[genes, ]
  out <- data.table(gene = genes, log2FC = tt$logFC, SE = tt$logFC / tt$t, ci95_low = tt$CI.L, ci95_high = tt$CI.R,
                    stat = tt$t, pvalue = tt$P.Value, AveExpr = tt$AveExpr)
  out[, padj := p.adjust(pvalue, "BH")]
  for (k in names(extra)) out[[k]] <- extra[[k]]
  out[]
}

# ================= FASTING: GSE107787 =================
e <- el[contrast_id == "FASTING_primary"]
R <- ready("GSE107787"); genes <- intersect(U$main, rownames(R))
y <- log2(R[genes, e$sample_id] + 1)
zt <- factor(e$zt, levels = paste0("ZT", c(0, 4, 8, 12, 16, 20)))
st <- factor(ifelse(e$arm == "case", "FAST", "FED"), levels = c("FED", "FAST"))
X <- model.matrix(~ zt + st); stopifnot(colnames(X)[ncol(X)] == "stFAST")
fit <- eBayes(lmFit(y, X), trend = TRUE)
png(file.path(FG, "FAST_GSE107787_meanvar.png"), 900, 400); par(mfrow = c(1, 2))
plotSA(fit, main = "limma-trend: sqrt(sigma) vs mean log2(RPKM+1)")
hist(rowMeans(y), 80, main = "mean log2(RPKM+1), universe genes", xlab = ""); dev.off()
det <- rowSums(R[genes, e$sample_id] >= 1) >= 18
tab <- limma_table(fit, "stFAST", genes, "FAST", list(detected = det, in_universe_sensitivity = genes %in% U$sens))
tab[, `:=`(dataset = "GSE107787", contrast = "FAST_vs_FED_ZTadjusted", model = "limma-trend log2(RPKM+1) ~ ZT + state")]
fwrite(tab, file.path(RE, "FAST_GSE107787_FASTvsFED_ZTadj.csv.gz"))
# interaction
Xi <- model.matrix(~ zt * st)
fi <- eBayes(lmFit(y, Xi), trend = TRUE)
icol <- grep(":stFAST", colnames(Xi))
ti <- topTable(fi, coef = icol, number = Inf, sort.by = "none")[genes, ]
fdr_i <- p.adjust(ti$P.Value, "BH")
deg <- tab$padj < 0.05 & abs(tab$log2FC) >= 0.38
# per-ZT fasting effects (descriptive)
base <- which(colnames(Xi) == "stFAST")
perzt <- sapply(levels(zt), function(z) {
  cv <- rep(0, ncol(Xi)); cv[base] <- 1
  if (z != "ZT0") cv[which(colnames(Xi) == paste0("zt", z, ":stFAST"))] <- 1
  cv
})
cf <- contrasts.fit(fi, perzt); cf <- eBayes(cf, trend = TRUE)
pz <- data.table(gene = genes, cf$coefficients[genes, ], interaction_F_p = ti$P.Value, interaction_fdr = fdr_i)
setnames(pz, levels(zt), paste0("fast_effect_", levels(zt)))
fwrite(pz, file.path(P3, "results/fasting/GSE107787_fasting_effect_by_ZT.csv.gz"))
info$fasting <- list(n_fast = sum(st == "FAST"), n_fed = sum(st == "FED"), n_genes = length(genes),
                     n_padj05 = sum(tab$padj < .05), n_deg = sum(deg), n_interaction_fdr05 = sum(fdr_i < .05),
                     frac_deg_with_interaction_fdr05 = mean(fdr_i[deg] < .05),
                     design = "~ ZT + state (ref FED)", interaction_design = "~ ZT * state")
# clean subset
ec <- e[in_clean_sensitivity == TRUE]
yc <- log2(R[genes, ec$sample_id] + 1)
ztc <- factor(ec$zt); stc <- factor(ifelse(ec$arm == "case", "FAST", "FED"), levels = c("FED", "FAST"))
fc <- eBayes(lmFit(yc, model.matrix(~ ztc + stc)), trend = TRUE)
tc <- limma_table(fc, "stcFAST", genes, "FASTclean", list(detected = det))
tc[, `:=`(dataset = "GSE107787", contrast = "FAST_vs_FED_ZT0_ZT4_muscleLow", model = "limma-trend ~ ZT + state")]
fwrite(tc, file.path(RE, "FASTclean_GSE107787_ZT0ZT4_muscleLow.csv.gz"))
info$fasting_clean <- list(samples = ec$sample_id, n_fast = sum(stc == "FAST"), n_fed = sum(stc == "FED"),
                           n_padj05 = sum(tc$padj < .05))
pw <- pathway_tests(y, X, c(rep(0, ncol(X) - 1), 1), "FAST_GSE107787")
fwrite(pw, file.path(P3, "results/fasting/FAST_GSE107787_pathways.csv"))

# ================= FASTING secondary: GSE225695 =================
e2 <- el[contrast_id == "FASTING_secondary"]
C2 <- ready("GSE225695"); g2 <- intersect(U$main, rownames(C2))
cts <- C2[g2, e2$sample_id]; storage.mode(cts) <- "integer"
cd <- data.frame(state = factor(ifelse(e2$arm == "case", "Fast", "Fed"), levels = c("Fed", "Fast")), row.names = e2$sample_id)
dds <- DESeq(DESeqDataSetFromMatrix(cts, cd, ~ state), quiet = TRUE)
stopifnot(resultsNames(dds)[2] == "state_Fast_vs_Fed")
t2 <- deseq_table(dds, coef_name = "state_Fast_vs_Fed", genes = g2, n_detect = 3)
t2[, `:=`(dataset = "GSE225695", contrast = "Fast_vs_Fed", model = "DESeq2 ~ state")]
fwrite(t2, file.path(RE, "FAST2_GSE225695_FastvsFed.csv.gz"))
info$fasting2 <- list(n = ncol(dds), n_padj05 = sum(t2$padj < .05, na.rm = TRUE))

# ================= CIRCADIAN: harmonic regression =================
harmonic <- function(acc, cid, label) {
  e <- el[contrast_id == cid]
  R <- ready(acc); genes <- intersect(U$main, rownames(R))
  y <- log2(R[genes, e$sample_id] + 1)
  z <- as.numeric(sub("ZT", "", e$zt))
  cs <- cos(2 * pi * z / 24); sn <- sin(2 * pi * z / 24)
  X <- cbind(Intercept = 1, cos24 = cs, sin24 = sn)
  fit <- eBayes(lmFit(y, X), trend = TRUE)
  F <- topTable(fit, coef = 2:3, number = Inf, sort.by = "none")[genes, ]
  a <- fit$coefficients[genes, "cos24"]; b <- fit$coefficients[genes, "sin24"]
  s2 <- fit$s2.post[match(genes, rownames(fit))]
  XtXi <- solve(crossprod(X))
  va <- s2 * XtXi[2, 2]; vb <- s2 * XtXi[3, 3]; cab <- s2 * XtXi[2, 3]
  amp <- sqrt(a^2 + b^2)
  phase <- (atan2(b, a) * 24 / (2 * pi)) %% 24
  ga <- -b / amp^2; gb <- a / amp^2
  phase_se_h <- sqrt(ga^2 * va + gb^2 * vb + 2 * ga * gb * cab) * 24 / (2 * pi)
  amp_se <- sqrt((a / amp)^2 * va + (b / amp)^2 * vb + 2 * (a / amp) * (b / amp) * cab)
  # categorical sensitivity
  fz <- eBayes(lmFit(y, model.matrix(~ factor(z))), trend = TRUE)
  Fz <- topTable(fz, coef = 2:6, number = Inf, sort.by = "none")[genes, ]
  rh <- data.table(gene = genes, cos_coef = a, sin_coef = b, amplitude_log2 = amp, amplitude_se = amp_se,
                   peak_phase_ZT = phase, phase_se_h = phase_se_h, F_harmonic = F$F, p_harmonic = F$P.Value,
                   fdr_harmonic = p.adjust(F$P.Value, "BH"), p_categorical = Fz$P.Value,
                   fdr_categorical = p.adjust(Fz$P.Value, "BH"), AveExpr = F$AveExpr,
                   detected = rowSums(R[genes, e$sample_id] >= 1) >= 9)
  fwrite(rh, file.path(P3, "results/circadian", paste0(label, "_rhythm.csv.gz")))
  # signed time contrasts T1 = f(6)-f(2), T2 = f(12)-f(0)
  f <- function(t) c(0, cos(2 * pi * t / 24), sin(2 * pi * t / 24))
  CM <- cbind(T1_ZT6_minus_ZT2 = f(6) - f(2), T2_ZT12_minus_ZT0 = f(12) - f(0)); rownames(CM) <- colnames(X)
  cf <- eBayes(contrasts.fit(fit, CM), trend = TRUE)
  for (k in colnames(CM)) {
    tt <- limma_table(cf, k, genes, k, list(detected = rh$detected, in_universe_sensitivity = genes %in% U$sens))
    tt[, `:=`(dataset = acc, contrast = k, model = "limma-trend log2(RPKM+1) ~ cos + sin (24 h); predicted difference")]
    fwrite(tt, file.path(RE, paste0(label, "_", k, ".csv.gz")))
  }
  png(file.path(FG, paste0(label, "_meanvar.png")), 600, 400); plotSA(fit, main = label); dev.off()
  list(n = ncol(y), n_genes = length(genes), n_rhythmic_fdr05 = sum(rh$fdr_harmonic < .05),
       n_rhythmic_categorical_fdr05 = sum(rh$fdr_categorical < .05),
       both = sum(rh$fdr_harmonic < .05 & rh$fdr_categorical < .05), sex = unique(e$group)[1])
}
info$circadian_primary <- harmonic("GSE117134", "CIRCADIAN_primary", "CIRC_GSE117134_WTfemale")
info$circadian_male_sens <- harmonic("GSE107787", "CIRCADIAN_sens_male", "CIRCsens_GSE107787_FEDmale")

# ================= PPARalpha: GSE282998 =================
PP <- ready("GSE282998"); gp <- intersect(U$main, rownames(PP))
sm <- fread(file.path(P2, "datasets/GSE282998/sample_metadata.csv"), colClasses = c(sample_id = "character"))
ew <- el[contrast_id == "PPARA_primary"]; ek <- el[contrast_id == "PPARA_specificity_KO"]
dose_of <- function(ids) as.integer(sub(" mg/kg/day", "", sm$dose[match(ids, sm$sample_id)]))
fitp <- function(ids, design_df, design) {
  cts <- PP[gp, ids]; storage.mode(cts) <- "integer"
  DESeq(DESeqDataSetFromMatrix(cts, design_df, design), quiet = TRUE)
}
dw <- data.frame(dose = factor(dose_of(ew$sample_id), levels = c(0, 5, 10, 20)),
                 muscle_high = factor(ew$muscle_high), row.names = ew$sample_id)
ddw <- fitp(ew$sample_id, dw, ~ dose)
stopifnot(all(c("dose_5_vs_0", "dose_10_vs_0", "dose_20_vs_0") %in% resultsNames(ddw)))
for (d in c(5, 10, 20)) {
  nm <- paste0("dose_", d, "_vs_0")
  tt <- deseq_table(ddw, coef_name = nm, genes = gp)
  tt <- add_ashr(tt, ddw, coef_name = nm)
  tt[, `:=`(dataset = "GSE282998", contrast = paste0("WT_GW7647_", d, "v0"), model = "DESeq2 ~ dose (B6129SF2/J)")]
  fwrite(tt, file.path(RE, paste0("PPARA_GSE282998_WT_", d, "v0.csv.gz")))
}
diag_plots(ddw, deseq_table(ddw, coef_name = "dose_20_vs_0", genes = gp), file.path(FG, "PPARA_WT_diagnostics.png"), "PPARA WT")
info$ppara_size_factors <- as.list(setNames(sizeFactors(ddw), paste0(colnames(ddw), "_dose", dw$dose)))
# trend
dwn <- data.frame(dose10 = dose_of(ew$sample_id) / 10, row.names = ew$sample_id)
ddt <- fitp(ew$sample_id, dwn, ~ dose10)
tt <- deseq_table(ddt, coef_name = "dose10", genes = gp)
tt[, `:=`(dataset = "GSE282998", contrast = "WT_GW7647_trend_per10mgkg", model = "DESeq2 ~ dose (numeric/10)")]
fwrite(tt, file.path(RE, "PPARA_GSE282998_WT_trend.csv.gz"))
# muscle sensitivity
ddm <- fitp(ew$sample_id, dw, ~ muscle_high + dose)
tt <- deseq_table(ddm, coef_name = "dose_20_vs_0", genes = gp)
tt[, `:=`(dataset = "GSE282998", contrast = "WT_GW7647_20v0_muscleAdjusted", model = "DESeq2 ~ muscle_high + dose")]
fwrite(tt, file.path(RE, "PPARA_GSE282998_WT_20v0_muscleAdj.csv.gz"))
# genotype x dose
ids <- c(ew$sample_id, ek$sample_id)
dg <- data.frame(genotype = factor(c(rep("WT", nrow(ew)), rep("KO", nrow(ek))), levels = c("WT", "KO")),
                 dose = factor(dose_of(ids), levels = c(0, 5, 10, 20)), row.names = ids)
ddg <- fitp(ids, dg, ~ genotype + genotype:dose)
rn <- resultsNames(ddg); info$ppara_genotype_coefs <- rn
stopifnot(all(c("genotypeWT.dose20", "genotypeKO.dose20") %in% rn))
wt <- deseq_table(ddg, coef_name = "genotypeWT.dose20", genes = gp)
ko <- deseq_table(ddg, coef_name = "genotypeKO.dose20", genes = gp)
it <- deseq_table(ddg, contrast = list("genotypeWT.dose20", "genotypeKO.dose20"), genes = gp)
for (z in list(list(ko, "KO_20v0"), list(it, "WTminusKO_20v0_interaction"))) {
  z[[1]][, `:=`(dataset = "GSE282998", contrast = z[[2]], model = "DESeq2 ~ genotype + genotype:dose (WT B6129SF2/J + Ppara-null)")]
  fwrite(z[[1]], file.path(RE, paste0("PPARA_GSE282998_", z[[2]], ".csv.gz")))
}
cls <- ifelse(wt$padj < .05 & it$padj < .05 & sign(it$log2FC) == sign(wt$log2FC), "PPARa_dependent",
              ifelse(wt$padj < .05, "agonist_associated_dependence_unresolved", "not_responsive"))
fwrite(data.table(gene = gp, WT_20v0_joint = wt$log2FC, WT_padj = wt$padj, KO_20v0 = ko$log2FC, KO_padj = ko$padj,
                  interaction = it$log2FC, interaction_se = it$lfcSE, interaction_padj = it$padj, class = cls),
       file.path(P3, "results/ppara/genotype_classification.csv.gz"))
info$ppara <- list(n_wt = nrow(ew), n_ko = nrow(ek), n_genes = length(gp), class_counts = as.list(table(cls)))
pw <- pathway_tests(log2(counts(ddw, normalized = TRUE) + 1), model.matrix(~ dw$dose), c(0, 0, 0, 1), "PPARA_WT_20v0")
fwrite(pw, file.path(P3, "results/ppara/PPARA_WT_20v0_pathways.csv"))
write_json(info, file.path(P3, "results/other_refs_run_info.json"), auto_unbox = TRUE, pretty = TRUE)
cat("done\n")
