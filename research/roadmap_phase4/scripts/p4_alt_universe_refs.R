# Variant V6: refit FAST, PPARA 20v0 and CIRC_T1 on the post-hoc alternative muscle-restricted universe
# (12,010 genes; research/roadmap_phase3/results/posthoc_audit/alt_universe.csv). MASH1 alt refit already exists there.
# Same models as Phase 3 d01_other_references.R. Writes ONLY to research/roadmap_phase4/results/alt_universe_refs/.
source("/mnt/c/Users/Soham/Desktop/microplastic/research/roadmap_phase3/scripts/common.R")
P4 <- file.path(REPO, "research/roadmap_phase4"); OUT <- file.path(P4, "results/alt_universe_refs")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
alt <- fread(file.path(P3, "results/posthoc_audit/alt_universe.csv"))$vm35_symbol
el <- eligible()
lt <- function(fit, coef, genes) {
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")[genes, ]
  data.table(gene = genes, log2FC = tt$logFC, SE = tt$logFC / tt$t, pvalue = tt$P.Value, padj = p.adjust(tt$P.Value, "BH"))
}
# fasting
e <- el[contrast_id == "FASTING_primary"]; R <- ready("GSE107787"); g <- intersect(alt, rownames(R))
y <- log2(R[g, e$sample_id] + 1)
X <- model.matrix(~ factor(e$zt) + factor(ifelse(e$arm == "case", "FAST", "FED"), levels = c("FED", "FAST")))
fwrite(lt(eBayes(lmFit(y, X), trend = TRUE), ncol(X), g), file.path(OUT, "FAST.csv.gz"))
# circadian T1
e <- el[contrast_id == "CIRCADIAN_primary"]; R <- ready("GSE117134"); g <- intersect(alt, rownames(R))
y <- log2(R[g, e$sample_id] + 1); z <- as.numeric(sub("ZT", "", e$zt))
X <- cbind(Intercept = 1, cos24 = cos(2 * pi * z / 24), sin24 = sin(2 * pi * z / 24))
f <- function(t) c(0, cos(2 * pi * t / 24), sin(2 * pi * t / 24))
CM <- cbind(T1 = f(6) - f(2)); rownames(CM) <- colnames(X)
fwrite(lt(eBayes(contrasts.fit(lmFit(y, X), CM), trend = TRUE), "T1", g), file.path(OUT, "CIRC_T1.csv.gz"))
# PPARa WT 20v0
PP <- ready("GSE282998"); g <- intersect(alt, rownames(PP))
sm <- fread(file.path(P2, "datasets/GSE282998/sample_metadata.csv"), colClasses = c(sample_id = "character"))
ew <- el[contrast_id == "PPARA_primary"]
dose <- factor(as.integer(sub(" mg/kg/day", "", sm$dose[match(ew$sample_id, sm$sample_id)])), levels = c(0, 5, 10, 20))
cts <- PP[g, ew$sample_id]; storage.mode(cts) <- "integer"
dds <- DESeq(DESeqDataSetFromMatrix(cts, data.frame(dose = dose, row.names = ew$sample_id), ~ dose), quiet = TRUE)
r <- results(dds, name = "dose_20_vs_0", independentFiltering = FALSE)[g, ]
fwrite(data.table(gene = g, log2FC = r$log2FoldChange, SE = r$lfcSE, pvalue = r$pvalue, padj = p.adjust(r$pvalue, "BH")),
       file.path(OUT, "PPARA_20v0.csv.gz"))
cat("alt refs done\n")
