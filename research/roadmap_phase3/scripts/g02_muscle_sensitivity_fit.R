# POST-HOC EXPLORATORY (step 2/3): refit MASH1 and MASH2 on the pre-registered alternative universe.
# Same samples, design (~ diet, ref chow), no offset, independent filtering off. Outputs go to
# results/posthoc_audit/ only; reference_effects/ is untouched.
source("/mnt/c/Users/Soham/Desktop/microplastic/research/roadmap_phase3/scripts/common.R")
OUT <- file.path(P3, "results/posthoc_audit")
alt <- fread(file.path(OUT, "alt_universe.csv"))$vm35_symbol
el <- eligible()
for (z in list(c("GSE316023", "MASH_primary", "MASH1"), c("GSE246328", "MASH_secondary", "MASH2"))) {
  e <- el[contrast_id == z[2]]
  cts <- ready(z[1]); genes <- intersect(alt, rownames(cts))
  cts <- cts[genes, e$sample_id]; storage.mode(cts) <- "integer"
  cd <- data.frame(diet = factor(ifelse(e$arm == "case", "GAN", "chow"), levels = c("chow", "GAN")), row.names = e$sample_id)
  dds <- DESeq(DESeqDataSetFromMatrix(cts, cd, ~ diet), quiet = TRUE)
  stopifnot(resultsNames(dds)[2] == "diet_GAN_vs_chow")
  tab <- deseq_table(dds, coef_name = "diet_GAN_vs_chow", genes = genes)
  fwrite(tab, file.path(OUT, paste0(z[3], "_altuniverse.csv.gz")))
  cat(z[3], length(genes), "genes\n")
}
