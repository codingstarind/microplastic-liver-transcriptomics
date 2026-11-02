# Exploratory data analysis: descriptive statistics, PCA and clustering

Scripts:

- `scripts/05_descriptive_eda.py` → figures 1–9 in `figures/eda/`
- `scripts/06_pca_clustering.py` → figures 10–15 in `figures/pca_clustering/`

Figures are saved as PNG (slides) and PDF (vector). Tables are in `results/tables/`.

**No hypothesis tests were run.** Statements about groups describe these 19 mice; they are not inferences. Group colours are the same in every figure: CD+Veh blue, CD+PE orange, MD+Veh green, MD+PE yellow. Triangles mark PE-gavaged mice and circles vehicle mice.

---

## Fig. 1: Library sizes and size factors (`fig01_library_sizes`)

1. **Axes.**
   - A: one bar per sample; y = total gene-assigned counts (millions).
   - B: x = library size / mean library size; y = median-of-ratios size factor; dotted line y = x.
2. **Why.** Raw counts scale with sequencing depth, so depth must be known before any comparison. B shows whether the robust size factor agrees with simple total-count scaling.
3. **Observation.**
   - Libraries range from 10.4 to 15.5 M.
   - **Every PE sample has a smaller library (10.4–12.8 M) than every vehicle sample (13.4–15.5 M).**
   - Size factors (0.78–1.26) track library size closely (Spearman ρ = 0.95), so no sample has a grossly distorted composition.
4. **Alternative explanations.**
   - PE libraries (S48–S57) were prepared or pooled after the vehicle libraries (S39–S47): a processing batch or pooling imbalance.
   - PE changes RNA yield or composition.
   - Chance (unlikely given a perfect split, but n is small).
   - This is a technical property, **not** a biological expression change.
5. **Next.** Check SRA run metadata (flowcell, lane, date). In Phase 2 models, check whether the PE effects are robust to adding library size as a covariate.

## Fig. 2: Distribution of counts (`fig02_count_distributions`)

1. **Axes.**
   - A: x = log10(mean count across 19 samples + 1) per gene; y = number of genes. The dashed line marks a mean of 10.
   - B: x = log2(count + 1) for genes with count > 0; y = density; one line per sample.
2. **Why.** To see the shape of expression before choosing filters and transforms.
3. **Observation.**
   - The distribution is extremely right-skewed and zero-inflated. 53% of the 55,836 genes are zero everywhere, a second large group has means below 10, and the expressed genes form a broad hump between ~10 and ~10⁴ counts with a long tail up to ~5 × 10⁵ (*Alb*).
   - Per-sample shapes (B) are very similar. No sample has an abnormal distribution.
4. **Alternative explanations.** The all-zero genes reflect the annotation (non-liver genes, small RNAs not captured by poly(A) selection, pseudogenes), not failed samples.
5. **Next.** This justifies the low-expression filter (≥ 10 counts in ≥ 4 samples; 13,777 genes kept, 99.86% of reads retained) and a log-type transform.

## Fig. 3: Per-sample boxplots before and after normalisation (`fig03_sample_boxplots_before_after_normalisation`)

1. **Axes.** x = sample; y = log2 expression of the 13,777 filtered genes.
   - A: log2(count + 1).
   - B: log2(normalised count + 1).
   - Boxes show median and IQR; whiskers 1.5 × IQR.
2. **Why.** The standard check that normalisation aligns the bulk of the distributions. Most genes should not change between conditions.
3. **Observation.**
   - Before normalisation, PE samples have visibly lower medians (the depth effect from Fig. 1).
   - After median-of-ratios normalisation all 19 medians and IQRs line up.
4. **Alternative explanations.** Aligned medians do not prove that normalisation is correct for every gene. If PE truly shifted many genes in one direction, median-based methods would partly absorb it.
5. **Next.** Accept median-of-ratios for exploration. Compare with TMM (edgeR) in Phase 2.

## Fig. 4: How concentrated are the reads? (`fig04_expression_dominance`)

1. **Axes.**
   - A: x = gene rank within sample (log scale); y = cumulative fraction of the library; one line per sample.
   - B: mean % of library for the 15 most abundant genes.
2. **Why.** The question was whether some genes are expressed far more strongly than others. Liver is known to be dominated by a few secreted-protein transcripts, which affects total-count normalisation.
3. **Observation.**
   - **The top gene (*Alb*) takes 3.3–6.7% of each library.** The top 10 genes take ~19%, the top 100 ~47% and the top 1,000 ~74%.
   - Abundant genes include *Alb*, *Apoe*, *Scd1*, *Trf*, *Serpina1b*, *Ahsg*, *Apoa1*, several major urinary proteins (*Mup3/7/19/20/22*) and mitochondrial transcripts (*mt-Co1*, *mt-Cytb*, *mt-Nd1/2/4/5*).
   - MD samples devote a smaller share to *Alb* (3.3–4.5%) than CD samples (4.7–6.7%).
4. **Alternative explanations.** Composition shifts like this (diet changing a few dominant genes) are exactly why total-count / CPM scaling can mislead and why median-of-ratios is preferred. The *Mup* genes are highly variable between individual male mice (hormonal, circadian), so they add noise.
5. **Next.** Keep median-of-ratios. Watch *Mup* genes in Phase 2 results.

## Fig. 5: Mean–variance relationship (`fig05_mean_variance`)

1. **Axes.**
   - A: x = mean normalised count; y = variance across the 19 samples (both log10). The dashed line is the Poisson expectation, variance = mean.
   - B: x = mean log2(normalised + 1); y = SD of log2 expression.
2. **Why.** To decide what transformation PCA, clustering and later tests need.
3. **Observation.**
   - On the count scale the variance rises steeply with the mean, and well above the Poisson line for expressed genes. This is **overdispersion**: biological variation between mice adds to sampling noise. It is the reason to use negative-binomial models (DESeq2, edgeR) rather than Poisson or *t*-tests on raw counts.
   - After log2, the SD is roughly flat for well-expressed genes, but **low-count genes (mean log2 < ~4) have inflated SD**, an artefact of the log at small values.
4. **Alternative explanations.** Some high-variance genes are real biology (diet-responsive genes) or contamination (muscle genes in a few samples), not noise.
5. **Next.** Use a variance-stabilising transform (DESeq2 VST) in Phase 2 and compare with log2(x + 1). Restrict PCA to well-expressed / most variable genes, as done here.

## Fig. 6: Per-sample quality indicators (`fig06_sample_qc_by_group`)

1. **Axes.** x = group; y = metric. Points are samples and the bar is the group median. Metrics: genes with ≥ 10 counts; fraction of counts from chrM genes; fraction of counts from the single top gene.
2. **Why.** To look for unusual samples and group-level technical differences.
3. **Observation.**
   - No sample is a gross failure. Every sample detects 12,176–13,878 genes, and the chrM fraction is 7.4–11.9%.
   - **PE groups have slightly fewer detected genes and a higher chrM fraction** (CD+PE median 0.116 vs CD+Veh 0.088).
   - MD groups have a lower top-gene (*Alb*) share.
4. **Alternative explanations.**
   - Detection counts partly track library depth, which is lower in PE samples (Fig. 1).
   - A higher mitochondrial fraction can come from RNA quality or extraction batch, or from real biology. The paper reports *mt-Nd6* up with PE, and mitochondrial changes are plausible under oxidative or toxic stress.
   - The design cannot separate these.
5. **Next.** Treat mitochondrial genes with caution when interpreting PE effects. Compare with the paper's claims.

## Fig. 7: Selected genes by group (`fig07_selected_genes_by_group`)

1. **Axes.** x = group; y = log2(normalised count + 1); points = mice; bar = group mean.
   - The 12 genes were chosen **before looking at the data**: genes highlighted by the paper plus textbook steatosis, inflammation and fibrosis markers.
2. **Why.** To see whether the groups show obvious differences in biologically interpretable genes, and to compare with the published claims. Descriptive only.
3. **Observations** (group-mean differences in log2 units; Δ = difference):
   - **Diet effects are large and consistent.** Under MD: *Cidea* +5.2, *Cd36* +2.4, *Plin4* +2.4, *Lcn2* +2.0, *Anxa2* +1.5, *Col1a1* +1.1.
   - **The lipogenic genes go down with both MD and PE**, and lowest in MD+PE:

     | Gene | PE in CD | PE in MD |
     |---|---|---|
     | *Fasn* | −1.6 | −1.5 |
     | *Pnpla3* | −3.3 | −2.3 |
     | *Acly* | −1.8 | −1.1 |
     | *Me1* | −1.2 | −0.8 |

     This agrees with the paper.
   - **PE raises *Plin4*** (+1.0 in CD, +1.1 in MD), ***Cyp4a14*** (a PPARα target: +2.6 in CD, +1.7 in MD) and ***Nr1i3*** (+0.4 / +0.6).
   - ***Ppara* itself changes little** (+0.15 / +0.06; overlapping points). That is weaker than the paper's emphasis implies.
   - ***Anxa2*** responds to diet, not PE, in bulk (−0.4 in CD, 0.0 in MD). The paper's *Anxa2* finding is cell- and zone-specific (Xenium), so a weak bulk signal is not a contradiction.
4. **Alternative explanations.**
   - **Lipogenic genes are among the most circadian and feeding-responsive in liver.** The clock-output gene *Dbp* is lower in CD+PE than CD+Veh (−1.2), so differences in time of sacrifice or recent food intake between groups could produce the same pattern.
   - Food intake was measured but is not reported per group in the data we have.
   - n = 4–5 per group; single mice can move a mean (e.g. one CD+PE outlier for *Ppara*).
5. **Next.** Check clock genes systematically (*Dbp*, *Per1–3*, *Nr1d1/2*, *Bmal1*) across groups. In Phase 2, test diet, PE and their interaction formally with a count model.

## Fig. 8: How much variation lines up with the design? (`fig08_variance_share_by_factor`)

1. **Axes.** x = η², the share of each gene's variance across 19 samples explained by grouping on diet, PE or the 4 groups; y = number of genes.
   - The blue outline uses the actual labels. The grey bars use the same labels shuffled once: what "no relationship" looks like with this sample size.
   - 13,164 genes with mean log2 > 3.
2. **Why.** A single summary of whether the groups differ systematically, without picking genes and without p-values.
3. **Observation.**

   | Factor | Genes with η² > 0.5 | Same, shuffled labels |
   |---|---|---|
   | Diet | 984 | 3 |
   | PE | 252 | 11 |
   | 4 groups | 2,871 | 134 |

   **Diet has a strong, widespread footprint; PE has a smaller but clearly non-random one.**
4. **Alternative explanations.**
   - The PE signal is not purely biological. PE is confounded with library order and depth (Fig. 1), so part of it may be technical.
   - With 4 groups and 19 samples, η² is biased upwards (it never reaches 0 even with shuffled labels), which is why the shuffled reference matters.
   - A single shuffle is a visual reference, not a test.
5. **Next.** Use a formal permutation or model-based approach in Phase 2, after the team agrees hypotheses.

## Fig. 9: Tissue and cell-type marker screen (`fig09_tissue_marker_screen`)

1. **Axes.** x = sample; y = summed CPM of each marker set, + 1, log scale. Sets: skeletal muscle, adipocyte (*Adipoq*, *Lep*, *Retn*, *Plin1*), erythrocyte, exocrine pancreas, cholangiocyte, Kupffer/macrophage, stellate/fibrogenic.
2. **Why.** Added after the PCA loadings (Fig. 11) showed PC1 was made of muscle genes. Liver pieces can carry diaphragm, fat or blood.
3. **Observation.**
   - **Skeletal-muscle reads:**

     | Sample | Muscle CPM |
     |---|---|
     | MD_1 | 5,190 |
     | CD_MP_1 | 539 |
     | MD_2 | 316 |
     | MD_3 | 161 |
     | CD_2 | 58 |
     | All others | ≤ 20 |

     MD_1 also has traces of adipocyte genes. Pancreas markers are negligible everywhere.
   - CD_MP_2 is high in Kupffer and erythrocyte markers.
   - MD_MP_2 and MD_MP_4 are high in cholangiocyte markers, and MD_MP_4 in stellate markers.
   - Kupffer and stellate markers are generally higher in MD groups.
4. **Alternative explanations.**
   - Muscle genes are not expressed in hepatocytes, so dissection carry-over is by far the most likely source.
   - The other differences could be real biology (MASH diet → inflammation, ductular reaction, early fibrosis) or simply which part of the lobe was sampled.
   - The first version of the adipocyte set included *Cfd* (adipsin). *Cfd* turned out to be ~100× higher in **every** MD sample while *Adipoq*/*Lep*/*Plin1* stayed at ~0, so it reflects diet-induced expression in liver, not fat contamination. It was removed from the set.
5. **Next.** The team should decide how to handle MD_1 (exclude, or add a muscle covariate). Cell-composition differences mean bulk "PE effects" could be composition effects. A deconvolution against a liver single-cell reference would be a Phase 2 option.

---

## PCA and clustering

### Fig. 10: Primary PCA (`fig10_pca_primary`)

1. **Axes.**
   - A: scree plot, % variance per PC with the cumulative line.
   - B: PC1 vs PC2 scores.
   - C: PC3 vs PC4 scores.
   - Axis labels give each PC's % of variance.
2. **Features.**
   - The **500 most variable genes** of log2(median-of-ratios normalised + 1), after the low-expression filter.
   - Genes are mean-centred but not scaled (the DESeq2 `plotPCA` convention).
   - Samples are observations and genes are features.
   - **Why reduce dimensions:** 19 samples × 13,777 genes cannot be inspected directly. PCA finds the few directions that carry most of the sample-to-sample variation, so structure (groups, outliers, batches) becomes visible in 2D.
3. **Variance explained:** PC1 30.3%, PC2 17.0%, PC3 8.5%, PC4 6.2%, PC5 4.8%. The first 5 PCs together explain 67%.
4. **Observations.**
   - **PC1 is not the experimental design.** It is driven by the muscle-contaminated samples: MD_1 at +66, MD_2 +33, MD_3 +27, CD_MP_1 +25, CD_2 +15. Its correlation with the muscle score is Spearman ρ = 0.79, and group explains only 37% of its spread (η²).
   - **PC2 is the design.** η² for group is 0.97. The groups are ordered without overlap:

     | Group | PC2 scores |
     |---|---|
     | CD+Veh | +23 to +29 |
     | CD+PE | +4 to +8 |
     | MD+Veh | −1 to −13 |
     | MD+PE | −18 to −24 |

     **On this axis PE moves samples in the same direction as the MASH diet, on both diets**, and MD+PE is the most extreme group.
   - PC3 is linked to library size (ρ = 0.66). On PC4, MD_MP_2 is extreme.
5. **Separation by diet, PE, both or neither?** Both, on PC2. Diet is the larger effect. PE separates within each diet (CD+Veh vs CD+PE does not overlap; MD+Veh vs MD+PE barely does).
6. **Outliers.** No sample fails QC, but these stand apart:
   - MD_1 (contamination);
   - CD_MP_1 and CD_2 (contamination);
   - MD_MP_2 (PC4; also high in cholangiocyte markers).

   In the PC1–4 outlier screen, these four have the largest distances to their group centroid (1.8–2.7 × the group median). MD_1 has the largest absolute distance (44.6).
7. **Can batch effects or normalisation explain it?**
   - Contamination clearly shapes PC1.
   - PC2 also correlates with library index (ρ = −0.72), because PE libraries come after vehicle libraries within each diet. However, the PC2 order (CD+Veh → CD+PE → MD+Veh → MD+PE) is **not** the library order (CD → MD → CD+PE → MD+PE). So library order alone does not explain PC2, but some technical contribution cannot be excluded.
   - The sensitivity analyses below show how robust PC2 is.

### Fig. 11: Genes driving PC1 and PC2 (`fig11_pca_top_loadings`)

1. **Axes.** Bars = the 6 most negative and 6 most positive gene loadings for PC1 (A) and PC2 (B).
2. **Why.** To find out what each axis means biologically.
3. **Observation.**
   - **PC1 (+):** *Myh1*, *Mb*, *Acta1*, *Ttn*, *Atp2a1*, *Myh2*, all skeletal-muscle genes. This is what flagged the contamination.
   - **PC2 (−, toward MD/PE):** *Cfd*, *Cidea*, *Cidec*, *Themis*, *H3c14*. Lipid-droplet and steatosis genes.
   - **PC2 (+, toward CD+Veh):** *Chrna4*, *Sqle*, *Fabp5*, *Rps6-ps4*, *Rps12-ps3*, *Gm11518*. Cholesterol synthesis (*Sqle*) is suppressed by dietary cholesterol.
4. **Alternative explanations.** The ribosomal pseudogenes (*Rps6-ps4*, *Rps12-ps3*) among the top PC2 loadings are suspicious. Pseudogene counts from Salmon can reflect multi-mapping and may be sensitive to library quality, which could tie them to the PE / library-order confound.
5. **Next.** Check whether PC2 is robust when pseudogenes are excluded. Run gene-set analysis of PC2 loadings in Phase 2.

### Fig. 12 and 13: PCA sensitivity analyses (`fig12_pca_all_filtered_genes`, `fig13_pca_muscle_genes_removed`)

Summarised from `results/tables/pca_variance_explained.csv` and `pca_pc_associations.csv`:

| Feature set | PC1 % | PC2 % | Which PC carries the design (η² group) | Muscle ρ with PC1 |
|---|---|---|---|---|
| Top 500, log2 normalised (primary) | 30.3 | 17.0 | PC2 (0.97) | 0.79 |
| Top 1,000 | 26.3 | 17.0 | PC2 (0.96) | – |
| Top 2,000 | 23.5 | 16.9 | PC2 (0.88) | – |
| All 13,777 filtered genes | 20.7 | 16.2 | PC1 (0.62) and PC2 (0.71) | 0.59 |
| Top 500, log2 CPM | 34.4 | 24.9 | PC2 (0.86) | 0.69 |
| **Top 500 after removing 135 muscle-tracking genes** | 20.8 | 17.9 | **PC1 (0.97; diet 0.80)** | 0.18 |
| Top 500 without MD_1 | 20.8 | 19.5 | **PC1 (0.96)** | 0.07 |

- **Robust:** the design-aligned axis (diet strongly, PE within diet) appears under every choice. Once contamination is removed (genes or sample), it becomes PC1.
- **Not robust:** the size of PC1 and its meaning depend entirely on how contamination is handled. Normalisation (CPM vs median-of-ratios) changes the percentages but not the picture.
- In Fig. 13 the four groups again form distinct regions on PC1/PC2. MD_1, MD_2 and MD_3 still sit high on PC2: residual contamination, or a real difference within MD+Veh.
- **Caveat:** the "muscle-tracking" list (|Spearman ρ| ≥ 0.7 with the muscle score) contains clear muscle genes but also some liver genes (*Bhmt*, *Atg14*) correlated by chance or through diet, since 3 of the 4 contaminated samples are MD+Veh. It is a crude filter for exploration only.

### Fig. 14: Sample correlation heatmap (`fig14_sample_correlation_heatmap`)

1. **Axes.** Rows and columns = samples, ordered by the clustering tree. Colour = Pearson correlation of log2 normalised expression over the 500 most variable genes. Side bars show group, diet (grey shade) and PE (grey shade).
2. **Distance and linkage.**
   - Distance = 1 − Pearson correlation. It compares the *shape* of expression profiles and ignores overall level, which suits normalised log data.
   - Linkage = average (UPGMA): the distance between two clusters is the mean of all pairwise distances. It is less sensitive to single outliers than single linkage and less forced into equal-sized blobs than Ward.
3. **Observation.**
   - The four muscle-contaminated samples (MD_1–3, CD_MP_1) form their own branch, joining the rest last.
   - The remaining samples split into a CD branch (CD+Veh and CD+PE sub-clusters) and an MD branch (MD+Veh and MD+PE mixed).
   - All between-sample correlations are high (mean correlation to other samples is 0.933–0.965, computed over all 13,777 filtered genes). The groups differ in a modest number of genes against a shared liver background.
4. **Alternative explanations.** The top branch reflects contamination, not treatment (see Fig. 15).

### Fig. 15: Dendrograms under different choices (`fig15_dendrograms_sensitivity`)

Four trees, each over the same 500 genes:

- 1 − Pearson, average linkage (primary);
- 1 − Spearman, average linkage;
- Euclidean distance, Ward linkage;
- 1 − Pearson, average linkage, after removing the muscle-tracking genes.

Leaf labels are coloured by group. **No number of clusters was imposed.** Agreement of tree cuts at k = 2–6 with the known labels is in `clustering_label_agreement.csv` (adjusted Rand index, ARI: 0 = chance, 1 = identical).

- **Primary tree:** at k = 2 it splits off the 4 contaminated samples (ARI with group 0.05). At k = 4: diet ARI 0.65, group 0.46.
- **Muscle genes removed:** the **k = 2 cut is exactly diet** (ARI = 1.00). At k = 3, MD_1–3 separate from the other MD samples.
- **Euclidean/Ward and Spearman:** group ARI peaks at 0.65 (k = 5). CD+Veh and CD+PE consistently form separate sub-clusters. MD+Veh and MD+PE are more intermixed.
- **PE alone is never a top-level split** (ARI with PE ≤ 0.37 for all trees and k). PE is a secondary effect nested within diet.

**Interpretation:** unsupervised clustering recovers diet cleanly once contamination is handled. PE is visible as a within-diet sub-structure, clearer under the control diet than the MASH diet.

---

## Summary of what normalisation or transformation is needed, and why

- **Library-size normalisation is essential.** Depth varies 1.5-fold and is confounded with PE. Median-of-ratios is appropriate and was applied.
- **A variance-stabilising transform is needed** for PCA and clustering, because count variance grows with the mean (Fig. 5). log2(x + 1) on filtered genes works for exploration; DESeq2's VST is the better choice and is planned.
- **For hypothesis tests (Phase 2):** use count models (negative binomial: DESeq2 or edgeR) on the original integer counts with size factors. Do not use *t*-tests on log values of low-count genes.
- **Contamination must be handled explicitly** before any gene-level conclusions: a covariate, exclusion, or at least sensitivity analysis.
