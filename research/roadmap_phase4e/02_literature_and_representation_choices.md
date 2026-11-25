# Literature and representation choices

**Sources** (verified from full text in Phases 4B/4D unless noted):
- singscore: Foroutan 2018, BMC Bioinformatics 19:404 (within-sample rank scores, stable at small n);
- GSVA: Hänzelmann 2013, BMC Bioinformatics 14:7 (depends on the sample population, so unsuitable for 8 v 8 cross-study use; not run, not installed);
- Hallmark: Liberzon 2015, Cell Systems 1:417;
- pathway vs gene features with cross-dataset testing: Lee 2008, PLoS Comput Biol 4:e1000217;
- cross-study validation: Bernau 2014, Bioinformatics 30:i105 (cited, not re-read);
- selection bias and leakage: Ambroise & McLachlan 2002, PNAS 99:6562; Varma & Simon 2006, BMC Bioinformatics 7:91 (cited);
- models: logistic regression (Cox 1958), ridge/L2 (Hoerl & Kennard 1970), CART/Gini (Breiman 1984), entropy splits (Quinlan 1986), random forest (Breiman 2001), naive Bayes at small n (Ng & Jordan 2002).

## Representations
All features are per-mouse, within-sample rank scores over the fixed common gene universe.

| Code | Features | Definition | Count |
|---|---|---|---|
| **R0** | signed MASH program score | singscore of up-50 / down-50 genes chosen by Welch t on log2CPM **within the training mice of the current split only** (\|mean difference\| ≥ 1, ranked by t) | 1 |
| **R1** | curated pathway scores | singscore (up-set form) of 12 frozen external sets (Phase 3 `gene_sets.csv`): LIPOGENESIS_CORE, CHOL_SYNTHESIS_CORE, SANDERSON_PPARA_TARGETS, and HALLMARK_ FATTY_ACID_METABOLISM, CHOLESTEROL_HOMEOSTASIS, INFLAMMATORY_RESPONSE, EPITHELIAL_MESENCHYMAL_TRANSITION, XENOBIOTIC_METABOLISM, BILE_ACID_METABOLISM, PEROXISOME, OXIDATIVE_PHOSPHORYLATION, MTORC1_SIGNALING. **No data-driven selection** | 12 |
| **R2** | R0 + R1 | compact combination (program score + pathways) | 13 |
| **PCA** (comparator) | PC1–PC3 | PCA of per-sample z-scored log2CPM, **fitted on training mice only**; test mice projected | 3 |

R2 deliberately has no interaction products: with 15–16 training mice, extra products would exceed the information available.
