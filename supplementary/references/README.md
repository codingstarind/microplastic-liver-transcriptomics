# References and dataset citations

## Report references

These are the references of the final report, in the report's numbering.

1. Djouina M. *et al.* (2023). Oral exposure to polyethylene microplastics induces inflammatory and metabolic changes and promotes fibrosis in mouse liver. *Ecotoxicology and Environmental Safety* 264:115417.
2. Jung W. *et al.* (2026). Spatial transcriptome mapping identifies Ppara–Anxa2 cross-talk in microplastic-induced hepatotoxicity. *Science Advances* 12:eaec8681. [doi:10.1126/sciadv.aec8681](https://doi.org/10.1126/sciadv.aec8681)
3. Li *et al.* (2025). Polystyrene microplastics induce liver fibrosis and lipid deposition in mice through three hub genes revealed by the RNA-seq. *Scientific Reports* 15:2583.
4. Saigal *et al.* (2026). AI digital pathology as a key tool providing in-depth understanding of the progression and regression of MASH and fibrosis in male mouse models. *Nature Communications* 17:6964.
5. Flensted-Jensen *et al.* (2024). Dietary intervention reverses molecular markers of hepatocellular senescence in the GAN diet-induced obese and biopsy-confirmed mouse model of NASH. *BMC Gastroenterology* 24:59.
6. Kinouchi K. *et al.* (2018). Fasting imparts a switch to alternative daily pathways in liver and muscle. *Cell Reports* 25:3299–3314.e6.
7. Heintz *et al.* (2025). Comparison of phenotypic and transcriptomic profiles between HFPO-DA and prototypical PPARα, PPARγ, and cytotoxic agents in wild-type and Ppara-null mouse livers. *Toxicological Sciences* 206:183–201.
8. Love M.I., Huber W., Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550.
9. Benjamini Y., Hochberg Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *JRSS B* 57:289–300.
10. Liberzon A. *et al.* (2015). The Molecular Signatures Database Hallmark gene set collection. *Cell Systems* 1:417–425.
11. Wu D. *et al.* (2010). ROAST: rotation gene set tests for complex microarray experiments. *Bioinformatics* 26:2176–2182.
12. Wu D., Smyth G.K. (2012). Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Research* 40:e133.
13. Foroutan M. *et al.* (2018). Single sample scoring of molecular phenotypes. *BMC Bioinformatics* 19:404.
14. Zou H., Hastie T. (2005). Regularization and variable selection via the elastic net. *JRSS B* 67:301–320.
15. Friedman J.H. (2001). Greedy function approximation: a gradient boosting machine. *Annals of Statistics* 29:1189–1232.
16. Prokhorenkova L. *et al.* (2018). CatBoost: unbiased boosting with categorical features. *NeurIPS* 31.
17. Holm S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics* 6:65–70.
18. Makdissi *et al.* (2026). Size-dependent plastic exposure disrupts macrophage function and tissue-specific metabolism. *Nature Metabolism* 8:2084–2101.

## Datasets

| Name in report | Accession | Role | Ref. |
|---|---|---|---|
| Main PE | [GSE294180](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294180) (BioProject PRJNA1248512) | 19 mouse livers, 2 × 2 diet × PE design; all main analyses | [2] |
| MASH-1 | [GSE316023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316023) | GAN diet vs chow; classifier training/testing; MASH reference signature | [4] |
| MASH-2 | [GSE246328](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246328) | GAN diet vs chow; classifier testing/training | [5] |
| FAST | [GSE107787](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107787) | Fasted vs fed reference signature | [6] |
| CLOCK | [GSE117134](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117134) | Circadian reference signature; source of the 20 co-expression modules | |
| PPARα | [GSE282998](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282998) | GW7647 vs control reference signature | [7] |
| Hong-PE | [GSE261893](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261893) | Independent PE-response transfer test | |

Other series audited in the research phases but not used for the report's results:
- **Phase 2:** GSE225695.
- **Phase 5:** GSE261725, GSE286116 and GSE262399 (Hong-lab PS, PP and PB series), GSE245069 and GSE277772.
- **Phase 5C:** PRJNA916285, which is incomplete and excluded from the report.

Details are in [`research/roadmap_phase2/`](../../research/roadmap_phase2/) and [`research/roadmap_phase5/dataset_provenance.csv`](../../research/roadmap_phase5/dataset_provenance.csv).

## Annotation and gene-set resources

| Resource | Version | Link |
|---|---|---|
| GENCODE mouse annotation | Release M35 (GRCm39) | https://www.gencodegenes.org/mouse/release_M35.html |
| MSigDB mouse collections, accessed with msigdbr 26.1.1 | 2026.1.Mm (50 Hallmark sets; SANDERSON_PPARA_TARGETS; WP_PPAR_SIGNALING_PATHWAY; REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA) | https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp |

## Further reading

- The full literature review is in [`research/literature/literature_review.md`](../../research/literature/literature_review.md), and its papers are in [`research/literature/papers/`](../../research/literature/papers/).
- Methods references from Phases 1 and 2 are in [`references/README.md`](../../references/README.md).
