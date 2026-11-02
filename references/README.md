# References and data sources

## Primary publication

Jung W, Abushukair H, Patil NY, Ampadu F, Firouzi M, Rus I, Choi J, Muthukrishnan SD,
Shukla S, Tarantini S, Csiszar A, Mehla K, Kim D, Yoon J-H, Kim D, Kim J, Jung J,
Klementieva O, Shah YM, Yoshihara E, Singh PK, Joshi AD, Oh TG.
**Spatial transcriptome mapping identifies Ppara-Anxa2 cross-talk in microplastic-induced
hepatotoxicity.** *Science Advances* 12(25): eaec8681 (2026).

- doi: [10.1126/sciadv.aec8681](https://doi.org/10.1126/sciadv.aec8681)
- PubMed 42308314; PMCID PMC13274612
- Licence: CC BY 4.0
- Full text read via the Europe PMC REST API (`/PMC13274612/fullTextXML`) on 2026-11-02

## Data

| Resource | Accession / link | Used for |
|---|---|---|
| GEO series | [GSE294180](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294180) | Gene-level counts and sample metadata (this project) |
| BioProject | PRJNA1248512 | — |
| SRA | SRX28316092 – SRX28316110 | Raw FASTQ files (not used yet) |
| Xenium spatial data | [Zenodo 19603085](https://zenodo.org/records/19603085) | Not used yet |
| GENCODE mouse release M35 | [gencode.vM35.annotation.gtf.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/) | Gene biotypes and chromosomes |
| PS-microplastic liver RNA-seq (used in the paper) | [GSE245069](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245069) | Possible Phase 2 comparison |

Checksums and download dates: `metadata/data_provenance.json`.

## Methods references

- Anders S, Huber W. Differential expression analysis for sequence count data. *Genome Biology* 11:R106 (2010). The median-of-ratios size factors.
- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550 (2014).
- Soneson C, Love MI, Robinson MD. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences (tximport). *F1000Research* 4:1521 (2015).
- Patro R *et al.* Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods* 14:417–419 (2017).
