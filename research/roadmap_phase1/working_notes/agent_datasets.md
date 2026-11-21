# Feasibility screen: public transcriptomic datasets related to GSE294180 (microplastics x MASH, mouse liver)

Screen date: 2026-11-21. Reference study: GSE294180 (Jung et al. 2026 Sci Adv, "Spatial transcriptome mapping identifies Ppara-Anxa2 cross-talk in microplastic-induced hepatotoxicity", doi 10.1126/sciadv.aec8681, PMID 42308314). The design is 2x2: control diet D09100304 vs MASH diet D09100310, crossed with vehicle vs PE 10-150 um at 2 mg/day by gavage for 8 weeks, male C57BL/6, n = 4/5/5/5.

## How the records were read

- GEO SOFT text was fetched for each series (`acc.cgi?acc=GSExxx&targ=self&form=text&view=brief`) and for its samples (`...&targ=gsm&form=text&view=brief`). All group sizes below were counted from the GSM characteristics and titles, not taken from the paper text.
- Supplementary file lists and sizes came from the anonymous FTP listing (`ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSExxxx/suppl/`) and from `filelist.txt`. The HTTPS directory listing returned 403, and the GEO download page returned a reCAPTCHA. That CAPTCHA was not bypassed.
- File headers were read by streaming each file through `gunzip | head` over HTTPS. No file was saved to disk, and nothing over 50 MB was opened beyond its first bytes.
- Caveat: streaming over plain FTP from this Windows machine corrupted the gzip streams. The GSE245069 file first looked corrupt, then passed `gzip -t` over HTTPS. Integrity checks were therefore redone over HTTPS.
- Linked papers were checked through PubMed E-utilities (efetch/esummary), Europe PMC REST and PMC. The PMC and Nature full text for the Nat Metab paper was blocked (reCAPTCHA or login redirect), so only its abstract was read.
- Searches run: NCBI esearch on db=gds with several queries (microplastic/nanoplastic + liver + Mus/Rattus; polymer names + liver + rodentia + HTS; microplastic + HFD/NAFLD/MASH; GAN/AMLN/Gubra; NASH + fructose + cholesterol; diet catalogue codes D09100310 etc.). Europe PMC was also searched.

---

## 1. GSE245069: PS-MPs (3 sizes), ICR mice, drinking water

- **URLs read:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245069 (SOFT self + gsm); ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245069/suppl/; https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245069/suppl/GSE245069_gene_expression.txt.gz (streamed header); PubMed 39833454; https://pmc.ncbi.nlm.nih.gov/articles/PMC11747203/ (methods).
- **Paper:** Li Q, Zhu K, et al. "Polystyrene microplastics induce liver fibrosis and lipid deposition in mice through three hub genes revealed by the RNA-seq." *Sci Rep* 2025;15:2583. doi:10.1038/s41598-025-86810-5. PMID 39833454, PMC11747203 (open access).
- **Organism:** Mus musculus, **ICR** (outbred), male. Age: GEO says "Six-week-old" but the paper methods say "Five-week-old" (**discrepancy**).
- **Tissue:** liver.
- **Polymer / size:** polystyrene at 80 nm, 0.5 um and 5 um (Baseline ChromTech, Tianjin). Surface modification is not stated.
- **Dose / route / duration:** 1 mg/L in drinking water, ad libitum, for 12 weeks.
- **Diet:** not specified ("free access to their diets"). There is no metabolic-stress factor.
- **Groups (from GSM):** H2O control n=4, 80 nm n=4, 0.5 um n=4, 5 um n=4 (16 samples). The paper had 10 mice per group, of which 4 were sequenced.
- **Controls:** plain water (no surfactant vehicle). One H2O control group is shared by all three size arms within this series.
- **Platform:** GPL28330, DNBSEQ-T7 (BGI). Assembly GRCm38.p6. The processing text mentions limma, which is unusual for RNA-seq; the quantifier is not stated (UNVERIFIED, BGI's pipeline is typically RSEM).
- **Processed data:** `GSE245069_gene_expression.txt.gz` (2,182,937 B; 20,925 genes). For each sample it has columns `fpkm_*`, `tpm_*` and `read_count_*`, plus group averages. The read_count columns are **mostly integer, but about 25k of roughly 335k values are non-integer**, which is consistent with RSEM expected counts. They would need rounding for DESeq2/edgeR. Library sizes are about 15.8-17.0 M per sample. A second file, `GSE245069_rna_expression_.txt.gz` (4.15 MB), holds transcript-level FPKM/TPM. Per-sample GSM supplementary files: NONE.
- **SRA:** yes (SRX22058298-SRX22058313), BioProject PRJNA1026899.
- **Overlap:** none detected. It is not a SuperSeries.
- **Released:** 23 Jan 2025.
- **Oddities:** age discrepancy (5 vs 6 weeks); limma for counts; diet unspecified; outbred strain.

## 2. GSE277772: PS-MPs, C57BL/6 male, oral (twice weekly)

- **URLs read:** GEO SOFT self + gsm; FTP suppl listing + filelist.txt; streamed https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8529nnn/GSM8529295/suppl/GSM8529295_MC_4.txt.gz; PubMed 39956511 abstract.
- **Paper:** Chiu HW, Chu CW, Huang CC, Chia ZC, Wang YL, Lee YH. "Polystyrene microplastics induce hepatic lipid metabolism and energy disorder by upregulating the NR4A1-AMPK signaling pathway." *Environ Pollut* 2025;369:125850. doi:10.1016/j.envpol.2025.125850. PMID 39956511. Paywalled; only the abstract was read.
- **Organism:** C57BL/6, male, 6 weeks old at start.
- **Tissue:** liver.
- **Polymer:** polystyrene. **Particle size UNVERIFIED**: it is not in the GEO record, and the abstract only cites "0.08-100 um" from earlier literature.
- **Dose / route / duration:** "orally dosed with 0.4 mg/day PS-MPs twice weekly for eight weeks". Whether "0.4 mg/day" means per dosing day is ambiguous. Gavage vs drinking is UNVERIFIED; GEO says "orally".
- **Diet:** not specified. There is no metabolic factor.
- **Groups (from GSM):** water control n=3 (MC_4, MC_5, MC_6) and PS-MP n=3 (MH_4, MH_5, MH_6). The sample labels start at _4, which suggests a subset of a larger cohort (UNVERIFIED).
- **Controls:** water, given orally on the same schedule.
- **Platform:** GPL24247, Illumina NovaSeq 6000. Processing: Trimmomatic, HISAT2, featureCounts. The processing text says "Assembly: GRCh38", which is a **template error**: the features are mouse Ensembl IDs (ENSMUSG).
- **Processed data:** only per-sample `GSM85292xx_M?_?.txt.gz` (inside `GSE277772_RAW.tar`, 2.0 MB). The values I read are **FPKM-like, non-integer** (e.g. 0.41397153; the column sum is about 401,576). The processing text claims a "matrix table with raw gene counts", but **no raw count matrix is deposited**. Raw counts would have to be regenerated from SRA.
- **SRA:** yes (SRX26170271-SRX26170276), PRJNA1164147.
- **Overlap:** none detected.
- **Released:** 22 Mar 2025.
- **Oddities:** n=3 per group; FPKM only; GRCh38 typo; sample numbering suggests a subset; dose wording.

## 3. GSE261893: PE-MP in chow, 71 weeks (Hong lab, Korea)

- **URLs read:** GEO SOFT self + gsm; FTP suppl listing + filelist.txt; streamed `GSE261893_readcount_Con_PE.txt.gz` and `GSM8153700_readcount_Con1.txt.gz`.
- **Paper:** **none linked** (no PubMed ID). Europe PMC searches for the accession, "71 weeks" and the authors found nothing. UNVERIFIED / likely unpublished.
- **Organism:** C57BL/6, male. Age: "reared over a period of 71 weeks" (age at start UNVERIFIED).
- **Tissue:** liver.
- **Polymer:** **polyethylene** (PE-MP). **Size and dose UNVERIFIED**: the record gives neither. The general summary only defines microplastics as 0.1 um-5 mm.
- **Route / duration:** dietary, "Chow containing PE-MP". The summary says over a year, and the titles say "71weeks". Exact exposure start is UNVERIFIED.
- **Diet:** normal chow, with no metabolic stress.
- **Groups (from GSM):** normal chow n=3 (CON1-3) and PE chow n=3 (PE-1 to PE-3).
- **Controls:** normal chow without plastic. **The same 3 control samples are reused across four series (see below).**
- **Platform:** GPL24247, NovaSeq 6000. GRCm39, "New Tuxedo" (HISAT2/StringTie/Ballgown).
- **Processed data:** `GSE261893_readcount_Con_PE.txt.gz` (430 kB), with columns gene_id (symbols) CON1 CON2 CON3 PE-1 PE-2 PE-3 and 56,674 genes. The values are **integer read counts**. Per-sample `GSM*_readcount_*.txt.gz` files are in the RAW.tar (1.2 MB). The processing text also mentions FPKM files, but none is listed.
- **SRA:** yes (SRX23990487-SRX23990492), PRJNA1089521.
- **Overlap: CONFIRMED shared controls.** The CON1/CON2/CON3 columns of GSE261893, GSE261725 (PS), GSE286116 (PP) and GSE262399 (PB) are **byte-identical** (diff of all rows). The per-sample control files in GSE261893 and GSE261725 also have identical sizes (195,807 / 196,147 / 198,416 B). They were deposited under different GSM/SRX IDs, in separate BioProjects, with no SuperSeries.
- **Released:** 17 Apr 2025 (submitted 19 Mar 2024).
- **Oddities:** title typo "expousre"; n=3; size/dose not reported; aged animals (about 71 weeks); controls shared and re-deposited as if independent.

## 4. GSE261725: PS-MP in chow, 71 weeks (same lab and design as GSE261893)

- **URLs read:** GEO SOFT self + gsm; FTP suppl listing + filelist.txt; streamed `GSE261725_readcount_Con_PS.txt.gz`.
- **Paper:** none linked. UNVERIFIED.
- **Organism:** C57BL/6 male (from titles and growth protocol; the GSM characteristics carry only tissue and treatment).
- **Tissue:** liver. **Polymer:** polystyrene. **Size and dose UNVERIFIED.** Route: chow. Duration: 71 weeks (as above).
- **Groups (from GSM):** normal chow n=3 and PS chow n=3.
- **Controls:** **the same 3 controls as GSE261893, GSE286116 and GSE262399** (identical counts).
- **Platform / processing:** as for GSE261893 (GPL24247, GRCm39, HISAT2/StringTie/Ballgown).
- **Processed data:** `GSE261725_readcount_Con_PS.txt.gz` (425 kB), with columns gene_id CON1-3 PS-1 to PS-3 and 56,674 genes; integer counts. Per-sample files are in RAW.tar.
- **SRA:** yes (SRX23963505-SRX23963510), PRJNA1088491.
- **Released:** 17 Apr 2025.
- **Oddities:** as GSE261893.

## 5. GSE291616: snRNA-seq, PS 50 nm vs 500 nm, weekly gavage, 12 weeks (Mass lab)

- **URLs read:** GEO SOFT self + gsm; FTP suppl listing + filelist.txt; streamed the head of `GSM8838590_Sample_1_PBS_features.tsv.gz` and the first 2 kB of `..._matrix.mtx.gz`; Europe PMC record for PMID 42811166 (abstract only). PMC (reCAPTCHA) and nature.com (login redirect) could not be read.
- **Paper:** Makdissi N, Viola MF, Steinheuer LM, ... Mass E. "Size-dependent plastic exposure disrupts macrophage function and tissue-specific metabolism." *Nat Metab* 2026. doi:10.1038/s42255-026-01615-8. PMID 42811166, PMC13690586 (full text not read).
- **Organism:** C57Bl6/J male. Age 8-10 weeks at start; GSM says "20-22 weeks old" at sacrifice. Strain and sex appear only in the series overall design, not in the GSM characteristics.
- **Tissue:** liver, snap-frozen, single-nucleus RNA-seq.
- **Polymer / size:** polystyrene, 50 nm and 500 nm.
- **Dose / route / duration:** 100 mg/kg, oral gavage, **weekly**, for 12 weeks.
- **Diet:** chow in this series (not stated explicitly; UNVERIFIED). The abstract describes a separate diet-induced-obesity experiment ("Under diet-induced obesity, microplastics exacerbate hepatic lipid accumulation"), but **no DIO transcriptomic series is linked** and none was found in GEO.
- **Groups (from GSM):** PBS n=2, 50 nm n=2, 500 nm n=2 (6 libraries). **n=2 per group is confirmed.** Pooling of nuclei across mice is UNVERIFIED (each GSM is "Sample N").
- **Controls:** PBS gavage, shared by both size arms.
- **Platform:** GPL30172, NextSeq 2000, 10x Genomics. Cell Ranger 7.1.0, refdata-gex-mm10-2020-A.
- **Processed data:** per-sample `barcodes.tsv.gz` / `features.tsv.gz` / `matrix.mtx.gz` (integer). The matrix header for Sample 1 reads **32,285 genes x 1,312,451 barcodes**. That looks like the **raw (unfiltered) Cell Ranger matrix**, so cell calling and QC must be redone. `GSE291616_RAW.tar` is **563 MB** (matrices are 80-95 MB each) and was **not downloaded**. No cell annotation or metadata file is deposited.
- **SRA:** yes (SRX27954032-SRX27954037), PRJNA1234520.
- **Overlap:** none detected. A related Xenium brain series from the same consortium is GSE291378 (not screened).
- **Released:** 9 Jun 2026 (last update 29 Sep 2026).
- **Oddities:** n=2; unfiltered matrices; sex/strain only in free text.

---

## Additional datasets found by search

### 6. GSE286116 (PP-MP) and GSE262399 (PB-MP, "polypropylene black"): same Hong lab 71-week chow series

- **URLs read:** GEO SOFT self + gsm for both; streamed `GSE286116_readcount_Con_PP.txt.gz` and `GSE262399_readcount_Con_PB.txt.gz`.
- Both are C57BL/6 male liver on GPL24247 NovaSeq 6000, chow with or without MP, 71 weeks, n=3+3, integer read counts, released 17 Apr 2025. BioProjects PRJNA1207961 and PRJNA1091649. SRA: SRX27284480-485 and SRX24056716-721.
- **Their CON1-3 columns are identical to GSE261725/GSE261893.** Across the four series there are 3 unique controls and 12 exposed animals (PE, PS, PP, PB with 3 each). That family can be analysed **only as one 15-sample study with 3 shared controls**, never as four independent replicates. Paper, size and dose are UNVERIFIED.

### 7. GSE303044: PET single gavage, Sprague-Dawley rats, liver/stomach/duodenum

- **URLs read:** GEO SOFT self + gsm; FTP suppl listing; streamed `GSE303044_liver.gene_counts.tsv.gz`; PubMed esummary 41853959.
- **Paper:** "A single oral exposure to polyethylene terephthalate microplastics causes mild metabolic and gastrointestinal disruption: dose and sex determinants." *Drug Chem Toxicol* 2026. doi:10.1080/01480545.2026.2643619. PMID 41853959.
- **Organism / timing:** Rattus norvegicus Sprague-Dawley, males and females, adult (3 months). **Single** oral gavage of PET at 5 or 50 mg/kg, with tissues collected about 18.5 h later. **Particle size UNVERIFIED.** The vehicle is called "Vehicle Control" (composition UNVERIFIED).
- **Liver groups (from GSM):** male vehicle 8, male low 8, male high 8; female vehicle 8, female low 8, female high 6. That is 46 liver samples out of 141 total, split across library batches 1-3.
- **Platform:** GPL36046 (Lexogen CORALL mRNA-seq V2, UMI dedup), STAR + Salmon, mRatBN7.2.
- **Processed data:** `GSE303044_liver.gene_counts.tsv.gz` (1.67 MB), 30,562 genes x 46 columns labelled L_01 to L_48 (L_43 and L_47 absent). The counts are **Salmon estimated counts and non-integer**. The mapping of L_xx to GSM is UNVERIFIED (it needs the GSM titles). Transcript-level counts are also deposited.
- **SRA:** BioProject PRJNA1293402. Released 7 Jul 2026.
- **Oddities:** acute (18.5 h) exposure; rat; batch partially unbalanced across groups.

### 8. GSE316023: GAN diet (D09100310) time course, C57BL/6J male, two housing temperatures (MASH reference candidate)

- **URLs read:** GEO SOFT self + gsm; FTP listing; streamed `GSE316023_COUNT_NCOMMS-25-40823-T.csv.gz`; PubMed esummary 42486853.
- **Paper:** "AI digital pathology as a key tool providing in-depth understanding of the progression and regression of MASH and fibrosis in male mouse models." *Nat Commun* 2026. doi:10.1038/s41467-026-73370-z.
- **Design:** male C57BL/6J (JAX), 7-8 weeks old at start. Diet: **GAN diet D09100310 (the same catalogue diet as the GSE294180 MASH arm)** vs normal chow (LabDiet 5001). Housing at 21 C (TA) or 29 C (TN).
- **Groups (from GSM):** GAN-TA at wk12/16/23/32 with n = 8/8/8/7; GAN-TN at the same weeks with n = 8/8/8/6; chow-TA wk23 n=8; chow-TN wk23 n=8 (77 samples). Chow controls exist **only at week 23**, so there is no 8-week GAN-vs-chow contrast. The **earliest is GAN wk12 vs chow wk23**, which is age-confounded. A clean diet contrast is GAN wk23 vs chow wk23 (n=8 vs 8 per housing).
- **Controls:** chow is **LabDiet 5001 grain-based chow**, not the purified D09100304 control used in GSE294180. That is a known source of a large "purified vs chow" signal.
- **Platform:** GPL24247 NovaSeq 6000.
- **Processed data:** `GSE316023_COUNT_...csv.gz` (3.3 MB) and an FPKM file. COUNT holds gene symbols x 77 columns (SD_wkXX_A/B_Rn and NC_wk23_A/B_Rn; A/B is presumably TA/TN, UNVERIFIED). The values are **non-integer estimated counts** (e.g. 1857.9999, 2451.3762). There are **malformed rows** (e.g. `"Gm26206\t1"` embedded in a gene name).
- **SRA:** PRJNA1400485. Released 2 Jun 2026.

### 9. GSE246328: GAN DIO-NASH with chow-reversal (Gubra), plus its overlap with GSE294629

- **URLs read:** GEO SOFT self + gsm for both; streamed both count files; PubMed esummary 38308212 and 40691365.
- **GSE246328 paper:** "Dietary intervention reverses molecular markers of hepatocellular senescence in the GAN diet-induced obese and biopsy-confirmed mouse model of NASH." *BMC Gastroenterol* 2024. doi:10.1186/s12876-024-03141-x.
- **GSE246328 groups (from GSM):** GAN vehicle at 8/16/24 wk with n = 9/8/7; chow-reversal at 8/16/24 wk with n = 7/8/8; chow vehicle at 24 wk n=8 (55 samples). The GAN feeding before the study weeks is long; per Gubra practice it is usually more than 30 weeks (UNVERIFIED for this series). C57BL/6J; sex is not stated in the GSM fields. Platform GPL19057, NextSeq 500, STAR.
- **GSE246328 data:** `GSE246328_GAN-DIO-NASH-Senescence_counts.txt.gz`, Ensembl IDs x 55 columns named by library ID (16174...), **integer counts**.
- **GSE294629 paper:** "Modulation of metabolic, inflammatory and fibrotic pathways by semaglutide in MASH." *Nat Med* 2025. doi:10.1038/s41591-025-03799-0.
- **GSE294629 design and data:** 77 samples, C57BL/6JRj male, GAN D09100310 for 34 weeks vs Altromin 1324 chow. Arms: vehicle, semaglutide and chow-reversal at 8/16/24 wk, plus chow vehicle 24 wk n=7. Platform GPL21626. Data: `GSE294629_processed_file.tsv.gz`, Salmon counts formatted as x.0.
- **OVERLAP CONFIRMED:** GSE294629 GSM titles name the fastq libraries, e.g. "16174.fastq.gz,...,Chow_Vehicle_24_w". These are the same library IDs as the GSE246328 columns (16174 = "Chow Vehicle 24 w, Rep 1 [16174]"). **GSE294629 re-deposits the GSE246328 samples plus semaglutide arms. Never use both.**
- **Reference value:** the chow comparator is again grain chow, and the GAN duration is very long. It is a mature MASH-fibrosis signature, not an 8-week one.

### 10. GSE256501: FPC diet, male and female C57BL/6J, 16 weeks (bulk liver; SubSeries of GSE256505)

- **URLs read:** GEO SOFT self + gsm; streamed `GSE256501_Bulk_liver_raw_counts.txt.gz`; PubMed esummary 39911943.
- **Paper:** "Sex differences in murine MASH induced by a fructose-palmitate-cholesterol-enriched diet." *JHEP Rep* 2025. doi:10.1016/j.jhepr.2024.101222.
- **Design:** FPC-NASH diet vs standard chow for 16 weeks.
- **Groups (from GSM):** male SC 3, male NASH 3, female SC 3, female NASH 3.
- **Platform:** GPL17021 HiSeq 2500.
- **Data:** `GSE256501_Bulk_liver_raw_counts.txt.gz`, 25,238 genes (symbols) x 12, integer counts. The columns are mouse IDs (F5464...), and their mapping to groups is UNVERIFIED (it needs GSM descriptions).
- **Related:** this is a SubSeries of SuperSeries GSE256505, alongside snRNA-seq (GSE256503) and isolated-cell bulk (GSE256502). Sample overlap between the SubSeries has not been checked.

### 11. Human MASLD cohorts (brief)

- **GSE135251** (URL read: SOFT self): 216 snap-frozen biopsies, NextSeq 500 (GPL18573), PMIDs 33268509 / 33762733 / 39566466. Supplementary data is RAW.tar of per-sample files; whether those are counts was not checked (UNVERIFIED).
- **GSE130970** (SOFT self): 78 biopsies (6 normal, 72 across the NAFLD spectrum), GPL16791, PMID 31467298. Supplementary data includes `GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz` and a TPM file.

### Searched but not retained

- GSE347380 (rat PS nanoplastic inhalation; lung-focused).
- GSE306246 (LDLR-/- MP atherosclerosis, n=2).
- GSE306011 (paternal MP, sperm sncRNA).
- Other non-liver MP series (placenta, lung, testis, ovary, heart, brain).
- GSE223193 (NASH time course, but Ldlr-/-).
- The many GAN/MASH knockout or drug series. They are usable only for their WT control-vs-diet arms, with small n.

No GEO series was found with a **microplastic x HFD/MASH diet factor in liver** other than GSE294180 itself. No other **PE + liver** series exists apart from the Hong 71-week chow series (GSE261893).

---

## Suitability judgements

| Dataset | (1) Replicate GSE294180 PE effect | (2) Independent MASH-diet reference | (3) Pool for meta-analysis / ML |
|---|---|---|---|
| GSE245069 | Weak. Different polymer (PS), sizes (nm to 5 um), outbred ICR, drinking water at 1 mg/L, 12 wk, no diet factor. Usable only as a loose "MP-in-liver" concordance check of direction. | No | Do not pool: strain (ICR), polymer, route, platform (DNBSEQ, RSEM-like counts). |
| GSE277772 | Weak. PS, size unknown, n=3 vs 3, FPKM only (raw counts need re-alignment from SRA). | No | Do not pool: FPKM only, n=3, unknown size. |
| GSE261893 (PE) | The only other PE liver set, but the design is very different: dietary PE for 71 wk, aged mice, chow only, size/dose unknown, n=3 vs 3. A null result would not refute GSE294180. Directional overlap check only. | No | Only as a "family" with GSE261725/286116/262399 (3 shared controls). Never treat as independent studies, or the controls are double-counted. |
| GSE261725 (PS) | As GSE261893, but PS. | No | Same shared-control problem. |
| GSE286116 / GSE262399 | As above (PP / black PP). | No | Same shared-control problem. |
| GSE291616 | Not comparable for bulk replication: snRNA-seq, PS 50/500 nm, n=2 per group, weekly 100 mg/kg, chow. Could provide pseudobulk Kupffer/hepatocyte signatures for a qualitative cell-type check. | No | Do not pool: single-nucleus vs bulk; n=2; raw unfiltered matrices (563 MB). |
| GSE303044 (rat PET) | Very weak. Rat, single dose, 18.5 h acute, PET. | No | Do not pool: species, acute exposure, Salmon non-integer counts, both sexes, batches. |
| GSE316023 (GAN) | n/a | **Best candidate.** Same diet (D09100310), male C57BL/6J, n=8 per group. Caveats: grain chow comparator (not D09100304), chow only at wk23, housing-temperature factor, non-integer counts with malformed rows. Use TA wk23 GAN vs chow (and/or the TN pair) as the reference contrast. | Not as raw pooling with GSE294180; use signature/rank-based comparison. |
| GSE246328 (GAN DIO-NASH) | n/a | Good secondary reference with integer counts (vehicle vs chow at 24 wk, n=7 vs 8). Advanced, long-duration MASH; grain chow. **Do not also use GSE294629 (same libraries).** | Signature-level only. |
| GSE256501 (FPC) | n/a | Small (n=3 per group) but integer counts and a fructose + cholesterol diet; adds a different MASH diet; males usable (3 vs 3). | Signature-level only. |
| GSE135251 / GSE130970 (human) | n/a | Human translational check of a diet signature, mapped by ortholog. | Separate-species validation, not pooling. |

### Main reasons NOT to pool the microplastic datasets with GSE294180

- **Polymer:** only the Hong 71-week series has PE; all others are PS, PP or PET.
- **Size:** GSE294180 uses 10-150 um. The others are 50 nm to 5 um, or not reported.
- **Route and dose:** GSE294180 is 2 mg/day by gavage. Others use drinking water, chow admixture, weekly gavage at 100 mg/kg, twice-weekly oral dosing, or a single gavage.
- **Duration:** 8 weeks here, against 12 weeks, 71 weeks, or 18.5 hours elsewhere.
- **Strain and species:** C57BL/6 here, against ICR or SD rat in two of the others.
- **Diet:** no MP dataset other than GSE294180 has a MASH or HFD factor. GSE294180 uses purified control and MASH diets; the others use chow or an unspecified diet.
- **Platform and quantification:** DNBSEQ with RSEM-like counts, HISAT2/featureCounts FPKM only, StringTie/Ballgown, Cell Ranger single-nucleus, and Salmon UMI.
- **Shared controls:** 4 Hong series share identical controls.
- **Sample size:** n=2-4 per group in most series.
- **Missing sex/strain fields:** several GSMs lack these.

Any cross-study use should be limited to **direction-of-effect or rank-based signature concordance** on matched genes, with study as a fixed effect. ML trained on a pooled MP label would mostly learn study, batch and platform.
