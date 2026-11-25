# Citation map: which reference supports which statement, and where

- **BibTeX:** `references.bib` (79 entries).
- **Verification status per entry:** `reference_verification.csv`.
  - 64 references were VERIFIED against Crossref: the DOI resolves, and the first author and year match.
  - 8 have a resolving DOI with a year that differs only by online-first vs issue date (Park 2025/2026, Sweeney 2016/2017, Haynes 2016/2017, Robinson 2009/2010, Stephens 2016/2017). Two have no author/year match because the citing text gave only a journal: GSE316023 paper = Saigal et al. 2026; GSE246328 paper = Flensted-Jensen et al. 2024. Use the Crossref year/first author.
  - 7 have no DOI (Cawley 2010, CART 1984, LightGBM 2017, CatBoost 2018, Ng & Jordan 2002, Holm 1979, scikit-learn 2011). Confirm these manually.
- **"Access"** records what the project actually read, taken from the phase reference files. We have **not** read the full text of every methods reference cited for an algorithm's properties.

## Introduction

| Statement to support | Keys | Access / caution |
|---|---|---|
| Microplastics reach liver; PE in human tissue | horvatits2022 (polymer list UNVERIFIED) | earlier review only; verify before use |
| Source study: oral PE causes hepatotoxicity and worsens MASH; PPARα–Anxa2 | jung2026 | FT+SI |
| Conflicting directions of lipogenic genes across plastic studies | jung2026, djouina2023, makdissi2026, okamura2023, xup2026, chiu2025 | Chiu: ABS only. Djouina fold changes read from bar charts |
| Diet × plastic designs rarely analyse interactions; reporting issues | ahn2025, okamura2023, fang2026, zhai2024, han2025 | FT |
| Small-n RNA-seq power | schurch2016 | FT |
| Batch / confounding cannot be removed when aligned with treatment | leek2010, zhang2020combatseq | FT |
| Murine MASH models vs human; pathway-level agreement | vacca2024, teufel2016 | FT / ABS |

## Materials and Methods

| Method | Keys |
|---|---|
| Quantification and import | patro2017 (Salmon), chen2018fastp, soneson2015 (tximport) |
| Differential expression | love2014 (DESeq2), robinson2010 (edgeR), law2014 (voom), ritchie2015 (limma), stephens2017 (ashr) |
| Multiple testing | benjamini1995, holm1979 |
| Gene sets | liberzon2015 (Hallmark); wu2010roast (fry); wu2012camera (camera) |
| Reference datasets | gse316023paper (Saigal 2026), gse246328paper (Flensted-Jensen 2024), kinouchi2018 (GSE107787), heintz2025 (GSE282998), li2025 (GSE245069), chiu2025 (GSE277772) |
| Precision weights / τ² | dersimonian1986 |
| Errors in variables | fuller1987 |
| Linear and regularised regression | hoerl1970, tibshirani1996, zou2005 |
| Trees and ensembles | breiman1984, breiman1996, breiman2001, friedman2001, ke2017, prokhorenkova2018; strobl2007/strobl2008 (importance caveats); freund1997 (considered, not run) |
| Nested and grouped CV | varma2006, cawley2010, krstajic2014, roberts2017; ambroise2002 (leakage) |
| Comparing models across CV splits | nadeau2003, bouckaert2004, bates2023 |
| Feature engineering rationale | sweeney2017, schubert2018, lee2008, vacca2024 |
| Per-mouse scoring | foroutan2018 (singscore); hanzelmann2013 (GSVA, considered, not run) |
| Classification | cox1958, ng2002, handyu2001, quinlan1986, raileanu2004, hand2006 |
| External / cross-study validation framing | bernau2014, steyerberg2016, riley2016, debray2015 |
| Robust regression | huber1964 |
| Software | pedregosa2011 (scikit-learn) |

## Results

Cite methods inline as above. The literature comparisons that belong in Results:
- **Ppara discrepancy:** jung2026 (Table S5, Fig. 2F legend) vs our estimate.
- **Lipogenic gene direction vs published PE:** djouina2023.
- **GSE245069 hub-gene / lipid claims:** li2025 (our study-level contrast only).

## Discussion

| Statement | Keys |
|---|---|
| Fasting/feeding and time of day drive SREBP and PPARα programmes; ZT-dependence | kinouchi2018 (reference data); xup2026 (ZT sampling in plastic work) |
| External validation needs independent cohorts; held-out genes ≠ held-out cohorts | bernau2014, steyerberg2016, riley2016, debray2015 |
| Single-experiment signatures are noisy; reproducibility improves with more studies | sweeney2017, schubert2018 |
| Small n limits interaction detection | schurch2016 |
| Microbiome/gut-liver mechanisms proposed elsewhere (context only; not tested here) | zhai2024, xu2024 |
| Confounded batch cannot be corrected | leek2010, zhang2020combatseq |

## Do not cite for these claims
- Do not cite jung2026 for **our** Ppara result. Cite it only for their reported values.
- Do not cite horvatits2022 for "PE is the dominant polymer in human liver". That claim is UNVERIFIED in our notes.
- Do not cite foroutan2018 for GSVA-type population scoring. We used singscore only.
- References in `roadmap_phase1/working_notes/` that the project never used for an analysis are **not** included. Examples: Maurya 2024 preprint, Lyu 2026, the PET rat study, Luo 2026, Roh 2024.
