# horvatits2022_human_cirrhosis (human context record; extracted by coordinator from PMC9386716)

## Citation
Horvatits T, Tamminga M, Liu B, Sebode M, Carambia A, Fischer L, Püschel K, Huber S, Fischer EK. "Microplastics detected in cirrhotic liver tissue." EBioMedicine 2022; 82:104147. DOI 10.1016/j.ebiom.2022.104147. PMID 35835713. PMCID PMC9386716. Open access. Corresponding: University Medical Center Hamburg-Eppendorf and Center for Earth System Research and Sustainability (CEN), University of Hamburg, Germany. (Full author list taken from the Europe PMC core record; the XML header shows only corresponding authors.)

## Access level read
Full text XML (Europe PMC). Supplement ("Further comments on methodology ... are given in the supplement") NOT read.

## Biological question
Do microplastics accumulate in human liver, and does cirrhosis favour that accumulation? Background: "two aspects remain unclear: whether MPs accumulate in peripheral organs, specifically in the liver, and if liver cirrhosis favours this process."

## Experimental design (groups, controls)
Proof-of-concept case series, Germany. "This proof-of-concept case series, conducted in Germany, Europe, analyzed tissue samples of 6 patients with liver cirrhosis and 5 individuals without underlying liver disease. A total of 17 samples (11 liver, 3 kidney and 3 spleen samples) were analyzed according to the final protocol." (Summary/Methods)
- Comparison groups: cirrhotic liver vs non-cirrhotic liver, kidney and spleen; plus **procedural blanks** ("blank samples ... which underwent the same processing procedures").
- Human observational study: no exposure assignment, no randomisation, no blinding reported.

## Species/tissue
Human. Liver (11 samples), kidney (3), spleen (3) from 11 individuals (6 cirrhosis, 5 without liver disease). Patient demographics, MELD, aetiology of cirrhosis: in Table 1 (not individually extracted here).

## Biological replicates per group
6 cirrhotic patients, 5 non-cirrhotic individuals; 17 tissue samples in total. Per-organ n is small (kidney 3, spleen 3). No power calculation.

## Polymer; particle size; detection method; characterisation; contamination control
- Detection window: **4–30 µm** ("A reliable method for detection of MP particles from 4 to 30 µm in human tissue was developed"; "Six different microplastic polymers ranging from 4 to 30 µm in size were detected").
- Method: chemical digestion of tissue → Nile red staining → fluorescence microscopy (quantification) → **Raman spectroscopy** (polymer identification). Anodisc and silver filters; particles "significantly distinguishable from the background by the emitted fluorescence were counted".
- Polymer abbreviations defined in the paper: PS (polystyrene), PET, PVC, PMMA, POM, PP. **Which six polymers were detected, and whether polyethylene was among them, is NOT VERIFIED here** (the per-polymer breakdown is in the Results/Table not extracted).
- Contamination control: procedural blanks processed identically; LOD = "mean + 3 SD" of blanks (LOD 2.4), LOQ = "mean + 10*SD" (7.9982). Degreasing with acetone, pre-cleaned filters.
- Endotoxin/other: not applicable (analytical study).

## Dose/route/duration
Not applicable — environmental human exposure, unquantified. "exact data on the dietary habits and living conditions are not available".

## Data types
Analytical chemistry/microscopy only. **No transcriptomics, proteomics or any molecular profiling.** Clinical metadata (Table 1).

## Statistics
"Kruskal-Wallis test (non-parametric rank sum test) followed by a post-hoc Dunn's test. For tests of the differences between the organs of the same patients (patient 1 to 3), a Wilcoxon test was also performed. All P values reported are two-sided and p<0.05 was considered as significant." LOD/LOQ as above.

## Main findings with effect directions and statistical evidence
- **Cirrhotic liver positive vs blanks:** "MP concentrations in liver samples derived from patients with cirrhosis were significantly higher than those from blank samples (which underwent the same processing procedures) (n = 6, p = 0.009, α = 0.05, Table 2)." ANALYTICAL DETECTION (Raman-confirmed).
- **Non-cirrhotic tissue negative:** "Considering the limit of detection, all liver, kidney and spleen samples from patients without underlying liver disease tested negative for MPs." / "no significant difference between MP concentrations could be detected" for healthy liver, kidney and spleen. ANALYTICAL NULL.
- **Six polymers, 4–30 µm**, detected in cirrhotic liver (Summary/Findings).
- **No causal or mechanistic claim.** Interpretation: "Future studies are needed to evaluate whether hepatic MP accumulation represents a potential cause in the pathogenesis of fibrosis, or a consequence of cirrhosis and portal hypertension."

## Pathway-theme directions
None — no molecular pathway data of any kind (no lipid, PPAR, inflammation, oxidative-stress, mitochondrial, fibrosis or circadian readouts).

## Author-acknowledged limitations (Discussion, "This study had several limitations")
1. Localisation: "because the methodology requires chemical digestion, we cannot determine exactly where in the liver MPs accumulated and whether MPs were deposited intracellularly, for example in Kupffer cells."
2. Detection floor: "this study sets a lower MP detection limit at 4 µm. Future approaches should consider increasing optical resolution ... down to 1 µm."
3. "the sample size was rather small. However, this is a proof-of-concept case series".
4. "exact data on the dietary habits and living conditions are not available", so no conclusion on alimentary source; iatrogenic sources considered but "none of the studied cirrhosis patients underwent liver surgery or radiological embolization procedures".
5. "Adjustment for potential confounders, such as age or sex was not possible due to the small sample size."
6. "a unified methodology to detect MPs in human tissue has not yet been established, making comparison with other scientific work difficult."
7. "future studies are needed to reproduce these findings in larger cohorts."

## Methodological points to verify independently
- Which six polymers were identified, and was **polyethylene** among them? This matters because Jung 2026 motivates its PE focus by asserting PE is the polymer most abundant in human liver — that claim should be traced to its own cited source rather than assumed to come from this paper.
- The primary statistical comparison is cirrhotic liver **vs procedural blanks**, not cirrhotic **vs** non-cirrhotic liver; the "higher in cirrhosis" conclusion should be checked against Table 2 for a direct group comparison.
- Cirrhosis vs control groups differ in age, sex, comorbidity and tissue-sampling route (explants/biopsy vs autopsy?) — unadjusted by the authors' own statement.
- Reverse causation is explicitly open (portal hypertension/altered barrier in cirrhosis may increase translocation rather than MPs causing cirrhosis).

## Data availability statement (verbatim) and accessions
"Data are available on request from the authors." (Data sharing statement.) **No repository accession of any kind.**

## Relevance to the team
Human-relevance anchor only. It establishes that microplastic particles are detectable in diseased human liver but not in healthy liver in a small case series, with no exposure quantification and no molecular data. It supports the motivation for RQ2 (does a diseased liver differ?) but cannot inform direction of gene expression, and it cannot be reanalysed.
