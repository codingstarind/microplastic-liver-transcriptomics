# Roadmap Phase 2: reference-dataset discovery, verification and preparation (2026-11-22)

## Main sources
- `research/roadmap_phase2/README.md`
- `reports/dataset_selection_report.md`, `gene_mapping_report.md`, `qc_report.md`, `external_microplastic_status.md`, `phase3_readiness_report.md`, `phase2_signoff.md`
- `metadata/*.csv`
- `qc/*.png`
- Scripts: `scripts/01_…08_*.py`, `gene_mapping.py`, `test_phase2_data.py`

## 1. Question
Which independent mouse-liver datasets can serve as MASH, fasting, time-of-day and PPARα reference axes? Can they be mapped onto the GSE294180 gene annotation without bias?

No differential expression, signature or model was run. **No GSE294180 PE effect was used** in this phase.

## 2. Why
Roadmap Phase 1 chose Track A, which needs references from independent mice. The fasting reference had not yet been located.

## 3. Selected references
| Axis | Dataset | Primary contrast (biological n) | Values | Key caveat |
|---|---|---|---|---|
| MASH 1° | GSE316023 (Nat Commun 2026) | GAN D09100310 vs LabDiet 5001 chow, 23 wk, 21 °C (8 v 8) | estimated counts, rounded | the control is chow, not purified diet; needs a length offset (603 genes > 1.2× implied-length change) |
| MASH 2° | GSE246328 (BMC Gastroenterol 2024) | GAN vehicle 24 w (≈ 58 wk GAN) vs chow (7 v 8) | integer STAR counts | advanced MASH; GSE294629 is the same animals and was rejected |
| Fasting 1° | GSE107787 (Kinouchi 2018, Cell Rep) | 24 h fasted vs fed, 6 ZT × 3 per state (18 v 18) | RPKM only | muscle carry-over aligned with FED at ZT8–20 |
| Fasting 2° | GSE225695 | fasted vs fed (3 v 3) | StringTie counts | sacrifice time unknown; direction check only |
| Circadian 1° | GSE117134 | WT female, ZT0–20 every 4 h (6 × 3) | RPKM only | female; diurnal sampling |
| PPARα 1° | GSE282998 (Heintz 2025, Toxicol Sci) | GW7647 5/10/20 mg/kg vs vehicle in B6129SF2/J (12 v 4); Ppara-null specificity arm | TempO-Seq probe counts (FFPE) | targeted assay; dose-group depth and muscle batch |
| PPARα 2° | GSE8295 | Wy14643 | GCRMA array | audited, not prepared |

## 4. Methods implemented
- **GEO SOFT parsing and ENA run tables.**
- **Checksummed downloads:** 33 files, about 150 MB.
- **Gene mapping to GENCODE vM35 symbols** (`gene_mapping.py`):
  - Ensembl version stripping;
  - exact, then case-insensitive symbol match, then a unique MGI synonym rescue (e.g. Arntl → Bmal1);
  - many-to-one handling: sum for Ensembl-keyed counts, highest-mean row for symbol-keyed data and probes.
- **QC** (`05_qc.py`):
  - depth ratio and Mann–Whitney test against treatment;
  - detected genes;
  - muscle score;
  - correlation-outlier rule (median − 3 × 1.4826 × MAD);
  - PCA per dataset;
  - overlap checks (shared GSMs, identical SRA read and base counts).
- **Fixed Phase 3 inputs** (`08_phase3_inputs.py`):
  - **muscle exclusion rule:** 309 genes (187 in the filtered set), defined as gastrocnemius median RPKM ≥ 10 and ≥ 20 × (clean liver + 1), taken from GSE107787's own muscle samples;
  - **gene universe:** main 11,888 genes; sensitivity 8,360;
  - **eligible-sample table:** 173 rows.
- **Tests:** 60 passed (`test_phase2_data.py`).

## 5. Key numbers
**Gene coverage:** GSE294180-filtered genes present in each reference: 13,242 / 13,552 / 12,566 / 12,600 / 12,245.

**Depth ratios (case/control):**
| Contrast | Ratio | Note |
|---|---|---|
| MASH 1° | 0.93 | p = 0.13 |
| MASH 2° | 1.02 | — |
| Fasting | 1.16 | p = 0.001 |
| PPARα | 1.68 | — |

**Muscle signal:**
- GSE107787: FED median 51.8 vs FAST 0.75 summed RPKM.
- GSE282998: 2 of 4 mice at 0 and 10 mg/kg (about 200–270 CPM).

**Hong-lab controls:** 3 control SRA runs appear identically in 4 series, so the 12 control GSMs are 3 mice.

## 6. Failures and limitations
- The AIBD GSE245069 analysis files were not found on this machine.
- GSE277772 deposits normalised values, not counts.
- Fasting and circadian references are RPKM only, so they need log-linear models, not negative-binomial models.
- None of the references addresses the PE ≡ depth confound.

## 7. Transition
With a fixed universe, muscle list and sample table in place, Phase 3 could estimate the reference effects. Its first step was the MASH positive control: our Diet_in_Veh effect vs the GSE316023 axis, as a go/no-go.
