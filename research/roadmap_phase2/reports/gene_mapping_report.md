# Gene mapping report (Phase 2)

**Date:** 2026-11-22.
- **Code:** `scripts/gene_mapping.py`, applied by `scripts/03_prepare_references.py`.
- **Per-identifier record:** `datasets/<GSE>/gene_mapping.csv.gz`. It has one row per original matrix row (`original_id`, `candidate`, `vm35_symbol`, `status`, `usable`, `kept_in_ready`, `collapse_note`).
- **Summary:** `metadata/gene_mapping_summary.csv`.

## 1. Target annotation

All references are mapped onto the **GSE294180 annotation: GENCODE vM35 gene symbols** (`metadata/gene_annotation_gencode_vM35.csv.gz`, 55,836 genes). That file also lists the Ensembl gene IDs behind each symbol; 228 symbols have 2–3 IDs. This keeps the gene universe identical to the one our PE effects were estimated on.

| Dataset | Deposited identifier | Source annotation | Original rows |
|---|---|---|---|
| GSE316023 | gene symbol | mm10 (release NOT STATED) | 48,795 |
| GSE246328 | Ensembl gene ID | GRCm38 Ensembl **v96** | 55,536 |
| GSE107787 | gene symbol | mm10 RefSeq (in-house, 2018) | 20,561 |
| GSE117134 | gene symbol | mm10 RefSeq (in-house) | 18,079 |
| GSE282998 | TempO-Seq probe `Symbol_probeNo` | TempO-Seq mouse whole-transcriptome panel | 30,146 probes |
| GSE225695 | Ensembl gene ID (+ symbol) | mm10 / StringTie | 22,565 |

## 2. Mapping rules, in order

1. **Original identifiers are preserved verbatim.** Duplicated identifiers get a `#row<n>` suffix in the original matrix so that every row stays addressable (GSE316023: 254 rows; for example `Rp1` appears twice).
2. **GSE316023 export artefact.**
   - 13,094 identifiers have the form `<symbol>\t1`. The row values are aligned (77 numeric fields), so the symbol before the tab is used (`tab_suffix_recovered`, 9,233 mapped).
   - 8 of these are dropped because a clean row with the same symbol exists (`tab_suffix_conflict_with_clean_row`).
   - 3,853 still do not map; almost all are `Gm` snoRNA/miRNA-like models.
   - **665 identifiers ending in `\t2` or `\t3` are flagged `malformed` and excluded.** The suffix meaning is unknown, possibly multi-locus aggregates. They are almost all Gm/RP pseudogene models.
3. **Ensembl IDs:** version stripped and matched to vM35 Ensembl IDs. Unmatched IDs (`ensembl_not_in_vM35`) are retired or pre-M35 models: GSE246328 4,263 (7.7%); GSE225695 1,000 (4.4%).
4. **Symbols:** exact vM35 match, then a unique case-insensitive match, then an **MGI synonym** (`symbol_mgi_synonym`).
   - The MGI source is `MRK_List1.rpt`, downloaded 2026-11-22 and checksummed in `metadata/download_manifest.csv`.
   - A synonym is used only if it points to exactly one current symbol that exists in vM35 and is not itself a vM35 symbol. Withdrawn symbols (`withdrawn, = X`) map to X.
   - Examples: `Arntl → Bmal1`, `Tceb1 → Eloc`, `Fam150a → Alkal1`, `Wars → Wars1`, `Atp5e → Atp5f1e`, `Hist2h2bb → H2bc18`.
   - Without this step, renamed genes (including the core clock gene Bmal1) would silently drop out.
   - Ambiguous synonyms (for example `Adss`) are **not** rescued.
5. **TempO-Seq probes:** the symbol part is mapped as in rule 4.
6. **Many-to-one resolution:**
   - **Ensembl-keyed count data** (GSE246328, GSE225695): rows mapping to the same vM35 symbol are **summed**, consistent with how a symbol-level count matrix aggregates several gene IDs (GSE246328: 60 rows merged).
   - **Symbol-keyed data and probes** (GSE316023, GSE107787, GSE117134, GSE282998): the row or probe with the **highest mean** across all samples in the dataset (GSE282998: within the GW7647 arm) is kept and the others are recorded as dropped. Summing is avoided here because duplicate symbols may be distinct loci. For TempO-Seq this is the rule the original authors used.
   - GSE282998: 8,735 probes dropped; the panel has several probes for about 6,300 genes.
7. **No values are changed** by mapping, except that the GSE316023 and GSE225695 estimated counts are rounded to integers in the ready matrix (section 4).

## 3. Gene counts at each step

"Filtered" = the 13,777 GSE294180 genes passing ≥ 10 counts in ≥ 4 samples.

"Detected in contrast" uses a reporting-only rule:
- counts: ≥ 10 in at least half (minimum 3) of the planned-contrast samples;
- RPKM: ≥ 1.

It is **not** a Phase 3 filter.

| Dataset | Total rows | Valid IDs | Mapped to vM35 | Unique vM35 genes (ready) | ∩ GSE294180 filtered (13,777) | Detected in contrast | Detected ∩ filtered |
|---|---|---|---|---|---|---|---|
| GSE316023 | 48,795 | 48,130 | 43,521 | 43,363 | 13,242 | 14,000 | 12,893 |
| GSE246328 | 55,536 | 55,536 | 51,273 | 51,213 | 13,552 | 13,204 | 12,632 |
| GSE107787 | 20,561 | 20,561 | 20,357 | 20,330 | 12,566 | 9,701 | 9,617 |
| GSE117134 | 18,079 | 18,079 | 17,895 | 17,880 | 12,600 | 10,359 | 10,137 |
| GSE282998 (GW7647) | 30,146 probes | 30,146 | 29,933 | 21,198 | 12,245 | 10,765 | 10,295 |
| GSE225695 | 22,565 | 22,565 | 21,565 | 21,561 | 12,395 | 11,696 | 11,497 |

**Candidate common universe** (`metadata/gene_universe_membership.csv.gz`, one row per GSE294180 filtered gene):
- present in all four primary references (GSE316023, GSE107787, GSE117134, GSE282998): **12,051** genes;
- detected in all four primary contrasts: **8,410** genes.

The large drop for the "detected" version comes from the RPKM ≥ 1 rule in the two RPKM datasets and from the targeted panel.

**The RPKM datasets cover fewer genes because their deposits contain RefSeq protein-coding and lncRNA models only, not because of mapping failure.** Unmapped symbols there are 1.0% (GSE107787) and 1.0% (GSE117134) after synonym rescue.

## 4. Data type in the ready matrices

| Dataset | Ready values | Comment |
|---|---|---|
| GSE316023 | **Integers rounded from deposited estimated counts** | About 6% of cells were non-integer (for example 10420.001), typical of EM-based estimated counts. The deposited FPKM file implies a gene length that is almost constant across samples (median CV 1.0%), so **no sample-specific length offset is needed** for a within-dataset comparison. Rounding follows the standard practice for estimated counts (as `DESeqDataSetFromTximport` does). TPM/FPKM were **not** rounded. |
| GSE246328 | Integer counts | As deposited |
| GSE107787, GSE117134 | **RPKM, unchanged** | No counts exist on GEO. They must not be rounded and fed into count models. Phase 3 should use log-scale linear models with an expression filter, or re-quantify from SRA if count-based methods are essential. |
| GSE282998 | Integer probe counts | Targeted assay. The gene universe is the probe panel, and genes absent from the panel are absent by **assay design**, not because they are unexpressed. |
| GSE225695 | Integers rounded from StringTie-derived counts | Secondary check only |

No dataset had missing values. Duplicated gene identifiers occur only in GSE316023 (254 rows), resolved as above.

## 5. Coverage of the programmes we care about (`metadata/key_gene_coverage.csv`)

Detection in the planned-contrast samples, for genes in the GSE294180 filtered set:
- **Lipogenesis/SREBP1, cholesterol/SREBP2 and PPARα-target genes** (Fasn, Acly, Acaca, Me1, Scd1, Hmgcr, Sqle, Ldlr, Cyp4a14, Acot1, Ehhadh, Cpt1a, Hmgcs2, Pdk4, Angptl4, Ppara and others) are detected in all six references. The exceptions are Fgf21 and Plin4 (not detected in GSE117134) and Pnpla3 (not detected in GSE225695).
- **Clock genes:** Bmal1, Clock, Cry1, Per3, Nr1d1/2, Dbp, Tef and Nfil3 are detected everywhere. Per1 and Per2 are not detected in the GSE282998 primary contrast.
- **MASH/fibrosis/inflammation:** Col1a1, Col3a1, Lgals3 and Cd68 are detected everywhere. Gpnmb, Trem2, Mmp12 and Ccl2 are detected only in the MASH references, which is expected because they are disease-induced.
- **Skeletal-muscle genes** must be excluded from all similarity work (see QC report).

## 6. Decisions to make before Phase 3 analysis (not made here)

1. **Universe choice:**
   - (a) the "present in all primary references" universe (12,051 genes), with each reference's own low-expression handling; or
   - (b) the "detected in all" universe (8,410).

   Recommendation: fix (a) for the main analysis and (b) as a sensitivity analysis, and state the choice before looking at PE-vs-reference results.
2. **A fixed muscle-gene exclusion list.** Suggested basis: the GSE294180 muscle-tracking genes identified in Phase 2 plus canonical sarcomere genes.
