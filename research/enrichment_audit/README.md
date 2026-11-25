# Enrichment audit: GSE294180 PE vs Veh (CD, MD) and diet x PE interaction

Computational outputs only. No interpretation.
- Existing files under `results/` and `analysis/` were read, not modified.
- Gene sets: MSigDB 2026.1.Mm (msigdbr 26.1.1), 50 Hallmark + 3 PPAR sets:
  - REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA
  - WP_PPAR_SIGNALING_PATHWAY
  - SANDERSON_PPARA_TARGETS

## Already available (reused, not rerun)
Source: `analysis/pathways/pathway_analysis.R` → `results/pathways/provisional_deposited/`. Covers all 53 sets × 4 contrasts.
- **limma fry:** directional and mixed.
- **camera.**
- **fgsea:** preranked on the DESeq2 Wald statistic.
- **Mean-z score:** fitted with lm.

Sensitivity variants (`results/pathways/`):
- covariate runs: `__cov_depth`, `__cov_muscle`, `__cov_clock_phase`, `__cov_quality_proxy`;
- `__drop_MD_1`, `__pilot_length_offset`;
- 19 leave-one-mouse-out refits (`loo/`).

## Added here (missing before)
| Item | Script | Output |
|---|---|---|
| Gene-symbol mapping audit (MSigDB vs deposited matrix; Ensembl-ID rescue via GENCODE vM35) | `scripts/ea01_mapping_and_sensitivity.R` | `results/gene_mapping_audit.csv`, `results/gene_mapping_unmatched_symbols.csv` |
| Reproduction check of the stored primary fry/camera | same | `results/reproduction_check_primary_fry_camera.csv` |
| Sensitivity: 309 muscle-exclusion genes removed | same | `results/sensitivity/new_sensitivity_fry_camera.csv`, `new_sensitivity_fgsea.csv` |
| Sensitivity: depth + muscle covariates jointly | same | same |
| Sensitivity: vehicle libraries binomially thinned to PE depth (seed 11), full DESeq2 refit | same | same, plus `results/sensitivity/depth_thinned_deseq2_*.csv.gz` and `run_info.json` |
| Alternative BH families (all 53 sets per contrast; 53 × 3 PE contrasts) | `scripts/ea02_tables.py` | `tables/T1` |
| Consolidated tables | `scripts/ea02_tables.py` | `tables/T1`–`T8` (below) |
| Pathway figure: CD \| MD \| interaction | `scripts/ea03_plot.py` | `figures/pathway_PE_vs_Veh_CD_MD_interaction.{png,pdf,svg}` |

## Tables
- `T1_enrichment_long_all_methods.csv`: set × contrast × method (statistic, direction, p, native FDR, two alternative FDR families, leading edge).
- `T2_enrichment_wide_CD_MD_interaction.csv`: one row per set, all methods, for PE_in_CD, PE_in_MD, Interaction and Diet_in_Veh.
- `T3_direction_consistency_across_methods.csv`: sign agreement between fry, camera, fgsea NES, mean shrunken log2FC, fraction of genes up, and mean-z.
- `T4_sensitivity_fry_camera_long.csv`: fry and camera in every variant, including the 19 LOO refits.
- `T5_sensitivity_summary_per_set.csv`: per set × contrast, fry call (FDR < 0.05) in each variant, plus agreement counts with primary and LOO.
- `T6_sensitivity_counts_by_variant.csv`: counts of significant sets per variant × contrast.
- `T7_sensitivity_fgsea_long.csv`: fgsea across variants.
- `T8_fdr_family_specification.csv`: tests and BH families.

## Rerun
```powershell
wsl.exe -d Ubuntu-24.04 -u root -e bash -lc "cd /mnt/c/Users/Soham/Desktop/microplastic/research/enrichment_audit/scripts && /root/mpliver-tools/mamba/envs/mpliver-r/bin/Rscript ea01_mapping_and_sensitivity.R"
```
```powershell
.venv\Scripts\python.exe research\enrichment_audit\scripts\ea02_tables.py
```
```powershell
.venv\Scripts\python.exe research\enrichment_audit\scripts\ea03_plot.py
```
Logs are in `logs/`.
