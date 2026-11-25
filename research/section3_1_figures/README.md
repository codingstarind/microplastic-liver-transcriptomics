# Section 3.1 figures: PCA and pathway enrichment

All figures come from `scripts/make_figures.py`. Run it from the repository root:

```powershell
.venv\Scripts\python.exe research\section3_1_figures\scripts\make_figures.py
```

The script takes about 10 seconds. No statistical method was changed. The enrichment values are read straight from the original output, and the PCA re-runs the original Phase 1 functions.

## Files

| File | Content |
|---|---|
| `figures/fig_pca.{pdf,png}` | Figure 1: PCA before and after muscle-gene exclusion |
| `figures/fig_enrichment.{pdf,png}` | Figure 2 (main paper): fry results for the 3 PPAR sets and for the 9 Hallmark sets with fry FDR < 0.25 in at least one diet |
| `figures/fig_enrichment_full_supplementary.{pdf,png}` | Supplementary figure: fry results for all 53 sets in both diets |
| `source_data_fig_*.csv` | The exact values plotted in each figure |
| `verification.json` | Checks of the plotted values against the original outputs |

PDFs are vector. PNGs are 300 dpi. Both figures use Arial at 7–8.5 pt and are 6.9 in wide, for full width in a two-column layout.

## Figure 1: PCA (source and transformations)

**Workflow.** This is the same workflow as `scripts/06_pca_clustering.py` (Phase 1), using the `mpliver` functions:

- **Input:** `data/processed/log2_normalized_filtered.csv.gz`, which holds log2(median-of-ratios normalised count + 1) for the 13,777 genes that pass the filter, in all 19 mice.
- **Muscle genes:** genes with |Spearman ρ| ≥ 0.7 against the skeletal-muscle marker score, log2(muscle-marker CPM + 1), via `marker_set_cpm` and `genes_tracking_score`.
  - This gives 135 genes.
  - The list is identical to `results/tables/muscle_tracking_genes.csv`.
- **Features:** the 500 most variable genes (`top_variable_genes`), taken before exclusion (a) and again after excluding the 135 muscle genes (b).
- **PCA:** `run_pca`, with genes mean-centred but not scaled.

**Verification.**

- Panel (a) scores match `results/tables/pca_primary_scores.csv` to within 5e-5. That file is rounded to 4 decimals.
- Explained variance matches `results/tables/pca_variance_explained.csv`:

| Panel | PC1 | PC2 |
|---|---|---|
| (a) before exclusion | 30.33% | 17.04% |
| (b) after exclusion | 20.84% | 17.94% |

**Note.** The transformation is the Phase 1 log2(normalised + 1), not VST.

## Figure 2: enrichment (source and tests)

**Source.** `results/pathways/provisional_deposited/limma_fry_camera.csv`, produced by `analysis/pathways/pathway_analysis.R`.

**Test.**

- limma `fry` on the VST matrix with the factorial design.
- Contrasts PE_in_CD and PE_in_MD.

**Columns used.**

- `fry_direction` sets the marker: ▲ up, ▼ down.
- `fry_fdr` sets the x position as −log10(FDR) multiplied by +1 for up or −1 for down. It also sets the fill: black for FDR < 0.05, grey for 0.05–0.10, white for ≥ 0.10.

**FDR families.** BH-FDR was computed in the original script separately for each family and contrast: 50 Hallmark sets and 3 PPAR sets. Nothing was recomputed.

**Thresholds.** Dashed lines mark FDR = 0.05 and dotted lines FDR = 0.10.

**camera.** camera results are not plotted. The `camera_*` columns were not used.

**Main-figure selection rule.** All 3 PPAR sets are shown, plus every Hallmark set with fry FDR < 0.25 in either diet. This is a display rule only, not a significance criterion. The supplementary figure shows all 53 sets.

**Verified values** (from `verification.json`):

| Set | Diet | fry FDR |
|---|---|---|
| PPARα targets (Sanderson) | CD | **0.0406** (the only FDR < 0.05) |
| HALLMARK_MTORC1_SIGNALING | CD | 0.0796 |
| HALLMARK_MTORC1_SIGNALING | MD | 0.0916 |
| HALLMARK_CHOLESTEROL_HOMEOSTASIS | MD | 0.0916 |
| HALLMARK_DNA_REPAIR | MD | 0.0916 |
| HALLMARK_GLYCOLYSIS | MD | 0.0916 |

The other five sets all have FDR between 0.05 and 0.10.

## Captions

**Figure 1.** Principal component analysis of liver transcriptomes from all 19 mice in GSE294180.

- **(a)** The 500 most variable genes of log2(median-of-ratios normalised count + 1).
- **(b)** The same procedure after excluding 135 genes whose expression tracked a skeletal-muscle marker score (|Spearman ρ| ≥ 0.7).

Genes were mean-centred, not scaled. Axis labels give the percentage of variance explained. Fill shows diet (dark: control diet, CD; light: metabolic-dysfunction diet, MD) and marker shape shows treatment (circle: vehicle; triangle: polyethylene, PE).

**Figure 2.** Directional gene-set tests (limma fry) of the PE versus vehicle effect under control diet (CD; a) and metabolic-dysfunction diet (MD; b).

- **Position:** the x-axis shows −log10(BH-adjusted FDR), signed by the direction of change in PE.
- **Markers:** ▲ means up and ▼ means down. Black means FDR < 0.05, grey means 0.05 ≤ FDR < 0.10, and white means FDR ≥ 0.10.
- **Threshold lines:** dashed at FDR = 0.05, dotted at FDR = 0.10.
- **FDR families:** FDR was controlled separately within the Hallmark family (50 sets) and the PPAR family (3 sets), for each contrast.
- **Sets shown:** all PPAR sets, and Hallmark sets with fry FDR < 0.25 in either diet. All 53 sets are in Supplementary Figure S1.
- **camera:** competitive camera results are not shown.

**Supplementary Figure S1.** As Figure 2, for all 50 Hallmark and 3 PPAR gene sets.
