# Phase 5 readiness (planning only; nothing executed)

## Frozen models / assays
- **Gene level:** F8-OLS (primary) and F8+C-Ridge (secondary), as recommended by Phase 4C. Block B (target-informed) and FULL (collinear) are **not** frozen.
- **Mouse level:** the Phase 4D signatures (`frozen_signature_manifest.csv`: MASH_main_n50, FAST_main_n50; PPARA_main_n50 secondary) with singscore on a TPM-proxy. No re-selection of genes.

## Independent microplastic datasets (Phase 2 status)

| Study | Independence | Usable data | Role |
|---|---|---|---|
| Hong lab GSE261893 (PE) + GSE261725/GSE286116/GSE262399 (PS/PP/PB) | One lab; **3 controls shared by all 4 series** (12 GSMs = 3 mice) | integer counts | **Primary external test:** PE vs control (3 v 3) as one study |
| GSE245069 (PS 80 nm / 0.5 µm / 5 µm, ICR) | Independent lab; one shared water control | estimated counts + FPKM | **Secondary:** treat the 3 sizes as one study with a shared control |
| GSE277772 (PS) | Independent lab | **no raw counts** (TMM/FPKM-like); needs SRA re-quantification | Not usable without re-quantification |
| AIBD GSE245069 analysis files | — | **not found locally** | — |

## Valid external evaluation
1. **Gene level:** fit the frozen models on all GSE294180 genes; predict each external study's plastic log2FC vector (estimated with DESeq2 in that study, its own controls). Score R² vs that study's own mean, Spearman ρ and sign agreement, against the null and single-axis baselines.
2. **Mouse level:** score the external mice with the frozen signatures. Test plastic vs control within study (Welch / exact permutation; n = 3 v 3 gives a minimum two-sided p of 0.1, so effect sizes and CIs are primary).

## What counts as failure (fixed before looking)
- Gene level: R² ≤ 0 or Spearman CI including 0 in the primary external study.
- Mouse level: fasting-score difference with sign opposite to GSE294180, or a CI centred near 0.
- "Plastic-specific" claims fail outright if both polymers and controls show the same pattern within the Hong lab.

## Limits on Phase 5
- n = 3 per arm; polymer, strain, diet, duration and lab all differ, and none of them are separable.
- PE is represented by one external lab.
- External studies also need their own depth/yield QC.
