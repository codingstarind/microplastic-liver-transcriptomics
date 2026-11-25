# Phase 2: hypothesis-driven analysis (in progress)

**Status labels used throughout:**

- **Deposited counts, pilot-verified under A1–A3**: the deposited GEO matrix without a tximport length offset, adopted as the working dataset after the re-quantification pilot (team decision 2026-11-02; see the quantification report §5.1). Output folders keep their original `provisional_deposited` names.
- **EXPLORATORY**: not a pre-specified confirmatory test.

## Reports

| Report | Contents |
|---|---|
| [analysis_plan.md](analysis_plan.md) | Choices fixed before results were seen, plus recorded deviations |
| [../quantification_validation/pilot_protocol.md](../quantification_validation/pilot_protocol.md) | Re-quantification pilot protocol and decision rules (pre-registered) |
| [../quantification_validation/quantification_validation_report.md](../quantification_validation/quantification_validation_report.md) | Count-format investigation and pilot results |
| [model_and_diagnostics.md](model_and_diagnostics.md) | Factorial DESeq2 model, contrasts, diagnostics |
| [pathway_report.md](pathway_report.md) | Hallmark/PPAR pathway responses (fry, camera, fgsea, mouse-level scores) |
| [pe_mash_similarity_report.md](pe_mash_similarity_report.md) | Direction A: PE vs diet concordance with a shared-reference correction |
| [interaction_report.md](interaction_report.md) | Direction B: diet × PE interaction |
| [robustness_report.md](robustness_report.md) | Direction C: robustness matrix and unresolvable confounders |
| [literature_comparison.md](literature_comparison.md) | Reproduced vs extended relative to Jung et al. 2026 |
| [interpretation.md](interpretation.md) | Scientific interpretation and decisions for the team |

## Reproducing Phase 2

Phase 2 uses R/Bioconductor and Salmon inside WSL (Ubuntu 24.04), installed with micromamba from `envs/quantification.yml` and `envs/r-analysis.yml`. Large files live outside the repository, in `D:\mpliver_work` (`/mnt/d/mpliver_work` in WSL).

1. Tools (inside WSL):

   ```bash
   micromamba create -f envs/quantification.yml
   ```

   ```bash
   micromamba create -f envs/r-analysis.yml
   ```

2. Factorial model on deposited counts (provisional):

   ```bash
   Rscript analysis/deseq2/run_deseq2_factorial.R $PWD deposited provisional_deposited
   ```

3. Pathways:

   ```bash
   Rscript analysis/pathways/pathway_analysis.R $PWD provisional_deposited
   ```

4. Concordance (Windows `.venv` or any Python with the repository installed):

   ```bash
   python analysis/concordance/pe_vs_diet_concordance.py data/processed/normalized_counts_filtered.csv.gz provisional_deposited
   ```

5. Interaction patterns:

   ```bash
   python analysis/interaction/interaction_patterns.py provisional_deposited
   ```

6. Sensitivity variants (WSL):

   ```bash
   bash analysis/robustness/run_sensitivity.sh provisional_deposited
   ```

7. Robustness summary:

   ```bash
   python analysis/robustness/summarise_robustness.py provisional_deposited
   ```

8. Figures:

   ```bash
   python analysis/figures/phase2_figures.py provisional_deposited provisional_deposited__cov_muscle
   ```

**Re-quantification:**

- `pipelines/quantification/quantify_run.sh` processes one SRA run (fastp → Salmon).
- `run_pilot_from_downloads.sh` processes browser-downloaded `.sra` files.
- `compare_tximport.R` evaluates the pre-specified metrics.
