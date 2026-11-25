# Phase 3: independent reference signatures, positive control, baselines

**Status:** complete for review (2026-11-22). Phase 4 not started. **Start with:** `reports/phase3_decision_report.md`.

## Contents

| Path | What |
|---|---|
| `analysis_plan.md` | Pre-specified models, contrasts and acceptance criteria. Locked before the first fit (SHA-256 in `results/plan_lock.json`; test-checked). |
| `reports/` | `count_offset_validation`, `mash_reference_report`, `fasting_reference_report`, `circadian_reference_report`, `ppara_reference_report`, `reference_comparison_report`, `baseline_comparison_report`, `phase3_decision_report` |
| `reference_effects/` | One gene-level table per reference contrast (all universe genes): `gene, log2FC, SE/lfcSE, ci95_low/high, stat, pvalue, padj (BH over the universe), baseMean/AveExpr, detected, dataset, contrast, model` |
| `results/` | `offset_audit/`, `mash/` (run_info, LOO, pathways, influence), `positive_control/`, `fasting/`, `circadian/`, `ppara/`, `ppara_discrepancy/`, `comparisons/` (pairwise, collinearity, baselines), `gene_sets.csv` |
| `figures/` | `mash_positive_control.png`, `reference_similarity_heatmap.png`, `baseline_comparisons.png`, `fasting_reference.png`, `ppara_reference.png`, `circadian_reference.png`, `offset_audit.png`, `diagnostics/` |
| `scripts/` | `a01` offset audit → `b01` MASH (R) → `c01` positive control → `d01` fasting, circadian, PPARα (R) → `e01` Ppara discrepancy (R) → `f01` comparisons and baselines |
| `tests/test_phase3.py` | Plan lock; reference-table schema and universe; BH correctness; contrast directions (GAN > chow on disease markers, fasting Hmgcs2↑/Fasn↓, GW7647 targets ↑ and KO-attenuated, Dbp/Bmal1 T2 signs); fixed sample lists; positive-control decision recomputation; deliverables |

## Review addendum (post-hoc, 2026-11-22)

`reports/phase3_review_addendum.md` and `reports/phase4_specification_memo.md`. These are post-hoc verifications only; no original output changed.

| Script | What it does | Outputs |
|---|---|---|
| `scripts/g01_exact_permutation_audit.py` | Exact enumeration of all 12,870 allocations | `results/posthoc_audit/` |
| `posthoc_muscle_rule_prereg.md`, then `scripts/g02_muscle_sensitivity_{list.py,fit.R,eval.py}` | Pre-registered, exploratory muscle sensitivity | `results/posthoc_audit/` |

## Inputs (unchanged)
- Phase 2 fixed universe (11,888 main / 8,360 sensitivity), muscle list (309), eligible samples, prepared matrices in `D:/mpliver_work/phase2_refs/`.
- GSE294180 effects from `results/deseq2/provisional_deposited/` and its `loo/` refits (read only).
- Nothing in Phase 1/2 or the original analyses was modified.

## Software
- **R** 4.4.3 in WSL Ubuntu-24.04, env `/root/mpliver-tools/mamba/envs/mpliver-r`: DESeq2 1.46.0, limma 3.62.1, edgeR 4.4.0, ashr 2.2.63, msigdbr (MSigDB 2026.1.Mm, as in Phase 2).
- **Python** 3.14 (pandas 3.0.5, numpy 2.5.1, scipy, matplotlib 3.11.1).
- **Seed** 20261122.

## Run

From `research/roadmap_phase3/scripts/` (R scripts via WSL):

```powershell
python a01_offset_audit.py
```
```powershell
wsl -d Ubuntu-24.04 -e /root/mpliver-tools/mamba/envs/mpliver-r/bin/Rscript /mnt/c/Users/Soham/Desktop/microplastic/research/roadmap_phase3/scripts/b01_mash_references.R
```
```powershell
python c01_positive_control.py
```
```powershell
wsl -d Ubuntu-24.04 -e /root/mpliver-tools/mamba/envs/mpliver-r/bin/Rscript /mnt/c/Users/Soham/Desktop/microplastic/research/roadmap_phase3/scripts/d01_other_references.R
```
```powershell
wsl -d Ubuntu-24.04 -e /root/mpliver-tools/mamba/envs/mpliver-r/bin/Rscript /mnt/c/Users/Soham/Desktop/microplastic/research/roadmap_phase3/scripts/e01_ppara_discrepancy.R
```
```powershell
python f01_compare_and_baselines.py
```
```powershell
python -m pytest ..\tests -q
```

Runtime is about 8 minutes in total.

## Multiple-testing families
- **Genes:** BH within each (dataset, contrast) over the universe genes reported; independent filtering off.
- **Pathways:** BH within each (contrast, family: Hallmark | PPAR + custom).
- **Hypergeometric overlaps and correlations** are descriptive; they are not corrected across pairs.

## Deviations from `analysis_plan.md` (the plan text itself is not edited)
1. **Offset audit check (ii):** the script used ρ ≥ 0.50 as the pass mark, but the plan gave no number. The result (0.497) does not change the decision, which was already set by the undocumented provenance.
2. **Positive-control c1 null:**
   - The pre-specified label-relabelling null includes partially correct labellings, so it cannot reject for strong effects.
   - Reported as FAIL as written. A balanced-split null is **proposed, not applied**.
3. **MASH LOO:** GSE246328 LOO was added for symmetry; GSE316023 TN and offset models were fitted without LOO (sensitivity models).
4. **Circadian outputs:** `results/other_refs_run_info.json` field `sex` is mislabelled (it holds a group name). The sexes are female (GSE117134) and male (GSE107787) as documented.
5. **Positive-control test genes:** Col1a1 was replaced by Col3a1/Trem2 in the direction test, because Col1a1 is in the fixed muscle exclusion list (finding reported, list not changed).
6. **R .gz reading:** `fread` cannot read .gz without R.utils, which is not installed. Files are streamed through `zcat` (no numerical effect).
