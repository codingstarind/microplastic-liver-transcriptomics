"""Read-only extraction of publication tables from saved research outputs (no model is fitted, nothing is modified).
Every row carries `source` = the repository-relative file it came from.
Run: python research/paper_evidence_audit/scripts/extract_tables.py"""
from pathlib import Path

import numpy as np
import pandas as pd

RS = Path(__file__).resolve().parents[2]
OUT = RS / "paper_evidence_audit" / "tables"
OUT.mkdir(exist_ok=True, parents=True)
rel = lambda p: str(p.relative_to(RS.parent)).replace("\\", "/")  # noqa: E731
master = []


def add(phase, analysis, quantity, value, source, ci_low=None, ci_high=None, n=None, context=""):
    master.append(dict(phase=phase, analysis=analysis, quantity=quantity, value=value, ci_low=ci_low, ci_high=ci_high, n=n,
                       context=context, source=source))


# ---------------- Phase 4 task A ----------------
f = RS / "roadmap_phase4/results/taskA/summary_by_run_model.csv"
d = pd.read_csv(f)
d = d[d.run.isin(["primary__PE_in_CD", "primary__PE_in_MD", "positive__Diet_in_Veh"])].copy()
d["model"] = d.model.fillna("null")
d["target"] = d.run.str.split("__").str[1]
T1 = d.pivot_table(index="model", columns="target", values=["r2_oos_mean", "rmse_mean", "spearman_mean"]).round(4)
T1.columns = [f"{a}|{b}" for a, b in T1.columns]
order = ["null", "ols_MASH", "ols_FAST", "ols_PPARA", "ols_CIRC", "ols_all", "ridge", "elasticnet", "tree", "rf", "hgb", "catboost"]
T1 = T1.reindex([o for o in order if o in T1.index])
T1["source"] = rel(f)
T1.to_csv(OUT / "T01_phase4_taskA_models.csv")
for m, r in T1.iterrows():
    for t in ("PE_in_CD", "PE_in_MD", "Diet_in_Veh"):
        add("4", f"taskA {m}", f"R2_oos {t}", r.get(f"r2_oos_mean|{t}"), rel(f), n="11758 genes; 15 outer splits")
        add("4", f"taskA {m}", f"Spearman {t}", r.get(f"spearman_mean|{t}"), rel(f), n="11758 genes; 15 outer splits")
f = RS / "roadmap_phase4/results/taskA/paired_comparisons.csv"
p = pd.read_csv(f)
p["model"] = p.model.fillna("null")
p["comparator"] = p.comparator.fillna("null")
p["source"] = rel(f)
p.to_csv(OUT / "T02_phase4_paired_comparisons.csv", index=False)

# ---------------- Phase 4B ----------------
f = RS / "roadmap_phase4b/results/metrics_main.csv"
m = pd.read_csv(f)
m = m[(m.stage == "main") & (~m.shuffle_b.astype(bool))]
T3 = m.groupby(["target", "rep", "model"])[["r2_oos", "rmse", "spearman", "n_features"]].mean().round(4).reset_index()
T3["n_splits"] = m.groupby(["target", "rep", "model"]).size().values
T3["source"] = rel(f)
T3.to_csv(OUT / "T03_phase4b_representation_x_model.csv", index=False)
for _, r in T3.iterrows():
    add("4B", f"{r.rep} {r.model}", f"R2_oos {r.target}", r.r2_oos, rel(f), n=f"{r.n_splits} splits", context="EXPLORATORY (features chosen after Phase 4)")
f = RS / "roadmap_phase4b/results/all_paired_comparisons.csv"
a = pd.read_csv(f)
a["source"] = rel(f)
a.to_csv(OUT / "T04_phase4b_paired_vs_F8_F4_OLS.csv", index=False)
f = RS / "roadmap_phase4b/results/extra_vs_tree_baselines.csv"
e = pd.read_csv(f)
e["source"] = rel(f)
e.to_csv(OUT / "T05_phase4b_vs_tree_baselines.csv", index=False)
for _, r in e.iterrows():
    add("4B", f"{r.rep}-{r.model} vs {r.vs_rep}-{r.vs_model}", f"dRMSE {r.target}", r.d_rmse, rel(f), r.ci_low, r.ci_high, f"{r.n_better}/{r.n_splits} splits better",
        "EXPLORATORY; NB-corrected CI")

f = RS / "roadmap_phase4b/feature_ablation_results.csv"
x = pd.read_csv(f)
x["source"] = rel(f)
x.to_csv(OUT / "T05b_phase4b_block_ablation.csv", index=False)
for _, r in x.iterrows():
    add("4B", f"ablation {r.model} {r.block_removed}", f"dRMSE (removed - full) {r.target}", r.d_rmse_removed_minus_full, rel(f), r.ci_low, r.ci_high,
        f"{r.n_splits_full_better}/15 splits full better", f"block_useful={r.block_useful}; EXPLORATORY")

# ---------------- Phase 4C ----------------
for name, fn in (("T06_phase4c_collinearity.csv", "collinearity_summary.csv"), ("T07_phase4c_fix_experiment.csv", "fix_experiment_summary.csv"),
                 ("T08_phase4c_target_quality.csv", "target_quality.csv")):
    f = RS / "roadmap_phase4c_diagnostics/tables" / fn
    x = pd.read_csv(f)
    x["source"] = rel(f)
    x.to_csv(OUT / name, index=False)
f = RS / "roadmap_phase4c_diagnostics/tables/collinearity_summary.csv"
for _, r in pd.read_csv(f).iterrows():
    add("4C", f"collinearity {r.matrix}", "max VIF", r.max_vif, rel(f), context=f"median condition number {r.cond_number_median:.2f}")
f = RS / "roadmap_phase4c_diagnostics/tables/fix_experiment_summary.csv"
for _, r in pd.read_csv(f).iterrows():
    add("4C", f"fix experiment {r.rep}-{r.model}{' (historical)' if r.historical else ''}", f"R2_oos {r.target}", r.r2_oos, rel(f),
        r.ci_low, r.ci_high, context=f"dRMSE vs F8-OLS {r.d_rmse_vs_F8_ols:.4f}; CI refers to dRMSE; descriptive (folds reused)")
f = RS / "roadmap_phase4c_diagnostics/tables/subgroup_errors.csv"
sg = pd.read_csv(f)
sg["source"] = rel(f)
sg.to_csv(OUT / "T09_phase4c_subgroup_errors.csv", index=False)

# ---------------- Phase 4D ----------------
f = RS / "roadmap_phase4d/results/contrast_estimates_primary.csv"
x = pd.read_csv(f)
x["source"] = rel(f)
pe = RS / "roadmap_phase4d/results/permutation_exact.csv"
if pe.exists():
    pp = pd.read_csv(pe)
    pp.to_csv(OUT / "T10b_phase4d_permutation_exact.csv", index=False)
x.to_csv(OUT / "T10_phase4d_contrasts.csv", index=False)
for _, r in x.iterrows():
    add("4D", f"singscore {r.signature}", r.contrast, r.estimate, rel(f), r.ci_low, r.ci_high, "19 mice (4/5/5/5)", f"p={r.p:.4g}; Holm={r.p_holm}; family={r.family}")

# ---------------- Phase 4E ----------------
f = RS / "roadmap_phase4e/results/external_metrics.csv"
x = pd.read_csv(f)
x["source"] = rel(f)
x.to_csv(OUT / "T11_phase4e_external_metrics_all.csv", index=False)
for _, r in x[x.task == "mash"].iterrows():
    add("4E", f"{r.split} {r.rep} {r.model}", "AUC", r.auc, rel(f), r.auc_ci_low, r.auc_ci_high, f"{r.n} test mice")
    add("4E", f"{r.split} {r.rep} {r.model}", "balanced accuracy", r.bal_acc, rel(f), n=f"{r.n} test mice")
    add("4E", f"{r.split} {r.rep} {r.model}", "Brier (misclassification rate for R0 rule)", r.brier, rel(f), n=f"{r.n} test mice")
f = RS / "roadmap_phase4e/results/gse294180_frozen_model_contrasts.csv"
x = pd.read_csv(f)
x["source"] = rel(f)
x.to_csv(OUT / "T12_phase4e_gse294180_contrasts.csv", index=False)
for _, r in x.iterrows():
    add("4E", f"frozen classifier {r.score} on GSE294180", r.contrast, r.estimate, rel(f), r.ci_low, r.ci_high, "19 mice", f"p={r.p:.4g}; Holm={r.p_holm}")
f = RS / "roadmap_phase4e/results/r2_feature_contributions_gse294180.csv"
x = pd.read_csv(f)
x["source"] = rel(f)
x.to_csv(OUT / "T13_phase4e_r2_contributions.csv", index=False)

# ---------------- Phase 5 ----------------
f = RS / "roadmap_phase5/gene_level_external_results.csv"
x = pd.read_csv(f)
x["source"] = rel(f)
x.to_csv(OUT / "T14_phase5_gene_level.csv", index=False)
for _, r in x[(x.model.isin(["F8_ols", "F8C_ridge", "FAST_only_ols"])) & (x.training_target == "PE_in_CD")].iterrows():
    add("5", f"{r.study} {r.contrast} {r.model}", "Spearman", r.spearman, rel(f), r.spearman_ci_low, r.spearman_ci_high, f"{r.n_genes} genes",
        f"gene-block CI (anti-conservative); perm p={r.get('perm_p_one_sided')}; R2w own mean={r.R2w_vs_own_mean:.3f}; verdict={r.verdict_locked_rule}")
f = RS / "roadmap_phase5/sample_level_external_results.csv"
x = pd.read_csv(f)
x["source"] = rel(f)
x.to_csv(OUT / "T15_phase5_sample_level.csv", index=False)
for _, r in x[x.score.isin(["FAST_singscore", "MASH_singscore", "PPARA_singscore", "R0_logistic_logodds", "R2_logistic_logodds"])].iterrows():
    add("5", f"{r.study} {r.contrast} {r.score}", "exposed - control", r.estimate, rel(f), r.ci_low, r.ci_high, f"{r.n_exposed} v {r.n_control}",
        f"exact p={r.p_perm_exact_two_sided:.3g} (min {r.p_perm_min_attainable:.3g}); g={r.hedges_g:.2f}; LOO {r.loo_n_same_sign}/{r.loo_n}; GSE294180 PE_in_CD={r.gse294180_PE_in_CD:.3f}")

# ---------------- Phase 5B ----------------
for name, fn in (("T16_phase5b_intervention_verdicts.csv", "intervention_verdicts.csv"), ("T17_phase5b_mash_classifier_bootstrap.csv", "mash_classifier_bootstrap.csv"),
                 ("T18_phase5b_gap_train_cv_external.csv", "gap_train_cv_external.csv"), ("T19_phase5b_reliability.csv", "gap_reliability_and_disattenuation.csv"),
                 ("T20_phase5b_target_distributions.csv", "gap_target_distributions.csv"), ("T21_phase5b_study_shift.csv", "gap_study_shift.csv"),
                 ("T22_phase5b_count_unit.csv", "count_unit_per_sample.csv"), ("T23_phase5b_equivalence_checks.csv", "equivalence_checks.csv")):
    f = RS / "roadmap_phase5b/tables" / fn
    if f.exists():
        x = pd.read_csv(f)
        x["source"] = rel(f)
        x.to_csv(OUT / name, index=False)
f = RS / "roadmap_phase5b/tables/intervention_verdicts.csv"
for _, r in pd.read_csv(f).iterrows():
    add("5B", r.candidate, "internal CV Spearman", r.cv_spearman, rel(f), n="15 outer splits", context=f"SD {r.cv_spearman_sd:.3f}; train {r.train_spearman:.3f}; CV R2 {r.cv_r2:.3f}")
    add("5B", r.candidate, "Hong PE Spearman (EXPLORATORY post-external)", r.rho_hong_pe, rel(f), r.rho_hong_lo, r.rho_hong_hi, "3 v 3 mice",
        f"CI = within-arm mouse bootstrap B=1000; P(delta>0)={r.P_delta_gt_0}; rank {r.rank_of_true}/{r.n_labellings}; verdict={r.verdict}")
f = RS / "roadmap_phase5b/tables/mash_classifier_bootstrap.csv"
for _, r in pd.read_csv(f).iterrows():
    add("5B", f"{r.split} {r.rep} {r.model}", "balanced accuracy (bootstrap)", r.bal_acc, rel(f), r.bal_acc_lo, r.bal_acc_hi, f"{r.n} test mice",
        f"AUC={r.auc}; frac boot AUC=1: {r.frac_boot_auc_eq_1}; {r.metric_prob}={r.value:.3f} [{r.value_lo:.3f},{r.value_hi:.3f}]")

M = pd.DataFrame(master)
M["value"] = pd.to_numeric(M.value, errors="coerce").round(5)
M.to_csv(RS / "paper_evidence_audit" / "results_master_extracted.csv", index=False)
print(len(M), "rows;", sorted(p.name for p in OUT.glob("*.csv")))
