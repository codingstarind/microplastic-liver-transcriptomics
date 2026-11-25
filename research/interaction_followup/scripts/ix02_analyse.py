"""Interaction follow-up, step 2: stability, null comparison and held-out full-vs-additive tests (reads ix01 outputs).
Interpreter: .venv or C:/Python314 (no catboost needed).

Criteria fixed before any ix01 result was inspected:
  Strength = absolute pure-interaction SD of the 2-D partial dependence (log2FC units, training subsample).
  Per split, a pair 'exceeds null' if its observed SD > all 10 interaction-free null replicates of that split.
  Reproducible pair (per model): exceeds null in >= 12/15 splits in CD AND in >= 12/15 splits in MD.
  Strongest reproducible pair: largest mean over diets of median (observed - null median) SD.
  Fold stability: Kendall's W of pair ranks over 15 splits; CD-MD repeatability: Spearman of median SD over 28 pairs.
  Held-out benefit (primary, 12 tests): per-gene gain = MSE_additive - MSE_full (repeat-averaged), groups
    {ALL, LIPO, PPARA} x {catboost, hgb} x {CD, MD}; module-cluster t (G-1 df); Holm. Secondary (BH): CHOL, OTHER.
Run: .venv/Scripts/python.exe research/interaction_followup/scripts/ix02_analyse.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

RS = Path(__file__).resolve().parents[2]
OUT = RS / "interaction_followup"
RES = OUT / "results"
SEED, B = 20261125, 2000
H = pd.read_csv(RES / "friedman_H_training_folds.csv.gz")
H["pair"] = H.f1 + " x " + H.f2
H["split"] = H["repeat"] * 5 + H["fold"]
H["pair_type"] = np.where(H.f1.str.endswith("_eff") & H.f2.str.endswith("_eff"), "eff x eff",
                          np.where(H.f1.str.endswith("_se") & H.f2.str.endswith("_se"), "se x se", "eff x se"))

# ---------------- per split: observed vs null
obs = H[H.kind == "observed"].set_index(["target", "model", "split", "pair"])
nul = H[H.kind == "null_additive"].groupby(["target", "model", "split", "pair"]).agg(
    null_max=("inter_sd", "max"), null_median=("inter_sd", "median"), null_H2_median=("H2", "median"))
S = obs.join(nul)
S["exceeds_null"] = S.inter_sd > S.null_max
S["excess"] = S.inter_sd - S.null_median
S = S.reset_index()
S["rank_in_split"] = S.groupby(["target", "model", "split"]).inter_sd.rank(ascending=False)
S.to_csv(RES / "pair_by_split_observed_vs_null.csv", index=False)

summ = S.groupby(["target", "model", "pair", "pair_type"]).agg(
    median_inter_sd=("inter_sd", "median"), q25=("inter_sd", lambda v: v.quantile(.25)), q75=("inter_sd", lambda v: v.quantile(.75)),
    median_null_sd=("null_median", "median"), median_excess=("excess", "median"), n_splits_exceed=("exceeds_null", "sum"),
    median_H2=("H2", "median"), median_null_H2=("null_H2_median", "median"), median_rank=("rank_in_split", "median"),
    n_top1=("rank_in_split", lambda v: int((v == 1).sum()))).reset_index()

# Kendall's W (fold stability) and CD-MD repeatability
stab = []
for (t, m), g in S.groupby(["target", "model"]):
    R = g.pivot(index="split", columns="pair", values="inter_sd").rank(axis=1)
    mR, n = R.shape
    Rs = R.sum(0)
    W = 12 * ((Rs - Rs.mean()) ** 2).sum() / (mR**2 * (n**3 - n))
    stab.append(dict(target=t, model=m, kendall_W=float(W), n_splits=mR, n_pairs=n))
rep = []
for m in ("catboost", "hgb"):
    a = summ[(summ.model == m) & (summ.target == "PE_in_CD")].set_index("pair")
    b = summ[(summ.model == m) & (summ.target == "PE_in_MD")].set_index("pair").loc[a.index]
    rep.append(dict(model=m, spearman_CD_MD_median_sd=float(stats.spearmanr(a.median_inter_sd, b.median_inter_sd).statistic),
                    spearman_CD_MD_excess=float(stats.spearmanr(a.median_excess, b.median_excess).statistic),
                    top_pair_CD=a.median_inter_sd.idxmax(), top_pair_MD=b.median_inter_sd.idxmax()))
for m in ("catboost", "hgb"):
    for t in ("PE_in_CD", "PE_in_MD"):
        g = S[(S.model == m) & (S.target == t)]
        cb = S[(S.model == ("hgb" if m == "catboost" else "catboost")) & (S.target == t)]
        x = g.groupby("pair").inter_sd.median()
        z = cb.groupby("pair").inter_sd.median().loc[x.index]
        rep.append(dict(model=f"{m} vs {'hgb' if m == 'catboost' else 'catboost'}", target=t,
                        spearman_between_models=float(stats.spearmanr(x, z).statistic)))
pd.DataFrame(stab).to_csv(RES / "fold_stability_kendallW.csv", index=False)
pd.DataFrame(rep).to_csv(RES / "repeatability_CD_MD_and_models.csv", index=False)

wide = summ.pivot_table(index=["model", "pair", "pair_type"], columns="target",
                        values=["median_inter_sd", "median_null_sd", "median_excess", "n_splits_exceed", "median_H2", "median_rank", "n_top1"])
wide.columns = [f"{a}_{b.replace('PE_in_', '')}" for a, b in wide.columns]
wide = wide.reset_index()
wide["reproducible"] = (wide.n_splits_exceed_CD >= 12) & (wide.n_splits_exceed_MD >= 12)
wide["mean_excess"] = (wide.median_excess_CD + wide.median_excess_MD) / 2
wide = wide.sort_values(["model", "mean_excess"], ascending=[True, False])
wide.to_csv(RES / "pair_summary_CD_MD.csv", index=False)
strongest = {}
for m, g in wide.groupby("model"):
    r = g[g.reproducible]
    strongest[m] = None if r.empty else r.iloc[0].pair

# ---------------- CatBoost SHAP interactions (descriptive; training genes)
SH = pd.read_csv(RES / "catboost_shap_interactions_training_folds.csv.gz")
SH["split"] = SH["repeat"] * 5 + SH["fold"]
inter = SH[SH.f2 != "MAIN"].copy()
inter["pair"] = inter.f1 + " x " + inter.f2
mains = SH[SH.f2 == "MAIN"].groupby(["target", "genes", "split"]).mean_abs_shap.sum().rename("sum_main_abs_shap")
inter = inter.join(mains, on=["target", "genes", "split"])
inter["share_of_main"] = inter.mean_abs_shap / inter.sum_main_abs_shap
shs = inter.groupby(["target", "genes", "pair"]).agg(median_mean_abs_shap=("mean_abs_shap", "median"),
                                                     median_share_of_main=("share_of_main", "median"),
                                                     median_n_genes=("n_genes", "median")).reset_index()
shs.to_csv(RES / "catboost_shap_interaction_summary.csv", index=False)

# ---------------- held-out full vs additive
P = pd.read_csv(RES / "heldout_full_vs_additive_predictions.csv.gz")
P["se_full"], P["se_add"], P["se0"] = (P.y - P.pred_full) ** 2, (P.y - P.pred_additive) ** 2, (P.y - P.null) ** 2
P["dir_full"] = (np.sign(P.pred_full) == np.sign(P.y)).astype(float)
P["dir_add"] = (np.sign(P.pred_additive) == np.sign(P.y)).astype(float)
Gq = P.groupby(["target", "model", "gene"]).agg(module=("module", "first"), mse_full=("se_full", "mean"), mse_add=("se_add", "mean"),
                                                 mse0=("se0", "mean"), dir_full=("dir_full", "mean"), dir_add=("dir_add", "mean")).reset_index()
gs = pd.read_csv(RS / "roadmap_phase3/results/gene_sets.csv")
SETS = {"LIPO": "LIPOGENESIS_CORE", "PPARA": "SANDERSON_PPARA_TARGETS", "CHOL": "CHOL_SYNTHESIS_CORE"}
genes = np.array(sorted(Gq.gene.unique()))
mem = {k: np.isin(genes, gs.gene[gs.set == v]) for k, v in SETS.items()}
GROUPS = {"ALL": np.ones(len(genes), bool), "LIPO": mem["LIPO"], "PPARA": mem["PPARA"], "CHOL": mem["CHOL"],
          "OTHER": ~(mem["LIPO"] | mem["PPARA"] | mem["CHOL"])}
mods = np.array(sorted(Gq.module.unique()))
BOOTW = np.random.default_rng(SEED).multinomial(len(mods), np.ones(len(mods)) / len(mods), size=B)


def cluster_t(v, m):
    u = np.unique(m)
    Gc, n, mean = len(u), len(v), v.mean()
    se = np.sqrt((Gc / (Gc - 1)) * sum(((v[m == g] - mean).sum()) ** 2 for g in u)) / n
    q = stats.t.ppf(.975, Gc - 1)
    return mean, mean - q * se, mean + q * se, float(2 * stats.t.sf(abs(mean / se), Gc - 1)), Gc


rows = []
for t in ("PE_in_CD", "PE_in_MD"):
    for m in ("catboost", "hgb"):
        d = Gq[(Gq.target == t) & (Gq.model == m)].set_index("gene").loc[genes]
        for g, mk in GROUPS.items():
            x = d[mk]
            gain = (x.mse_add - x.mse_full).values
            mean, lo, hi, p, Gc = cluster_t(gain, x.module.values)
            Sm = x.groupby("module")[["mse_full", "mse_add", "mse0"]].sum().reindex(mods, fill_value=0).values
            bs = BOOTW @ Sm
            with np.errstate(invalid="ignore", divide="ignore"):
                bd = (bs[:, 1] - bs[:, 0]) / bs[:, 2]
            r2f, r2a = 1 - x.mse_full.sum() / x.mse0.sum(), 1 - x.mse_add.sum() / x.mse0.sum()
            rows.append(dict(diet=t, model=m, group=g, n_genes=int(mk.sum()), n_modules=Gc, R2_full=r2f, R2_additive=r2a,
                             dR2_full_minus_additive=r2f - r2a, dR2_ci_low=np.nanquantile(bd, .025), dR2_ci_high=np.nanquantile(bd, .975),
                             mean_gain_mse=mean, gain_ci_low=lo, gain_ci_high=hi, p_cluster_t=p,
                             n_genes_full_better=int((gain > 0).sum()), dir_full=float(x.dir_full.mean()), dir_additive=float(x.dir_add.mean())))
HO = pd.DataFrame(rows)


def holm(p):
    p = np.asarray(p, float)
    o, adj, run = np.argsort(p), np.empty(len(p)), 0
    for r_, j in enumerate(o):
        run = max(run, (len(p) - r_) * p[j])
        adj[j] = min(1, run)
    return adj


def bh(p):
    p = np.asarray(p, float)
    n, o = len(p), np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


pm = HO.group.isin(["ALL", "LIPO", "PPARA"])
HO.loc[pm, "p_holm_primary12"] = holm(HO.loc[pm, "p_cluster_t"])
HO.loc[~pm, "p_BH_secondary"] = bh(HO.loc[~pm, "p_cluster_t"])
HO.to_csv(RES / "heldout_full_vs_additive_by_group.csv", index=False)

info = dict(strongest_reproducible_pair=strongest,
            n_reproducible_pairs={m: int(g.reproducible.sum()) for m, g in wide.groupby("model")},
            primary_heldout_tests_holm_lt_0_05_gain_pos=int(((HO.p_holm_primary12 < .05) & (HO.mean_gain_mse > 0)).sum()))
json.dump(info, open(RES / "run_info_ix02.json", "w"), indent=2)
pd.set_option("display.width", 260)
print(json.dumps(info, indent=1))
print(pd.DataFrame(stab).round(3).to_string(index=False))
print(pd.DataFrame(rep).round(3).to_string(index=False))
print(wide.round(4).groupby("model").head(8).to_string(index=False))
print(HO.round(4).to_string(index=False))
