"""Locked Stage 3-4 runs (03_ and 04_ documents): nested grouped CV of comparators + interventions on the Phase 4 outer
folds, final refits, exploratory external evaluation, exhaustive relabelling, paired mouse bootstrap, MASH-classifier
bootstrap. All train-dependent steps (RepBuilder, scaler, tau2, hyperparameters) are fitted on training genes only."""
import itertools
import json
import sys
import time

import joblib
import numpy as np
import pandas as pd
from b00_common import N_JOBS, P5, PB, REPO, RES, RS, SEED, TAB, WORK5, dl_tau2, json_dump, mod_diff, sha256
from c00_common import RAW_COLS, RepBuilder, outer_splits, raw_table, target
from joblib import Parallel, delayed
from scipy import stats
from sklearn.base import clone
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import HuberRegressor, LinearRegression, Ridge
from sklearn.metrics import make_scorer
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

B_BOOT = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
IT = RES / "intervention_targets"
MOD = RES / "models"
MOD.mkdir(parents=True, exist_ok=True)

# ---------------- data (module level, read-only in workers) ----------------
X0 = raw_table()
g294 = pd.read_csv(REPO / "data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz", index_col=0)
sm = pd.read_csv(REPO / "metadata/sample_metadata.csv").set_index("counts_column").loc[g294.columns]
ARM = sm.group.values  # CD+Veh, CD+PE, MD+Veh, MD+PE
keep294 = (g294 >= 10).sum(axis=1) >= 4
G = X0.index.intersection(g294.index[keep294])
HX = pd.read_csv(WORK5 / "HONG_expression.csv.gz", index_col=0)
man = pd.read_csv(P5 / "results/external_sample_manifest.csv")
mh = man[man.study == "HONG"].set_index("mouse_id")
HX = HX[mh.index]
HGRP = mh.group.values
keepH = (HX >= 10).sum(axis=1) >= 3
GH = X0.index.intersection(HX.index[keepH])


def cpm(M):
    return M / M.sum(axis=0) * 1e6


def arith(C, a, b):
    """MLE surrogate: log2 ratio of group-mean CPM (+1)."""
    return np.log2(C[:, a].mean(1) + 1) - np.log2(C[:, b].mean(1) + 1)


C294 = cpm(g294.loc[G]).values
L294 = np.log2(C294 + 1)
CH = cpm(HX.loc[GH]).values  # CPM over all prefiltered Hong genes
LH = np.log2(CH + 1)


def source_targets(C, L, arm):
    """All candidate training targets from one (possibly resampled) GSE294180 matrix."""
    cd_pe, cd_v, md_pe, md_v = (arm == "CD+PE"), (arm == "CD+Veh"), (arm == "MD+PE"), (arm == "MD+Veh")
    d_cd, s_cd = mod_diff(L, cd_pe, cd_v)
    d_md, s_md = mod_diff(L, md_pe, md_v)
    out = {"arith_CD": (arith(C, cd_pe, cd_v), s_cd),
           "arith_avg": ((arith(C, cd_pe, cd_v) + arith(C, md_pe, md_v)) / 2, np.sqrt(s_cd**2 + s_md**2) / 2),
           "robust_CD": (d_cd, s_cd), "robust_avg": ((d_cd + d_md) / 2, np.sqrt(s_cd**2 + s_md**2) / 2)}
    return {k: pd.DataFrame({"y": v[0], "s": v[1]}, index=G) for k, v in out.items()}


def normal_scores(y):
    r = stats.rankdata(y)
    return stats.norm.ppf((r - 0.375) / (len(y) + 0.25))


# real-data targets
T_MLE = target("PE_in_CD")[["y", "s"]]
avg = pd.read_csv(IT / "gse294180_PE_avg_over_diets.csv.gz").set_index("gene").rename(columns={"log2FC": "y", "lfcSE": "s"})[["y", "s"]]
ST = source_targets(C294, L294, ARM)
TARGETS = {"MLE_CD": T_MLE, "MLE_avg": avg, "robust_CD": ST["robust_CD"], "robust_avg": ST["robust_avg"]}
spearman_scorer = make_scorer(lambda y, p: stats.spearmanr(y, p).statistic)

CANDS = {  # id: (target key, bootstrap proxy key, rep, estimator, grid, scoring, scale_ok)
    "B0_F8_OLS": ("MLE_CD", "arith_CD", "F8", LinearRegression(), {}, None, True),
    "B1_F8C_Ridge": ("MLE_CD", "arith_CD", "F8+C", Ridge(alpha=1000), {}, None, True),
    "B2_F8_HGB": ("MLE_CD", "arith_CD", "F8", HistGradientBoostingRegressor(max_iter=300, max_leaf_nodes=8, early_stopping=False, random_state=SEED),
                  {"model__learning_rate": [0.03, 0.1], "model__min_samples_leaf": [50, 200]}, "neg_mean_squared_error", True),
    "I2a_robustCD_Ridge": ("robust_CD", "robust_CD", "F8", Ridge(), {"model__alpha": [0.1, 1, 10, 100, 1000]}, "neg_mean_squared_error", True),
    "I2b_MLEavg_Ridge": ("MLE_avg", "arith_avg", "F8", Ridge(), {"model__alpha": [0.1, 1, 10, 100, 1000]}, "neg_mean_squared_error", True),
    "I2c_robustavg_Ridge": ("robust_avg", "robust_avg", "F8", Ridge(), {"model__alpha": [0.1, 1, 10, 100, 1000]}, "neg_mean_squared_error", True),
    "I3_Huber_normalscore": ("MLE_CD", "arith_CD", "F8", HuberRegressor(max_iter=500),
                             {"model__epsilon": [1.1, 1.35, 1.75], "model__alpha": [1e-4, 1e-2]}, spearman_scorer, False),
}


def design(cid, tdf, genes):
    tk, _, rep, est, grid, scoring, _ = CANDS[cid]
    j = X0.loc[genes].join(tdf, how="inner").dropna(subset=["y", "s"])
    j = j[j.s > 0]
    y = normal_scores(j.y.values) if cid.startswith("I3") else j.y.values
    return j, y


def pipe(cid, params=None):
    _, _, rep, est, _, _, _ = CANDS[cid]
    p = Pipeline([("rep", RepBuilder(rep)), ("scale", StandardScaler()), ("model", clone(est))])
    if params:
        p.set_params(**params)
    return p


def fit(cid, j, y, tr, params=None, grp=None, tune=False):
    w = None if cid.startswith("I3") else 1 / (j.s.values**2 + dl_tau2(j.y.values[tr], j.s.values[tr]))
    p = pipe(cid, params)
    _, _, _, _, grid, scoring, _ = CANDS[cid]
    kw = {} if w is None else {"model__sample_weight": w[tr]}
    if tune and grid:
        gs = GridSearchCV(p, grid, cv=GroupKFold(4), scoring=scoring, n_jobs=1)
        gs.fit(j[RAW_COLS].values[tr], y[tr], groups=grp[tr], **kw)
        return gs.best_estimator_, gs.best_params_
    p.fit(j[RAW_COLS].values[tr], y[tr], **kw)
    return p, params or {}


COMMON = G.intersection(T_MLE.dropna().index[T_MLE.dropna().s > 0])
for k in TARGETS:
    COMMON = COMMON.intersection(TARGETS[k].dropna().index)


def cv_job(cid, r_i, k_i):
    j, y = design(cid, TARGETS[CANDS[cid][0]], COMMON)
    yref = T_MLE.loc[j.index, "y"].values  # common internal outcome
    for r, k, tr, te, grp in outer_splits(j.index):
        if (r, k) == (r_i, k_i):
            break
    est, bp = fit(cid, j, y, tr, grp=grp, tune=True)
    p = est.predict(j[RAW_COLS].values[te])
    null = np.average(yref[tr])
    out = dict(candidate=cid, repeat=r_i, fold=k_i, n_train=int(tr.sum()), n_test=int(te.sum()), best_params=json.dumps(bp),
               spearman_vs_PE_in_CD=stats.spearmanr(yref[te], p).statistic, spearman_vs_own_target=stats.spearmanr(y[te], p).statistic,
               train_spearman_vs_PE_in_CD=stats.spearmanr(yref[tr], est.predict(j[RAW_COLS].values[tr])).statistic)
    if CANDS[cid][6]:
        out["r2_vs_train_mean_PE_in_CD"] = 1 - np.sum((yref[te] - p) ** 2) / np.sum((yref[te] - null) ** 2)
    return out


t0 = time.time()
jobs = [(c, r, k) for c in CANDS for r in range(3) for k in range(5)]
cv = pd.DataFrame(Parallel(n_jobs=N_JOBS, backend="loky")(delayed(cv_job)(*x) for x in jobs))
cv.to_csv(TAB / "cv_fold_metrics.csv", index=False)
print(f"nested CV done {time.time() - t0:.0f}s", flush=True)

# ---------------- final models (modal inner choice) ----------------
final, hp_rows = {}, []
for cid in CANDS:
    choices = cv[cv.candidate == cid].best_params.value_counts()
    bp = json.loads(choices.index[0])
    hp_rows.append(dict(candidate=cid, modal_params=choices.index[0], counts=json.dumps(choices.to_dict())))
    j, y = design(cid, TARGETS[CANDS[cid][0]], COMMON)
    est, _ = fit(cid, j, y, np.ones(len(j), bool), params=bp)
    joblib.dump(est, MOD / f"{cid}.joblib")
    final[cid] = (est, bp)
pd.DataFrame(hp_rows).to_csv(TAB / "hyperparameter_choices.csv", index=False)
PRED = pd.DataFrame({cid: est.predict(X0[RAW_COLS].values) for cid, (est, _) in final.items()}, index=X0.index)
PRED.to_csv(RES / "candidate_predictions_all_genes.csv.gz")
f5 = pd.read_csv(P5 / "frozen/gene_level_frozen_predictions.csv.gz", index_col=0)
repro_b0 = dict(max_abs_pred_diff_vs_phase5_frozen=float(np.max(np.abs(PRED["B0_F8_OLS"] - f5["F8_ols__PE_in_CD"]))),
                spearman_vs_phase5_frozen=float(stats.spearmanr(PRED["B0_F8_OLS"], f5["F8_ols__PE_in_CD"]).statistic),
                n_common_genes=len(COMMON), n_phase5_training_genes=int(X0.index.intersection(T_MLE.dropna().index).size))

# ---------------- exploratory external evaluation ----------------
DE = P5 / "results/external_de"
ext = {"HONG_PE": pd.read_csv(DE / "HONG__PE_vs_CON.csv.gz").set_index("gene"),
       "GSE245069_PSpooled": pd.read_csv(DE / "GSE245069__PSpooled_vs_CON.csv.gz").set_index("gene")}
for c in ("PS", "PP", "PB"):
    ext[f"HONG_{c}"] = pd.read_csv(DE / f"HONG__{c}_vs_CON.csv.gz").set_index("gene")
hu = pd.read_csv(IT / "hong_unit_harmonised_targets.csv.gz")
for f in (0.5, 0.472):
    ext[f"HONG_PE_unit{f}"] = hu[(hu.unit_factor == f) & (hu.contrast == "PE_vs_CON")].set_index("gene")
# control-free gene vector: PE minus mean of other polymers (no shared-control term)
cf = np.log2(CH[:, HGRP == "PE"].mean(1) + 1) - np.mean([np.log2(CH[:, HGRP == g].mean(1) + 1) for g in ("PS", "PP", "PB")], axis=0)
CF = pd.Series(cf, index=GH)
ev_rows = []
for cid in CANDS:
    p = PRED[cid]
    for name, t in ext.items():
        t = t[np.isfinite(t.log2FC) & np.isfinite(t.lfcSE) & (t.lfcSE > 0)]
        g = p.index.intersection(t.index)
        y, s = t.loc[g, "log2FC"].values, t.loc[g, "lfcSE"].values
        w = 1 / (s**2 + dl_tau2(y, s))
        row = dict(candidate=cid, outcome=name, n_genes=len(g), spearman=stats.spearmanr(y, p[g]).statistic,
                   sign_agree=float(np.mean(np.sign(y) == np.sign(p[g]))))
        if CANDS[cid][6]:
            row["R2w_vs_own_mean"] = 1 - np.sum(w * (y - p[g]) ** 2) / np.sum(w * (y - np.average(y, weights=w)) ** 2)
        ev_rows.append(row)
    g = p.index.intersection(CF.index)
    ev_rows.append(dict(candidate=cid, outcome="HONG_controlfree_PE_minus_otherpolymers", n_genes=len(g), spearman=stats.spearmanr(CF[g], p[g]).statistic))
EV = pd.DataFrame(ev_rows)
# exhaustive relabelling (arith proxy over the 6 PE/CON mice)
six = np.where(np.isin(HGRP, ["PE", "CON"]))[0]
lab = []
gp = PRED.index.intersection(GH)
gi = np.array([GH.get_loc(x) for x in gp])
for cid in CANDS:
    pr = stats.rankdata(PRED.loc[gp, cid])
    vals = []
    for e in itertools.combinations(six, 3):
        a = np.zeros(len(HGRP), bool)
        a[list(e)] = True
        b = np.zeros(len(HGRP), bool)
        b[[i for i in six if i not in e]] = True
        vals.append((set(e) == set(np.where(HGRP == "PE")[0]), stats.spearmanr(pr, arith(CH[gi], a, b)).statistic))
    obs = [v for t_, v in vals if t_][0]
    lab.append(dict(candidate=cid, true_labelling_rho_proxy=obs, rank_of_true=int(sum(v >= obs - 1e-12 for _, v in vals)), n_labellings=len(vals),
                    p_one_sided=sum(v >= obs - 1e-12 for _, v in vals) / len(vals), null_q95=float(np.quantile([v for _, v in vals], .95))))
LABR = pd.DataFrame(lab)
EV.to_csv(TAB / "external_candidates.csv", index=False)
LABR.to_csv(TAB / "labelling_ranks.csv", index=False)
print("external + labellings done", flush=True)

# ---------------- paired mouse bootstrap ----------------
FINAL_PARAMS = {cid: bp for cid, (_, bp) in final.items()}


def boot(i):
    rng = np.random.default_rng(SEED + i)
    idx = np.concatenate([rng.choice(np.where(ARM == a)[0], (ARM == a).sum()) for a in ("CD+Veh", "CD+PE", "MD+Veh", "MD+PE")])
    C, L, arm = C294[:, idx], L294[:, idx], ARM[idx]
    st = source_targets(C, L, arm)
    hid = np.concatenate([rng.choice(np.where(HGRP == g)[0], (HGRP == g).sum()) for g in ("CON", "PE", "PS", "PP", "PB")])
    Ch, hg = CH[:, hid], HGRP[hid]
    yH = pd.Series(arith(Ch, hg == "PE", hg == "CON"), index=GH)
    cfb = pd.Series(np.log2(Ch[:, hg == "PE"].mean(1) + 1) - np.mean([np.log2(Ch[:, hg == g].mean(1) + 1) for g in ("PS", "PP", "PB")], axis=0), index=GH)
    out = {"rep": i}
    for cid in CANDS:
        tdf = st[CANDS[cid][1]]
        j, y = design(cid, tdf, COMMON)
        est, _ = fit(cid, j, y, np.ones(len(j), bool), params=FINAL_PARAMS[cid])
        p = pd.Series(est.predict(X0.loc[GH.intersection(X0.index), RAW_COLS].values), index=GH.intersection(X0.index))
        out[f"{cid}|rho_hong_pe"] = stats.spearmanr(yH[p.index], p).statistic
        out[f"{cid}|rho_controlfree"] = stats.spearmanr(cfb[p.index], p).statistic
    return out


t0 = time.time()
BS = pd.DataFrame(Parallel(n_jobs=N_JOBS, backend="loky", batch_size=4)(delayed(boot)(i) for i in range(B_BOOT)))
BS.to_csv(TAB / "bootstrap_transfer_replicates.csv.gz", index=False)
print(f"bootstrap {B_BOOT} done {time.time() - t0:.0f}s", flush=True)
bsum = []
for cid in CANDS:
    a, c = BS[f"{cid}|rho_hong_pe"], BS[f"{cid}|rho_controlfree"]
    d = a - BS["B0_F8_OLS|rho_hong_pe"]
    bsum.append(dict(candidate=cid, rho_hong_median=a.median(), rho_hong_lo=a.quantile(.025), rho_hong_hi=a.quantile(.975),
                     delta_vs_B0_median=d.median(), delta_lo=d.quantile(.025), delta_hi=d.quantile(.975), P_delta_gt_0=float((d > 0).mean()),
                     rho_controlfree_median=c.median(), rho_controlfree_lo=c.quantile(.025), rho_controlfree_hi=c.quantile(.975),
                     P_controlfree_gt_0=float((c > 0).mean())))
BSUM = pd.DataFrame(bsum)
BSUM.to_csv(TAB / "bootstrap_summary.csv", index=False)

# ---------------- verdicts (locked criteria A1-A4) ----------------
cvs = cv.groupby("candidate").agg(cv_spearman=("spearman_vs_PE_in_CD", "mean"), cv_spearman_sd=("spearman_vs_PE_in_CD", "std"),
                                  train_spearman=("train_spearman_vs_PE_in_CD", "mean"), cv_r2=("r2_vs_train_mean_PE_in_CD", "mean")).reset_index()
hp = EV[EV.outcome == "HONG_PE"].set_index("candidate").spearman
b0 = hp["B0_F8_OLS"]
V = cvs.merge(BSUM, on="candidate").merge(LABR, on="candidate")
V["rho_hong_pe"] = V.candidate.map(hp)
V["delta_rho_hong_vs_B0"] = V.rho_hong_pe - b0
V["A1"] = V.delta_rho_hong_vs_B0 >= 0.05
V["A2"] = V.P_delta_gt_0 >= 0.90
V["A3"] = V.cv_spearman >= V.loc[V.candidate == "B0_F8_OLS", "cv_spearman"].iloc[0] - 0.02
V["A4"] = (V.rank_of_true == 1) | (V.P_controlfree_gt_0 >= 0.90)
V["verdict"] = np.where(V.candidate.str.startswith("B"), "comparator",
                        np.where(V.A1 & V.A2 & V.A3 & V.A4, "exploratory external improvement WITH PE-specificity",
                                 np.where(V.A1 & V.A2 & V.A3, "exploratory external rank improvement (not PE-specific)",
                                          np.where(V.A3 & ~V.A1, "no external improvement", "fails internal or bootstrap criterion"))))
V.to_csv(TAB / "intervention_verdicts.csv", index=False)

# ---------------- MASH classifier: stratified mouse bootstrap of saved Phase 4E predictions ----------------
from sklearn.metrics import balanced_accuracy_score, brier_score_loss, roc_auc_score  # noqa: E402

ep = pd.read_csv(RS / "roadmap_phase4e/results/external_predictions.csv")
mrows = []
sel = [("R0", "R0_rule"), ("R1", "logistic_default"), ("R2", "logistic_default"), ("R2", "l2_logistic"), ("R0", "l2_logistic"), ("PCA", "l2_logistic")]
rng = np.random.default_rng(SEED)
for sp in ("S1_train316TA_test246", "S2_train246_test316TA"):
    for rep, mo in sel:
        d = ep[(ep.split == sp) & (ep.rep == rep) & (ep.model == mo)]
        y, p, s = d.y.values, d.p.values, d.score.values
        pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
        au, ba, br = [], [], []
        for _ in range(2000):
            ii = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
            au.append(roc_auc_score(y[ii], s[ii]))
            ba.append(balanced_accuracy_score(y[ii], (p[ii] >= .5).astype(int)))
            br.append(brier_score_loss(y[ii], p[ii]))
        hard = mo == "R0_rule"
        mrows.append(dict(split=sp, rep=rep, model=mo, n=len(y), auc=roc_auc_score(y, s), auc_lo=np.quantile(au, .025), frac_boot_auc_eq_1=float(np.mean(np.isclose(au, 1))),
                          bal_acc=balanced_accuracy_score(y, (p >= .5).astype(int)), bal_acc_lo=np.quantile(ba, .025), bal_acc_hi=np.quantile(ba, .975),
                          metric_prob="misclassification rate (hard 0/1 rule)" if hard else "Brier", value=brier_score_loss(y, p),
                          value_lo=np.quantile(br, .025), value_hi=np.quantile(br, .975)))
pd.DataFrame(mrows).to_csv(TAB / "mash_classifier_bootstrap.csv", index=False)

json_dump(dict(B_boot=B_BOOT, n_common_genes=len(COMMON), n_hong_genes=len(GH), b0_reproduces_phase5_frozen_max_abs=repro_b0,
               seed=SEED, workers=N_JOBS, model_hashes={f.name: sha256(f) for f in sorted(MOD.glob("*.joblib"))}), RES / "b05_run_info.json")
with pd.option_context("display.width", 250, "display.max_columns", 40):
    print(V[["candidate", "train_spearman", "cv_spearman", "cv_r2", "rho_hong_pe", "delta_rho_hong_vs_B0", "P_delta_gt_0", "rank_of_true",
             "rho_controlfree_median", "P_controlfree_gt_0", "A1", "A2", "A3", "A4", "verdict"]].round(3).to_string(index=False))
    print(EV.pivot(index="candidate", columns="outcome", values="spearman").round(3).to_string())
    print(pd.DataFrame(hp_rows).to_string(index=False))
    print(pd.DataFrame(mrows).round(3).to_string(index=False))
