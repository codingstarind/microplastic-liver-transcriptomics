"""Interaction follow-up, step 1 (exploratory; 2026-11-25). Interpreter: C:/Python314/python.exe (has catboost 1.2.10;
same scikit-learn 1.9.1 / numpy 2.5.1 as the project .venv, which lacks catboost).

For each diet (PE_in_CD, PE_in_MD) x Phase 4 outer split (20 modules, 5 folds x 3 repeats) x model (catboost, hgb):
 1. Refit the Phase 4 model on the training genes with that split's STORED training-only tuned hyperparameters
    (metrics_primary_positive.csv best_params) and DL weights, exactly as roadmap_phase4/scripts/p4_run.py; check the
    held-out predictions against the saved Phase 4 predictions.
 2. Training-fold diagnostics (training genes only, fixed random subsample n = 200):
    Friedman H^2 and absolute pure-interaction SD for all 28 feature pairs (partial dependence on the subsample);
    CatBoost only: mean |SHAP interaction value| (native ShapInteractionValues) on a training subsample n = 1000,
    and on lipogenesis / PPARa-target / cholesterol-synthesis genes that are in the training fold.
 3. Interaction-free null (R = 10 per split): fit an additive HGB (interaction_cst='no_interactions', same params as the
    split's HGB) on the training genes; y* = additive fit + permuted training residuals; refit the same model class
    with the same params on y*; recompute H^2 / interaction SD on the same subsample. Any interaction here is spurious.
 4. Held-out additive comparator: HGB with interaction_cst='no_interactions' (same params) and CatBoost with depth = 1
    (same learning rate / l2 / iterations); predictions on the held-out genes.
No hyperparameter is chosen here; nothing uses held-out genes except prediction.
Run: C:/Python314/python.exe research/interaction_followup/scripts/ix01_fit_and_diagnose.py
"""
import itertools
import json
import os
import sys
import time
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"
RS = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RS / "roadmap_phase4/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from sklearn.base import clone  # noqa: E402

from p4_data import dl_tau2, feature_names, feature_table, outer_splits, target  # noqa: E402
from p4_models import make  # noqa: E402

OUT = RS / "interaction_followup"
RES = OUT / "results"
RES.mkdir(parents=True, exist_ok=True)
TARGETS = ("PE_in_CD", "PE_in_MD")
MODELS = ("catboost", "hgb")
N_H, N_SHAP, R_NULL, SEED = 200, 1000, 10, 20261125
FEATS = feature_names("primary")
PAIRS = list(itertools.combinations(range(len(FEATS)), 2))
GS = pd.read_csv(RS / "roadmap_phase3/results/gene_sets.csv")
SETS = {"LIPO": "LIPOGENESIS_CORE", "PPARA": "SANDERSON_PPARA_TARGETS", "CHOL": "CHOL_SYNTHESIS_CORE"}
METR = pd.read_csv(RS / "roadmap_phase4/results/taskA/metrics_primary_positive.csv")
_C = {}


def data(t):
    if t not in _C:
        X = feature_table("primary")
        j = X.join(target(t), how="inner").dropna(subset=["y", "s"])
        _C[t] = j[j["s"] > 0]
    return _C[t]


def build(model, params):
    est, _ = make(model, FEATS)
    est.set_params(**params)
    return est


def pd_table(est, Xs):
    """Centered 1-D and 2-D partial dependence on the subsample points; returns H2 and interaction SD per pair."""
    n, p = Xs.shape

    def pdv(cols):
        Z = np.repeat(Xs[None, :, :], n, axis=0)  # (i, l, p)
        for c in cols:
            Z[:, :, c] = Xs[:, c][:, None]
        f = est.predict(Z.reshape(-1, p)).reshape(n, n).mean(1)
        return f - f.mean()
    one = {c: pdv([c]) for c in range(p)}
    rows = []
    for a, b in PAIRS:
        jk = pdv([a, b])
        inter = jk - one[a] - one[b]
        rows.append(dict(f1=FEATS[a], f2=FEATS[b], H2=float((inter**2).sum() / max((jk**2).sum(), 1e-300)),
                         inter_sd=float(np.sqrt((inter**2).mean())), joint_sd=float(np.sqrt((jk**2).mean()))))
    main = {FEATS[c]: float(np.sqrt((one[c] ** 2).mean())) for c in range(p)}
    return rows, main


def shap_int(est, X, idx_sets):
    from catboost import Pool
    m = est.named_steps["model"]
    out = {}
    for name, ix in idx_sets.items():
        if len(ix) == 0:
            continue
        v = m.get_feature_importance(Pool(X[ix]), type="ShapInteractionValues")  # (n, p+1, p+1)
        v = v[:, :-1, :-1]
        for a, b in PAIRS:
            out[(name, FEATS[a], FEATS[b])] = (float(np.abs(v[:, a, b] + v[:, b, a]).mean()), len(ix))
        for a in range(len(FEATS)):
            out[(name, FEATS[a], "MAIN")] = (float(np.abs(v[:, a, a]).mean()), len(ix))
    return out


def task(t, model, r_i, k_i):
    t0 = time.time()
    j = data(t)
    for r, k, tr, te, grp in outer_splits(j.index, "module"):
        if (r, k) == (r_i, k_i):
            break
    X, y, s = j[FEATS].values, j["y"].values, j["s"].values
    w = 1.0 / (s**2 + dl_tau2(y[tr], s[tr]))
    mrow = METR[(METR.run == f"primary__{t}") & (METR.model == model) & (METR.repeat == r_i) & (METR.fold == k_i)].iloc[0]
    params = json.loads(mrow.best_params)
    hrow = METR[(METR.run == f"primary__{t}") & (METR.model == "hgb") & (METR.repeat == r_i) & (METR.fold == k_i)].iloc[0]
    hparams = json.loads(hrow.best_params)
    est = build(model, params).fit(X[tr], y[tr], model__sample_weight=w[tr])
    pred = est.predict(X[te])
    trI = np.where(tr)[0]
    rng = np.random.default_rng(SEED + 1000 * r_i + 10 * k_i + (0 if t == "PE_in_CD" else 5))
    sub = rng.choice(trI, N_H, replace=False)
    Hrows, main = pd_table(est, X[sub])
    meta = dict(target=t, model=model, repeat=r_i, fold=k_i)
    H = [dict(**meta, kind="observed", rep_null=-1, **h) for h in Hrows]
    M = [dict(**meta, kind="observed", rep_null=-1, feature=f, main_sd=v) for f, v in main.items()]
    # CatBoost SHAP interactions on training genes
    SH = []
    if model == "catboost":
        genes_tr = j.index.values[trI]
        idx = {"ALL_train_sub": rng.choice(trI, N_SHAP, replace=False)}
        for g, sname in SETS.items():
            idx[f"{g}_train"] = trI[np.isin(genes_tr, GS.gene[GS.set == sname])]
        for (name, f1, f2), (val, n) in shap_int(est, X, idx).items():
            SH.append(dict(**meta, genes=name, n_genes=n, f1=f1, f2=f2, mean_abs_shap=val))
    # additive comparator fits (held-out predictions)
    add_h = build("hgb", hparams)
    add_h.set_params(model__interaction_cst="no_interactions")
    add_h.fit(X[tr], y[tr], model__sample_weight=w[tr])
    if model == "catboost":
        add_m = build("catboost", {**params, "model__depth": 1})
        add_m.fit(X[tr], y[tr], model__sample_weight=w[tr])
    else:
        add_m = add_h
    preds = pd.DataFrame({"gene": j.index[te], "y": y[te], "target": t, "model": model, "repeat": r_i, "fold": k_i,
                          "module": grp[te], "pred_full": pred, "pred_additive": add_m.predict(X[te]),
                          "null": np.average(y[tr], weights=w[tr])})
    # interaction-free null
    fa = add_h.predict(X[tr])
    res = y[tr] - fa
    for q in range(R_NULL):
        ys = fa + res[rng.permutation(len(res))]
        e0 = clone(est).fit(X[tr], ys, model__sample_weight=w[tr])
        h0, m0 = pd_table(e0, X[sub])
        H += [dict(**meta, kind="null_additive", rep_null=q, **h) for h in h0]
        M += [dict(**meta, kind="null_additive", rep_null=q, feature=f, main_sd=v) for f, v in m0.items()]
    return dict(H=pd.DataFrame(H), M=pd.DataFrame(M), SH=pd.DataFrame(SH), P=preds,
                info=dict(**meta, params=json.dumps(params), seconds=round(time.time() - t0, 1)))


if __name__ == "__main__":
    jobs = [(t, m, r, k) for t in TARGETS for m in MODELS for r in range(3) for k in range(5)]
    jobs.sort(key=lambda a: a[1] != "catboost")
    out = Parallel(n_jobs=int(os.environ.get("IX_JOBS", 14)), backend="loky", verbose=5, batch_size=1)(delayed(task)(*a) for a in jobs)
    H = pd.concat([o["H"] for o in out], ignore_index=True)
    M = pd.concat([o["M"] for o in out], ignore_index=True)
    SH = pd.concat([o["SH"] for o in out], ignore_index=True)
    P = pd.concat([o["P"] for o in out], ignore_index=True)
    H.to_csv(RES / "friedman_H_training_folds.csv.gz", index=False)
    M.to_csv(RES / "main_effect_sd_training_folds.csv.gz", index=False)
    SH.to_csv(RES / "catboost_shap_interactions_training_folds.csv.gz", index=False)
    P.to_csv(RES / "heldout_full_vs_additive_predictions.csv.gz", index=False)
    pd.DataFrame([o["info"] for o in out]).to_csv(RES / "fit_info.csv", index=False)
    # refit check vs saved Phase 4 predictions
    chk = []
    for t in TARGETS:
        sv = pd.read_csv(RS / f"roadmap_phase4/results/taskA/predictions/primary__{t}.csv.gz")
        for m in MODELS:
            a = sv[sv.model == m].set_index(["gene", "repeat", "fold"]).pred
            b = P[(P.target == t) & (P.model == m)].set_index(["gene", "repeat", "fold"]).pred_full
            d = (a.loc[b.index] - b).abs()
            chk.append(dict(target=t, model=m, n=len(d), max_abs_diff=float(d.max()), mean_abs_diff=float(d.mean()),
                            corr=float(np.corrcoef(a.loc[b.index], b)[0, 1])))
    C = pd.DataFrame(chk)
    C.to_csv(RES / "refit_check_vs_saved_phase4.csv", index=False)
    print(C.to_string(index=False))
