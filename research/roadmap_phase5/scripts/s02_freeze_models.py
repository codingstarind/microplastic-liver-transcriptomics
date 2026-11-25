"""Step 2: freeze + serialise every model used in Phase 5 and verify it reproduces earlier-phase outputs.

Gene level  : F8-OLS and F8+C-Ridge(alpha=1000) fitted on ALL GSE294180 universe genes (targets PE_in_CD, PE_in_MD),
              plus single-axis comparators (FAST_eff, MASH_eff). Predictions for every feature-table gene are saved.
Mouse level : Phase 4E R0 / R2 L2-logistic (C=10) refitted exactly as e03 on the 31 MASH mice, exported as explicit
              parameters (gene lists, scaler, coefficients) and re-applied by a pure function that must reproduce
              results/gse294180_frozen_model_scores.csv; Phase 4D signatures hash-checked against the manifest.
No external data are read in this script."""
import sys

import joblib
import numpy as np
import pandas as pd
from p5_common import FRZ, REPO, RS, json_dump, sha256, singscore_signed, singscore_up, versions
from scipy import stats
from sklearn.linear_model import LinearRegression, LogisticRegression, Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

FRZ.mkdir(exist_ok=True)
sys.path.insert(0, str(RS / "roadmap_phase4c_diagnostics" / "scripts"))
from c00_common import RAW_COLS, REPS, RepBuilder, dl_tau2, raw_table, target  # noqa: E402

man = {"versions": versions(), "files": {}}

# ---------------- gene level ----------------
X0 = raw_table()
gene_models = [("F8_ols", "F8", LinearRegression()), ("F8C_ridge", "F8+C", Ridge(alpha=1000)),
               ("FAST_only_ols", None, LinearRegression()), ("MASH_only_ols", None, LinearRegression())]
preds, coefs = {}, []
for tgt in ("PE_in_CD", "PE_in_MD"):
    j = X0.join(target(tgt), how="inner").dropna(subset=["y", "s"])
    j = j[j.s > 0]
    y, s = j.y.values, j.s.values
    tau2 = dl_tau2(y, s)
    w = 1 / (s**2 + tau2)
    for name, rep, est in gene_models:
        if rep is None:
            col = name.split("_")[0] + "_eff"
            pipe = Pipeline([("scale", StandardScaler()), ("model", est)])
            pipe.fit(j[[col]].values, y, model__sample_weight=w)
            p_all = pipe.predict(X0[[col]].values)
            feats = [col]
        else:
            pipe = Pipeline([("rep", RepBuilder(rep)), ("scale", StandardScaler()), ("model", est)])
            pipe.fit(j[RAW_COLS].values, y, model__sample_weight=w)
            p_all = pipe.predict(X0[RAW_COLS].values)
            feats = REPS[rep]
        key = f"{name}__{tgt}"
        f = FRZ / f"gene_{key}.joblib"
        joblib.dump(pipe, f)
        man["files"][f.name] = sha256(f)
        preds[key] = pd.Series(p_all, index=X0.index)
        fit_r2 = 1 - np.sum(w * (y - pipe.predict(j[RAW_COLS].values if rep else j[[col]].values)) ** 2) / np.sum(w * (y - np.average(y, weights=w)) ** 2)
        coefs += [dict(model=key, feature=ft, coef_per_SD=c) for ft, c in zip(feats, pipe.named_steps["model"].coef_)]
        coefs.append(dict(model=key, feature="(intercept)", coef_per_SD=float(pipe.named_steps["model"].intercept_)))
        coefs.append(dict(model=key, feature="(in-sample weighted R2)", coef_per_SD=float(fit_r2)))
        coefs.append(dict(model=key, feature="(training weighted mean y)", coef_per_SD=float(np.average(y, weights=w))))
        coefs.append(dict(model=key, feature="(tau2_DL)", coef_per_SD=float(tau2)))
        coefs.append(dict(model=key, feature="(n_genes)", coef_per_SD=len(j)))
P = pd.DataFrame(preds)
# observed GSE294180 PE effects as a direct comparator ("does the original PE vector replicate?")
for tgt in ("PE_in_CD", "PE_in_MD"):
    t = pd.read_csv(REPO / f"results/deseq2/provisional_deposited/results_{tgt}.csv.gz").set_index("gene")
    P[f"observed_GSE294180_ashr__{tgt}"] = t["log2FC_ashr"].reindex(P.index)
P.index.name = "gene"
P.to_csv(FRZ / "gene_level_frozen_predictions.csv.gz")
pd.DataFrame(coefs).to_csv(FRZ / "gene_level_frozen_coefficients.csv", index=False)
man["files"]["gene_level_frozen_predictions.csv.gz"] = sha256(FRZ / "gene_level_frozen_predictions.csv.gz")
# sanity: F8 OLS coefficients equal the Phase 4C per-fold sign pattern
c4 = pd.read_csv(RS / "roadmap_phase4c_diagnostics/tables/fix_experiment_coefficients.csv")
c4 = c4[(c4.rep == "F8") & (c4.model == "ols")].groupby(["target", "feature"]).coef_per_SD.median()
C = pd.DataFrame(coefs)
sign_check = {}
for tgt in ("PE_in_CD", "PE_in_MD"):
    cf = C[C.model == f"F8_ols__{tgt}"].set_index("feature").coef_per_SD
    sign_check[tgt] = bool(all(np.sign(cf[f]) == np.sign(c4[(tgt, f)]) for f in REPS["F8"] if abs(c4[(tgt, f)]) > 1e-3))
man["F8_full_fit_signs_match_phase4c_fold_median"] = sign_check

# ---------------- mouse level: Phase 4E R0/R2 ----------------
sys.path.insert(0, str(RS / "roadmap_phase4e" / "scripts"))
import e01_mash_state_models as E  # noqa: E402

E.build_ctx("mash")
T, Xg = E.CTX["mash"]["T"], E.CTX["mash"]["X"]
genes4e = list(E.CTX["mash"]["genes"])
train = ((T.study == "GSE246328") | ((T.study == "GSE316023") & (T.contrast == "MASH_primary"))).values


class BothStudySignature(E.FeatureBuilder):  # exact copy of e03 logic
    def fit(self, idx, y):
        idx = np.asarray(idx).ravel().astype(int)
        st = T.study.values[idx]
        ts, ds = [], []
        for s_ in ("GSE316023", "GSE246328"):
            jj = idx[st == s_]
            yy = np.asarray(y)[st == s_]
            a, b = Xg[jj][yy == 1], Xg[jj][yy == 0]
            ts.append(np.nan_to_num(stats.ttest_ind(a, b, equal_var=False).statistic))
            ds.append(a.mean(0) - b.mean(0))
        (t1, t2), (d1, d2) = ts, ds
        upc = np.where((d1 >= 1) & (d2 >= 1))[0]
        dnc = np.where((d1 <= -1) & (d2 <= -1))[0]
        self.up_ = upc[np.argsort(-np.minimum(t1[upc], t2[upc]))][:50]
        self.dn_ = dnc[np.argsort(np.maximum(t1[dnc], t2[dnc]))][:50]
        return self


saved = pd.read_csv(RS / "roadmap_phase4e/results/gse294180_frozen_model_scores.csv", index_col=0)
tr = np.where(train)[0]
te = np.where((T.study == "GSE294180").values)[0]
frozen4e, repro = {}, {}
for rep in ("R0", "R2"):
    pipe = Pipeline([("feat", BothStudySignature(rep, "mash")), ("sc", StandardScaler()), ("clf", LogisticRegression(C=10, max_iter=5000))])
    pipe.fit(tr.reshape(-1, 1), T.y.values[tr])
    s_pipe = pipe.decision_function(te.reshape(-1, 1))
    fb = pipe.named_steps["feat"]
    par = dict(rep=rep, C=10, background_genes=genes4e, up=[genes4e[i] for i in fb.up_], dn=[genes4e[i] for i in fb.dn_],
               sets={s: [genes4e[i] for i in E.CTX["mash"]["set_idx"][s]] for s in E.SETS} if rep == "R2" else {},
               feature_names=(["R0_signature"] if rep in ("R0", "R2") else []) + ([f"R1_{s}" for s in E.SETS] if rep == "R2" else []),
               scaler_mean=pipe.named_steps["sc"].mean_.tolist(), scaler_scale=pipe.named_steps["sc"].scale_.tolist(),
               coef=pipe.named_steps["clf"].coef_.ravel().tolist(), intercept=float(pipe.named_steps["clf"].intercept_[0]),
               train_mice=T.key.values[tr].tolist(), train_y=T.y.values[tr].tolist())
    frozen4e[rep] = par
    repro[rep] = dict(max_abs_diff_vs_phase4e=float(np.max(np.abs(s_pipe - saved.loc[T.sample_id.values[te], f"{rep}_l2_logistic"].values))))


def apply_4e(par, counts: pd.DataFrame, background=None):
    """Pure re-application: counts (genes x samples, any monotone per-sample scale) -> (log-odds, feature matrix)."""
    bg = [g for g in (background or par["background_genes"]) if g in counts.index]
    R = np.apply_along_axis(stats.rankdata, 1, counts.loc[bg].T.values.astype(float))
    gi = {g: i for i, g in enumerate(bg)}
    up = np.array([gi[g] for g in par["up"] if g in gi])
    dn = np.array([gi[g] for g in par["dn"] if g in gi])
    F = [singscore_signed(R, up, dn)]
    for s in par["sets"]:
        F.append(singscore_up(R, np.array([gi[g] for g in par["sets"][s] if g in gi])))
    F = np.column_stack(F)
    Z = (F - np.array(par["scaler_mean"])) / np.array(par["scaler_scale"])
    return Z @ np.array(par["coef"]) + par["intercept"], F


g294 = pd.read_csv(REPO / "data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz", index_col=0)
for rep, par in frozen4e.items():
    s_pure, _ = apply_4e(par, g294)
    repro[rep]["pure_function_max_abs_diff_vs_phase4e"] = float(np.max(np.abs(s_pure - saved.loc[g294.columns, f"{rep}_l2_logistic"].values)))
    f = FRZ / f"phase4e_{rep}_l2_logistic_params.json"
    json_dump(par, f)
    man["files"][f.name] = sha256(f)
man["phase4e_reproduction"] = repro
assert all(v["pure_function_max_abs_diff_vs_phase4e"] < 1e-6 for v in repro.values()), repro

# ---------------- Phase 4D signatures + background ----------------
sm = pd.read_csv(RS / "roadmap_phase4d/frozen_signature_manifest.csv").set_index("file")
for fn in ("MASH_main_n50.csv", "FAST_main_n50.csv", "PPARA_main_n50.csv"):
    assert sha256(RS / "roadmap_phase4d/signatures" / fn) == sm.loc[fn, "sha256"], fn
filt = pd.read_csv(REPO / "data/processed/gene_filter.csv").query("passes_filter").gene.tolist()
length = pd.read_csv(RS / "roadmap_phase4d/results/cache/gene_lengths.csv", index_col=0)["length"]
bg4d = [g for g in filt if g in length.index]
pd.Series(bg4d, name="gene").to_csv(FRZ / "phase4d_background_genes.csv", index=False)
man["phase4d_background_n"] = len(bg4d)
man["phase4d_signature_hashes_verified"] = True
man["files"]["phase4d_background_genes.csv"] = sha256(FRZ / "phase4d_background_genes.csv")
json_dump(man, FRZ / "frozen_manifest.json")
print(pd.DataFrame(coefs).query("feature.str.startswith('(')").pivot(index="model", columns="feature", values="coef_per_SD").round(4).to_string())
print("sign check", sign_check, "\n4E repro", repro, "\n4D background", len(bg4d))
