"""Step 5 (Q2): frozen mouse-level scores on independent plastic-exposed mice (locked_plan.md section 5).
Phase 4D singscores (TPM-proxy over frozen background n external genes) and Phase 4E R0/R2 L2-logistic log-odds
(frozen parameters, logCPM ranks over the 4E background n external genes). Within-study contrasts only:
Welch CI, Hedges g, exact permutation over all relabellings, leave-one-mouse-out signs, technical correlates,
R2 feature decomposition, and a coverage control that re-scores GSE294180 with each external gene background."""
import itertools
import json

import numpy as np
import pandas as pd
from p5_common import FRZ, P5, REPO, RES, RS, SETS12, WORK, holm, load_prepared, singscore_signed, singscore_up, welch_ci
from scipy import stats

length = pd.read_csv(RS / "roadmap_phase4d/results/cache/gene_lengths.csv", index_col=0)["length"]
BG4D = pd.read_csv(FRZ / "phase4d_background_genes.csv").gene.tolist()
SIGS = {k: pd.read_csv(RS / f"roadmap_phase4d/signatures/{k}_main_n50.csv") for k in ("FAST", "MASH", "PPARA")}
PAR = {r: json.loads((FRZ / f"phase4e_{r}_l2_logistic_params.json").read_text()) for r in ("R0", "R2")}
QC = pd.read_csv(RES / "external_sample_qc.csv").set_index("mouse_id")
CHOL_LIKE = ["R1_CHOL_SYNTHESIS_CORE", "R1_HALLMARK_CHOLESTEROL_HOMEOSTASIS", "R1_LIPOGENESIS_CORE", "R1_SANDERSON_PPARA_TARGETS",
             "R1_HALLMARK_MTORC1_SIGNALING"]
DISEASE_LIKE = ["R0_signature", "R1_HALLMARK_INFLAMMATORY_RESPONSE", "R1_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]


def score_4d(X: pd.DataFrame, kind: str, background=None):
    """X genes x samples. kind 'counts' -> TPM-proxy (vM35 median transcript length); 'norm' -> rank as deposited."""
    bg = [g for g in (background or BG4D) if g in X.index]
    M = X.loc[bg].astype(float)
    if kind == "counts":
        M = M.div(length[bg], axis=0)
        M = M / M.sum() * 1e6
    R = np.apply_along_axis(stats.rankdata, 1, np.log2(M.T.values + 1))
    gi = {g: i for i, g in enumerate(bg)}
    out, cov = {}, {}
    for k, s in SIGS.items():
        up = np.array([gi[g] for g in s[s.direction == "up"].gene if g in gi])
        dn = np.array([gi[g] for g in s[s.direction == "down"].gene if g in gi])
        out[k] = singscore_signed(R, up, dn)
        cov[k] = dict(n_up=len(up), n_down=len(dn), n_missing=100 - len(up) - len(dn))
    return pd.DataFrame(out, index=X.columns), cov, len(bg)


def score_4e(X: pd.DataFrame, rep: str, background=None):
    par = PAR[rep]
    bg = [g for g in (background or par["background_genes"]) if g in X.index]
    R = np.apply_along_axis(stats.rankdata, 1, X.loc[bg].T.values.astype(float))
    gi = {g: i for i, g in enumerate(bg)}
    F = [singscore_signed(R, np.array([gi[g] for g in par["up"] if g in gi]), np.array([gi[g] for g in par["dn"] if g in gi]))]
    cov = {"R0_signature": sum(g in gi for g in par["up"] + par["dn"])}
    for s, gl in par["sets"].items():
        idx = np.array([gi[g] for g in gl if g in gi])
        F.append(singscore_up(R, idx))
        cov[f"R1_{s}"] = len(idx)
    F = np.column_stack(F)
    Z = (F - np.array(par["scaler_mean"])) / np.array(par["scaler_scale"])
    lo = Z @ np.array(par["coef"]) + par["intercept"]
    return pd.Series(lo, index=X.columns), pd.DataFrame(F, index=X.columns, columns=par["feature_names"]), cov, len(bg)


# ---------------- score every external mouse ----------------
rows, feats, coverage = [], {}, []
for st, kind in (("HONG", "counts"), ("GSE245069", "counts"), ("GSE277772", "norm")):
    X, S = load_prepared(st)
    d4, cov4, nb4 = score_4d(X, kind)
    df = d4.rename(columns=lambda c: f"{c}_singscore")
    for rep in ("R0", "R2"):
        lo, F, cov, nb = score_4e(X, rep)
        df[f"{rep}_logistic_logodds"] = lo
        if rep == "R2":
            feats[st] = F
            coverage.append(dict(study=st, score="Phase4E_R0R2", background_genes=nb, background_frozen=len(PAR[rep]["background_genes"]), **cov))
    if st == "GSE245069":  # sensitivity: deposited TPM instead of the TPM-proxy
        T = pd.read_csv(WORK / "GSE245069_tpm.csv.gz", index_col=0)[X.columns]
        dt, _, _ = score_4d(T, "norm")
        for c in dt:
            df[f"{c}_singscore_depositedTPM"] = dt[c]
    for k, c in cov4.items():
        coverage.append(dict(study=st, score=f"Phase4D_{k}", background_genes=nb4, background_frozen=len(BG4D), **c))
    df = df.join(S[["group", "polymer", "size", "role"]]).join(QC[["counted_total", "detected_genes", "muscle_cpm_sum", "sra_read_pairs"]])
    df.insert(0, "study", st)
    rows.append(df)
SC = pd.concat(rows)
SC.index.name = "mouse_id"
SC.to_csv(P5 / "per_sample_scores.csv")
pd.DataFrame(coverage).to_csv(RES / "score_gene_coverage.csv", index=False)

# ---------------- GSE294180 reference (same functions, full background) ----------------
g294 = pd.read_csv(REPO / "data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz", index_col=0)
sm = pd.read_csv(REPO / "metadata/sample_metadata.csv").set_index("counts_column").loc[g294.columns]


def g294_pe_effects(background4d=None, background4e=None):
    d4, _, _ = score_4d(g294, "counts", background4d)
    out = {}
    for k in d4:
        v = d4[k]
        out[f"{k}_singscore"] = {c: v[(sm.diet == dd) & (sm.microplastic == "PE")].mean() - v[(sm.diet == dd) & (sm.microplastic == "Veh")].mean()
                                 for c, dd in (("PE_in_CD", "CD"), ("PE_in_MD", "MD"))}
    for rep in ("R0", "R2"):
        v, _, _, _ = score_4e(g294, rep, background4e)
        out[f"{rep}_logistic_logodds"] = {c: v[(sm.diet == dd) & (sm.microplastic == "PE")].mean() - v[(sm.diet == dd) & (sm.microplastic == "Veh")].mean()
                                          for c, dd in (("PE_in_CD", "CD"), ("PE_in_MD", "MD"))}
    return out


REF = g294_pe_effects()
# sanity: reproduces Phase 4D / 4E saved estimates
p4d = pd.read_csv(RS / "roadmap_phase4d/results/contrast_estimates_primary.csv")
chk = {k: abs(REF[f"{k}_singscore"]["PE_in_CD"] - p4d[(p4d.score == f"{k}|tpm_proxy|primary") & (p4d.contrast == "PE_in_CD")].estimate.iloc[0])
       for k in SIGS}
assert max(chk.values()) < 1e-9, chk
cov_ctrl = []
for st in ("HONG", "GSE245069", "GSE277772"):
    X, _ = load_prepared(st)
    r = g294_pe_effects([g for g in BG4D if g in X.index], [g for g in PAR["R0"]["background_genes"] if g in X.index])
    for sc in REF:
        for c in ("PE_in_CD", "PE_in_MD"):
            cov_ctrl.append(dict(external_background=st, score=sc, contrast=c, gse294180_full=REF[sc][c], gse294180_restricted=r[sc][c],
                                 rel_change=(r[sc][c] - REF[sc][c]) / abs(REF[sc][c])))
CC = pd.DataFrame(cov_ctrl)
CC.to_csv(RES / "coverage_control_gse294180_rescored.csv", index=False)

# ---------------- within-study contrasts ----------------
CONTR = [("HONG", ["PE"], "PE_vs_CON", "PRIMARY"),
         ("HONG", ["PS"], "PS_vs_CON", "exploratory specificity"), ("HONG", ["PP"], "PP_vs_CON", "exploratory specificity"),
         ("HONG", ["PB"], "PB_vs_CON", "exploratory specificity"),
         ("HONG", ["PS", "PP", "PB"], "nonPE_polymers_vs_CON", "exploratory specificity (pooled other polymers)"),
         ("GSE245069", ["PS_80nm", "PS_0.5um", "PS_5um"], "PSpooled_vs_CON", "SECONDARY cross-polymer"),
         ("GSE245069", ["PS_80nm"], "PS80nm_vs_CON", "exploratory size"), ("GSE245069", ["PS_0.5um"], "PS05um_vs_CON", "exploratory size"),
         ("GSE245069", ["PS_5um"], "PS5um_vs_CON", "exploratory size"),
         ("GSE277772", ["PS"], "PS_vs_CON", "exploratory (normalised values)")]
SCORES = ["FAST_singscore", "MASH_singscore", "PPARA_singscore", "R0_logistic_logodds", "R2_logistic_logodds"]
p4e = pd.read_csv(RS / "roadmap_phase4e/results/r2_feature_contributions_gse294180.csv").set_index("feature")
coef_sd = pd.Series(np.array(PAR["R2"]["coef"]) / np.array(PAR["R2"]["scaler_scale"]), index=PAR["R2"]["feature_names"])
res, decomp = [], []
for st, ex, name, role in CONTR:
    D = SC[(SC.study == st) & SC.group.isin(["CON"] + ex)]
    yv = D.group.isin(ex).values
    for sc in SCORES + ([f"{k}_singscore_depositedTPM" for k in ("FAST", "MASH")] if st == "GSE245069" else []):
        v = D[sc].values
        w = welch_ci(v[yv], v[~yv])
        obs = w["estimate"]
        null = []
        n, k = len(v), int((~yv).sum())
        for ctrl in itertools.combinations(range(n), k):
            m = np.ones(n, bool)
            m[list(ctrl)] = False
            null.append(v[m].mean() - v[~m].mean())
        null = np.array(null)
        loo_signs = []
        for i in range(n):
            mm = np.arange(n) != i
            vv, yy = v[mm], yv[mm]
            loo_signs.append(np.sign(vv[yy].mean() - vv[~yy].mean()))
        base = sc.replace("_depositedTPM", "")
        ref_cd = REF[base]["PE_in_CD"]
        tech = {f"rho_{c}": stats.spearmanr(v, D[c].values).statistic for c in ("counted_total", "detected_genes", "muscle_cpm_sum")}
        r = dict(study=st, contrast=name, role=role, score=sc, n_exposed=int(yv.sum()), n_control=int((~yv).sum()), **w,
                 perm_n_labellings=len(null), p_perm_exact_two_sided=float(np.mean(np.abs(null) >= abs(obs) - 1e-12)),
                 p_perm_min_attainable=float(2 / len(null)) if len(null) > 1 else np.nan,
                 loo_n_same_sign=int(np.sum(np.array(loo_signs) == np.sign(obs))), loo_n=n,
                 gse294180_PE_in_CD=ref_cd, gse294180_PE_in_MD=REF[base]["PE_in_MD"],
                 same_sign_as_gse294180_PE_in_CD=bool(np.sign(obs) == np.sign(ref_cd)),
                 ratio_to_gse294180_PE_in_CD=obs / ref_cd, **tech)
        res.append(r)
    # R2 decomposition: contribution_j = coef_j/SD_j * (mean feature diff)
    F = feats[st].loc[D.index]
    diff = F[yv].mean() - F[~yv].mean()
    contrib = coef_sd * diff
    tot = contrib.sum()
    for f_, c_ in contrib.items():
        decomp.append(dict(study=st, contrast=name, feature=f_, feature_diff=diff[f_], coef_per_unit=coef_sd[f_], contribution_logodds=c_,
                           gse294180_contrib_PE_in_CD=p4e.loc[f_, "contrib_PE_in_CD"], gse294180_contrib_PE_in_MD=p4e.loc[f_, "contrib_PE_in_MD"]))
    chol = contrib[CHOL_LIKE].sum()
    dis = contrib[DISEASE_LIKE].sum()
    decomp.append(dict(study=st, contrast=name, feature="(TOTAL)", contribution_logodds=tot,
                       share_abs_chol_fasting_like=float(np.abs(contrib[CHOL_LIKE]).sum() / np.abs(contrib).sum()),
                       share_abs_disease_like=float(np.abs(contrib[DISEASE_LIKE]).sum() / np.abs(contrib).sum()),
                       sum_chol_fasting_like=chol, sum_disease_like=dis))
R = pd.DataFrame(res)
R["p_welch_holm_secondary3"] = np.nan
for (st, name), d in R.groupby(["study", "contrast"]):
    i = d.index[d.score.isin(["MASH_singscore", "R0_logistic_logodds", "R2_logistic_logodds"])]
    R.loc[i, "p_welch_holm_secondary3"] = holm(R.loc[i, "p_welch"].values)


def q2_verdict(r):
    if r.score != "FAST_singscore" or r.contrast != "PE_vs_CON":
        return ""
    if (not r.same_sign_as_gse294180_PE_in_CD) or abs(r.estimate / r.se) < 1:
        return "FAIL"
    if r.loo_n_same_sign == r.loo_n and r.p_perm_exact_two_sided <= r.p_perm_min_attainable + 1e-12:
        return "DIRECTIONALLY CONSISTENT"
    return "INCONCLUSIVE"


R["primary_verdict_locked_rule"] = R.apply(q2_verdict, axis=1)
R.to_csv(P5 / "sample_level_external_results.csv", index=False)
DC = pd.DataFrame(decomp)
DC.to_csv(RES / "r2_feature_contributions_external.csv", index=False)
with pd.option_context("display.width", 300, "display.max_rows", 300):
    print(R[["study", "contrast", "score", "estimate", "ci_low", "ci_high", "hedges_g", "p_perm_exact_two_sided", "p_perm_min_attainable",
             "loo_n_same_sign", "loo_n", "gse294180_PE_in_CD", "same_sign_as_gse294180_PE_in_CD", "rho_counted_total", "primary_verdict_locked_rule"]].round(4).to_string(index=False))
    print(DC[DC.feature == "(TOTAL)"].round(3).to_string(index=False))
    print(CC.groupby(["external_background"]).rel_change.apply(lambda x: x.abs().max()).round(4))
    print(pd.DataFrame(coverage).to_string(index=False))
