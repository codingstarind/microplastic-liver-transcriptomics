"""Pathway follow-up, step 1 (exploratory, post hoc): 8 prespecified pathway scores in the 19 GSE294180 mice.

Reads existing data only; writes to research/pathway_followup/{results,figures}.
Primary score: singscore (up-set) on within-sample ranks over the frozen Phase 4E background (11,883 genes; the same
representation as the R1 classifier features). FAST = frozen Phase 4D signed FAST singscore (per_mouse_scores.csv).
Secondary score: mean z of the VST matrix over set genes (as analysis/pathways/pathway_analysis.R Method 1); FAST secondary
= mean z(up genes) - mean z(down genes).
Model per score: OLS score ~ diet * PE (19 mice); contrasts PE_in_CD, PE_in_MD, Interaction.
Families: (a) new family = the 8 scores within each contrast (BH and Holm); (b) original fry q-values copied unchanged.
Run: .venv/Scripts/python.exe research/pathway_followup/scripts/pf01_pathways.py
"""
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
RS = ROOT / "research"
OUT = RS / "pathway_followup"
RES, FIG = OUT / "results", OUT / "figures"
RES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(exist_ok=True)
sys.path.insert(0, str(RS / "roadmap_phase5" / "scripts"))
from p5_common import load_prepared, singscore_signed, singscore_up, welch_ci  # noqa: E402

PATHS = {"PPARA_targets": "SANDERSON_PPARA_TARGETS", "mTORC1": "HALLMARK_MTORC1_SIGNALING",
         "Chol_homeostasis": "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "Chol_synthesis": "CHOL_SYNTHESIS_CORE",
         "Glycolysis": "HALLMARK_GLYCOLYSIS", "DNA_repair": "HALLMARK_DNA_REPAIR", "Lipogenesis": "LIPOGENESIS_CORE"}
ORDER = list(PATHS) + ["FAST"]
CONTR = {"PE_in_CD": [0, 0, 1, 0], "PE_in_MD": [0, 0, 1, 1], "Interaction": [0, 0, 0, 1], "Diet_in_Veh": [0, 1, 0, 0]}
GROUPS = ["CD+Veh", "CD+PE", "MD+Veh", "MD+PE"]

gs = pd.read_csv(RS / "roadmap_phase3/results/gene_sets.csv")
SETG = {k: list(dict.fromkeys(gs.gene[gs.set == v])) for k, v in PATHS.items()}
fast_sig = pd.read_csv(RS / "roadmap_phase4d/signatures/FAST_main_n50.csv")
FAST_UP, FAST_DN = list(fast_sig.gene[fast_sig.direction == "up"]), list(fast_sig.gene[fast_sig.direction == "down"])
par = json.loads((RS / "roadmap_phase5/frozen/phase4e_R2_l2_logistic_params.json").read_text())
BG = par["background_genes"]
musc = set(pd.read_csv(RS / "roadmap_phase2/metadata/muscle_exclusion_genes.csv").vm35_symbol)

cnt = pd.read_csv(ROOT / "data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz", index_col=0)
meta = pd.read_csv(ROOT / "metadata/sample_metadata.csv").set_index("counts_column").loc[cnt.columns]
p4d = pd.read_csv(RS / "roadmap_phase4d/results/per_mouse_scores.csv", index_col=0).loc[cnt.columns]
vst = pd.read_csv(ROOT / "results/deseq2/provisional_deposited/vst_matrix.csv.gz", index_col=0)[cnt.columns]


def ss_scores(X, sets, bg=BG):
    """singscore up per set on within-sample ranks over background bg. X genes x samples."""
    bgx = [g for g in bg if g in X.index]
    R = np.apply_along_axis(stats.rankdata, 1, X.loc[bgx].T.values.astype(float))
    gi = {g: i for i, g in enumerate(bgx)}
    out, n = {}, {}
    for k, gl in sets.items():
        idx = np.array([gi[g] for g in gl if g in gi])
        out[k], n[k] = singscore_up(R, idx), len(idx)
    return pd.DataFrame(out, index=X.columns), n


def mz_scores(V, sets):
    Z = V.sub(V.mean(1), axis=0).div(V.std(1, ddof=1), axis=0)
    out = {k: Z.loc[[g for g in gl if g in Z.index]].mean() for k, gl in sets.items()}
    out["FAST"] = Z.loc[[g for g in FAST_UP if g in Z.index]].mean() - Z.loc[[g for g in FAST_DN if g in Z.index]].mean()
    return pd.DataFrame(out)


MD = (meta.diet == "MD").astype(int).values
PE = (meta.microplastic == "PE").astype(int).values
GRP = meta.group.values


def ols(y, md=MD, pe=PE, cov=None):
    cols = [np.ones(len(y)), md, pe, md * pe] + ([] if cov is None else [cov])
    X = np.column_stack(cols)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    df = len(y) - X.shape[1]
    s2 = r @ r / df
    V = s2 * np.linalg.inv(X.T @ X)
    q = stats.t.ppf(.975, df)
    out = {}
    for k, l in CONTR.items():
        l = np.array(l + [0] * (X.shape[1] - 4), float)
        e, se = l @ b, np.sqrt(l @ V @ l)
        out[k] = dict(estimate=e, se=se, ci_low=e - q * se, ci_high=e + q * se, t=e / se, df=df,
                      p=2 * stats.t.sf(abs(e / se), df), std_effect=e / np.sqrt(s2))
    return out


def exact_perm(y, diet):
    """Exact within-diet relabelling p for the PE-vs-Veh mean difference (two-sided, all C(n, n_PE) labellings)."""
    m = meta.diet.values == diet
    yy, pe = y[m], PE[m]
    obs = yy[pe == 1].mean() - yy[pe == 0].mean()
    n, k = len(yy), int(pe.sum())
    d = np.array([yy[list(c)].mean() - np.delete(yy, list(c)).mean() for c in itertools.combinations(range(n), k)])
    return float(np.mean(np.abs(d) >= abs(obs) - 1e-12)), len(d)


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)[::-1]
    a = np.minimum(1, np.minimum.accumulate(p[o] * n / np.arange(n, 0, -1))); out = np.empty(n); out[o] = a
    return out


def holm(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); out = np.empty(n); run = 0
    for i, j in enumerate(o):
        run = max(run, (n - i) * p[j]); out[j] = min(1, run)
    return out


def contrast_table(S, label):
    rows = []
    for k in S.columns:
        y = S[k].values
        for c, d in ols(y).items():
            rows.append(dict(score_type=label, pathway=k, contrast=c, **d))
    T = pd.DataFrame(rows)
    for c in CONTR:
        i = T.contrast == c
        T.loc[i, "p_BH_8"] = bh(T.loc[i, "p"])
        T.loc[i, "p_Holm_8"] = holm(T.loc[i, "p"])
    return T


# ---------------------------------------------------------------- scores
S_ss, n_ss = ss_scores(cnt, SETG)
S_ss["FAST"] = p4d["FAST"]
S_ss = S_ss[ORDER]
S_mz = mz_scores(vst, SETG)[ORDER]
scores = pd.concat([S_ss.add_prefix("ss_"), S_mz.add_prefix("mz_"), meta[["group", "diet", "microplastic"]],
                    p4d[["depth", "quality_proxy", "muscle", "clock_phase"]]], axis=1)
scores.to_csv(RES / "per_mouse_pathway_scores.csv")

T = pd.concat([contrast_table(S_ss, "singscore_primary"), contrast_table(S_mz, "meanz_VST_secondary")], ignore_index=True)
for i, r in T[T.contrast.isin(["PE_in_CD", "PE_in_MD"])].iterrows():
    S = S_ss if r.score_type == "singscore_primary" else S_mz
    T.loc[i, "p_exact_perm"], T.loc[i, "n_labellings"] = exact_perm(S[r.pathway].values, r.contrast[-2:])
fry = pd.read_csv(RS / "enrichment_audit/tables/T2_enrichment_wide_CD_MD_interaction.csv").set_index("pathway")
for i, r in T.iterrows():
    nm = PATHS.get(r.pathway)
    if nm in fry.index and r.contrast in ("PE_in_CD", "PE_in_MD", "Interaction"):
        T.loc[i, "orig_fry_direction"] = fry.loc[nm, f"{r.contrast}__fry_direction"]
        T.loc[i, "orig_fry_q_family"] = fry.loc[nm, f"{r.contrast}__fry_fdr"]
        T.loc[i, "orig_fry_family"] = "PPAR (3 sets)" if nm.startswith("SANDERSON") else ("Hallmark (50 sets)" if nm.startswith("HALLMARK") else "")
for col in ("p_BH_8", "orig_fry_q_family"):
    T[f"{col}_lt_0.05"] = T[col] < 0.05
    T[f"{col}_lt_0.10"] = T[col] < 0.10
T.to_csv(RES / "pathway_contrasts.csv", index=False)

# ---------------------------------------------------------------- gene-set overlap
allsets = {**SETG, "FAST_up": FAST_UP, "FAST_down": FAST_DN}
bgs = set(BG)
ov = []
for a, b in itertools.product(allsets, allsets):
    A, B = set(allsets[a]) & bgs, set(allsets[b]) & bgs
    ov.append(dict(set_a=a, set_b=b, n_a=len(A), n_b=len(B), n_shared=len(A & B),
                   jaccard=len(A & B) / max(1, len(A | B)), frac_of_a_in_b=len(A & B) / max(1, len(A))))
ov = pd.DataFrame(ov)
ov.to_csv(RES / "gene_set_overlap.csv", index=False)

# unique-gene scores: each set minus genes in any other of the 8 sets (incl. FAST up/down)
uniq = {k: [g for g in gl if not any(g in allsets[o] for o in allsets if o != k)] for k, gl in SETG.items()}
S_u, n_u = ss_scores(cnt, uniq)
S_u["FAST"] = p4d["FAST"]
Tu = contrast_table(S_u[ORDER], "singscore_unique_genes")
Tu["n_genes_unique"] = Tu.pathway.map(n_u)
Tu["n_genes_full"] = Tu.pathway.map(n_ss)
Tu.to_csv(RES / "pathway_contrasts_unique_genes.csv", index=False)

# ---------------------------------------------------------------- correlations between scores
resid = S_ss.groupby(GRP).transform(lambda x: x - x.mean())
cor = []
for a, b in itertools.combinations(ORDER, 2):
    for lab, D in (("all_19_mice", S_ss), ("within_group_residuals", resid), ("unique_genes_all_19", S_u[ORDER]),
                   ("unique_genes_within_group", S_u[ORDER].groupby(GRP).transform(lambda x: x - x.mean()))):
        r = stats.spearmanr(D[a], D[b])
        rp = stats.pearsonr(D[a], D[b])
        cor.append(dict(basis=lab, a=a, b=b, spearman=r.statistic, spearman_p=r.pvalue, pearson=rp.statistic, pearson_p=rp.pvalue))
cor = pd.DataFrame(cor)
cor["spearman_p_BH_within_basis"] = cor.groupby("basis").spearman_p.transform(bh)
cor.to_csv(RES / "pathway_score_correlations.csv", index=False)

# PE-effect vector agreement: per-mouse PE deviation (score - mean of same-diet vehicle), standardized by vehicle SD
dev = []
for d in ("CD", "MD"):
    veh = (meta.diet == d) & (meta.microplastic == "Veh")
    pe = (meta.diet == d) & (meta.microplastic == "PE")
    for k in ORDER:
        mu, sd = S_ss.loc[veh, k].mean(), S_ss.loc[veh, k].std(ddof=1)
        vmin, vmax = S_ss.loc[veh, k].min(), S_ss.loc[veh, k].max()
        for m in S_ss.index[pe]:
            v = S_ss.loc[m, k]
            dev.append(dict(diet=d, mouse=m, pathway=k, score=v, z_vs_vehicle=(v - mu) / sd,
                            below_vehicle_min=v < vmin, above_vehicle_max=v > vmax))
dev = pd.DataFrame(dev)
dev.to_csv(RES / "per_PE_mouse_deviation.csv", index=False)
dsum = dev.groupby(["diet", "pathway"]).agg(n_PE=("mouse", "size"), n_z_pos=("z_vs_vehicle", lambda z: int((z > 0).sum())),
                                           n_z_neg=("z_vs_vehicle", lambda z: int((z < 0).sum())),
                                           n_below_vehicle_min=("below_vehicle_min", "sum"), n_above_vehicle_max=("above_vehicle_max", "sum"),
                                           z_min=("z_vs_vehicle", "min"), z_max=("z_vs_vehicle", "max")).reset_index()

# ---------------------------------------------------------------- robustness: LOO, Cook's D, covariates, muscle genes
loo = []
for k in ORDER:
    y = S_ss[k].values
    full = ols(y)
    for i in range(19):
        m = np.arange(19) != i
        o = ols(y[m], MD[m], PE[m])
        for c in ("PE_in_CD", "PE_in_MD", "Interaction"):
            loo.append(dict(pathway=k, dropped=S_ss.index[i], contrast=c, estimate=o[c]["estimate"], p=o[c]["p"],
                            estimate_full=full[c]["estimate"]))
loo = pd.DataFrame(loo)
loo_s = loo.groupby(["pathway", "contrast"]).apply(lambda d: pd.Series(dict(
    estimate_full=d.estimate_full.iloc[0], loo_min=d.estimate.min(), loo_max=d.estimate.max(),
    n_same_sign=int((np.sign(d.estimate) == np.sign(d.estimate_full)).sum()), loo_max_p=d.p.max(),
    most_influential=d.dropped.iloc[int(np.argmax(np.abs(d.estimate - d.estimate_full)))])), include_groups=False).reset_index()
loo_s.to_csv(RES / "loo_summary.csv", index=False)

cook = []
Xd = np.column_stack([np.ones(19), MD, PE, MD * PE])
H = Xd @ np.linalg.inv(Xd.T @ Xd) @ Xd.T
h = np.diag(H)
for k in ORDER:
    y = S_ss[k].values
    r = y - H @ y
    s2 = r @ r / (19 - 4)
    D = r**2 / (4 * s2) * h / (1 - h) ** 2
    for i, m in enumerate(S_ss.index):
        cook.append(dict(pathway=k, mouse=m, group=GRP[i], cooks_d=D[i], leverage=h[i], flag_gt_4_over_n=D[i] > 4 / 19))
pd.DataFrame(cook).to_csv(RES / "cooks_distance.csv", index=False)

sens = []
for k in ORDER:
    y = S_ss[k].values
    for cv in ("depth", "quality_proxy", "muscle", "clock_phase", "FAST"):
        if cv == "FAST" and k == "FAST":
            continue
        z = scores[cv if cv != "FAST" else "ss_FAST"].values
        o = ols(y, cov=(z - z.mean()) / z.std())
        for c in ("PE_in_CD", "PE_in_MD", "Interaction"):
            sens.append(dict(pathway=k, adjustment=f"cov_{cv}", contrast=c, **{x: o[c][x] for x in ("estimate", "ci_low", "ci_high", "p")}))
S_nm, _ = ss_scores(cnt, {k: [g for g in gl if g not in musc] for k, gl in SETG.items()}, bg=[g for g in BG if g not in musc])
for k in PATHS:
    o = ols(S_nm[k].values)
    for c in ("PE_in_CD", "PE_in_MD", "Interaction"):
        sens.append(dict(pathway=k, adjustment="muscle_genes_excluded", contrast=c, **{x: o[c][x] for x in ("estimate", "ci_low", "ci_high", "p")}))
sens = pd.DataFrame(sens)
base = T[T.score_type == "singscore_primary"].set_index(["pathway", "contrast"])
sens["estimate_primary"] = [base.loc[(a, b), "estimate"] for a, b in zip(sens.pathway, sens.contrast)]
sens["p_primary"] = [base.loc[(a, b), "p"] for a, b in zip(sens.pathway, sens.contrast)]
sens["same_sign_as_primary"] = np.sign(sens.estimate) == np.sign(sens.estimate_primary)
sens.to_csv(RES / "covariate_and_muscle_sensitivity.csv", index=False)
tech = pd.DataFrame([dict(pathway=k, covariate=c, spearman_all19=stats.spearmanr(S_ss[k], scores[c]).statistic,
                          spearman_within_group=stats.spearmanr(resid[k], scores[c] - scores.groupby("group")[c].transform("mean")).statistic)
                     for k in ORDER for c in ("depth", "quality_proxy", "muscle", "clock_phase")])
tech.to_csv(RES / "technical_correlations.csv", index=False)

# ---------------------------------------------------------------- independent plastic studies (existing Phase 5 data)
ext, ext_scores = [], []
for st in ("HONG", "GSE245069"):
    X, sm = load_prepared(st)
    E, n_e = ss_scores(X, SETG)
    bgx = [g for g in par["background_genes"] if g in X.index]
    p5 = pd.read_csv(RS / "roadmap_phase5/per_sample_scores.csv", index_col=0)
    E["FAST"] = p5.loc[E.index, "FAST_singscore"]
    E = E[ORDER]
    E["study"], E["group"] = st, sm.loc[E.index, "group"]
    ext_scores.append(E)
    con = E.group == "CON"
    for g in sorted(set(E.group) - {"CON"}) + (["PS_pooled"] if st == "GSE245069" else []):
        ex = E.group.str.startswith("PS") if g == "PS_pooled" else E.group == g
        for k in ORDER:
            a, b = E.loc[ex, k].values, E.loc[con, k].values
            w = welch_ci(a, b)
            y = np.r_[a, b]
            obs = a.mean() - b.mean()
            ds = [y[list(c)].mean() - np.delete(y, list(c)).mean() for c in itertools.combinations(range(len(y)), len(a))]
            ext.append(dict(study=st, exposure=g, pathway=k, n_exposed=len(a), n_control=len(b), **w,
                            p_exact_perm=float(np.mean(np.abs(ds) >= abs(obs) - 1e-12)), n_labellings=len(ds),
                            n_set_genes_in_background=n_e.get(k, np.nan), background_genes=len(bgx)))
ext = pd.DataFrame(ext)
g294_est = base.reset_index()
for i, r in ext.iterrows():
    for c in ("PE_in_CD", "PE_in_MD"):
        ext.loc[i, f"gse294180_{c}"] = g294_est[(g294_est.pathway == r.pathway) & (g294_est.contrast == c)].estimate.iloc[0]
ext["same_sign_as_gse294180_CD"] = np.sign(ext.estimate) == np.sign(ext.gse294180_PE_in_CD)
ext["same_sign_as_gse294180_MD"] = np.sign(ext.estimate) == np.sign(ext.gse294180_PE_in_MD)
ext.to_csv(RES / "external_studies_pathway_contrasts.csv", index=False)
pd.concat(ext_scores).to_csv(RES / "external_per_mouse_pathway_scores.csv")

json.dump(dict(n_genes_in_background=len(BG), n_set_genes_singscore=n_ss, n_set_genes_unique=n_u,
               fast_signature=dict(n_up=len(FAST_UP), n_down=len(FAST_DN), source="roadmap_phase4d/signatures/FAST_main_n50.csv"),
               gene_sets_source="roadmap_phase3/results/gene_sets.csv"), open(RES / "run_info_pf01.json", "w"), indent=2)
dsum.to_csv(RES / "per_PE_mouse_deviation_summary.csv", index=False)

# ---------------------------------------------------------------- figure: all 19 mice, 8 scores
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COL = {"CD+Veh": "#86b6ef", "CD+PE": "#1c5cab", "MD+Veh": "#f2a38a", "MD+PE": "#c0392b"}
MK = {"CD+Veh": "o", "CD+PE": "o", "MD+Veh": "s", "MD+PE": "s"}
TITLE = {"PPARA_targets": "PPARα targets (Sanderson)", "mTORC1": "mTORC1 signalling (Hallmark)",
         "Chol_homeostasis": "Cholesterol homeostasis (Hallmark)", "Chol_synthesis": "Cholesterol synthesis (core)",
         "Glycolysis": "Glycolysis (Hallmark)", "DNA_repair": "DNA repair (Hallmark)", "Lipogenesis": "Lipogenesis (core)",
         "FAST": "Fasting signature (Phase 4D, signed)"}
plt.rcParams.update({"font.size": 7.5, "font.family": "DejaVu Sans", "axes.edgecolor": "#52514e"})
fig, axes = plt.subplots(2, 4, figsize=(11, 5.6))
rng = np.random.default_rng(1)
tp = T[T.score_type == "singscore_primary"].set_index(["pathway", "contrast"])
for ax, k in zip(axes.ravel(), ORDER):
    for j, g in enumerate(GROUPS):
        v = S_ss.loc[GRP == g, k]
        x = j + rng.uniform(-.12, .12, len(v))
        ax.scatter(x, v, s=22, color=COL[g], marker=MK[g], edgecolor="white", lw=.6, zorder=3)
        ax.hlines(v.mean(), j - .28, j + .28, color="#0b0b0b", lw=1.2, zorder=4)
        for m, xx, yy in zip(v.index, x, v):
            pass
    ax.set_xticks(range(4)); ax.set_xticklabels(GROUPS, fontsize=7)
    ax.set_title(TITLE[k], fontsize=8, loc="left")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ylo, yhi = ax.get_ylim()
    txt = []
    for c, lab in (("PE_in_CD", "CD"), ("PE_in_MD", "MD")):
        r = tp.loc[(k, c)]
        txt.append(f"{lab}: Δ={r.estimate:+.3f} [{r.ci_low:+.3f}, {r.ci_high:+.3f}]  q8={r.p_BH_8:.3f}" if r.p_BH_8 >= .001 else f"{lab}: Δ={r.estimate:+.3f} [{r.ci_low:+.3f}, {r.ci_high:+.3f}]  q8<0.001")
    ax.text(0.0, -0.30, "\n".join(txt), transform=ax.transAxes, fontsize=6.2, color="#52514e", va="top")
    ax.set_ylabel("singscore" if k != "FAST" else "signed singscore", fontsize=7)
fig.suptitle("Per-mouse pathway scores, GSE294180 (n = 4/5/5/5); bar = group mean; Δ = PE − Veh within diet (OLS, 95% CI); "
             "q8 = BH across the 8 scores within contrast", fontsize=8.5, x=0.01, ha="left")
fig.tight_layout(rect=(0, 0.02, 1, 0.95), h_pad=4.2)
for ext_ in ("png", "pdf", "svg"):
    fig.savefig(FIG / f"fig1_pathway_scores_19_mice.{ext_}", dpi=300)
print(T[(T.score_type == "singscore_primary") & T.contrast.isin(["PE_in_CD", "PE_in_MD", "Interaction"])]
      [["pathway", "contrast", "estimate", "ci_low", "ci_high", "std_effect", "p", "p_exact_perm", "p_BH_8", "orig_fry_q_family"]].round(4).to_string(index=False))
