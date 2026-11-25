"""Step 2: per-mouse singscore + 2x2 factorial contrasts + all locked controls (analysis_plan.md).
Fast path: inputs, gene lengths and rank matrices are cached in results/cache/; the 2x2 model is saturated, so
contrasts are exact differences of group means (identical to OLS), computed vectorised for permutations and null sets."""
import gzip
import itertools
import json
import os
import re
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

P4D = Path(__file__).resolve().parents[1]
REPO = P4D.parents[1]
RS = P4D.parent
RES, FIG, SIG = P4D / "results", P4D / "figures", P4D / "signatures"
CACHE = RES / "cache"
CACHE.mkdir(parents=True, exist_ok=True)
SEED = 20261122

# ---------------- inputs (cached) ----------------
sm = pd.read_csv(REPO / "metadata/sample_metadata.csv").set_index("counts_column")
cnt = pd.read_csv(REPO / "data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz", index_col=0)
cols = list(cnt.columns)
sm = sm.loc[cols]
MD = (sm.diet == "MD").values.astype(int)
PE = (sm.microplastic == "PE").values.astype(int)
filt = pd.read_csv(REPO / "data/processed/gene_filter.csv").query("passes_filter").gene.tolist()
lf = CACHE / "gene_lengths.csv"
if not lf.exists():
    tx, gn = {}, {}
    with gzip.open(REPO / "data/raw/gencode.vM35.annotation.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "exon":
                continue
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            tx[t] = tx.get(t, 0) + int(f[4]) - int(f[3]) + 1
            gn[t] = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
    pd.Series(tx).groupby(pd.Series(gn)).median().rename("length").to_csv(lf)
length = pd.read_csv(lf, index_col=0)["length"]
bg = [g for g in filt if g in length.index]
C = cnt.loc[bg].astype(float)
tpm = C.div(length[bg], axis=0)
tpm = tpm / tpm.sum() * 1e6
inputs = {"tpm_proxy": np.log2(tpm + 1), "logcpm": np.log2(C / C.sum() * 1e6 + 1)}
vst = pd.read_csv(REPO / "results/deseq2/provisional_deposited/vst_matrix.csv.gz", index_col=0)
inputs["vst"] = vst.loc[[g for g in bg if g in vst.index], cols]
RANKS = {}
for k, M in inputs.items():
    f = CACHE / f"ranks_{k}.npy"
    if not f.exists():
        np.save(f, M.rank(axis=0, method="average").values)
    RANKS[k] = (np.load(f), list(M.index))


def singscore(rank_mat, genes, up, down):
    """Vectorised singscore total score. up/down: lists of gene-index arrays (one per set) or a single array."""
    N = rank_mat.shape[0]

    def part(idx, reverse):
        r = rank_mat[idx]  # (..., n, 19)
        if reverse:
            r = N + 1 - r
        n = idx.shape[-1]
        smin, smax = (n + 1) / 2, (2 * N - n + 1) / 2
        return (r.mean(axis=-2) - smin) / (smax - smin) - 0.5
    return part(up, False) + part(down, True)


def contrasts(S):
    """S: (..., 19). Exact saturated-model contrasts = group-mean differences."""
    def m(d, p):
        return S[..., (MD == d) & (PE == p)].mean(axis=-1)
    return {"PE_in_CD": m(0, 1) - m(0, 0), "PE_in_MD": m(1, 1) - m(1, 0),
            "Interaction": (m(1, 1) - m(1, 0)) - (m(0, 1) - m(0, 0)), "Diet_in_Veh": m(1, 0) - m(0, 0)}


def ols_table(S, label):
    X = np.column_stack([np.ones(19), MD, PE, MD * PE])
    beta, *_ = np.linalg.lstsq(X, S, rcond=None)
    resid = S - X @ beta
    s2 = (resid**2).sum() / 15
    cov = s2 * np.linalg.inv(X.T @ X)
    L = {"PE_in_CD": [0, 0, 1, 0], "PE_in_MD": [0, 0, 1, 1], "Interaction": [0, 0, 0, 1], "Diet_in_Veh": [0, 1, 0, 0]}
    rows = []
    for k, l in L.items():
        l = np.array(l, float)
        est, se = l @ beta, np.sqrt(l @ cov @ l)
        t = est / se
        q = stats.t.ppf(.975, 15)
        rows.append(dict(score=label, contrast=k, estimate=est, se=se, ci_low=est - q * se, ci_high=est + q * se, t=t,
                         p=2 * stats.t.sf(abs(t), 15), std_effect=est / np.sqrt(s2)))
    return rows


def load_sig(tag, genes):
    s = pd.read_csv(SIG / f"{tag}.csv")
    gi = {g: i for i, g in enumerate(genes)}
    up = np.array([gi[g] for g in s[s.direction == "up"].gene if g in gi])
    dn = np.array([gi[g] for g in s[s.direction == "down"].gene if g in gi])
    return up, dn, len(s) - len(up) - len(dn), set(s.gene)


chol = set(pd.read_csv(RS / "roadmap_phase3/results/gene_sets.csv").query("set=='CHOL_SYNTHESIS_CORE'").gene)
mash_g, fast_g = set(pd.read_csv(SIG / "MASH_main_n50.csv").gene), set(pd.read_csv(SIG / "FAST_main_n50.csv").gene)
shared = mash_g & fast_g


def variant_sets(sig, genes):
    """(label, up, down, n_missing) for primary + prespecified sensitivities."""
    out = []
    base = f"{sig}_main_n50"
    up, dn, miss, gs = load_sig(base, genes)
    out.append(("primary", up, dn, miss))
    for tag, lab in ((f"{sig}_main_n25", "size_n25"), (f"{sig}_alt_n50", "alt_muscle_rule")):
        u, d, m_, _ = load_sig(tag, genes)
        out.append((lab, u, d, m_))
    gi = {g: i for i, g in enumerate(genes)}
    s = pd.read_csv(SIG / f"{base}.csv")
    if sig == "MASH":
        s2 = s[~s.gene.isin(chol)]
        out.append(("minus_CHOL_SYNTHESIS_CORE", np.array([gi[g] for g in s2[s2.direction == "up"].gene if g in gi]),
                    np.array([gi[g] for g in s2[s2.direction == "down"].gene if g in gi]), 0))
    if sig in ("MASH", "FAST"):
        s3 = s[~s.gene.isin(shared)]
        out.append(("minus_MASH_FAST_shared", np.array([gi[g] for g in s3[s3.direction == "up"].gene if g in gi]),
                    np.array([gi[g] for g in s3[s3.direction == "down"].gene if g in gi]), 0))
    return out


# ---------------- primary scores, contrasts, sensitivities ----------------
scores, table, sens = {}, [], []
for inp, (R, genes) in RANKS.items():
    for sig in ("MASH", "FAST", "PPARA"):
        for lab, up, dn, miss in variant_sets(sig, genes):
            if inp != "tpm_proxy" and lab != "primary":
                continue
            S = singscore(R, genes, up, dn)
            key = f"{sig}|{inp}|{lab}"
            scores[key] = S
            rows = ols_table(S, key)
            for r in rows:
                r.update(signature=sig, input=inp, variant=lab, n_up=len(up), n_down=len(dn), n_missing=miss)
            (table if (inp == "tpm_proxy" and lab == "primary") else sens).extend(rows)
T = pd.DataFrame(table)
fam = T.signature.isin(["MASH", "FAST"]) & T.contrast.isin(["PE_in_CD", "PE_in_MD", "Interaction"])


def holm(p):
    o = np.argsort(p)
    adj = np.empty(len(p))
    run = 0
    for i, j in enumerate(o):
        run = max(run, (len(p) - i) * p[j])
        adj[j] = min(1, run)
    return adj


T["family"] = np.where(fam, "primary(Holm,6)", np.where(T.signature.eq("PPARA") & T.contrast.ne("Diet_in_Veh"),
                                                           "secondary_PPARA(Holm,3)", np.where(T.signature.eq("MASH") & T.contrast.eq("Diet_in_Veh"),
                                                                                              "positive_control(single)", "descriptive")))
T["p_holm"] = np.nan
for f_ in ("primary(Holm,6)", "secondary_PPARA(Holm,3)"):
    i = T.family == f_
    T.loc[i, "p_holm"] = holm(T.loc[i, "p"].values)
pc = T[(T.signature == "MASH") & (T.contrast == "Diet_in_Veh")].iloc[0]
pos_pass = bool(pc.estimate > 0 and pc.ci_low > 0)

# ---------------- LOO ----------------
loo = []
for key in [k for k in scores if k.endswith("|tpm_proxy|primary")]:
    S = scores[key]
    full = {r["contrast"]: r["estimate"] for r in ols_table(S, key)}
    for i, c in enumerate(cols):
        keep = np.arange(19) != i
        Sx, MDx, PEx = S[keep], MD[keep], PE[keep]
        X = np.column_stack([np.ones(18), MDx, PEx, MDx * PEx])
        b, *_ = np.linalg.lstsq(X, Sx, rcond=None)
        e = {"PE_in_CD": b[2], "PE_in_MD": b[2] + b[3], "Interaction": b[3], "Diet_in_Veh": b[1]}
        for k, v in e.items():
            loo.append(dict(score=key.split("|")[0], dropped=c, contrast=k, estimate=v, full=full[k], sign_same=np.sign(v) == np.sign(full[k])))
LOO = pd.DataFrame(loo)

# ---------------- exact within-diet permutations ----------------
cd, mdi = np.where(MD == 0)[0], np.where(MD == 1)[0]
perm_cd = [np.isin(cd, c).astype(int) for c in itertools.combinations(cd, 5)]
perm_md = [np.isin(mdi, c).astype(int) for c in itertools.combinations(mdi, 5)]
perm = []
for key in [k for k in scores if k.endswith("|tpm_proxy|primary")]:
    S = scores[key]
    obs = contrasts(S)
    Scd, Smd = S[cd], S[mdi]
    Pcd, Pmd = np.array(perm_cd), np.array(perm_md)
    dcd = (Pcd * Scd).sum(1) / 5 - ((1 - Pcd) * Scd).sum(1) / 4
    dmd = (Pmd * Smd).sum(1) / 5 - ((1 - Pmd) * Smd).sum(1) / 5
    inter = (dmd[None, :] - dcd[:, None]).ravel()
    for k, null in (("PE_in_CD", dcd), ("PE_in_MD", dmd), ("Interaction", inter)):
        perm.append(dict(score=key.split("|")[0], contrast=k, observed=obs[k], n_perm=len(null),
                         p_exact_two_sided=float((np.abs(null) >= abs(obs[k]) - 1e-12).mean())))
PERM = pd.DataFrame(perm)

# ---------------- abundance-matched random signed gene sets ----------------
R, genes = RANKS["tpm_proxy"]
mean_ab = inputs["tpm_proxy"].mean(axis=1).values
dec = pd.qcut(mean_ab, 10, labels=False)
rng = np.random.default_rng(SEED)
nullrows = []
for sig in ("MASH", "FAST", "PPARA"):
    up, dn, _, gs = load_sig(f"{sig}_main_n50", genes)
    sig_idx = set(up) | set(dn)
    pool = {d: np.array([i for i in np.where(dec == d)[0] if i not in sig_idx]) for d in range(10)}
    draw = lambda idx: np.stack([np.array([rng.choice(pool[dec[i]]) for i in idx]) for _ in range(1000)])  # noqa: E731
    U, D = draw(up), draw(dn)
    Snull = singscore(R, genes, U, D)  # (1000, 19)
    cn = contrasts(Snull)
    obs = contrasts(scores[f"{sig}|tpm_proxy|primary"])
    for k in cn:
        nullrows.append(dict(score=sig, contrast=k, observed=obs[k], null_mean=cn[k].mean(), null_sd=cn[k].std(),
                             null_q025=np.quantile(cn[k], .025), null_q975=np.quantile(cn[k], .975),
                             percentile_of_observed=float((cn[k] < obs[k]).mean()),
                             p_two_sided_vs_null=float((np.abs(cn[k] - cn[k].mean()) >= abs(obs[k] - cn[k].mean())).mean())))
NULL = pd.DataFrame(nullrows)

# ---------------- QC associations ----------------
cov = pd.read_csv(REPO / "results/deseq2/covariates.csv").set_index("counts_column").loc[cols]
sc = pd.DataFrame({k.split("|")[0]: v for k, v in scores.items() if k.endswith("|tpm_proxy|primary")}, index=cols)
sc = sc.join(sm[["group", "diet", "microplastic"]]).join(cov)
sc.to_csv(RES / "per_mouse_scores.csv")
qc = []
for s_ in ("MASH", "FAST", "PPARA"):
    for c_ in ("depth", "quality_proxy", "muscle", "clock_phase"):
        r, p = stats.spearmanr(sc[s_], sc[c_])
        qc.append(dict(score=s_, covariate=c_, spearman_all19=r, p=p))
QC = pd.DataFrame(qc)
sep = dict(depth_PE_min=float(cov.depth[PE == 1].min()), depth_PE_max=float(cov.depth[PE == 1].max()),
           depth_Veh_min=float(cov.depth[PE == 0].min()), depth_Veh_max=float(cov.depth[PE == 0].max()),
           quality_PE_max=float(cov.quality_proxy[PE == 1].max()), quality_Veh_min=float(cov.quality_proxy[PE == 0].min()))

# ---------------- write ----------------
T.to_csv(RES / "contrast_estimates_primary.csv", index=False)
pd.DataFrame(sens).to_csv(RES / "contrast_estimates_sensitivity.csv", index=False)
LOO.to_csv(RES / "loo_contrasts.csv", index=False)
PERM.to_csv(RES / "permutation_exact.csv", index=False)
NULL.to_csv(RES / "null_gene_sets.csv", index=False)
QC.to_csv(RES / "qc_associations.csv", index=False)
loo_s = LOO.groupby(["score", "contrast"]).agg(min=("estimate", "min"), max=("estimate", "max"), n_sign_same=("sign_same", "sum")).reset_index()
loo_s.to_csv(RES / "loo_summary.csv", index=False)
json.dump(dict(positive_control_MASH_Diet_in_Veh=dict(estimate=pc.estimate, ci=[pc.ci_low, pc.ci_high], p=pc.p, PASS=pos_pass),
               background_genes=len(bg), depth_yield_separation=sep, parallelism="none needed: vectorised numpy; BLAS threads=1",
               seed=SEED), open(RES / "run_summary.json", "w"), indent=2, default=float)

# ---------------- figures ----------------
order = ["CD+Veh", "CD+PE", "MD+Veh", "MD+PE"]
fig, ax = plt.subplots(1, 3, figsize=(14, 4))
for i, s_ in enumerate(("MASH", "FAST", "PPARA")):
    for j, g in enumerate(order):
        v = sc[sc.group == g][s_]
        ax[i].scatter(np.full(len(v), j) + np.random.default_rng(j).uniform(-.1, .1, len(v)), v, s=30)
        ax[i].hlines(v.mean(), j - .25, j + .25, color="k")
    ax[i].set_xticks(range(4))
    ax[i].set_xticklabels(order)
    ax[i].set_title(f"{s_} singscore per mouse ({'secondary' if s_ == 'PPARA' else 'primary'})", fontsize=9)
fig.tight_layout()
fig.savefig(FIG / "fig1_scores_by_group.png", dpi=110)
fig, ax = plt.subplots(3, 4, figsize=(15, 10))
colr = {"CD+Veh": "tab:blue", "CD+PE": "tab:cyan", "MD+Veh": "tab:red", "MD+PE": "tab:orange"}
for i, s_ in enumerate(("MASH", "FAST", "PPARA")):
    for j, c_ in enumerate(("depth", "quality_proxy", "muscle", "clock_phase")):
        for g in order:
            d = sc[sc.group == g]
            ax[i, j].scatter(d[c_], d[s_], c=colr[g], label=g, s=25)
        ax[i, j].set_xlabel(c_)
        ax[i, j].set_ylabel(s_)
        ax[i, j].set_title(f"rho={QC.query('score==@s_ and covariate==@c_').spearman_all19.iloc[0]:.2f}", fontsize=8)
ax[0, 0].legend(fontsize=7)
fig.tight_layout()
fig.savefig(FIG / "fig2_scores_vs_technical_covariates.png", dpi=100)
fig, ax = plt.subplots(1, 3, figsize=(14, 4))
for i, s_ in enumerate(("MASH", "FAST", "PPARA")):
    d = T[(T.signature == s_)]
    nn = NULL[NULL.score == s_].set_index("contrast")
    x = np.arange(len(d))
    ax[i].errorbar(x, d.estimate, yerr=[d.estimate - d.ci_low, d.ci_high - d.estimate], fmt="o", label="observed (95% CI)")
    ax[i].fill_between(x, nn.loc[d.contrast, "null_q025"], nn.loc[d.contrast, "null_q975"], color="grey", alpha=.3,
                       label="random matched gene sets 95%")
    ax[i].axhline(0, c="k", lw=.6)
    ax[i].set_xticks(x)
    ax[i].set_xticklabels(d.contrast, rotation=20, fontsize=8)
    ax[i].set_title(s_, fontsize=9)
    ax[i].legend(fontsize=7)
fig.tight_layout()
fig.savefig(FIG / "fig3_contrasts_vs_null_sets.png", dpi=110)
S2 = pd.DataFrame(sens)
fig, ax = plt.subplots(figsize=(11, 4.5))
lab = S2[S2.contrast.isin(["PE_in_CD", "PE_in_MD"])].assign(name=lambda d: d.signature + "|" + d.input + "|" + d.variant + "|" + d.contrast)
allr = pd.concat([T[T.contrast.isin(["PE_in_CD", "PE_in_MD"])].assign(name=lambda d: d.signature + "|primary|" + d.contrast), lab])
allr = allr[allr.signature != "PPARA"]
ax.errorbar(range(len(allr)), allr.std_effect, fmt="o")
ax.axhline(0, c="k", lw=.6)
ax.set_xticks(range(len(allr)))
ax.set_xticklabels(allr.name, rotation=80, fontsize=6)
ax.set_ylabel("standardised effect (estimate / residual SD)")
ax.set_title("Method / input / signature sensitivity (PE contrasts)", fontsize=9)
fig.tight_layout()
fig.savefig(FIG / "fig4_sensitivity.png", dpi=110)

with pd.option_context("display.width", 230):
    print(T[["signature", "contrast", "estimate", "ci_low", "ci_high", "std_effect", "p", "p_holm", "family"]].round(4).to_string(index=False))
    print("positive control PASS:", pos_pass)
    print(PERM.round(4).to_string(index=False))
    print(NULL[["score", "contrast", "observed", "null_mean", "null_q025", "null_q975", "percentile_of_observed", "p_two_sided_vs_null"]].round(4).to_string(index=False))
    print(loo_s.round(4).to_string(index=False))
    print(QC.round(3).to_string(index=False))
    print(S2[S2.contrast.isin(["PE_in_CD", "PE_in_MD", "Diet_in_Veh"])][["signature", "input", "variant", "contrast", "estimate", "ci_low", "ci_high", "p"]].round(4).to_string(index=False))
    print(json.dumps(sep))
