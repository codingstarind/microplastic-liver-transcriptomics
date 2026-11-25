"""Phase 5C step 4 (prespecified DESCRIPTIVE biology diagnostics; not confirmatory).
(a) Marker panels fixed before scoring (fasting/PPARa, lipogenesis, cholesterol synthesis, clock, inflammation,
    fibrosis/stellate, cell-composition/zonation, stress, technical). Per panel: mean within-gene z-scored log2CPM per
    mouse -> exposed - control difference, exact permutation p over all labellings, Holm across panels.
(b) Pathway panel (Phase 4E SETS12): Z23 median log2FC and camera, concordance (Spearman over 12 sets) with
    GSE294180 PE_in_CD, the fasting and MASH references, Hong PE and GSE245069 PS (Phase 5 values, unchanged).
(c) Cross-study gene-effect correlations of Z23 log2FC with other studies' log2FC, with the Z23 relabelling null.
(d) Technical balance: depth, mapping, muscle carry-over, sex genes."""
import itertools

import numpy as np
import pandas as pd
from c00_common import P5, RES, RS, TAB, holm, load_counts, log2cpm
from scipy import stats

PANELS = {
    "fasting_PPARa": ["Pck1", "G6pc1", "Hmgcs2", "Cpt1a", "Angptl4", "Fgf21", "Cyp4a14", "Acot1", "Pdk4", "Ehhadh"],
    "lipogenesis_SREBP1_ChREBP": ["Srebf1", "Fasn", "Acaca", "Scd1", "Elovl6", "Me1", "Gck", "Pklr", "Acly", "Insig1"],
    "cholesterol_synthesis_SREBP2": ["Hmgcr", "Hmgcs1", "Sqle", "Ldlr", "Cyp51", "Fdps", "Idi1", "Dhcr7", "Pcsk9"],
    "clock_core": ["Arntl", "Npas2", "Clock", "Nr1d1", "Nr1d2", "Dbp", "Per1", "Per2", "Per3", "Cry1", "Tef"],
    "inflammation_myeloid": ["Ptprc", "Cd68", "Adgre1", "Clec4f", "Ccl2", "Tnf", "Il1b", "Lcn2", "Saa1", "Saa2", "Cxcl1"],
    "fibrosis_stellate": ["Col1a1", "Col1a2", "Col3a1", "Timp1", "Acta2", "Lrat", "Dcn", "Tgfb1", "Lox"],
    "zonation_pericentral": ["Glul", "Cyp2e1", "Oat", "Cyp1a2", "Axin2", "Lgr5"],
    "zonation_periportal": ["Cyp2f2", "Hal", "Sds", "Pck1", "Ass1", "Arg1"],
    "cholangiocyte_endothelial": ["Krt19", "Epcam", "Sox9", "Kdr", "Pecam1", "Stab2"],
    "oxidative_xenobiotic_stress": ["Hmox1", "Nqo1", "Gsta1", "Gclc", "Mt1", "Mt2", "Atf3", "Gadd45a", "Hspa1a", "Cyp1a1"],
}
CONTR = [("MPs_vs_PBS", ["PBS"], ["MPs"]), ("EGCGMPs_vs_EGCG", ["EGCG"], ["EGCGMPs"]), ("EGCG_vs_PBS", ["PBS"], ["EGCG"])]
X, S = load_counts()
L = log2cpm(X)
mus = set(pd.read_csv(RS / "roadmap_phase2/metadata/muscle_exclusion_genes.csv").vm35_symbol)

# (d) technical
cpm = X / X.sum() * 1e6
T = pd.DataFrame(dict(counted_total=X.sum(), detected_genes=(X > 0).sum(), muscle_cpm_sum=cpm.loc[[g for g in mus if g in X.index]].sum(),
                      Xist_cpm=cpm.loc["Xist"] if "Xist" in X.index else np.nan, Ddx3y_cpm=cpm.loc["Ddx3y"] if "Ddx3y" in X.index else np.nan,
                      top100_share=X.apply(lambda c: np.sort(c.values)[::-1][:100].sum() / c.sum()))).join(S[["group", "salmon_percent_mapped", "read_pairs"]])
Cc = np.corrcoef(L.loc[(X > 10).all(1)].T.values)
np.fill_diagonal(Cc, np.nan)
T["max_corr_other_mouse"] = np.nanmax(Cc, axis=1)
T["identical_count_column_exists"] = [bool(any((X[c] == X[o]).all() for o in X.columns if o != c)) for c in X.columns]
T.to_csv(TAB / "z23_technical_qc.csv")


def perm_all(v, y):
    k = int(y.sum())
    obs = v[y].mean() - v[~y].mean()
    null = []
    for ex in itertools.combinations(range(len(v)), k):
        m = np.zeros(len(v), bool)
        m[list(ex)] = True
        null.append(v[m].mean() - v[~m].mean())
    return obs, float(np.mean(np.abs(null) >= abs(obs) - 1e-12))


rows, grows = [], []
for name, ctrl, expo in CONTR:
    m = S[S.group.isin(ctrl + expo)]
    if (m.group.isin(expo)).sum() < 3 or (m.group.isin(ctrl)).sum() < 3:
        continue
    y = m.group.isin(expo).values
    Lm = L[m.index]
    Z = Lm.sub(Lm.mean(1), axis=0).div(Lm.std(1, ddof=1).replace(0, np.nan), axis=0)
    for p, genes in PANELS.items():
        g = [x for x in genes if x in Z.index and X.loc[x, m.index].mean() >= 5]
        sc = Z.loc[g].mean().values
        obs, pp = perm_all(sc, y)
        rows.append(dict(contrast=name, panel=p, n_genes=len(g), genes_missing_or_low=";".join(sorted(set(genes) - set(g))),
                         z_diff_exposed_minus_control=obs, p_perm_two_sided=pp))
        for x in g:
            o, pg = perm_all(Lm.loc[x].values, y)
            grows.append(dict(contrast=name, panel=p, gene=x, log2cpm_diff=o, p_perm_two_sided=pg, mean_count=float(X.loc[x, m.index].mean())))
    for c in ("counted_total", "detected_genes", "muscle_cpm_sum", "top100_share", "salmon_percent_mapped"):
        o, pp = perm_all(T.loc[m.index, c].values.astype(float), y)
        rows.append(dict(contrast=name, panel=f"TECH_{c}", z_diff_exposed_minus_control=o, p_perm_two_sided=pp))
R = pd.DataFrame(rows)
R["p_holm_bio_panels"] = np.nan
for c, d in R.groupby("contrast"):
    i = d.index[~d.panel.str.startswith("TECH_")]
    R.loc[i, "p_holm_bio_panels"] = holm(R.loc[i, "p_perm_two_sided"].values)
R.to_csv(TAB / "z23_marker_panels.csv", index=False)
pd.DataFrame(grows).to_csv(TAB / "z23_marker_genes.csv", index=False)

# (b) pathway panel concordance
gs = pd.read_csv(RS / "roadmap_phase3/results/gene_sets.csv")
pc = pd.read_csv(P5 / "pathway_concordance.csv")
pc = pc[pc.panel == "primary12"] if "panel" in pc else pc
cam = pd.read_csv(RES / "de/Z23__camera.csv")
prow = []
for name in ("MPs_vs_PBS", "EGCGMPs_vs_EGCG", "MPmain_additive", "EGCG_vs_PBS"):
    f = RES / f"de/Z23__{name}.csv.gz"
    if not f.exists():
        continue
    t = pd.read_csv(f).set_index("gene")
    for s_ in pc.set.unique():
        g = [x for x in gs[gs.set == s_].gene if x in t.index]
        cm = cam[(cam.contrast == name) & (cam.set == s_)]
        prow.append(dict(contrast=name, set=s_, n_genes=len(g), median_log2FC_Z23=float(t.loc[g, "log2FC"].median()),
                         camera_direction=cm.direction.iloc[0] if len(cm) else None, camera_p=float(cm.p.iloc[0]) if len(cm) else np.nan))
PW = pd.DataFrame(prow)
ref = pc[(pc.study == "HONG") & (pc.contrast == "PE_vs_CON")].set_index("set")
refcols = {"GSE294180_PE_in_CD": "median_log2FC_GSE294180_PE_in_CD", "REF_FAST": "median_log2FC_REF_FAST", "REF_MASH1": "median_log2FC_REF_MASH1",
           "REF_PPARA": "median_log2FC_REF_PPARA"}
for k, c in refcols.items():
    PW[k] = PW.set.map(ref[c])
for st, ct in (("HONG", "PE_vs_CON"), ("GSE245069", "PSpooled_vs_CON")):
    PW[f"{st}_{ct}"] = PW.set.map(pc[(pc.study == st) & (pc.contrast == ct)].set_index("set").median_log2FC_external)
PW.to_csv(TAB / "z23_pathway_panel.csv", index=False)
conc = []
for name, d in PW.groupby("contrast"):
    for k in list(refcols) + ["HONG_PE_vs_CON", "GSE245069_PSpooled_vs_CON"]:
        conc.append(dict(contrast=name, reference=k, n_sets=len(d), spearman_over_12_sets=stats.spearmanr(d.median_log2FC_Z23, d[k]).statistic))
pd.DataFrame(conc).to_csv(TAB / "z23_pathway_concordance.csv", index=False)

# (c) cross-study gene-effect correlations with Z23 relabelling null (proxy), primary contrast only + secondary
others = {"GSE294180_PE_in_CD_MLE": None, "HONG_PE_vs_CON": P5 / "results/external_de/HONG__PE_vs_CON.csv.gz",
          "HONG_PS_vs_CON": P5 / "results/external_de/HONG__PS_vs_CON.csv.gz",
          "GSE245069_PSpooled_vs_CON": P5 / "results/external_de/GSE245069__PSpooled_vs_CON.csv.gz"}
g294 = pd.read_csv(RS / "roadmap_phase5b/results/intervention_targets/gse294180_PE_in_CD_recomputed.csv.gz").set_index("gene") \
    if (RS / "roadmap_phase5b/results/intervention_targets/gse294180_PE_in_CD_recomputed.csv.gz").exists() else None
xs = []
for name, ctrl, expo in CONTR:
    f = RES / f"de/Z23__{name}.csv.gz"
    if not f.exists():
        continue
    t = pd.read_csv(f).set_index("gene")
    m = S[S.group.isin(ctrl + expo)]
    y = m.group.isin(expo).values
    for k, fp in others.items():
        if fp is None:
            if g294 is None:
                continue
            o = g294
        else:
            o = pd.read_csv(fp).set_index("gene")
        genes = t.index.intersection(o.index)
        ov = stats.rankdata(o.loc[genes, "log2FC"].values)
        C = (X.loc[genes, m.index] / X[m.index].sum() * 1e6).values
        null = []
        for ex in itertools.combinations(range(len(y)), int(y.sum())):
            v = np.zeros(len(y), bool)
            v[list(ex)] = True
            null.append(stats.spearmanr(ov, np.log2(C[:, v].mean(1) + 1) - np.log2(C[:, ~v].mean(1) + 1)).statistic)
        null = np.array(null)
        obs = stats.spearmanr(ov, np.log2(C[:, y].mean(1) + 1) - np.log2(C[:, ~y].mean(1) + 1)).statistic
        xs.append(dict(contrast=name, other=k, n_genes=len(genes), spearman_mle=stats.spearmanr(t.loc[genes, "log2FC"], o.loc[genes, "log2FC"]).statistic,
                       spearman_proxy=obs, rank_true=int((null > obs + 1e-12).sum() + 1), n_labellings=len(null),
                       null_q95=float(np.quantile(null, .95)), null_q05=float(np.quantile(null, .05))))
pd.DataFrame(xs).to_csv(TAB / "z23_cross_study_gene_correlations.csv", index=False)
with pd.option_context("display.width", 250):
    print(R.round(3).to_string(index=False))
    print(pd.DataFrame(conc).round(3).to_string(index=False))
    print(pd.DataFrame(xs).round(3).to_string(index=False))
    print(T.round(3).to_string())
