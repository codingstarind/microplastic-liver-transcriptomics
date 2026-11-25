"""Step 6: fixed-panel pathway concordance (locked_plan.md section 6) + validation diagnostics.

Pathways: per (study contrast x set) median gene log2FC (MLE) over set genes tested, camera direction/p (from s03), and the
same median for GSE294180 PE_in_CD / PE_in_MD and the Phase 3 reference axes (FAST, MASH1, MASH2, PPARA). Primary panel =
the 12 Phase 4E R1 sets; 50 Hallmarks exploratory. Concordance across sets is descriptive (sets overlap; not independent).
Diagnostics: shared-control structure in HONG (cross-polymer gene-effect correlation vs the 0.5 expected under the null,
control-group rank on each score), depth/detection balance, muscle carry-over vs GSE294180."""
import itertools

import numpy as np
import pandas as pd
from p5_common import P5, REPO, RES, RS, SETS12, load_prepared
from scipy import stats

DE = RES / "external_de"
RE3 = RS / "roadmap_phase3/reference_effects"
gs = pd.read_csv(RS / "roadmap_phase3/results/gene_sets.csv")
SETS = {s: set(d.gene) for s, d in gs.groupby("set")}
HALL = sorted(s for s in SETS if s.startswith("HALLMARK_"))
PANEL = SETS12 + [s for s in HALL if s not in SETS12]


def med_by_set(fc: pd.Series):
    return {s: (fc.reindex(list(SETS[s])).dropna().median(), int(fc.reindex(list(SETS[s])).notna().sum())) for s in PANEL}


vectors = {}
for tgt in ("PE_in_CD", "PE_in_MD"):
    t = pd.read_csv(REPO / f"results/deseq2/provisional_deposited/results_{tgt}.csv.gz").set_index("gene")
    vectors[f"GSE294180_{tgt}"] = t.log2FC.dropna()
for lab, f in (("REF_FAST", "FAST_GSE107787_FASTvsFED_ZTadj"), ("REF_MASH1", "MASH1_GSE316023_GANvsChow_TA"),
               ("REF_MASH2", "MASH2_GSE246328_GANvsChow_24w"), ("REF_PPARA", "PPARA_GSE282998_WT_20v0")):
    vectors[lab] = pd.read_csv(RE3 / f"{f}.csv.gz").set_index("gene").log2FC.dropna()
EXT = [("HONG", "PE_vs_CON"), ("HONG", "PS_vs_CON"), ("HONG", "PP_vs_CON"), ("HONG", "PB_vs_CON"),
       ("GSE245069", "PSpooled_vs_CON"), ("GSE245069", "PS80nm_vs_CON"), ("GSE245069", "PS05um_vs_CON"), ("GSE245069", "PS5um_vs_CON"),
       ("GSE277772", "PS_vs_CON")]
cam = pd.concat([pd.read_csv(DE / f"{s}__camera.csv") for s in ("HONG", "GSE245069", "GSE277772")])
g294cam = pd.read_csv(REPO / "results/pathways/provisional_deposited/limma_fry_camera.csv")
rows = []
for st, c in EXT:
    fc = pd.read_csv(DE / f"{st}__{c}.csv.gz").set_index("gene").log2FC.dropna()
    m = med_by_set(fc)
    ref = {k: med_by_set(v) for k, v in vectors.items()} if not rows else ref  # noqa: F821
    for s in PANEL:
        cr = cam[(cam.study == st) & (cam.contrast == c) & (cam.set == s)]
        g1 = g294cam[(g294cam.pathway == s) & (g294cam.contrast == "PE_in_CD")]
        rows.append(dict(study=st, contrast=c, set=s, panel="primary12" if s in SETS12 else "hallmark_exploratory",
                         n_genes_external=m[s][1], n_genes_set=len(SETS[s]), median_log2FC_external=m[s][0],
                         camera_direction=cr.camera_direction.iloc[0] if len(cr) else None, camera_p=cr.camera_p.iloc[0] if len(cr) else np.nan,
                         camera_fdr=cr.camera_fdr.iloc[0] if len(cr) else np.nan,
                         **{f"median_log2FC_{k}": ref[k][s][0] for k in vectors}, n_genes_GSE294180=ref["GSE294180_PE_in_CD"][s][1],
                         GSE294180_PE_in_CD_camera_direction=g1.camera_direction.iloc[0] if len(g1) else None,
                         GSE294180_PE_in_CD_camera_p=g1.camera_pvalue.iloc[0] if len(g1) else np.nan))
PC = pd.DataFrame(rows)
PC.to_csv(P5 / "pathway_concordance.csv", index=False)
summ = []
for (st, c, panel), d in PC.assign(panel=PC.panel.where(PC.panel == "primary12", "hallmark_all")).groupby(["study", "contrast", "panel"]):
    if panel == "hallmark_all":
        d = PC[(PC.study == st) & (PC.contrast == c) & PC.set.str.startswith("HALLMARK_")]
    r = dict(study=st, contrast=c, panel=panel, n_sets=len(d))
    for k in vectors:
        a, b = d.median_log2FC_external.values, d[f"median_log2FC_{k}"].values
        r[f"spearman_vs_{k}"] = stats.spearmanr(a, b).statistic
        r[f"sign_agree_vs_{k}"] = float(np.mean(np.sign(a) == np.sign(b)))
    summ.append(r)
PS_ = pd.DataFrame(summ)
PS_.to_csv(RES / "pathway_concordance_summary.csv", index=False)

# ---------------- diagnostics ----------------
diag = []
# (1) HONG cross-polymer correlation of gene effects (shared control => ~0.5 expected under the null)
fcs = {c: pd.read_csv(DE / f"HONG__{c}_vs_CON.csv.gz").set_index("gene").log2FC for c in ("PE", "PS", "PP", "PB")}
X, S = load_prepared("HONG")
L = np.log2(X / X.sum() * 1e6 + 1)
L = L[(X >= 10).sum(1) >= 3]
gm = {g: L[S.index[S.group == g]].mean(1) for g in ("CON", "PE", "PS", "PP", "PB")}
for a, b in itertools.combinations(("PE", "PS", "PP", "PB"), 2):
    g = fcs[a].dropna().index.intersection(fcs[b].dropna().index)
    # same correlation computed WITHOUT the shared control: (a - mean of other polymers) vs (b - mean of other polymers) is
    # still shared; instead use the direct polymer-vs-polymer difference as the control-free comparison
    diag.append(dict(check="HONG cross-polymer gene-effect correlation (DESeq2 log2FC vs shared CON)", pair=f"{a}-CON vs {b}-CON",
                     value=stats.spearmanr(fcs[a][g], fcs[b][g]).statistic, null_expectation="~0.5 from the shared control alone"))
# (2) is the shared control group extreme on each score? rank of CON mean among the 5 group means
SC = pd.read_csv(P5 / "per_sample_scores.csv", index_col=0)
H = SC[SC.study == "HONG"]
for sc in ("FAST_singscore", "MASH_singscore", "PPARA_singscore", "R0_logistic_logodds", "R2_logistic_logodds"):
    m = H.groupby("group")[sc].mean().sort_values()
    diag.append(dict(check="HONG: position of shared CON group among 5 group means (1 = lowest)", pair=sc,
                     value=int(list(m.index).index("CON") + 1), null_expectation="each position has probability 1/5 under exchangeability"))
    # polymer-vs-polymer (control-free) FAST contrast: PE vs other polymers
    if sc == "FAST_singscore":
        pe, oth = H[H.group == "PE"][sc], H[H.group.isin(["PS", "PP", "PB"])][sc]
        diag.append(dict(check="HONG control-free contrast: PE minus other polymers (FAST)", pair=sc, value=pe.mean() - oth.mean(),
                         null_expectation="> 0 would indicate PE-specific excess over other polymers"))
# (3) depth / detection balance and muscle carry-over vs GSE294180
Q = pd.read_csv(RES / "external_sample_qc.csv")
for st, d in Q.groupby("study"):
    for col in ("counted_total", "sra_read_pairs", "detected_genes", "muscle_cpm_sum"):
        ex, co = d[d.group != "CON"][col], d[d.group == "CON"][col]
        diag.append(dict(check=f"{st}: exposed/control ratio of mean {col}", pair="all exposed vs CON",
                         value=ex.mean() / co.mean(), null_expectation="1.0 = balanced"))
g294 = pd.read_csv(REPO / "data/raw/GSE294180_MP_raw_counts_gene_level.csv.gz", index_col=0)
mus = [g for g in pd.read_csv(RS / "roadmap_phase2/metadata/muscle_exclusion_genes.csv").vm35_symbol if g in g294.index]
mcpm = (g294.loc[mus] / g294.sum() * 1e6).sum()
diag.append(dict(check="GSE294180 muscle_cpm_sum (median of 19 mice) for comparison", pair="reference", value=float(mcpm.median()),
                 null_expectation="external values in external_sample_qc.csv"))
DG = pd.DataFrame(diag)
DG.to_csv(RES / "validation_diagnostics_table.csv", index=False)
with pd.option_context("display.width", 250, "display.max_rows", 200):
    print(PS_[PS_.panel == "primary12"].round(3).to_string(index=False))
    print(PC[(PC.panel == "primary12") & (PC.contrast.isin(["PE_vs_CON", "PSpooled_vs_CON"]))][
        ["study", "contrast", "set", "n_genes_external", "median_log2FC_external", "camera_direction", "camera_p",
         "median_log2FC_GSE294180_PE_in_CD", "median_log2FC_REF_FAST", "median_log2FC_REF_MASH1"]].round(3).to_string(index=False))
    print(DG.round(3).to_string(index=False))
