"""Enrichment audit, step 2: consolidate existing + new enrichment outputs into complete numerical tables.

Reads only (never writes) results/pathways/** and results/deseq2/**; writes to research/enrichment_audit/tables/.
Scope: 50 MSigDB Hallmark sets + 3 PPAR sets; contrasts PE_in_CD, PE_in_MD, Interaction (Diet_in_Veh kept as
the design's positive-control contrast). Run: .venv/Scripts/python.exe research/enrichment_audit/scripts/ea02_tables.py
"""
from pathlib import Path
import numpy as np
import pandas as pd


def _bh_core(p):
    """Benjamini-Hochberg adjusted p-values (identical to R p.adjust(method = 'BH'))."""
    n = len(p)
    order = np.argsort(p)[::-1]
    ranked = p[order] * n / np.arange(n, 0, -1)
    adj = np.minimum(1, np.minimum.accumulate(ranked))
    out = np.empty(n)
    out[order] = adj
    return out

REPO = Path(__file__).resolve().parents[3]
PW = REPO / "results" / "pathways"
EA = REPO / "research" / "enrichment_audit"
TAB = EA / "tables"
TAB.mkdir(parents=True, exist_ok=True)
FAMS = ["Hallmark", "PPAR"]
PE_CONTRASTS = ["PE_in_CD", "PE_in_MD", "Interaction"]
ALL_CONTRASTS = PE_CONTRASTS + ["Diet_in_Veh"]
ALPHA = 0.05


def bh(p):
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        out[ok] = _bh_core(p[ok])
    return out


def sign_of(x):
    return np.sign(x).astype(int)


def dir_to_sign(d):
    return d.map({"Up": 1, "Down": -1}).astype(float)


# ---------------------------------------------------------------- primary results (existing, read-only)
prim = PW / "provisional_deposited"
fry = pd.read_csv(prim / "limma_fry_camera.csv")
fg = pd.read_csv(prim / "fgsea_results.csv.gz")
fg = fg[fg.family.isin(FAMS)].copy()
mz = pd.read_csv(prim / "pathway_score_effects.csv.gz")
mz = mz[mz.family.isin(FAMS)].copy()
for name, d in [("fry", fry), ("fgsea", fg), ("meanz", mz)]:
    n = d.groupby("contrast").pathway.nunique()
    assert (n == 53).all(), (name, n)

# ---------------------------------------------------------------- T1 long table with three FDR families
rows = []
for _, r in fry.iterrows():
    base = dict(family=r.family, pathway=r.pathway, contrast=r.contrast, n_genes=r.n_genes,
                mean_gene_log2FC_ashr=r.mean_gene_log2FC_ashr, frac_genes_up=r.frac_genes_up)
    rows.append({**base, "method": "fry_directional", "statistic": np.nan, "direction": r.fry_direction,
                 "pvalue": r.fry_pvalue, "fdr_native": r.fry_fdr})
    rows.append({**base, "method": "fry_mixed", "statistic": np.nan, "direction": "Mixed (non-directional)",
                 "pvalue": r.fry_pvalue_mixed, "fdr_native": r.fry_fdr_mixed})
    rows.append({**base, "method": "camera", "statistic": np.nan, "direction": r.camera_direction,
                 "pvalue": r.camera_pvalue, "fdr_native": r.camera_fdr})
for _, r in fg.iterrows():
    rows.append(dict(family=r.family, pathway=r.pathway, contrast=r.contrast, n_genes=r["size"],
                     method="fgsea_preranked_wald", statistic=r.NES, direction="Up" if r.NES > 0 else "Down",
                     pvalue=r.pval, fdr_native=r.padj, leading_edge_first15=r.leadingEdge))
for _, r in mz.iterrows():
    rows.append(dict(family=r.family, pathway=r.pathway, contrast=r.contrast, n_genes=r.n_genes,
                     method="meanz_score_lm", statistic=r.estimate, se=r.se, ci95_low=r.ci95_low,
                     ci95_high=r.ci95_high, direction="Up" if r.estimate > 0 else "Down",
                     pvalue=r.pvalue, fdr_native=r.padj))
long = pd.DataFrame(rows)
# Recompute the native family FDR as a check (BH within method x gene-set family x contrast)
long["fdr_native_recomputed"] = long.groupby(["method", "family", "contrast"]).pvalue.transform(bh)
long["fdr_all53_within_contrast"] = long.groupby(["method", "contrast"]).pvalue.transform(bh)
sub = long.contrast.isin(PE_CONTRASTS)
long["fdr_all53_x3contrasts"] = np.nan
long.loc[sub, "fdr_all53_x3contrasts"] = long[sub].groupby("method").pvalue.transform(bh)
long["fdr_native_check_maxdiff"] = (long.fdr_native - long.fdr_native_recomputed).abs()
cols = ["family", "pathway", "contrast", "method", "n_genes", "statistic", "se", "ci95_low", "ci95_high",
        "direction", "pvalue", "fdr_native", "fdr_all53_within_contrast", "fdr_all53_x3contrasts",
        "mean_gene_log2FC_ashr", "frac_genes_up", "leading_edge_first15"]
long = long.sort_values(["contrast", "method", "family", "pathway"])
long[cols].to_csv(TAB / "T1_enrichment_long_all_methods.csv", index=False)
fdr_check = long.groupby(["method"]).fdr_native_check_maxdiff.max()

# ---------------------------------------------------------------- T2 wide: CD | MD | Interaction (+ Diet_in_Veh)
fry_i = fry.set_index(["family", "pathway", "contrast"])
fg_i = fg.set_index(["family", "pathway", "contrast"])
mz_i = mz.set_index(["family", "pathway", "contrast"])
fdr53 = long.set_index(["method", "family", "pathway", "contrast"]).fdr_all53_within_contrast
wide = []
for (fam, pw), _ in fry.groupby(["family", "pathway"]):
    row = {"family": fam, "pathway": pw, "n_genes_tested": fry_i.loc[(fam, pw, "PE_in_CD"), "n_genes"]}
    for k in ALL_CONTRASTS:
        f, g, m = fry_i.loc[(fam, pw, k)], fg_i.loc[(fam, pw, k)], mz_i.loc[(fam, pw, k)]
        row.update({
            f"{k}__mean_log2FC_ashr": f.mean_gene_log2FC_ashr, f"{k}__frac_genes_up": f.frac_genes_up,
            f"{k}__fry_direction": f.fry_direction, f"{k}__fry_p": f.fry_pvalue, f"{k}__fry_fdr": f.fry_fdr,
            f"{k}__fry_fdr_all53": fdr53.loc[("fry_directional", fam, pw, k)],
            f"{k}__fry_mixed_p": f.fry_pvalue_mixed, f"{k}__fry_mixed_fdr": f.fry_fdr_mixed,
            f"{k}__camera_direction": f.camera_direction, f"{k}__camera_p": f.camera_pvalue,
            f"{k}__camera_fdr": f.camera_fdr,
            f"{k}__fgsea_NES": g.NES, f"{k}__fgsea_p": g.pval, f"{k}__fgsea_padj": g.padj,
            f"{k}__meanz_est": m.estimate, f"{k}__meanz_ci95_low": m.ci95_low,
            f"{k}__meanz_ci95_high": m.ci95_high, f"{k}__meanz_p": m.pvalue, f"{k}__meanz_fdr": m.padj,
        })
    wide.append(row)
wide = pd.DataFrame(wide).sort_values(["family", "pathway"])
wide.to_csv(TAB / "T2_enrichment_wide_CD_MD_interaction.csv", index=False)

# ---------------------------------------------------------------- T3 direction consistency across methods
dc = []
for (fam, pw, k), f in fry_i.iterrows():
    g, m = fg_i.loc[(fam, pw, k)], mz_i.loc[(fam, pw, k)]
    s = {"fry": dir_to_sign(pd.Series([f.fry_direction]))[0],
         "camera": dir_to_sign(pd.Series([f.camera_direction]))[0],
         "fgsea_NES": np.sign(g.NES), "mean_log2FC_ashr": np.sign(f.mean_gene_log2FC_ashr),
         "frac_up_minus_half": np.sign(f.frac_genes_up - 0.5), "meanz": np.sign(m.estimate)}
    vals = np.array(list(s.values()))
    dc.append({"family": fam, "pathway": pw, "contrast": k, **{f"sign_{a}": b for a, b in s.items()},
               "n_methods_agree_with_fry": int((vals == s["fry"]).sum() - 1), "n_other_methods": len(vals) - 1,
               "all_signs_agree": bool((vals == vals[0]).all()), "fry_fdr": f.fry_fdr,
               "fry_significant_and_sign_disagreement": bool(f.fry_fdr < ALPHA and not (vals == vals[0]).all())})
dc = pd.DataFrame(dc).sort_values(["contrast", "family", "pathway"])
dc.to_csv(TAB / "T3_direction_consistency_across_methods.csv", index=False)

# ---------------------------------------------------------------- T4 / T5 sensitivity
variants = {
    "primary": prim,
    "E_cov_depth": PW / "provisional_deposited__cov_depth",
    "E_cov_muscle": PW / "provisional_deposited__cov_muscle",
    "E_cov_clock_phase": PW / "provisional_deposited__cov_clock_phase",
    "E_cov_quality_proxy": PW / "provisional_deposited__cov_quality_proxy",
    "E_drop_MD_1": PW / "provisional_deposited__drop_MD_1",
    "E_pilot_length_offset": PW / "provisional_deposited__pilot_length_offset",
}
sens = []
for v, d in variants.items():
    x = pd.read_csv(d / "limma_fry_camera.csv"); x["variant"] = v; x["source"] = "existing:" + str(d.relative_to(REPO))
    sens.append(x)
new = pd.read_csv(EA / "results" / "sensitivity" / "new_sensitivity_fry_camera.csv")
new["source"] = "new:research/enrichment_audit/results/sensitivity/new_sensitivity_fry_camera.csv"
sens.append(new)
loo_dirs = sorted((PW / "loo").glob("provisional_deposited__without_*"))
for d in loo_dirs:
    x = pd.read_csv(d / "limma_fry_camera.csv"); x["variant"] = "LOO_" + d.name.split("without_")[1]
    x["source"] = "existing:" + str(d.relative_to(REPO)); sens.append(x)
sens = pd.concat(sens, ignore_index=True)
sens = sens[["variant", "family", "pathway", "contrast", "n_genes", "fry_direction", "fry_pvalue", "fry_fdr",
             "fry_pvalue_mixed", "fry_fdr_mixed", "camera_direction", "camera_pvalue", "camera_fdr",
             "mean_gene_log2FC_ashr", "source"]]
sens.to_csv(TAB / "T4_sensitivity_fry_camera_long.csv", index=False)


def call(d, q):
    return np.where(q < ALPHA, np.where(d == "Up", "Up", "Down"), "ns")


sens["fry_call"] = call(sens.fry_direction, sens.fry_fdr)
sens["camera_call"] = call(sens.camera_direction, sens.camera_fdr)
key = ["family", "pathway", "contrast"]
pv = sens.pivot_table(index=key, columns="variant", values="fry_call", aggfunc="first")
pv_p = sens.pivot_table(index=key, columns="variant", values="fry_pvalue", aggfunc="first")
pv_dir = sens.pivot_table(index=key, columns="variant", values="fry_direction", aggfunc="first")
non_loo = [c for c in pv.columns if not c.startswith("LOO_") and c != "primary"]
loo = [c for c in pv.columns if c.startswith("LOO_")]
summ = pd.DataFrame(index=pv.index)
summ["primary_fry_call"] = pv["primary"]
summ["primary_fry_p"] = pv_p["primary"]
for c in non_loo:
    summ[f"{c}__call"] = pv[c]
    summ[f"{c}__p"] = pv_p[c]
summ["n_variants"] = len(non_loo)
summ["n_variants_same_call"] = (pv[non_loo].values == pv[["primary"]].values).sum(1)
summ["n_variants_same_direction"] = (pv_dir[non_loo].values == pv_dir[["primary"]].values).sum(1)
summ["n_LOO"] = len(loo)
summ["n_LOO_same_call"] = (pv[loo].values == pv[["primary"]].values).sum(1)
summ["n_LOO_same_direction"] = (pv_dir[loo].values == pv_dir[["primary"]].values).sum(1)
summ["LOO_max_fry_p"] = pv_p[loo].max(1)
summ = summ.reset_index().sort_values(["contrast", "family", "pathway"])
summ.to_csv(TAB / "T5_sensitivity_summary_per_set.csv", index=False)

vs = (sens.assign(sig=sens.fry_fdr < ALPHA)
      .groupby(["variant", "contrast"]).agg(n_sets=("pathway", "size"), n_fry_fdr05=("sig", "sum"),
                                            n_up=("fry_call", lambda s: (s == "Up").sum()),
                                            n_down=("fry_call", lambda s: (s == "Down").sum())).reset_index())
agree = (sens.merge(sens[sens.variant == "primary"][key + ["fry_call"]], on=key, suffixes=("", "_primary"))
         .assign(same=lambda d: d.fry_call == d.fry_call_primary)
         .groupby(["variant", "contrast"]).same.sum().rename("n_calls_same_as_primary").reset_index())
vs = vs.merge(agree, on=["variant", "contrast"])
vs.to_csv(TAB / "T6_sensitivity_counts_by_variant.csv", index=False)

# fgsea sensitivity (new variants only have fgsea for muscle-excluded and depth-thinned; existing for the rest)
fgs = []
for v, d in variants.items():
    x = pd.read_csv(d / "fgsea_results.csv.gz"); x = x[x.family.isin(FAMS)]; x["variant"] = v; fgs.append(x)
xn = pd.read_csv(EA / "results" / "sensitivity" / "new_sensitivity_fgsea.csv"); fgs.append(xn)
fgs = pd.concat(fgs, ignore_index=True)[["variant", "family", "pathway", "contrast", "size", "NES", "pval", "padj"]]
fgs.to_csv(TAB / "T7_sensitivity_fgsea_long.csv", index=False)

# ---------------------------------------------------------------- T8 FDR family specification
spec = pd.DataFrame([
    ("fry_directional / fry_mixed", "limma::fry rotation test on VST matrix; design ~ diet*microplastic (+ covariate in variants); contrast vector per estimand; mouse = unit; self-contained",
     "BH within gene-set family (Hallmark: 50; PPAR: 3) x contrast", "fdr_native"),
    ("camera", "limma::camera (competitive; inter-gene correlation estimated) on same VST matrix and design",
     "BH within gene-set family x contrast", "fdr_native"),
    ("fgsea_preranked_wald", "fgseaMultilevel on DESeq2 Wald statistic of the contrast (gene = unit; ignores inter-gene correlation); minSize 15, maxSize 500, eps 0, seed 20261102",
     "BH within gene-set family x contrast (fgsea padj)", "fdr_native"),
    ("meanz_score_lm", "per-mouse mean z of VST over set genes (>=10 genes); lm(score ~ diet*microplastic); t-test of contrast",
     "BH within gene-set family x contrast", "fdr_native"),
    ("all methods", "alternative family", "BH across all 53 sets (Hallmark + PPAR) within method x contrast", "fdr_all53_within_contrast"),
    ("all methods", "alternative family", "BH across 53 sets x 3 contrasts (PE_in_CD, PE_in_MD, Interaction) within method; Diet_in_Veh excluded", "fdr_all53_x3contrasts"),
], columns=["method", "test", "BH family", "column in T1"])
spec.to_csv(TAB / "T8_fdr_family_specification.csv", index=False)

# ---------------------------------------------------------------- QC log
with open(EA / "logs" / "ea02_checks.txt", "w") as fh:
    fh.write(f"primary rows fry={len(fry)} fgsea={len(fg)} meanz={len(mz)}\n")
    fh.write(f"max |stored FDR - recomputed BH within family x contrast| by method:\n{fdr_check.to_string()}\n")
    fh.write(f"LOO refits found: {len(loo_dirs)}\n")
    fh.write(f"fry-significant (FDR<0.05) set x contrast with any cross-method sign disagreement: "
             f"{int(dc.fry_significant_and_sign_disagreement.sum())}\n")
    fh.write(vs.to_string() + "\n")
print(open(EA / "logs" / "ea02_checks.txt").read())
