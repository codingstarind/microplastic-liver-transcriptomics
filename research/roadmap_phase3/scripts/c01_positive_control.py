"""Stage C/D: MASH positive control and MASH reproducibility (analysis_plan.md sections 4-5).

Target = GSE294180 Diet_in_Veh (existing DESeq2 results). References from b01_mash_references.R.
"""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from p3lib import P3, REPO, compare, read_gz_csv, ref, relabel_null, target, universe  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

OUT = P3 / "results" / "positive_control"
OUT.mkdir(parents=True, exist_ok=True)
U, US = universe("main"), universe("sensitivity")
M1 = "MASH1_GSE316023_GANvsChow_TA"
M1o = M1 + "_lengthOffset"
M1tn = "MASH1sens_GSE316023_GANvsChow_TN"
M2 = "MASH2_GSE246328_GANvsChow_24w"
CORE = ["LIPOGENESIS_CORE", "CHOL_SYNTHESIS_CORE", "SANDERSON_PPARA_TARGETS", "HALLMARK_FATTY_ACID_METABOLISM",
        "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]


def pathway_agreement(ref_label, tgt_label="GSE294180_Diet_in_Veh"):
    rp = pd.read_csv(P3 / "results" / "mash" / f"{ref_label}_pathways.csv").set_index("set")
    tp = pd.read_csv(P3 / "results" / "mash" / "GSE294180_pathways.csv")
    tp = tp[tp.contrast == tgt_label].set_index("set")
    rows = []
    for s in sorted(set(rp.index) & set(tp.index)):
        rows.append(dict(set=s, core=s in CORE, ref_dir=rp.loc[s, "fry_direction"], ref_fry_fdr=rp.loc[s, "fry_fdr"],
                         ref_camera_dir=rp.loc[s, "camera_direction"], ref_camera_fdr=rp.loc[s, "camera_fdr"],
                         tgt_dir=tp.loc[s, "fry_direction"], tgt_fry_fdr=tp.loc[s, "fry_fdr"],
                         tgt_camera_dir=tp.loc[s, "camera_direction"], tgt_camera_fdr=tp.loc[s, "camera_fdr"]))
    df = pd.DataFrame(rows)
    df["both_sig"] = (df.ref_fry_fdr < 0.05) & (df.tgt_fry_fdr < 0.05)
    df["agree"] = df.ref_dir == df.tgt_dir
    return df


def main():
    tgt = target("Diet_in_Veh")
    r1, r1o, r1tn, r2 = ref(M1), ref(M1o), ref(M1tn), ref(M2)
    res = {"criteria": "analysis_plan.md section 4 (locked, see results/plan_lock.json)"}
    res["M1_vs_target"] = compare(tgt, r1, U, "GSE294180_Diet_in_Veh", M1)
    # relabelling null, mean-difference on log2 normalised counts (same estimator for observed)
    ln = read_gz_csv(P3 / "results" / "mash" / f"{M1}_log2norm.csv.gz", index_col=0)
    ln = ln.loc[ln.index.intersection(list(U))]
    case_cols = [c for c in ln.columns if c.startswith("SD_")]
    obs_md = ln[case_cols].mean(1) - ln.drop(columns=case_cols).mean(1)
    g = obs_md.index.intersection(tgt["log2FC"].dropna().index)
    obs_rho_md = spearmanr(obs_md[g], tgt.loc[g, "log2FC"]).statistic
    null = relabel_null(ln, 8, tgt["log2FC"], n_perm=1000)
    res["M5_relabel_null"] = dict(n_perm=1000, observed_rho_meandiff=float(obs_rho_md),
                                  null_q99=float(np.quantile(null, .99)), null_max=float(null.max()),
                                  null_median=float(np.median(null)),
                                  perm_p=float((1 + (null >= obs_rho_md).sum()) / 1001))
    # LOO reference
    lo = read_gz_csv(P3 / "results" / "mash" / f"{M1}_loo_log2FC.csv.gz", index_col=0)
    gg = lo.index.intersection(tgt["log2FC"].dropna().index).intersection(list(U))
    loo_ref = {c: float(spearmanr(lo.loc[gg, c], tgt.loc[gg, "log2FC"], nan_policy="omit").statistic) for c in lo.columns}
    # LOO target (existing refits)
    loo_t = {}
    for d in sorted((REPO / "results" / "deseq2" / "loo").glob("provisional_deposited__without_*")):
        t2 = pd.read_csv(d / "results_Diet_in_Veh.csv.gz").set_index("gene")
        gq = t2.index.intersection(r1.index).intersection(list(U))
        loo_t[d.name.split("__")[-1]] = float(spearmanr(t2.loc[gq, "log2FC"], r1.loc[gq, "log2FC"], nan_policy="omit").statistic)
    res["M5_loo"] = dict(reference_min=min(loo_ref.values()), reference_max=max(loo_ref.values()),
                         target_min=min(loo_t.values()), target_max=max(loo_t.values()),
                         reference=loo_ref, target=loo_t)
    res["M5_offset_model_vs_target"] = compare(tgt, r1o, U, "GSE294180_Diet_in_Veh", M1o)
    res["M5_sensitivity_universe"] = compare(tgt, r1, US, "GSE294180_Diet_in_Veh", M1 + "@8360")
    res["TN_vs_target"] = compare(tgt, r1tn, U, "GSE294180_Diet_in_Veh", M1tn)
    pw = pathway_agreement(M1)
    pw.to_csv(OUT / "pathway_agreement_MASH1_vs_Diet_in_Veh.csv", index=False)
    core = pw[pw.core]
    cb = core[core.both_sig]
    res["M4_core_pathways"] = dict(n_core=int(len(core)), n_both_sig=int(len(cb)), n_agree=int(cb.agree.sum()),
                                   table=core[["set", "ref_dir", "ref_fry_fdr", "tgt_dir", "tgt_fry_fdr", "agree"]].to_dict("records"))
    allsig = pw[pw.both_sig]
    res["all_sets_both_sig_agree"] = f"{int(allsig.agree.sum())}/{len(allsig)}"

    # ---- decision (plan section 4) ----
    m = res["M1_vs_target"]
    c1 = m["spearman"] >= 0.30 and res["M5_relabel_null"]["observed_rho_meandiff"] > res["M5_relabel_null"]["null_q99"]
    c2 = m["sign_agree_sig_both"] >= 0.70
    c3 = res["M4_core_pathways"]["n_both_sig"] >= 2 and res["M4_core_pathways"]["n_agree"] >= 2 / 3 * res["M4_core_pathways"]["n_both_sig"]
    c4 = (res["M5_loo"]["reference_min"] >= 0.20 and res["M5_loo"]["target_min"] >= 0.20
          and res["M5_offset_model_vs_target"]["spearman"] >= 0.20)
    res["criteria_results"] = dict(c1_rho_and_null=c1, c2_sign_agreement=c2, c3_pathways=c3, c4_loo_offset=c4)
    res["decision"] = "PASS" if all([c1, c2, c3, c4]) else ("CONDITIONAL" if c1 and c2 else "FAIL")

    # ---- Stage D: reproducibility ----
    rep = [compare(r1, r2, U, M1, M2), compare(tgt, r2, U, "GSE294180_Diet_in_Veh", M2),
           compare(tgt, r1, U, "GSE294180_Diet_in_Veh", M1), compare(r1, r1tn, U, M1, M1tn),
           compare(r1, r1o, U, M1, M1o)]
    ln2 = read_gz_csv(P3 / "results" / "mash" / f"{M2}_log2norm.csv.gz", index_col=0)
    ln2 = ln2.loc[ln2.index.intersection(list(U))]
    null2 = relabel_null(ln2, 7, tgt["log2FC"], n_perm=1000)
    null12 = relabel_null(ln2, 7, r1["log2FC"], n_perm=1000)
    res["stageD"] = dict(pairs=rep, null_M2_vs_target_q99=float(np.quantile(null2, .99)),
                         null_M2_vs_M1_q99=float(np.quantile(null12, .99)))
    pw2 = pathway_agreement(M2)
    pw2.to_csv(OUT / "pathway_agreement_MASH2_vs_Diet_in_Veh.csv", index=False)
    p1 = pd.read_csv(P3 / "results" / "mash" / f"{M1}_pathways.csv").set_index("set")
    p2 = pd.read_csv(P3 / "results" / "mash" / f"{M2}_pathways.csv").set_index("set")
    j = p1.join(p2, lsuffix="_M1", rsuffix="_M2")
    j["both_sig"] = (j.fry_fdr_M1 < .05) & (j.fry_fdr_M2 < .05)
    j["agree"] = j.fry_direction_M1 == j.fry_direction_M2
    j.to_csv(OUT / "pathway_agreement_MASH1_vs_MASH2.csv")
    res["stageD"]["pathways_M1_M2_both_sig_agree"] = f"{int(j[j.both_sig].agree.sum())}/{int(j.both_sig.sum())}"
    res["stageD"]["pathways_M2_target_both_sig_agree"] = f"{int(pw2[pw2.both_sig].agree.sum())}/{int(pw2.both_sig.sum())}"
    # offset-sensitive genes
    a, b = pd.read_csv(P3 / "reference_effects" / f"{M1}.csv.gz").set_index("gene"), pd.read_csv(
        P3 / "reference_effects" / f"{M1o}.csv.gz").set_index("gene")
    sa = (a.padj < .05) * np.sign(a.log2FC)
    sb = (b.padj < .05) * np.sign(b.log2FC)
    sens = (sa != sb)
    pd.DataFrame({"no_offset_log2FC": a.log2FC, "offset_log2FC": b.log2FC, "no_offset_padj": a.padj,
                  "offset_padj": b.padj, "offset_sensitive": sens}).to_csv(OUT / "offset_sensitive_genes.csv.gz")
    res["offset_sensitive_genes"] = int(sens.sum())
    res["offset_sign_flips_among_sig_either"] = int(((a.padj < .05) & (b.padj < .05) & (np.sign(a.log2FC) != np.sign(b.log2FC))).sum())
    (OUT / "positive_control.json").write_text(json.dumps(res, indent=2, default=float))

    # ---- figure ----
    fig, ax = plt.subplots(1, 4, figsize=(18, 4.3))
    for k, (rr, lab) in enumerate([(r1, "GSE316023 GAN vs chow (TA, wk23)"), (r2, "GSE246328 GAN vs chow (58 wk)")]):
        gq = sorted(set(rr.index) & set(tgt.index) & U)
        ax[k].scatter(rr.loc[gq, "log2FC"].clip(-6, 6), tgt.loc[gq, "log2FC"].clip(-6, 6), s=2, alpha=.25)
        rho = spearmanr(rr.loc[gq, "log2FC"], tgt.loc[gq, "log2FC"], nan_policy="omit").statistic
        ax[k].axhline(0, c="k", lw=.5)
        ax[k].axvline(0, c="k", lw=.5)
        ax[k].set_xlabel(lab + " log2FC")
        ax[k].set_ylabel("GSE294180 Diet_in_Veh log2FC")
        ax[k].set_title(f"Spearman {rho:.2f}", fontsize=9)
    gq = sorted(set(r1.index) & set(r2.index) & U)
    ax[2].scatter(r1.loc[gq, "log2FC"].clip(-6, 6), r2.loc[gq, "log2FC"].clip(-6, 6), s=2, alpha=.25, c="tab:green")
    ax[2].set_xlabel("GSE316023 log2FC")
    ax[2].set_ylabel("GSE246328 log2FC")
    ax[2].set_title(f"MASH1 vs MASH2 Spearman {rep[0]['spearman']:.2f}", fontsize=9)
    ax[3].hist(null, bins=40, color="grey")
    ax[3].axvline(res["M5_relabel_null"]["observed_rho_meandiff"], c="red")
    ax[3].set_title("Relabelling null (GSE316023 mice), red = observed", fontsize=9)
    ax[3].set_xlabel("Spearman with Diet_in_Veh")
    fig.tight_layout()
    fig.savefig(P3 / "figures" / "mash_positive_control.png", dpi=110)
    print(json.dumps({k: v for k, v in res.items() if k not in ("M5_loo",)}, indent=1, default=float)[:6000])
    print("LOO ref", res["M5_loo"]["reference_min"], res["M5_loo"]["reference_max"],
          "LOO tgt", res["M5_loo"]["target_min"], res["M5_loo"]["target_max"])


if __name__ == "__main__":
    main()
