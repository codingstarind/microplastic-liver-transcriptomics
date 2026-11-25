"""Step 05: basic QC and comparability checks for each prepared reference.

Per dataset, on the analysis-ready matrix:
  - samples per contrast group; biological vs technical units (SRA runs per GSM)
  - library size (matrix column sum; counts data only) and sequenced reads (ENA)
  - assignment yield proxy = counted reads / sequenced reads (counts data only)
  - detected genes per sample; skeletal-muscle carry-over score
  - within-group correlation outlier flag (median - 3*MAD on log scale)
  - case-vs-control difference in technical properties (Mann-Whitney, descriptive)
  - one PCA + technical-property panel per dataset
Plus cross-study sample-overlap checks. No differential expression.
"""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from p2_common import DATA_ROOT, REPO_P2  # noqa: E402
from scipy.stats import mannwhitneyu  # noqa: E402

QC = REPO_P2 / "qc"
QC.mkdir(exist_ok=True)
MUSCLE = ["Acta1", "Ckm", "Myh1", "Myh4", "Mylpf", "Tnnt3", "Myl1", "Actn3"]
CONTRASTS = {  # dataset -> (primary contrast prefix, kind)
    "GSE316023": ("MASH_primary", "counts"),
    "GSE246328": ("MASH_secondary", "counts"),
    "GSE107787": ("FASTING_primary", "rpkm"),
    "GSE117134": ("CIRCADIAN_primary", "rpkm"),
    "GSE282998": ("PPARA_primary", "counts"),
    "GSE225695": ("FASTING_secondary", "counts"),
}


def logexpr(x: pd.DataFrame, kind: str) -> pd.DataFrame:
    if kind == "counts":
        return np.log2(x / x.sum() * 1e6 + 1)
    return np.log2(x + 1)


def qc_dataset(acc: str, runs: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    prefix, kind = CONTRASTS[acc]
    meta = pd.read_csv(REPO_P2 / "datasets" / acc / "sample_metadata.csv", dtype={"sample_id": str})
    meta = meta[meta["in_ready_matrix"]].copy()
    x = pd.read_csv(DATA_ROOT / acc / "prepared" / "expression_ready.csv.gz", index_col=0)
    x.columns = x.columns.astype(str)
    assert list(x.columns) == list(meta["sample_id"])
    r = runs[runs["dataset"] == acc].groupby("gsm").agg(sra_runs=("run_accession", "size"),
                                                         sequenced_reads=("read_count", "sum"))
    meta = meta.join(r, on="gsm")
    le = logexpr(x, kind)
    m = meta.set_index("sample_id")
    m["detected_genes"] = ((x >= (10 if kind == "counts" else 1)).sum()).values
    mus = [g for g in MUSCLE if g in x.index]
    norm = x / x.sum() * 1e6 if kind == "counts" else x
    m["muscle_score"] = norm.loc[mus].sum().values  # CPM (counts) or summed RPKM
    if kind == "counts":
        m["counted_reads"] = x.sum().values
        m["yield_counted_per_read"] = m["counted_reads"] / m["sequenced_reads"]
    # within-group correlation outliers on detected genes
    expressed = le.index[(le > 1).sum(axis=1) >= max(3, le.shape[1] // 4)]
    cor = le.loc[expressed].corr(method="pearson")
    m["mean_r_to_group"] = np.nan
    for _g, idx in m.groupby("group").groups.items():
        if len(idx) > 1:
            sub = cor.loc[idx, idx]
            m.loc[idx, "mean_r_to_group"] = (sub.sum() - 1) / (len(idx) - 1)
    med = m["mean_r_to_group"].median()
    mad = (m["mean_r_to_group"] - med).abs().median()
    m["corr_outlier_flag"] = m["mean_r_to_group"] < med - 3 * 1.4826 * mad

    # treatment-aligned technical differences in the planned contrast(s)
    tests = []
    roles = m["contrast_role"].astype(str)
    for cname in sorted({rr.split(":")[0] for rr in roles if ":" in rr and "timecourse" not in rr}):
        case = m.index[roles.str.contains(f"{cname}:case")]
        ctrl = m.index[roles.str.contains(f"{cname}:control")]
        for col in ["sequenced_reads", "counted_reads", "yield_counted_per_read", "detected_genes", "muscle_score"]:
            if col not in m or m.loc[case, col].isna().all():
                continue
            a, b = m.loc[case, col].astype(float), m.loc[ctrl, col].astype(float)
            p = mannwhitneyu(a, b).pvalue if len(a) > 1 and len(b) > 1 else np.nan
            tests.append(dict(dataset=acc, contrast=cname, n_case=len(case), n_control=len(ctrl), metric=col,
                              median_case=a.median(), median_control=b.median(),
                              ratio_case_to_control=a.median() / b.median() if b.median() else np.nan,
                              mannwhitney_p=p, complete_separation=bool(a.min() > b.max() or a.max() < b.min())))
    # PCA on contrast samples (time courses: all ready samples)
    sel = m.index[roles.str.startswith(prefix) | roles.str.contains(prefix)]
    z = le.loc[expressed, sel]
    top = z.var(axis=1).sort_values(ascending=False).index[:500]
    zc = z.loc[top].sub(z.loc[top].mean(axis=1), axis=0).T
    u, s, _ = np.linalg.svd(zc.values, full_matrices=False)
    pcs = u[:, :2] * s[:2]
    ve = s**2 / (s**2).sum()
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    groups = m.loc[sel, "group"]
    cmap = plt.get_cmap("tab20")
    for i, g in enumerate(sorted(groups.unique())):
        k = (groups == g).values
        ax[0].scatter(pcs[k, 0], pcs[k, 1], label=g, s=30, color=cmap(i % 20))
    ax[0].set_xlabel(f"PC1 ({ve[0]:.0%})")
    ax[0].set_ylabel(f"PC2 ({ve[1]:.0%})")
    ax[0].set_title(f"{acc}: PCA, top-500 variable genes, {prefix} samples", fontsize=9)
    ax[0].legend(fontsize=6, ncol=2)
    tech = "counted_reads" if kind == "counts" else "sequenced_reads"
    vals = m.loc[sel, tech] / 1e6
    order = sorted(groups.unique())
    for i, g in enumerate(order):
        v = vals[groups == g]
        ax[1].scatter(np.full(len(v), i) + np.random.default_rng(1).uniform(-.12, .12, len(v)), v, s=20,
                      color=cmap(i % 20))
    ax[1].set_xticks(range(len(order)))
    ax[1].set_xticklabels(order, rotation=45, ha="right", fontsize=7)
    ax[1].set_ylabel(("counted" if kind == "counts" else "sequenced") + " reads (M)")
    ax[1].set_title("Technical depth by group", fontsize=9)
    fig.tight_layout()
    fig.savefig(QC / f"{acc}_qc.png", dpi=110)
    plt.close(fig)
    m["pc1_contrast"] = np.nan
    m.loc[sel, "pc1_contrast"] = pcs[:, 0]
    m.loc[sel, "pc2_contrast"] = pcs[:, 1]
    m["dataset"] = acc
    keep = ["dataset", "gsm", "group", "contrast_role", "sra_runs", "sequenced_reads", "counted_reads",
            "yield_counted_per_read", "detected_genes", "muscle_score", "mean_r_to_group", "corr_outlier_flag",
            "pc1_contrast", "pc2_contrast"]
    return m.reset_index()[["sample_id"] + [k for k in keep if k in m]], tests


def overlap_checks(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    # 1. GSM accession reuse across all fetched series
    allm = []
    for f in (REPO_P2 / "metadata" / "geo_raw").glob("*_samples_raw.csv"):
        d = pd.read_csv(f, usecols=["gsm", "title"])
        d["series"] = f.name.split("_")[0]
        allm.append(d)
    allm = pd.concat(allm)
    dup = allm[allm["gsm"].duplicated(keep=False)]
    sup = {"GSE225697"}  # SuperSeries contains its SubSeries' GSMs by design
    dup_ns = dup[~dup["series"].isin(sup)]
    rows.append(dict(check="GSM shared between series (excluding SuperSeries)", result=int(dup_ns["gsm"].duplicated().sum()),
                     detail="GSE225697 SuperSeries GSMs excluded by design"))
    # 2. GSE246328 vs GSE294629 library IDs
    a = pd.read_csv(REPO_P2 / "metadata" / "geo_raw" / "GSE246328_samples_raw.csv")["title"].str.extract(r"\[(\d+)\]")[0]
    b = pd.read_csv(REPO_P2 / "metadata" / "geo_raw" / "GSE294629_samples_raw.csv")["title"].str.extract(r"^(\d+)\.fastq")[0]
    shared = set(a) & set(b)
    rows.append(dict(check="GSE246328 library IDs also in GSE294629", result=len(shared),
                     detail=f"{len(shared)}/{len(set(a))} GSE246328 libraries re-deposited in GSE294629 ({len(set(b))} libraries)"))
    # 3. identical SRA read/base counts across series (same FASTQ re-uploaded)
    rr = runs.dropna(subset=["read_count"]).copy()
    rr["sig"] = rr["read_count"].astype(int).astype(str) + "_" + rr["base_count"].astype(int).astype(str)
    g = rr.groupby("sig")["dataset"].nunique()
    multi = g[g > 1]
    det = rr[rr["sig"].isin(multi.index)].groupby("sig").apply(
        lambda t: "|".join(f"{d}:{s}" for d, s in zip(t["dataset"], t["gsm"])), include_groups=False)
    rows.append(dict(check="SRA runs with identical read+base counts in >1 series", result=len(multi),
                     detail=" ; ".join(det.tolist())))
    # 4. byte-identical control count columns, Hong-lab series
    hong = {}
    for acc, f in [("GSE261893", "GSE261893_readcount_Con_PE.txt.gz"), ("GSE261725", "GSE261725_readcount_Con_PS.txt.gz"),
                   ("GSE286116", "GSE286116_readcount_Con_PP.txt.gz"), ("GSE262399", "GSE262399_readcount_Con_PB.txt.gz")]:
        hong[acc] = pd.read_csv(DATA_ROOT / acc / f, sep="\t", index_col=0)[["CON1", "CON2", "CON3"]]
    ref = hong["GSE261893"]
    same = [k for k, v in hong.items() if k != "GSE261893" and v.reindex(ref.index).equals(ref)]
    rows.append(dict(check="Hong-lab control count columns identical to GSE261893", result=len(same), detail=",".join(same)))
    # 5. selected references: any sample used in two reference contrasts from different series?
    ms = pd.read_csv(REPO_P2 / "metadata" / "sample_metadata_master.csv")
    rows.append(dict(check="Duplicate GSM in sample_metadata_master", result=int(ms["gsm"].duplicated().sum()), detail=""))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    runs = pd.read_csv(REPO_P2 / "metadata" / "run_metadata.csv")
    allm, alltests = [], []
    for acc in CONTRASTS:
        m, t = qc_dataset(acc, runs)
        allm.append(m)
        alltests += t
        print(acc, "flags:", int(m["corr_outlier_flag"].sum()))
    pd.concat(allm).to_csv(QC / "qc_sample_metrics.csv", index=False)
    pd.DataFrame(alltests).to_csv(QC / "qc_technical_vs_treatment.csv", index=False)
    ov = overlap_checks(runs)
    ov.to_csv(QC / "overlap_checks.csv", index=False)
    print(ov.to_string())
    (QC / "qc_run_info.json").write_text(json.dumps({"muscle_genes": MUSCLE, "pca_genes": 500,
                                                       "outlier_rule": "mean within-group r < median - 3*1.4826*MAD"}))
