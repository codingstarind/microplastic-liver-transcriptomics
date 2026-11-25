"""Direction B, gene level: what do the diet x PE interaction genes look like?

For every gene, the PE effect under each diet (MLE log2FC with 95% Wald CI) is classified
with the same rules as the pathway categories (reports/phase2/analysis_plan.md). Genes with
interaction FDR < 0.05 are summarised by pattern, and checked against technical gene classes
(muscle-tracking genes: MD+Veh holds 3 of the 4 contaminated samples).

Usage: python analysis/interaction/interaction_patterns.py <deseq2_label>
"""

import json
import sys

import numpy as np
import pandas as pd

from mpliver.config import REPO_ROOT, TABLES_DIR

DESEQ = REPO_ROOT / "results" / "deseq2"


def excl(df, prefix):
    return (df[f"{prefix}_low"] > 0) | (df[f"{prefix}_high"] < 0)


def classify(row) -> str:
    cd, md, it = row["cd_on"], row["md_on"], row["int_on"]
    if not cd and not md:
        return "Interaction without a clear PE effect in either diet"
    if cd and md and np.sign(row["PE_in_CD"]) != np.sign(row["PE_in_MD"]):
        return "Opposing"
    if cd and md:
        if not it:
            return "Similar in both diets"
        return ("Amplified under MD" if np.sign(row["Interaction"]) == np.sign(row["PE_in_CD"])
                else "Attenuated under MD")
    return ("CD-only" if cd else "MD-only") + (", interaction supported" if it else
                                               ", interaction not supported")


def main() -> int:
    label = sys.argv[1]
    out = REPO_ROOT / "results" / "interaction" / label
    out.mkdir(parents=True, exist_ok=True)
    tabs = {}
    for k in ["PE_in_CD", "PE_in_MD", "Interaction", "Diet_in_Veh"]:
        t = pd.read_csv(DESEQ / label / f"results_{k}.csv.gz").set_index("gene")
        tabs[k] = t
    df = pd.DataFrame({"baseMean": tabs["Interaction"]["baseMean"]})
    for k, t in tabs.items():
        df[k] = t["log2FC"]
        df[f"{k}_low"] = t["ci95_low"]
        df[f"{k}_high"] = t["ci95_high"]
        df[f"{k}_padj"] = t["padj"]
        df[f"{k}_ashr"] = t["log2FC_ashr"]
    df = df.dropna(subset=["Interaction_padj"])
    df["cd_on"], df["md_on"], df["int_on"] = (excl(df, "PE_in_CD"), excl(df, "PE_in_MD"),
                                              excl(df, "Interaction"))
    sig = df[df["Interaction_padj"] < 0.05].copy()
    sig["pattern"] = sig.apply(classify, axis=1)
    muscle = set(pd.read_csv(TABLES_DIR / "muscle_tracking_genes.csv").iloc[:, 0])
    sig["muscle_tracking"] = sig.index.isin(muscle)
    sig["mitochondrial"] = sig.index.str.startswith("mt-")
    cols = ["baseMean", "Interaction", "Interaction_low", "Interaction_high", "Interaction_padj",
            "PE_in_CD", "PE_in_CD_low", "PE_in_CD_high", "PE_in_MD", "PE_in_MD_low",
            "PE_in_MD_high", "Diet_in_Veh", "pattern", "muscle_tracking", "mitochondrial"]
    sig.sort_values("Interaction_padj")[cols].round(4).to_csv(out / "interaction_fdr_genes.csv")

    # Genes with a clear PE effect in both diets: is there evidence that they differ?
    both = df[df["cd_on"] & df["md_on"] & (np.sign(df["PE_in_CD"]) == np.sign(df["PE_in_MD"]))]
    summary = {
        "label": label,
        "n_genes_tested_interaction": int(len(df)),
        "n_interaction_fdr_lt_0.05": int(len(sig)),
        "patterns": sig["pattern"].value_counts().to_dict(),
        "n_muscle_tracking_among_interaction_genes": int(sig["muscle_tracking"].sum()),
        "frac_muscle_tracking_among_tested": float(df.index.isin(muscle).mean()),
        "median_abs_interaction_log2FC_fdr_genes": float(sig["Interaction"].abs().median()),
        "median_interaction_ci_width_all_genes":
            float((df["Interaction_high"] - df["Interaction_low"]).median()),
        "n_genes_PE_effect_same_sign_both_diets": int(len(both)),
        "of_which_interaction_fdr_lt_0.05": int((both["Interaction_padj"] < 0.05).sum()),
        "of_which_interaction_ci_excludes_0": int(both["int_on"].sum()),
        "n_genes_PE_effect_opposite_sign_both_diets_ci":
            int((df["cd_on"] & df["md_on"] &
                 (np.sign(df["PE_in_CD"]) != np.sign(df["PE_in_MD"]))).sum()),
    }
    (out / "interaction_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(sig.sort_values("Interaction_padj")[["Interaction", "PE_in_CD", "PE_in_MD",
                                               "pattern", "muscle_tracking"]].head(25)
          .round(2).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
