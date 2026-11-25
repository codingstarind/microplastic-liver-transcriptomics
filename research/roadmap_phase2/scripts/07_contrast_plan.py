"""Step 07: write metadata/reference_contrast_plan.csv.

Eligible sample IDs are pulled from the prepared sample metadata so the plan
cannot drift from the data. This is a PLAN: nothing here is estimated.
"""
import pandas as pd
from p2_common import REPO_P2

PLAN = [
    dict(contrast_id="MASH_primary", dataset="GSE316023", priority="primary",
         biological_question="Independent MASH-diet (GAN D09100310) liver response",
         comparison="GAN 23 wk vs LabDiet 5001 chow 23 wk, both at 21 C (TA)",
         design="~ diet (two groups, one factor)",
         proposed_method="DESeq2 Wald on rounded estimated counts (sensitivity: edgeR-QL; limma-voom); report unshrunk log2FC + SE and ashr-shrunk log2FC for every gene in the fixed universe",
         expected_output="per-gene log2FC, SE, p, FDR; program-level scores (same gene sets as Phase 2)",
         key_assumptions="independent mice; NB GLM; estimated counts treated as counts (length ~constant across samples, median CV 1%)",
         confounders="chow vs purified control diet (chow comparator inflates diet differences); 23 wk vs our 8 wk; 5 h fast before necropsy (our data: NOT STATED); clock time not stated"),
    dict(contrast_id="MASH_sens_TN", dataset="GSE316023", priority="sensitivity",
         biological_question="Is the MASH axis stable to housing temperature?",
         comparison="GAN 23 wk vs chow 23 wk at 29 C (TN)", design="~ diet",
         proposed_method="as MASH_primary; compare log2FC agreement with TA axis",
         expected_output="per-gene log2FC; TA-vs-TN axis correlation",
         key_assumptions="as MASH_primary", confounders="thermoneutral housing exacerbates MASH; not our housing condition"),
    dict(contrast_id="MASH_secondary", dataset="GSE246328", priority="secondary",
         biological_question="Second independent MASH axis (advanced, long-duration)",
         comparison="GAN vehicle 24 w (58 wk GAN total) vs chow vehicle 24 w (58 wk chow)",
         design="~ diet", proposed_method="DESeq2 Wald on integer counts; same outputs as MASH_primary",
         expected_output="per-gene log2FC/SE; agreement with MASH_primary (Track-B style reproducibility benchmark)",
         key_assumptions="independent mice; do NOT also use GSE294629 (same libraries)",
         confounders="58 wk vs our 8 wk; prior liver biopsy and single housing in GAN mice; Altromin grain chow; 8 w and 16 w GAN groups lack age-matched chow controls and are not eligible"),
    dict(contrast_id="FASTING_primary", dataset="GSE107787", priority="primary",
         biological_question="Fasting-vs-fed liver response at matched clock times",
         comparison="24 h fasted vs ad libitum fed liver, 3 vs 3 at each of ZT0/4/8/12/16/20",
         design="~ ZT (6-level factor) + state; heterogeneity check ~ ZT * state",
         proposed_method="limma-trend on log2(RPKM + 1) for genes passing a pre-specified expression filter (RPKM only; NB count models not applicable); muscle-marker genes excluded + muscle-score covariate sensitivity",
         expected_output="ZT-adjusted fasting log2FC + moderated SE per gene; per-ZT fasting effects (descriptive)",
         key_assumptions="state and ZT are crossed and balanced, so the fasting effect is estimable separately from time of day; RPKM log-scale approximately homoscedastic after trend",
         confounders="RPKM only (no depth weighting); muscle-gene carry-over in FED livers at ZT8-20 (treatment-aligned technical artefact); 24 h fast is far more extreme than any plausible PE feeding change; single housing"),
    dict(contrast_id="FASTING_secondary", dataset="GSE225695", priority="secondary (check only)",
         biological_question="Count-based replication of the fasting axis",
         comparison="24 h fasted vs ad libitum fed (3 vs 3); refed excluded", design="~ state",
         proposed_method="DESeq2 on StringTie-derived counts; used only to check direction agreement with FASTING_primary",
         expected_output="per-gene log2FC/SE; correlation with FASTING_primary",
         key_assumptions="independent mice", confounders="sacrifice clock time NOT STATED (fasting may be confounded with time of day); n = 3"),
    dict(contrast_id="CIRCADIAN_primary", dataset="GSE117134", priority="primary",
         biological_question="Time-of-day (diurnal) liver expression axis, independent of fasting mice",
         comparison="WT liver across ZT0/4/8/12/16/20 (3 per ZT), ad libitum",
         design="~ cos(2*pi*ZT/24) + sin(2*pi*ZT/24) (harmonic regression); alternative ~ ZT factor",
         proposed_method="limma-trend on log2(RPKM + 1); per-gene amplitude/phase; derive reference axes as predicted differences between pre-specified ZT pairs (e.g. the plausible sacrifice-time window of our data, which is UNKNOWN)",
         expected_output="per-gene rhythm coefficients, amplitude, phase, F-test p; ZT-pair contrasts",
         key_assumptions="24 h period; 4 h sampling resolves phase to ~ +/-2 h",
         confounders="female mice (ours male); RPKM only; diurnal (light-dark) not free-running, so clock and feeding rhythms are combined"),
    dict(contrast_id="CIRCADIAN_sens_male", dataset="GSE107787", priority="sensitivity",
         biological_question="Male diurnal axis check",
         comparison="FED liver across 6 ZTs (3 per ZT)", design="harmonic regression as CIRCADIAN_primary",
         proposed_method="as CIRCADIAN_primary",
         expected_output="rhythm coefficients; agreement with GSE117134",
         key_assumptions="as CIRCADIAN_primary",
         confounders="shares FED samples with FASTING_primary control arm (not independent of the fasting axis); FED muscle artefact"),
    dict(contrast_id="PPARA_primary", dataset="GSE282998", priority="primary",
         biological_question="Pharmacological PPARa activation liver axis",
         comparison="B6129SF2/J WT: GW7647 5/10/20 mg/kg/day (4 each) vs vehicle (4), 5 d oral",
         design="~ dose (factor); primary summary = 20 vs 0 and a dose-trend test; batch-matched sensitivity 10 vs 0",
         proposed_method="DESeq2 on TempO-Seq probe counts (one probe per gene); muscle-gene exclusion + muscle-score sensitivity",
         expected_output="per-gene log2FC/SE per dose; dose-trend coefficient",
         key_assumptions="independent mice; targeted assay so the gene universe is the probe panel",
         confounders="targeted FFPE assay vs bulk polyA RNA-seq; pharmacological agonist at high dose (5 d) is far stronger than any physiological PPARa activation; muscle signal in 2/4 mice of the 0 and 10 mg/kg groups; B6;129 mixed background; fasting/time not stated"),
    dict(contrast_id="PPARA_specificity_KO", dataset="GSE282998", priority="negative control",
         biological_question="Which GW7647 responses require PPARa?",
         comparison="Ppara-null: GW7647 5/10/20 vs vehicle (4 each)", design="~ dose",
         proposed_method="as PPARA_primary; genes responding in WT but not KO define the PPARa-dependent subset",
         expected_output="per-gene log2FC/SE in KO", key_assumptions="KO strain is the matched background for B6129SF2/J",
         confounders="KO and WT are different stocks; do not pool genotypes in one fit without a genotype term"),
    dict(contrast_id="PPARA_sens_outbred", dataset="GSE282998", priority="sensitivity",
         biological_question="Strain robustness of the PPARa axis",
         comparison="J:ARC(S) outbred: GW7647 vs vehicle", design="~ dose", proposed_method="as PPARA_primary",
         expected_output="per-gene log2FC/SE", key_assumptions="as PPARA_primary", confounders="outbred stock; same muscle pattern"),
]

if __name__ == "__main__":
    ms = pd.read_csv(REPO_P2 / "metadata" / "sample_metadata_master.csv", dtype={"sample_id": str})
    rows = []
    for p in PLAN:
        roles = ms["contrast_role"].astype(str).str.split(";")
        sel = ms[(ms["dataset"] == p["dataset"]) & roles.apply(lambda r, cid=p["contrast_id"]: any(x.startswith(cid + ":") for x in r))]
        rr = sel["contrast_role"].astype(str)
        case = sel[rr.str.contains(p["contrast_id"] + ":case")]
        ctrl = sel[rr.str.contains(p["contrast_id"] + ":control")]
        tc = sel[rr.str.contains(p["contrast_id"] + ":timecourse")]
        rows.append(dict(**p, n_case=len(case), n_control=len(ctrl), n_timecourse=len(tc),
                         case_gsm=";".join(case["gsm"]), control_gsm=";".join(ctrl["gsm"]),
                         timecourse_gsm=";".join(tc["gsm"]),
                         case_sample_ids=";".join(case["sample_id"]), control_sample_ids=";".join(ctrl["sample_id"]),
                         timecourse_sample_ids=";".join(tc["sample_id"])))
        assert len(sel) > 0, p["contrast_id"]
    out = pd.DataFrame(rows)
    out.to_csv(REPO_P2 / "metadata" / "reference_contrast_plan.csv", index=False)
    print(out[["contrast_id", "dataset", "n_case", "n_control", "n_timecourse"]].to_string())
