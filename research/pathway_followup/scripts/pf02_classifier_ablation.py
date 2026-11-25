"""Pathway follow-up, step 2 (exploratory, post hoc): leave-feature(s)-out ablation of the Phase 4E R2 l2-logistic
MASH-state classifier. Reuses roadmap_phase4e/scripts/e01_mash_state_models.py unchanged (imported, read-only):
same samples, gene universe, R0/R1 features, tuning grid (C in 0.01-10, 4x3 repeated stratified CV, Brier), study-held-out
splits S1 (train GSE316023 -> test GSE246328) and S2 (reverse), and the e03 final model (both MASH studies, R0 from
genes concordant in both). GSE294180 mice are never used for fitting or tuning.
For each ablation: S1/S2 external AUC, balanced accuracy, Brier, log loss, calibration intercept/slope; then the frozen
final model's GSE294180 decision function and contrasts (Holm over PE_in_CD, PE_in_MD, Interaction as in e03).
Scale handling: contrasts are also expressed (a) per residual SD in GSE294180, (b) as a fraction of the model's own
training MASH-minus-control mean log-odds separation, (c) as a fraction of the model's GSE294180 Diet_in_Veh contrast.
Run: .venv/Scripts/python.exe research/pathway_followup/scripts/pf02_classifier_ablation.py
"""
import json
import os
import sys
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "research/roadmap_phase4e/scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy import stats  # noqa: E402

from pf02_lib import ABLATIONS, FEATS, run_ablation  # noqa: E402

OUT = ROOT / "research/pathway_followup"
RES, FIG = OUT / "results", OUT / "figures"

if __name__ == "__main__":
    jobs = list(ABLATIONS.items())
    out = Parallel(n_jobs=int(os.environ.get("PF_JOBS", 12)), backend="loky", verbose=5)(delayed(run_ablation)(k, v) for k, v in jobs)
    M = pd.concat([o["metrics"] for o in out], ignore_index=True)
    C = pd.concat([o["contrasts"] for o in out], ignore_index=True)
    S = pd.concat([o["scores"] for o in out], axis=1)
    L = pd.concat([o["loo"] for o in out], ignore_index=True)
    CV = pd.concat([o["cov"] for o in out], ignore_index=True)
    INFO = {o["name"]: o["info"] for o in out}
    M.to_csv(RES / "ablation_external_metrics.csv", index=False)
    C.to_csv(RES / "ablation_gse294180_contrasts.csv", index=False)
    S.to_csv(RES / "ablation_gse294180_scores.csv")
    L.to_csv(RES / "ablation_loo.csv", index=False)
    CV.to_csv(RES / "ablation_covariate_adjusted.csv", index=False)
    json.dump(INFO, open(RES / "ablation_model_info.json", "w"), indent=2, default=str)

    # reproduction check against stored Phase 4E outputs
    ref_m = pd.read_csv(ROOT / "research/roadmap_phase4e/results/external_metrics.csv")
    ref_m = ref_m[(ref_m.rep == "R2") & (ref_m.model == "l2_logistic")].set_index("split")
    ref_c = pd.read_csv(ROOT / "research/roadmap_phase4e/results/gse294180_frozen_model_contrasts.csv")
    ref_c = ref_c[ref_c.score == "R2_l2_logistic"].set_index("contrast")
    full_m = M[M.ablation == "full_R2"].set_index("split")
    full_c = C[C.ablation == "full_R2"].set_index("contrast")
    rep = [dict(item=f"{sp} {m}", stored=ref_m.loc[sp, m], reproduced=full_m.loc[sp, m]) for sp in ("S1_train316TA_test246", "S2_train246_test316TA")
           for m in ("auc", "bal_acc", "brier")] + \
          [dict(item=f"{c} {m}", stored=ref_c.loc[c, m], reproduced=full_c.loc[c, m]) for c in ("PE_in_CD", "PE_in_MD", "Interaction")
           for m in ("estimate", "p_holm")]
    rep = pd.DataFrame(rep)
    rep["abs_diff"] = (rep.stored - rep.reproduced).abs()
    rep.to_csv(RES / "ablation_reproduction_check.csv", index=False)
    print(rep.to_string(index=False))

    # paired mouse bootstrap: difference in scaled PE contrast, ablation minus full (frozen scores; resample mice within group)
    sm = pd.read_csv(ROOT / "metadata/sample_metadata.csv").set_index("counts_column").loc[S.index]
    grp = sm.group.values
    MD, PE = (sm.diet == "MD").astype(int).values, (sm.microplastic == "PE").astype(int).values
    sep = {k: v["train_separation"] for k, v in INFO.items()}

    def pe_contr(y, md, pe):
        X = np.column_stack([np.ones(len(y)), md, pe, md * pe])
        b = np.linalg.lstsq(X, y, rcond=None)[0]
        return b[2], b[2] + b[3]

    rng = np.random.default_rng(20261125)
    idx_g = {g: np.where(grp == g)[0] for g in np.unique(grp)}
    B = 2000
    boots = {k: np.zeros((B, 2)) for k in S.columns}
    for bi in range(B):
        ii = np.concatenate([rng.choice(v, len(v)) for v in idx_g.values()])
        for k in S.columns:
            boots[k][bi] = np.array(pe_contr(S[k].values[ii], MD[ii], PE[ii])) / sep[k]
    rows = []
    for k in S.columns:
        pt = np.array(pe_contr(S[k].values, MD, PE)) / sep[k]
        ptf = np.array(pe_contr(S["full_R2"].values, MD, PE)) / sep["full_R2"]
        for j, c in enumerate(("PE_in_CD", "PE_in_MD")):
            d = boots[k][:, j] - boots["full_R2"][:, j]
            rows.append(dict(ablation=k, contrast=c, scaled_effect=pt[j], scaled_ci_low=np.quantile(boots[k][:, j], .025),
                             scaled_ci_high=np.quantile(boots[k][:, j], .975), diff_vs_full=pt[j] - ptf[j],
                             diff_ci_low=np.quantile(d, .025), diff_ci_high=np.quantile(d, .975),
                             retained_fraction_of_full=pt[j] / ptf[j]))
    BT = pd.DataFrame(rows)
    BT.to_csv(RES / "ablation_scaled_effects_bootstrap.csv", index=False)
    print(BT.round(3).to_string(index=False))
