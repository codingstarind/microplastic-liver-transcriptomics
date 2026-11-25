"""Interaction follow-up, step 3: the single pair that met the reproducibility rule in ix02 (HGB: FAST_eff x CIRC_se).
 (a) Held-out test: HGB allowed ONLY this interaction (interaction_cst = [[FAST_eff, CIRC_se], other features singly]),
     same stored per-split params, vs the additive HGB (no interactions). Family: {ALL, LIPO, PPARA} x {CD, MD} = 6 tests,
     module-cluster t, Holm. CHOL, OTHER reported with BH.
 (b) Pure-interaction surface of the full Phase 4 HGB refits on a fixed 15 x 15 quantile grid (2.5-97.5%) of the two
     features: PD2 - PD1(FAST_eff) - PD1(CIRC_se) (centred), background = 300 random training genes per split.
     Reported per split; figure uses the median over 15 splits.
Run: .venv/Scripts/python.exe research/interaction_followup/scripts/ix03_pair_followup.py
"""
import json
import os
import sys
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"
RS = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RS / "roadmap_phase4/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from p4_data import dl_tau2, feature_names, feature_table, outer_splits, target  # noqa: E402
from p4_models import make  # noqa: E402

OUT = RS / "interaction_followup"
RES = OUT / "results"
FEATS = feature_names("primary")
A, Bf = "FAST_eff", "CIRC_se"
ia, ib = FEATS.index(A), FEATS.index(Bf)
METR = pd.read_csv(RS / "roadmap_phase4/results/taskA/metrics_primary_positive.csv")
NG, NBG, SEED = 15, 300, 20261125
X0 = feature_table("primary")
GA = np.quantile(X0[A], np.linspace(.025, .975, NG))
GB = np.quantile(X0[Bf], np.linspace(.025, .975, NG))


def task(t, r_i, k_i):
    X = feature_table("primary")
    j = X.join(target(t), how="inner").dropna(subset=["y", "s"])
    j = j[j["s"] > 0]
    for r, k, tr, te, grp in outer_splits(j.index, "module"):
        if (r, k) == (r_i, k_i):
            break
    Xv, y, s = j[FEATS].values, j["y"].values, j["s"].values
    w = 1.0 / (s**2 + dl_tau2(y[tr], s[tr]))
    params = json.loads(METR[(METR.run == f"primary__{t}") & (METR.model == "hgb") & (METR.repeat == r_i) & (METR.fold == k_i)].iloc[0].best_params)

    def fit(cst):
        est, _ = make("hgb", FEATS)
        est.set_params(**params)
        if cst is not None:
            est.set_params(model__interaction_cst=cst)
        return est.fit(Xv[tr], y[tr], model__sample_weight=w[tr])
    full = fit(None)
    pair = fit([[ia, ib]] + [[c] for c in range(len(FEATS)) if c not in (ia, ib)])
    add = fit("no_interactions")
    P = pd.DataFrame({"gene": j.index[te], "y": y[te], "target": t, "repeat": r_i, "fold": k_i, "module": grp[te],
                      "pred_full": full.predict(Xv[te]), "pred_pair": pair.predict(Xv[te]), "pred_additive": add.predict(Xv[te]),
                      "null": np.average(y[tr], weights=w[tr])})
    rng = np.random.default_rng(SEED + 10 * r_i + k_i + (0 if t == "PE_in_CD" else 100))
    bg = Xv[rng.choice(np.where(tr)[0], NBG, replace=False)]

    def pdf(cols_vals):
        Z = bg.copy()
        for c, v in cols_vals:
            Z[:, c] = v
        return full.predict(Z).mean()
    pa = np.array([pdf([(ia, a)]) for a in GA])
    pb = np.array([pdf([(ib, b)]) for b in GB])
    p2 = np.array([[pdf([(ia, a), (ib, b)]) for b in GB] for a in GA])
    inter = p2 - pa[:, None] - pb[None, :]
    inter -= inter.mean()
    surf = pd.DataFrame([dict(target=t, repeat=r_i, fold=k_i, i=i, jx=jj, FAST_eff=GA[i], CIRC_se=GB[jj], pure_interaction=inter[i, jj],
                              joint_pd=p2[i, jj] - p2.mean()) for i in range(NG) for jj in range(NG)])
    return P, surf


if __name__ == "__main__":
    jobs = [(t, r, k) for t in ("PE_in_CD", "PE_in_MD") for r in range(3) for k in range(5)]
    out = Parallel(n_jobs=14, backend="loky", verbose=3)(delayed(task)(*a) for a in jobs)
    P = pd.concat([o[0] for o in out], ignore_index=True)
    SU = pd.concat([o[1] for o in out], ignore_index=True)
    P.to_csv(RES / "pair_only_heldout_predictions.csv.gz", index=False)
    SU.to_csv(RES / "pair_interaction_surface_by_split.csv.gz", index=False)
    # check full refit equals saved Phase 4 HGB
    mx = 0.0
    for t in ("PE_in_CD", "PE_in_MD"):
        sv = pd.read_csv(RS / f"roadmap_phase4/results/taskA/predictions/primary__{t}.csv.gz")
        sv = sv[sv.model == "hgb"].set_index(["gene", "repeat", "fold"]).pred
        b = P[P.target == t].set_index(["gene", "repeat", "fold"]).pred_full
        mx = max(mx, float((sv.loc[b.index] - b).abs().max()))
    json.dump({"pair": [A, Bf], "max_abs_diff_full_vs_saved_phase4_hgb": mx, "grid_FAST_eff": GA.tolist(), "grid_CIRC_se": GB.tolist()},
              open(RES / "run_info_ix03.json", "w"), indent=2)
    print("max diff vs saved", mx)
