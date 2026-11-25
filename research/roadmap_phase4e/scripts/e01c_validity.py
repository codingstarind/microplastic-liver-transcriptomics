"""Predeclared validity rule applied to saved external metrics (both primary splits:
AUC >= 0.80, bootstrap CI lower > 0.5, Brier < 0.25). Recomputed here because e01 stopped at its permutation step."""
import pandas as pd
from e01_mash_state_models import RES
M = pd.read_csv(RES / "external_metrics.csv")
prim = M[M.split.isin(["S1_train316TA_test246", "S2_train246_test316TA"])]
v = prim.groupby(["rep", "model"]).apply(lambda d: bool(len(d) == 2 and ((d.auc >= .8) & (d.auc_ci_low.fillna(0) > .5) & (d.brier < .25)).all()),
                                         include_groups=False).rename("valid").reset_index()
v.to_csv(RES / "validity.csv", index=False)
print(v.to_string(index=False))
