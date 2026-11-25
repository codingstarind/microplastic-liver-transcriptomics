"""B2 (demonstration): fasting-state classifier trained on the GSE107787 muscle-clean subset (ZT0/ZT4; 4 FAST v 6 FED),
R0 signature built within the training mice only; tested once on GSE225695 (3 v 3). No tuning comparison."""
import json

import numpy as np
import pandas as pd
from e01_mash_state_models import CTX, RES, build_ctx, fit_eval, metrics

build_ctx("fast")
T = CTX["fast"]["T"]
tr = ((T.study == "GSE107787") & T.clean).values
te = (T.study == "GSE225695").values
rows, preds = [], []
for model in ("R0_rule", "logistic_default"):
    _, _, p, s = fit_eval("fast", "R0", model, tr, te, "fast")
    y = T.y.values[te]
    m = metrics(y, p, s)
    rows.append(dict(task="fasting_demo", train="GSE107787 clean (n=%d, %d FAST)" % (tr.sum(), T.y.values[tr].sum()),
                     test="GSE225695 (n=%d)" % te.sum(), model=model, rep="R0", **m))
    preds.append(pd.DataFrame({"model": model, "key": T.key.values[te], "y": y, "p": p, "score": s}))
pd.DataFrame(rows).to_csv(RES / "fasting_demo_metrics.csv", index=False)
pd.concat(preds).to_csv(RES / "fasting_demo_predictions.csv", index=False)
json.dump({"n_genes": len(CTX["fast"]["genes"]), "note": "demonstration only; external n=6; sacrifice time unknown"},
          open(RES / "fasting_demo_info.json", "w"), indent=2)
print(pd.DataFrame(rows).round(3).to_string(index=False))
