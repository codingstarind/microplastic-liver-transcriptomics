"""Rerun determinism: snapshot the full b05 outputs, rerun b05 with B=50 (same seeds), compare CV folds, final
predictions, external metrics, labellings and the first 50 bootstrap replicates, then restore the full-run outputs."""
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
from b00_common import PB, RES, TAB

FILES = [TAB / "cv_fold_metrics.csv", TAB / "external_candidates.csv", TAB / "labelling_ranks.csv", TAB / "hyperparameter_choices.csv",
         TAB / "bootstrap_transfer_replicates.csv.gz", TAB / "bootstrap_summary.csv", TAB / "intervention_verdicts.csv",
         TAB / "mash_classifier_bootstrap.csv", RES / "candidate_predictions_all_genes.csv.gz", RES / "b05_run_info.json"]
snap = RES / "determinism_snapshot"
snap.mkdir(exist_ok=True)
for f in FILES:
    shutil.copy2(f, snap / f.name)
for m in (RES / "models").glob("*.joblib"):
    shutil.copy2(m, snap / m.name)
subprocess.run([sys.executable, str(PB / "scripts/b05_interventions.py"), "50"], check=True, cwd=PB / "scripts", capture_output=True)
rows = []


def cmp(name, a, b):
    num = a.select_dtypes("number").columns
    same_cols = list(a.columns) == list(b.columns) and len(a) == len(b)
    md = float(np.nanmax(np.abs(a[num].values - b[num].values))) if same_cols and len(num) else np.nan
    rows.append(dict(output=name, identical_shape=same_cols, max_abs_numeric_diff=md,
                     non_numeric_equal=bool(a.drop(columns=num).astype(str).equals(b.drop(columns=num).astype(str))) if same_cols else False))


for f in ("cv_fold_metrics.csv", "external_candidates.csv", "labelling_ranks.csv", "hyperparameter_choices.csv", "mash_classifier_bootstrap.csv"):
    cmp(f, pd.read_csv(snap / f), pd.read_csv(TAB / f))
cmp("candidate_predictions_all_genes.csv.gz", pd.read_csv(snap / "candidate_predictions_all_genes.csv.gz"), pd.read_csv(RES / "candidate_predictions_all_genes.csv.gz"))
a = pd.read_csv(snap / "bootstrap_transfer_replicates.csv.gz").iloc[:50]
b = pd.read_csv(TAB / "bootstrap_transfer_replicates.csv.gz")
cmp("bootstrap replicates 0-49", a.reset_index(drop=True), b.reset_index(drop=True))
D = pd.DataFrame(rows)
for f in FILES:  # restore the full run
    shutil.copy2(snap / f.name, f)
for m in snap.glob("*.joblib"):
    shutil.copy2(m, RES / "models" / m.name)
D.to_csv(TAB / "determinism_check.csv", index=False)
print(D.to_string(index=False))
