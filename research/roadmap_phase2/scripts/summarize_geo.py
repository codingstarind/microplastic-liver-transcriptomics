"""Print a compact design summary of a fetched GEO series (audit helper)."""
import sys

import pandas as pd
from p2_common import REPO_P2

for acc in sys.argv[1:]:
    raw = REPO_P2 / "metadata" / "geo_raw"
    print("=" * 30, acc)
    for line in open(raw / f"{acc}_series.txt", encoding="utf-8"):
        if line.split(":")[0] in ("title", "summary", "overall_design", "pubmed_id", "supplementary_file", "type"):
            print(line.strip()[:2500])
    d = pd.read_csv(raw / f"{acc}_samples_raw.csv")
    ch = [c for c in d.columns if c.startswith("ch_")]
    print("n samples:", len(d), "| instrument:", d.get("instrument_model", pd.Series()).unique()[:3],
          "| strategy:", d.get("library_strategy", pd.Series()).unique()[:3])
    print("data_processing:", str(d["data_processing"].iloc[0])[:900] if "data_processing" in d else "")
    print("suppl per-sample:", d.filter(like="supplementary_file").iloc[0].tolist()[:2])
    cols = ["gsm", "title", "source_name_ch1"] + ch
    with pd.option_context("display.width", 250, "display.max_rows", 300, "display.max_colwidth", 45):
        print(d[cols].to_string())
