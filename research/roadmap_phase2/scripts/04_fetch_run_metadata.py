"""Step 04: per-run sequencing metadata (read counts, layout, instrument)
from the ENA portal API, keyed by GSM (sample_alias). Small text only.
Output: metadata/run_metadata.csv (one row per SRA run; runs per GSM are
technical sequencing units, not biological replicates).
"""
import io
import urllib.request

import pandas as pd
from p2_common import REPO_P2, USER_AGENT

PROJECTS = {
    "GSE316023": "PRJNA1400485", "GSE246328": "PRJNA1032545", "GSE107787": "PRJNA421476",
    "GSE117134": "PRJNA481257", "GSE282998": "PRJNA1191425", "GSE225695": "PRJNA937184",
    "GSE245069": "PRJNA1026899", "GSE277772": "PRJNA1164147", "GSE261893": "PRJNA1089521",
    "GSE261725": "PRJNA1088491", "GSE286116": "PRJNA1207961", "GSE262399": "PRJNA1091649",
}
FIELDS = "run_accession,sample_alias,experiment_alias,read_count,base_count,library_layout,instrument_model,first_public"

if __name__ == "__main__":
    out = REPO_P2 / "metadata" / "run_metadata.csv"
    have = pd.read_csv(out) if out.exists() else pd.DataFrame(columns=["dataset"])
    frames = [have]
    for acc, prj in PROJECTS.items():
        if acc in set(have["dataset"]):
            print("present", acc)
            continue
        url = (f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={prj}&result=read_run"
               f"&fields={FIELDS}&format=tsv")
        txt = urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": USER_AGENT}),
                                     timeout=120).read().decode()
        df = pd.read_csv(io.StringIO(txt), sep="\t")
        df.insert(0, "bioproject", prj)
        df.insert(0, "dataset", acc)
        frames.append(df)
        print(acc, len(df), "runs")
    allr = pd.concat(frames, ignore_index=True)
    # ENA leaves sample_alias empty for some runs; the GSM is the experiment_alias prefix
    allr["gsm"] = allr["experiment_alias"].astype(str).str.extract(r"^(GSM\d+)")[0]
    allr.to_csv(out, index=False)
