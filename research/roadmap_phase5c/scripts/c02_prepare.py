"""Phase 5C step 2: sample manifest + sequencing QC for PRJNA916285 (Z23). Label-free except the group column,
which comes from the deposited BioSample titles (PBSi / MPsi / EGCGi / EGCG_MPsi). One row per biological mouse.
Writes results/z23_sample_manifest.csv and tables/z23_sequencing_qc.csv."""
import json

import pandas as pd
from c00_common import P5C, RES, TAB, WORK

dm = pd.read_csv(P5C / "manifests/download_manifest.csv")
ena = pd.read_csv(P5C / "manifests/ena_filereport_PRJNA916285.tsv", sep="\t").set_index("run_accession")
rows = []
for _, r in dm.iterrows():
    lab = r.label
    q = WORK / "requant" / lab
    meta = json.loads((q / "salmon/aux_info/meta_info.json").read_text())
    fp = json.loads((q / "fastp.json").read_text())["summary"]
    grp = lab.split("_")[0]
    rows.append(dict(
        mouse_id=f"Z23_{lab}", label=lab, group=grp, run=r.run, sample_accession=r.sample_accession, sample_title=r.sample_title,
        mp_exposed=grp in ("MPs", "EGCGMPs"), egcg=grp in ("EGCG", "EGCGMPs"),
        polymer="PS" if grp in ("MPs", "EGCGMPs") else "none", size="5 um (paper abstract)", route="oral (paper abstract)",
        dose="NOT VERIFIED (paper paywalled; not in SRA/BioSample)", duration="NOT VERIFIED", strain="C57BL/6", sex="male",
        age="10 weeks (BioSample 'age')", tissue="Liver", diet="NOT STATED",
        read_pairs=int(r.pairs_dumped), read_count_match_ena=bool(r.read_count_match),
        reads_after_fastp=int(fp["after_filtering"]["total_reads"] // 2), q30_rate_r1=fp["before_filtering"]["q30_rate"],
        salmon_percent_mapped=meta["percent_mapped"], salmon_num_mapped=meta["num_mapped"],
        library_types=";".join(meta.get("library_types", [])) if isinstance(meta.get("library_types"), list) else str(meta.get("library_types"))))
M = pd.DataFrame(rows).sort_values(["egcg", "mp_exposed", "label"])
M.to_csv(RES / "z23_sample_manifest.csv", index=False)
M.drop(columns=["dose", "duration"]).to_csv(TAB / "z23_sequencing_qc.csv", index=False)
print(M[["mouse_id", "run", "read_pairs", "read_count_match_ena", "reads_after_fastp", "salmon_percent_mapped"]].to_string(index=False))
