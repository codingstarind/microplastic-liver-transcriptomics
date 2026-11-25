"""Lock the Phase 5 plan and every frozen INPUT before any external effect/score is computed."""
import datetime as dt

from p5_common import DROOT, P5, REPO, RES, RS, json_dump, sha256, versions

inputs = [P5 / "locked_plan.md"]
inputs += sorted((RS / "roadmap_phase4d" / "signatures").glob("*_main_n50.csv"))
inputs += [RS / "roadmap_phase4d" / "frozen_signature_manifest.csv", RS / "roadmap_phase4d" / "results" / "cache" / "gene_lengths.csv",
           RS / "roadmap_phase3" / "results" / "gene_sets.csv", RS / "roadmap_phase2" / "metadata" / "phase3_gene_universe.csv",
           RS / "roadmap_phase4e" / "scripts" / "e01_mash_state_models.py", RS / "roadmap_phase4e" / "scripts" / "e03_apply_to_gse294180.py",
           RS / "roadmap_phase4e" / "results" / "gse294180_frozen_model_scores.csv",
           RS / "roadmap_phase4b" / "scripts" / "p4b_features.py", RS / "roadmap_phase4" / "scripts" / "p4_data.py",
           RS / "roadmap_phase4c_diagnostics" / "tables" / "fix_experiment_summary.csv",
           REPO / "results" / "deseq2" / "provisional_deposited" / "results_PE_in_CD.csv.gz",
           REPO / "results" / "deseq2" / "provisional_deposited" / "results_PE_in_MD.csv.gz",
           DROOT / "GSE261893" / "GSE261893_readcount_Con_PE.txt.gz", DROOT / "GSE245069" / "GSE245069_gene_expression.txt.gz"]
lock = dict(locked_at=dt.datetime.now().isoformat(timespec="seconds"), seed=20261122, versions=versions(),
            sha256={str(p.relative_to(RS) if RS in p.parents else p): sha256(p) for p in inputs})
RES.mkdir(exist_ok=True)
json_dump(lock, RES / "plan_lock.json")
print(f"locked {len(inputs)} files")
