"""Hash-lock the Phase 5B protocol (03, 04) and the intervention code BEFORE any intervention is scored."""
import datetime as dt

from b00_common import P5, PB, RES, json_dump, sha256

files = [PB / "03_resampling_and_nested_cv_protocol.md", PB / "04_targeted_interventions.md",
         PB / "scripts/b04_intervention_targets.R", PB / "scripts/b05_interventions.py", PB / "scripts/b00_common.py",
         P5 / "results/plan_lock.json", P5 / "frozen/frozen_manifest.json"]
RES.mkdir(exist_ok=True)
json_dump(dict(locked_at=dt.datetime.now().isoformat(timespec="seconds"),
               note="Phase 5 external outcomes were already seen; all 5B external comparisons are post-external exploratory.",
               sha256={f.name: sha256(f) for f in files}), RES / "protocol_lock.json")
print("locked", [f.name for f in files])
