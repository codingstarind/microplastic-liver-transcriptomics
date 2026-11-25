# paper_evidence_audit: manuscript evidence package (P0 + Roadmap Phases 1–5B)

**Prepared:** 2026-11-23.

**Scope:**
- Evidence collection and methodological documentation only. **No manuscript text.**
- No research file was modified, no model was retrained, and no experiment was run.
- Phase 5C is excluded.

## Contents

| File | What it is |
|---|---|
| `research_progression.md` | (A) Chronology and reasoning, P0 → Phase 5B |
| `phase_summaries/` | (B) One document per phase: P0, P1, P2, P3, P4, P4B, P4C, P4D, P4E, P5, P5B |
| `methods_inventory.csv` | (C) 56 rows: 50 implemented methods (M01–M50) + 6 considered or rejected (N01–N06), with script, inputs, settings, objective, outputs, literature |
| `results_master.csv` | (D) 701 numbers with phase, CI, n, context and **exact source path**. `source_type` distinguishes rows extracted programmatically from result tables vs hand-transcribed from reports |
| `tables/T*.csv` | Machine-readable copies of the main comparison tables, each row with a `source` column |
| `model_comparisons.md` | (E) Publication-style comparison and ablation tables, plus which comparisons are valid |
| `enrichment_and_signatures.md` | (F) Pathway, enrichment and signature evidence with robustness classes |
| `figure_catalog.csv`, `figure_recommendations.md` | (G) Inventory of every figure (visually inspected) and the proposed main/supplementary set |
| `references.bib`, `reference_verification.csv`, `citation_map.md` | (H) 79 references: 72 DOIs resolved on Crossref; 7 without DOI. Map from each claim to its citation |
| `limitations_and_corrections.md` | (I) Limitations A1–A13, corrections C1–C15, failed experiments, unresolved issues |
| `paper_outline.md` | (J) Section skeleton with phase mapping |
| `evidence_index.md` | (K) Claim → number → file → script → figure → citation |
| `scripts/` | Read-only generators: `extract_tables.py`, `make_results_master.py`, `make_methods_inventory.py`, `make_figure_catalog.py`, `make_candidate_references.py` + `verify_references.py` (Crossref; generic user agent). Also raw figure inspection notes (`figure_notes_raw.md`) |

## Regenerate

From the repository root, with the project `.venv`:

```powershell
.venv\Scripts\python.exe research\paper_evidence_audit\scripts\extract_tables.py
```
```powershell
.venv\Scripts\python.exe research\paper_evidence_audit\scripts\make_results_master.py
```
```powershell
.venv\Scripts\python.exe research\paper_evidence_audit\scripts\verify_references.py
```

## Caveats
- **Hand-transcribed rows** (P0–P3 report numbers) should be re-checked against the underlying CSVs before final citation.
- **Course methods:** the BT3041 slides have not been read. No statement here claims a method was taught in the course.
- **Pre-roadmap numbering:** the course "Phase 1 / Phase 2" work is called **P0** here, to avoid confusion with the roadmap phases.
