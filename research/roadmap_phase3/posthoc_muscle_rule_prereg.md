# Post-hoc exploratory muscle-exclusion sensitivity: pre-registration

**Written 2026-11-22, before computing the list or any result with it.** The SHA-256 is recorded in `results/posthoc_audit/muscle_prereg_lock.json` before the run. **Status: EXPLORATORY / POST-HOC.**
- It does not replace the locked 309-gene exclusion or any primary result.
- It is run once, and its output is reported whatever it shows.

**Why:** the Phase 2 rule (muscle median RPKM ≥ 10 and ≥ 20 × clean-liver median + 1) also excluded genes with real hepatic expression and regulation (e.g. Pdk4, Plin4, Col1a1/2, Lpl, Cpt1b). For a gene that the liver itself expresses, carry-over contamination is a minor share of its signal.

**Rule** (tissue expression only; independent of PE and of every Phase 3 outcome; same GSE107787 source data as Phase 2: 36 gastrocnemius samples and the 21 livers with canonical muscle score < 2 RPKM):

A gene is **muscle-restricted** if all three hold:
1. median muscle RPKM ≥ 10;
2. median muscle RPKM ≥ 100 × (median clean-liver RPKM + 1);
3. median clean-liver RPKM < 1, i.e. essentially not expressed in uncontaminated liver.

The 8 canonical sarcomere markers (Acta1, Ckm, Myh1, Myh4, Mylpf, Tnnt3, Myl1, Actn3) are always excluded.

**Alternative universe:** GSE294180 filtered genes present in all four primary references (Phase 2 `in_all_primary_present`), minus the muscle-restricted list.

**Single run, scope limited to the MASH positive control:**
- refit MASH1 (GSE316023, no offset) and MASH2 (GSE246328) with DESeq2 `~ diet` on the alternative universe, using the same samples and settings;
- report M1 (Spearman) and M2 (sign agreement among genes significant in both) vs Diet_in_Veh, and MASH1 vs MASH2;
- report which previously excluded genes are re-admitted.

No thresholds will be changed after seeing the output.
