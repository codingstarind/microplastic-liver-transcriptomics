# Representation go / no-go: per-mouse signature scores

**Decision: CONDITIONAL GO.**

| Check | Result |
|---|---|
| (a) 19-mouse matrix and metadata | **PASS.** 19 unique columns map 1:1 to GSM/group; 2×2 = 4/5/5/5 |
| (b) Expression input suitable for within-sample ranking | **CONDITIONAL.** singscore expects TPM/RPKM-like abundance after low-count filtering (Foroutan 2018, Methods). Deposited values are integer gene counts *without* a length offset (Phase 2 assumptions A1–A3). **Primary input:** TPM-proxy = count / GENCODE vM35 median transcript length, rescaled per sample, on the 13,777 filtered genes. A gene-constant length factor shifts every mouse's ranks equally, so between-mouse score differences are unaffected **except** for sample-specific isoform or length biases. Those cannot be corrected here and could be aligned with PE (library quality). **Sensitivities:** log-CPM (no length) and the existing VST. |
| (c) Signatures frozen from independent references without PE data | **PASS.** Built by `s01_build_signatures.py` from Phase 3 reference effects only; hashes in `frozen_signature_manifest.csv`, written before any mouse was scored |
| (d) IDs and universe | **PASS.** All references are mapped to vM35 symbols (the GSE294180 annotation); universe = Phase 2 main 11,888 (alt 12,010 as sensitivity); 50 up + 50 down available for every signature (candidates MASH 742 / 136, FAST 285 / 536, PPARα 584 / 584) |
| (e) MASH reference consistency | **CONDITIONAL.** Genes are required to agree in direction in BOTH MASH studies (padj < 0.05, \|log2FC\| ≥ 1). The **Phase 3 positive control formally FAILED c1 and stays failed.** The MASH-down list is dominated by cholesterol-synthesis genes (shared with fasting); a prespecified sensitivity removes CHOL_SYNTHESIS_CORE genes |
| (f) Fasting reference validity | **CONDITIONAL.** ZT-adjusted 24 h fast; genes must also agree in the muscle-clean ZT0/ZT4 subset. The fasting effect is time-dependent (Phase 3), so the score is a ZT-averaged fasting direction. GSE225695 is not used for selection |

A per-mouse score is a **different estimand** from gene-level R². It is an additional concordance assay. It does **not** correct the Phase 3 failure, and it does **not** remove the PE ≡ depth/yield alignment.
