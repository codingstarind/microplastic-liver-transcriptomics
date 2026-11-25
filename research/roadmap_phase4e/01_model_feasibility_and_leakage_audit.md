# Feasibility and leakage audit (before any fit)

| Item | Finding | Decision |
|---|---|---|
| MASH study 1: GSE316023 | wk23, 21 °C (TA): GAN 8 (GSM9443427–434) vs LabDiet 5001 chow 8 (GSM9443443–450); male C57BL/6J; 23 wk on diet; estimated counts (rounded). TN pair (8 v 8) as sensitivity training data | eligible |
| MASH study 2: GSE246328 | Vehicle 24 w (≈ 58 wk GAN; 7) vs Chow Vehicle 24 w (58 wk chow; 8); C57BL/6J(Rj) male; integer counts. GSE294629 re-deposits the same mice and is **not used** | eligible |
| Overlap | no shared GSM, SRA run or library between the two studies (Phase 2 overlap checks) | independent |
| Control diet | both studies use grain chow; GSE294180 uses a purified control | known mismatch; transfer to GSE294180 is checked descriptively, not assumed |
| **Phase 4D signature leakage** | **MASH_main_n50 required concordance in BOTH GSE316023 and GSE246328** | **cannot be used for study-held-out evaluation.** Signatures are rebuilt **inside each training study only** |
| Fasting: GSE107787 | 18 v 18 across 6 ZT; 12 of 18 FED livers carry muscle signal; muscle-free clean subset is ZT0/ZT4, 4 FAST v 6 FED; RPKM only | label partly confounded with muscle (FED, ZT8–20) and ZT |
| Fasting: GSE225695 | 3 v 3; sacrifice time unknown; counts | the only external fasting test; n = 6 |
| Fasting classifier | training on all 36 would let the model use muscle/ZT artefacts; training on the clean 10 is tiny; external n = 6 | **CONDITIONAL, demonstration only:** train on the clean subset only, one baseline + one regularised model, single direction; no model comparison |
| Expression comparability | studies differ in platform, counting tool, length handling and read length; no cross-study normalisation is valid with study ≡ batch | features use **within-sample ranks** (singscore) over a **fixed common gene universe**, so no between-study scaling is fitted. Ranks are invariant to CPM vs counts vs VST; length bias is platform-specific and may shift thresholds (calibration), so AUC is primary |
| GSE294180 length correction | Phase 4D used a median-transcript-length TPM-proxy | all datasets here use log-CPM ranks for consistency (no assumed TPM) |
| PE classifier on 19 mice | not built (forbidden; PE ≡ depth) | — |

## GO / NO-GO matrix

| Task | Decision |
|---|---|
| B1: MASH-state classifier, study-held-out (2 directions) | **GO** (2 studies only: two point estimates, wide CIs) |
| B2: fasting-state classifier | **CONDITIONAL:** clean-subset training, 1 external study of 6 mice, descriptive |
| C: apply a frozen MASH model to 19 mice | **only if B1 passes the predeclared validity criteria** (plan) |
| Fasting model applied to 19 mice | **NO-GO** (validation too weak); the Phase 4D FAST singscore remains the fasting-direction result |
