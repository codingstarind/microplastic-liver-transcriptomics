# Validation diagnostics (`results/validation_diagnostics_table.csv`)

| Issue | Evidence | Consequence |
|---|---|---|
| **Shared controls (Hong)** | CON1–3 byte-identical across 4 series; identical SRA read counts; 24 GSMs = 15 mice. Cross-polymer gene-effect ρ 0.57–0.77 (≈ 0.5 expected from the shared control alone) | One study, not four. Polymer contrasts are correlated; specificity is untestable beyond "PE vs other polymers" |
| **Atypical controls?** | Hong CON is the lowest of 5 groups on FAST and R2, and the highest on MASH and R0. Control-free PE − other polymers on FAST = +0.010 | The external FAST shift may be partly control-driven |
| **Depth/yield** | Exposed/control ratios: Hong 0.95–0.98, GSE245069 0.97–1.0, GSE277772 0.97–1.0 | Unlike GSE294180 (PE ≡ lower depth), the external contrasts are not depth-confounded |
| **Muscle carry-over** | External muscle CPM 570–1,030 vs GSE294180 median 1,076; exposed/control ratio 1.06–1.18 | Comparable; no strong imbalance |
| **Sex** | Xist ≈ 0 and Ddx3y > 0 in all external mice | All male, as annotated |
| **Gene coverage** | 4E background 11,853–11,883 of 11,883; R0 100/100 genes; 4D background 12,357–13,574 of 13,777; ≤ 10 signature genes missing | Coverage control: GSE294180 estimates change ≤ 9% (`coverage_control_gse294180_rescored.csv`) |
| **Noise ceiling** | Signal share: Hong PE 0.47, PS 0.20, PP 0.61, PB 0.52; GSE245069 pooled 0.65 | Hong PS is borderline uninterpretable; the others are noisy but usable |
| **Gene-wise vs mouse-wise inference** | F8 Hong PE gene-block CI [0.14, 0.28] excludes 0, but the mouse-label permutation p is 0.25 | Genes share the same 6 mice; the mouse permutation is the valid null |
| **Fasting/clock** | Sacrifice time / feeding state not reported in any external study | Fasting-like shifts cannot be separated from husbandry timing |
| **Domain shift** | Chow vs purified diet, 71 vs 8 wk, food vs gavage, ICR vs C57BL/6, different pipelines | Transfer is a hard test; failure is not evidence of no effect |
| **Small n** | 3 v 3 (min exact p 0.10); 4 v 4 (0.029); 12 v 4 (0.0011) | Effect sizes and CIs are primary |
