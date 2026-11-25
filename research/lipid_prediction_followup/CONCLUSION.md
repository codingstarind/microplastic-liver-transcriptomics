# Short conclusion: what the gene-level regression adds (exploratory)

**Scope:** F8-OLS and F8+C-Ridge are the models used to answer the question. FULL-OLS is reported but not used as evidence, because it contains membership indicators for these same three gene sets. All numbers are in `EVIDENCE_AUDIT.md` and `TABLE_results_concise.csv`.

## What was tested
- The models predict each gene's PE log2FC from that gene's responses in independent reference datasets: MASH, fasting, the PPARα agonist and circadian time. F8+C adds four response-shape terms.
- The prediction rule is learned across the genome on training modules and applied to genes in held-out co-expression modules.
- The question is whether this one rule predicts lipid-set genes better than other genes. No lipid-specific model was fitted.

## Findings

**1. Cholesterol synthesis (21 genes, 6 modules): clearly better predicted than other genes, in both diets and both models.**

| Model | Diet | Held-out R² in set | Held-out R² other genes | Gap vs genes matched on null error [95% CI] |
|---|---|---|---|---|
| F8 | CD | 0.65 | 0.05 | +0.58 [0.50, 0.80] |
| F8 | MD | 0.60 | 0.06 | +0.56 [0.39, 0.68] |
| F8+C | CD | 0.74 | 0.06 | +0.64 [0.60, 0.79] |
| F8+C | MD | 0.52 | 0.08 | +0.47 [0.35, 0.59] |

- Every gene in the set was better predicted than by the null.
- The gap holds within modules (module-adjusted p ≤ 0.0006).
- The gap holds against other genes with similar-sized changes.
- The gap holds when every gene in the other pathway sets is excluded from the reference.

**2. Lipogenesis (12 genes, 5 modules): better than other genes, with weaker and less stable support.**
- Module-adjusted gains are positive for both models in both diets (p 0.0002–0.026).
- Against genes matched on null error:
  - F8 in CD gives p = 0.054.
  - Most R² intervals have lower bounds close to 0.
- With only 5 modules, this result rests on little independent information.

**3. PPARα targets (14 genes, 9 modules): no reliable advantage. This is a negative result.**
- Module-adjusted differences from other genes are small, and every interval includes 0 (p 0.20–0.50).
- One exception points the other way: for F8 in MD, the module-bootstrap R² difference excludes 0 (+0.30 [0.25, 0.87]).
- F8+C in MD is indistinguishable from other genes (R² 0.07 vs 0.08).

**4. Adding response-shape terms (F8+C vs F8) does not consistently help the lipid genes.**
- **CD:** improves them more than other genes (+0.10 [0.02, 0.19] in squared error).
- **MD:** predicts them worse than F8 does (−0.12 [−0.19, −0.04]).

**5. FULL scores highest on these sets**, for example R² 0.89 for cholesterol synthesis in CD. It uses these sets' membership as features, so this is not independent evidence.

## What this contributes biologically

**What the regression shows:**
- The PE-associated changes in cholesterol-synthesis genes, and to a lesser extent lipogenesis genes, closely follow how those same genes respond in the independent reference states.
- A rule learned from other gene modules explains about 50–75% of these genes' held-out variation, compared with about 5–8% for genes outside the three sets.
- This holds in both diets.
- It provides model-based support, independent of gene-set testing, that the PE sterol and lipogenic changes share structure with known reference programmes and are not idiosyncratic to these genes.
- It does not support the same statement for the PPARα target set.

**What it does not show:**
- Which reference programme drives the fit. The per-gene contributions of each reference axis were not decomposed.
- Mechanism, causality, or that PE worsens MASH.
- That CD and MD differ. The CD and MD intervals overlap for the main comparisons.

**Limits:**
- These are exploratory reuses of splits already seen in Phases 4–4C.
- Only 5–13 modules carry the lipid-gene comparisons.
- Lipid genes have larger PE changes than other genes. The comparison against genes matched on null error addresses this only approximately; lipogenesis matching in CD is imperfect (median null error 2.52 vs 1.90).
