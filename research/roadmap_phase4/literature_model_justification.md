# Model choices: literature justification

Citations are to the original methodological sources or standard references. They were chosen for documented algorithm properties, not for performance in our data.

## Task A: predicting PE gene effects from reference effects

| Model | Why it could help here | Limitation in this setting | Source |
|---|---|---|---|
| Null (mean) | Defines R²_oos = 0; any useful predictor must beat it | — | — |
| Single-axis OLS | The regression form of the Phase 3 correlation baselines; tests each state separately | Ignores the other axes | — |
| OLS (all features) | The direct "multi-axis" model from the Phase 3 memo; coefficients are interpretable | Coefficients are attenuated by errors in the reference effects (errors in variables) and are not causal | Fuller 1987, *Measurement Error Models* |
| Ridge | Shrinks correlated coefficients (effect and SE features are correlated within an axis); lower variance | Shrinkage biases coefficients toward 0 | Hoerl & Kennard 1970, *Technometrics* 12:55 |
| Lasso / Elastic Net | Tests whether a sparse subset of axes suffices; Elastic Net handles correlated features better than Lasso | Selection is unstable under collinearity; only 8 features, so the gain from sparsity is small | Tibshirani 1996, *JRSS-B* 58:267; Zou & Hastie 2005, *JRSS-B* 67:301 |
| Shallow decision tree (CART) | Captures thresholds and interactions, e.g. "PE follows fasting only where the fasting effect is large" | High variance; piecewise-constant predictions | Breiman, Friedman, Olshen & Stone 1984, *CART* |
| Random forest / bagged trees | Averages many trees to reduce variance; nonlinear; bagging (`max_features = 1.0`) is a special case | Impurity importance is biased toward high-variance features, and permutation importance is diluted under correlation | Breiman 1996, *Mach Learn* 24:123; Breiman 2001, *Mach Learn* 45:5; Strobl et al. 2007, *BMC Bioinf* 8:25; Strobl et al. 2008, *BMC Bioinf* 9:307 |
| Histogram gradient boosting | Stage-wise additive trees fit structured nonlinearity with shrinkage; efficient on ~10⁴ rows | Can overfit noisy targets without regularisation; many hyperparameters | Friedman 2001, *Ann Stat* 29:1189; Ke et al. 2017 (LightGBM histogram algorithm), *NeurIPS* |
| CatBoost | A separate boosting implementation (ordered boosting, oblivious trees) that reduces prediction shift; a check that boosting results do not depend on one implementation | Mostly redundant with HGB when features are continuous; included as a single cross-implementation check | Prokhorenkova et al. 2018, *NeurIPS* |
| AdaBoost.R2 (**not run**) | Reweights hard examples | Sensitive to heavy-tailed noise and outliers, which log2FC targets have; adds nothing beyond gradient boosting here | Freund & Schapire 1997, *JCSS* 55:119; Drucker 1997, *ICML* |
| Neural networks (**not run**) | — | 8 features with weak signal; no structural reason to expect a gain; unnecessary complexity | — |

**Not used as features:** gene identity, target significance, gene-set labels. These would leak target information or give the model a lookup table (Ambroise & McLachlan 2002, *PNAS* 99:6562).

## Validation design
- **Nested CV:** hyperparameters chosen in an inner loop, separate from the outer evaluation, because using the same folds for tuning and evaluation is optimistically biased (Varma & Simon 2006, *BMC Bioinf* 7:91; Cawley & Talbot 2010, *JMLR* 11:2079; Krstajic et al. 2014, *J Cheminform* 6:10, on repeated nested CV).
- **Grouped (blocked) holdout:** correlated observations split across folds inflate CV scores; blocking by structure is recommended (Roberts et al. 2017, *Ecography* 40:913). Here genes are blocked by co-expression module, with chromosome blocks as a sensitivity analysis.
- **Comparing models across CV splits:** folds overlap in their training data, so naive paired t-tests are anti-conservative. The Nadeau–Bengio correction is used (Nadeau & Bengio 2003, *Mach Learn* 52:239; Bouckaert & Frank 2004). CV estimates have their own uncertainty that is not captured by fold SD (Bates, Hastie & Tibshirani 2023, *JASA*).
- **Held-out genes are not held-out cohorts.** Cross-study validation is a different and harder standard (Bernau et al. 2014, *Bioinformatics* 30:i105). It is not claimed.
- **Errors in variables:** each reference effect is an estimate with SE; regression on noisy predictors attenuates the coefficients (Fuller 1987). SE features let tree models down-weight unreliable reference effects, but this does not correct the bias.

## Task B (demonstration)
- **Regularised logistic regression:** the standard high-dimensional baseline.
- **Gaussian naive Bayes:** reaches its asymptotic error faster with little data, but assumes conditional independence (Ng & Jordan 2002, *NeurIPS*; Hand & Yu 2001, *Int Stat Rev* 69:385). It is used **only** for the class label, never for continuous targets.
- **Gini vs entropy trees:** they usually give near-identical splits (Raileanu & Stoffel 2004, *Ann Math AI* 41:77). Both are shown once, as requested, and then treated as redundant.
- **Random forest and boosting:** as in task A.
- **Small-n warning:** with 32 training mice and a near-perfectly separable class, accuracy differences between classifiers are not interpretable (Hand 2006, *Stat Sci* 21:1 on the "illusion of progress").
- **Redundancy pruning:** bagged trees ≈ RF with max_features = 1 (not run separately in B); CatBoost ≈ HGB for B (not run).
