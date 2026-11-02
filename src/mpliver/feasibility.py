"""Planning calculations for the Phase 2 methodology proposal.

Everything here is a rough, descriptive planning aid: no function tests genes
for differential expression. The approximations and their limits are spelled
out in each docstring and in reports/methodology_proposal.md.
"""

import numpy as np
import pandas as pd
from scipy.stats import norm

CELLS = ["CD+Veh", "CD+PE", "MD+Veh", "MD+PE"]

# Contrast weights over the four cells, in CELLS order.
CONTRASTS = {
    "PE within CD (CD+PE vs CD+Veh)": [-1, 1, 0, 0],
    "PE within MD (MD+PE vs MD+Veh)": [0, 0, -1, 1],
    "Diet within Veh (MD+Veh vs CD+Veh)": [-1, 0, 1, 0],
    "PE main effect (average over diets)": [-0.5, 0.5, -0.5, 0.5],
    "Diet x PE interaction (difference of PE effects)": [1, -1, -1, 1],
}


def design_matrix(diet_md: np.ndarray, pe: np.ndarray, index) -> pd.DataFrame:
    """Treatment-coded 2 x 2 factorial design (reference: CD, vehicle)."""
    diet_md = np.asarray(diet_md, dtype=float)
    pe = np.asarray(pe, dtype=float)
    return pd.DataFrame({"intercept": 1.0, "diet_MD": diet_md, "PE": pe,
                         "diet_MD:PE": diet_md * pe}, index=index)


def r2_on_design(X: pd.DataFrame, covariate) -> dict:
    """How much of a technical covariate the design explains (R^2), and the implied VIF.

    A high R^2 means adding the covariate to the model would inflate the
    standard errors of the design coefficients by roughly sqrt(VIF).
    """
    y = np.asarray(covariate, dtype=float)
    beta, *_ = np.linalg.lstsq(X.to_numpy(), y, rcond=None)
    resid = y - X.to_numpy() @ beta
    r2 = 1 - resid.var() / y.var()
    return {"r2_explained_by_design": float(r2),
            "vif": float(1 / (1 - r2)) if r2 < 1 else float("inf")}


def moment_dispersion(norm_counts: pd.DataFrame, groups: pd.Series,
                      size_factors: pd.Series) -> pd.Series:
    """Pooled within-group method-of-moments NB dispersion per gene.

    For normalised counts y = k / s, Var(y) ~ mu * mean(1/s) + phi * mu^2, so
    phi ~ (var - mu * mean(1/s)) / mu^2, pooled over groups by degrees of freedom.

    Limitations (why this is only a planning number): with 3-4 residual degrees
    of freedom per group each gene's estimate is very noisy; negative estimates
    are clipped to a small floor (upward bias for low-dispersion genes); it is
    not the DESeq2 estimator (no Cox-Reid adjustment, no shrinkage to a trend);
    and contamination or outlying mice inflate it.
    """
    weighted, df_total = 0.0, 0
    for g in CELLS:
        cols = groups.index[groups == g]
        sub = norm_counts[cols]
        mu = sub.mean(axis=1)
        poisson_part = mu * (1 / size_factors[cols]).mean()
        phi = (sub.var(axis=1) - poisson_part) / mu ** 2
        weighted = weighted + phi * (len(cols) - 1)
        df_total += len(cols) - 1
    return (weighted / df_total).clip(lower=1e-4)


def minimum_detectable_lfc(mean_count: float, dispersion: float, n_cells,
                           contrast, alpha: float, power: float = 0.8) -> float:
    """Smallest |log2FC| detectable by a two-sided Wald test (normal approximation).

    Var(log mu_hat) for a cell of n mice ~ (1/mu + phi) / n. Ignores
    dispersion-estimation error, shrinkage, the t-like behaviour of small
    samples, and the fact that the genome-wide FDR threshold is data-dependent;
    ``alpha`` is a fixed stand-in for that threshold.
    """
    var_cell = (1 / mean_count + dispersion) / np.asarray(n_cells, dtype=float)
    se = np.sqrt(np.sum(np.asarray(contrast, dtype=float) ** 2 * var_cell))
    return float((norm.ppf(1 - alpha / 2) + norm.ppf(power)) * se / np.log(2))


def leave_one_mouse_out(groups: pd.Series):
    """Yield (dropped mouse, remaining labels) for each of the mice in turn.

    Used instead of a bootstrap: resampling 4-5 mice with replacement duplicates
    animals, which shrinks within-group variance and biases dispersion downwards.
    """
    for mouse in groups.index:
        yield mouse, groups.drop(mouse)


def jackknife_interval(estimates: np.ndarray, full_estimate: float, z: float = 1.96):
    """Approximate interval from leave-one-out estimates (jackknife standard error)."""
    est = np.asarray(estimates, dtype=float)
    n = len(est)
    se = np.sqrt((n - 1) / n * np.sum((est - est.mean()) ** 2))
    return full_estimate - z * se, full_estimate + z * se
