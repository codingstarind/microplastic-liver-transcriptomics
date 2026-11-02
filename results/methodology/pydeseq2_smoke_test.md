# pydeseq2 smoke test (pydeseq2 0.5.4)

## 1. Synthetic data, real group sizes (4/5/5/5)

Design columns: `['Intercept', 'diet[T.MD]', 'microplastic[T.PE]', 'diet[T.MD]:microplastic[T.PE]']`

| Planted effect (log2) | Estimated coefficient | Mean estimate on planted genes (log2) | Mean on null genes |
|---|---|---|---|
| 1.5 | `diet[T.MD]` | 1.51 | -0.01 |
| 1.0 | `microplastic[T.PE]` | 1.00 | -0.01 |
| 1.5 | `diet[T.MD]:microplastic[T.PE]` | 1.48 | -0.02 |

## 2. Real data: size factors only (no model fit, no tests)

Maximum relative difference between pydeseq2 and Phase 1 size factors: **0.00e+00** across 19 samples (13,777 filtered genes).

