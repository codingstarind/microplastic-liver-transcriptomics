"""Model factory for task A (grids = model_search_spaces.csv). Every model is a Pipeline so scaling is
fitted on training rows only; tuned models are wrapped in GridSearchCV with GroupKFold(4) inner CV."""
from __future__ import annotations

from sklearn.compose import ColumnTransformer
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import ElasticNet, LinearRegression, Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeRegressor

SEED = 20261122
BASELINES = ["null", "ols_MASH", "ols_FAST", "ols_PPARA", "ols_CIRC", "ols_all"]
TUNED = ["ridge", "elasticnet", "tree", "rf", "hgb", "catboost"]
REDUCED = BASELINES + ["ridge", "rf", "hgb"]
LOO_SET = BASELINES + ["ridge", "hgb"]
COMPLEXITY = {"null": 0, "ols_MASH": 1, "ols_FAST": 1, "ols_PPARA": 1, "ols_CIRC": 1, "ols_all": 2, "ridge": 2,
              "elasticnet": 2, "tree": 3, "rf": 4, "hgb": 5, "catboost": 5}


def _select(cols, feats):
    return ColumnTransformer([("sel", "passthrough", [feats.index(c) for c in cols])], remainder="drop")


def make(name: str, feats: list[str]):
    """Return (estimator, param_grid or None)."""
    sc = ("scale", StandardScaler())
    if name == "null":
        return Pipeline([("model", DummyRegressor(strategy="mean"))]), None
    if name.startswith("ols_") and name != "ols_all":
        axis = name.split("_")[1]
        return Pipeline([("sel", _select([f"{axis}_eff"], feats)), sc, ("model", LinearRegression())]), None
    if name == "ols_all":
        return Pipeline([sc, ("model", LinearRegression())]), None
    if name == "ridge":
        return Pipeline([sc, ("model", Ridge())]), {"model__alpha": [0.01, 0.1, 1, 10, 100, 1000]}
    if name == "elasticnet":
        return Pipeline([sc, ("model", ElasticNet(max_iter=20000, random_state=SEED))]), {
            "model__alpha": [1e-4, 1e-3, 1e-2, 1e-1], "model__l1_ratio": [0.2, 0.5, 0.8, 1.0]}
    if name == "tree":
        return Pipeline([("model", DecisionTreeRegressor(random_state=SEED))]), {
            "model__max_depth": [2, 3, 4, 6], "model__min_samples_leaf": [50, 200]}
    if name == "rf":
        return Pipeline([("model", RandomForestRegressor(n_estimators=300, n_jobs=1, random_state=SEED))]), {
            "model__max_features": [0.33, 1.0], "model__min_samples_leaf": [20, 100], "model__max_depth": [None, 8]}
    if name == "hgb":
        return Pipeline([("model", HistGradientBoostingRegressor(max_iter=300, early_stopping=False, random_state=SEED))]), {
            "model__learning_rate": [0.03, 0.1], "model__max_leaf_nodes": [8, 31], "model__min_samples_leaf": [50, 200]}
    if name == "catboost":
        from catboost import CatBoostRegressor
        return Pipeline([("model", CatBoostRegressor(iterations=300, random_seed=SEED, thread_count=1, verbose=0,
                                                     allow_writing_files=False))]), {
            "model__depth": [4, 6], "model__learning_rate": [0.03, 0.1], "model__l2_leaf_reg": [3, 10]}
    raise ValueError(name)
