"""Predict N-degron-stability-index with the upstream N-FIVE model.

The model predicts Protein Stability Index (PSI). Lower values indicate a
more destabilizing N-terminus in the E. coli reporter assay.
"""

from functools import lru_cache
from math import isfinite
from pathlib import Path
from pickle import UnpicklingError

import numpy as np

DEFAULT_MODEL_PATH = (
    Path(__file__).resolve().parents[2]
    / "vendor/N-terminal-cluster-stability/N-FIVE/N-FIVE WT model.pkl"
)
N_DEGRON_STABILITY_GROUP = "n-degron-stability-index"
WINDOW_LENGTH = 5
MET_EXCISION_PREFIX = "MSG"
# The model's columns use this published order at each position.
_RESIDUE_COLUMNS = {aa: index for index, aa in enumerate("ARNDCEQGHILKMFPSTWYV")}
_FEATURE_COUNT = WINDOW_LENGTH * len(_RESIDUE_COLUMNS)


@lru_cache(maxsize=1)
def _load_model(path: Path):
    if not path.is_file():
        raise ValueError(
            f"N-FIVE model not found at {path}. Initialize the submodule with "
            "'git submodule update --init vendor/N-terminal-cluster-stability' "
            "or supply its model file with --nfive-model."
        )
    try:
        import joblib
        from xgboost import XGBRegressor

        model = joblib.load(path)
    except ImportError as exc:
        raise ValueError("N-FIVE requires the optional dependencies from 'protein-property-metrics[nfive]'") from exc
    except (EOFError, KeyError, UnpicklingError) as exc:
        raise ValueError(f"could not read the N-FIVE model at {path}: {exc}") from exc
    if not isinstance(model, XGBRegressor) or model.n_features_in_ != _FEATURE_COUNT:
        raise ValueError(f"expected an N-FIVE regressor with {_FEATURE_COUNT} features at {path}")
    # Parallelism is controlled by the CLI's worker processes.
    model.set_params(n_jobs=1, device="cpu")
    return model


def _encode(windows: list[str]) -> np.ndarray:
    residue_count = len(_RESIDUE_COLUMNS)
    encoded = np.zeros((len(windows), _FEATURE_COUNT), dtype=np.float32)
    for row, window in enumerate(windows):
        for position, residue in enumerate(window):
            encoded[row, position * residue_count + _RESIDUE_COLUMNS[residue]] = 1.0
    return encoded


def describe_n_terminus(sequence: str, *, model_path: str | Path | None = None) -> dict[str, str | float | None]:
    """Score a normalized construct before and, for MSG, after initiator Met excision."""
    candidates = {
        "": sequence[:WINDOW_LENGTH],
        "_met_excised": sequence[1:1 + WINDOW_LENGTH] if sequence.startswith(MET_EXCISION_PREFIX) else "",
    }
    windows = {suffix: window for suffix, window in candidates.items() if len(window) == WINDOW_LENGTH}
    scores = {}
    if windows:
        path = DEFAULT_MODEL_PATH if model_path is None else Path(model_path).expanduser().resolve()
        predictions = _load_model(path).predict(_encode(list(windows.values())))
        scores = dict(zip(windows, map(float, predictions), strict=True))
        if not all(map(isfinite, scores.values())):
            raise ValueError("N-FIVE returned a non-finite Protein Stability Index")

    metrics = {}
    for suffix in candidates:
        metrics[f"n_degron_sequence{suffix}"] = windows.get(suffix)
        metrics[f"n_degron_stability_index{suffix}"] = scores.get(suffix)
    return metrics
