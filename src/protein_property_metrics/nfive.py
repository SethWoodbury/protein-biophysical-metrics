"""Predict N-degron-stability-index with the upstream N-FIVE model.

The model predicts Protein Stability Index (PSI). Lower values indicate a
more destabilizing N-terminus in the E. coli reporter assay.
"""

from collections import OrderedDict
from collections.abc import Collection, Iterable, Iterator, Mapping
from functools import lru_cache
from hashlib import sha256
from itertools import product
from pathlib import Path
from pickle import UnpicklingError
from tempfile import TemporaryDirectory
from threading import Lock
from zipfile import BadZipFile

import numpy as np

DEFAULT_MODEL_PATH = (
    Path(__file__).resolve().parents[2]
    / "vendor/N-terminal-cluster-stability/N-FIVE/N-FIVE WT model.pkl"
)
N_DEGRON_STABILITY_GROUP = "n-degron-stability-index"
WINDOW_LENGTH = 5
MET_EXCISION_PREFIX = "MSG"
# The model's columns use this published order at each position.
_ALPHABET = "ARNDCEQGHILKMFPSTWYV"
_RESIDUE_COLUMNS = {aa: index for index, aa in enumerate(_ALPHABET)}
_FEATURE_COUNT = WINDOW_LENGTH * len(_ALPHABET)
_SCORE_SUFFIXES = ("", "_met_excised")
_LOOKUP_PREFIXES = (MET_EXCISION_PREFIX, MET_EXCISION_PREFIX[1:])
_LOOKUP_VERSION = 1
# Bound sparse live score entries and the dense encoding matrix per process.
_PREDICTION_BATCH_ROWS = 4096
_DIRECT_CACHE_CAPACITY = 16_384
_DIRECT_CACHE: OrderedDict[tuple[str, str], float] = OrderedDict()
_DIRECT_CACHE_LOCK = Lock()
_FileStamp = tuple[int, int, int, int, int]


def _file_stamp(path: Path) -> _FileStamp:
    stat = path.stat()
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


def _check_file_stamp(path: Path, stamp: _FileStamp) -> None:
    if _file_stamp(path) != stamp:
        raise ValueError(f"N-FIVE input changed while it was being read: {path}")


@lru_cache(maxsize=4)
def _model_digest(path: Path, stamp: _FileStamp) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    _check_file_stamp(path, stamp)
    return digest.hexdigest()


def _model_source(model_path: str | Path | None) -> tuple[Path, _FileStamp, str]:
    path = DEFAULT_MODEL_PATH if model_path is None else Path(model_path).expanduser().resolve()
    if not path.is_file():
        raise ValueError(
            f"N-FIVE model not found at {path}. Initialize the submodule with "
            "'git submodule update --init vendor/N-terminal-cluster-stability' "
            "or supply its model file with --nfive-model."
        )
    stamp = _file_stamp(path)
    return path, stamp, _model_digest(path, stamp)


@lru_cache(maxsize=1)
def _load_model(path: Path, stamp: _FileStamp):
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
    _check_file_stamp(path, stamp)
    # Score each batch on one CPU thread, independently of CLI workers.
    model.set_params(n_jobs=1, device="cpu")
    return model


def _encode(windows: list[str]) -> np.ndarray:
    residue_count = len(_RESIDUE_COLUMNS)
    encoded = np.zeros((len(windows), _FEATURE_COUNT), dtype=np.float32)
    for row, window in enumerate(windows):
        for position, residue in enumerate(window):
            encoded[row, position * residue_count + _RESIDUE_COLUMNS[residue]] = 1.0
    return encoded


def nfive_windows(sequence: str) -> dict[str, str]:
    """Return the scored windows of a normalized construct in output column order."""
    candidates = {
        "": sequence[:WINDOW_LENGTH],
        "_met_excised": sequence[1:1 + WINDOW_LENGTH] if sequence.startswith(MET_EXCISION_PREFIX) else "",
    }
    return {suffix: window for suffix, window in candidates.items() if len(window) == WINDOW_LENGTH}


def _predict_windows(path: Path, stamp: _FileStamp, windows: list[str]) -> Iterator[tuple[str, float]]:
    model = _load_model(path, stamp)
    for start in range(0, len(windows), _PREDICTION_BATCH_ROWS):
        batch = windows[start:start + _PREDICTION_BATCH_ROWS]
        predictions = np.asarray(model.predict(_encode(batch)))
        if predictions.shape != (len(batch),):
            raise ValueError("N-FIVE returned an unexpected number of scores")
        if not np.isfinite(predictions).all():
            raise ValueError("N-FIVE returned a non-finite Protein Stability Index")
        yield from zip(batch, map(float, predictions), strict=True)


def _lookup_shape(prefix: str) -> tuple[int, ...]:
    return (len(_ALPHABET),) * (WINDOW_LENGTH - len(prefix))


@lru_cache(maxsize=2)
def _load_lookup(path: Path, stamp: _FileStamp, model_hash: str) -> dict[str, np.ndarray]:
    try:
        with np.load(path, allow_pickle=False) as archive:
            if archive["version"].item() != _LOOKUP_VERSION:
                raise ValueError("unsupported N-FIVE lookup version")
            if archive["alphabet"].item() != _ALPHABET:
                raise ValueError("N-FIVE lookup alphabet does not match the model encoding")
            if archive["model_sha256"].item() != model_hash:
                raise ValueError("N-FIVE lookup model SHA256 does not match the selected model")
            tables = {prefix: archive[prefix.lower()] for prefix in _LOOKUP_PREFIXES}
        for prefix, scores in tables.items():
            if scores.shape != _lookup_shape(prefix) or scores.dtype != np.float32 or not np.isfinite(scores).all():
                raise ValueError(f"invalid float32 N-FIVE lookup table for {prefix}")
        _check_file_stamp(path, stamp)
    except (BadZipFile, EOFError, KeyError, OSError, TypeError, ValueError) as exc:
        raise ValueError(f"could not read N-FIVE lookup {path}: {exc}") from exc
    return tables


def _lookup_score(window: str, tables: Mapping[str, np.ndarray]) -> float | None:
    for prefix, values in tables.items():
        if window.startswith(prefix):
            indices = tuple(_RESIDUE_COLUMNS[aa] for aa in window[len(prefix):])
            return float(values[indices])
    return None


def prepare_nfive_scores(
    sequences: Iterable[str],
    *,
    model_path: str | Path | None = None,
    lookup_path: str | Path | None = None,
) -> dict[str, float]:
    """Score unique complete windows across normalized constructs in CPU batches."""
    windows = dict.fromkeys(window for sequence in sequences for window in nfive_windows(sequence).values())
    return _score_windows(windows, model_path, lookup_path)


def _score_windows(
    windows: Collection[str],
    model_path: str | Path | None,
    lookup_path: str | Path | None,
) -> dict[str, float]:
    if not windows:
        return {}

    path, stamp, model_hash = _model_source(model_path)
    tables = {}
    if lookup_path is not None:
        lookup = Path(lookup_path).expanduser().resolve()
        tables = _load_lookup(lookup, _file_stamp(lookup), model_hash)

    scores = {}
    missing = []
    with _DIRECT_CACHE_LOCK:
        for window in windows:
            score = _lookup_score(window, tables)
            if score is None:
                key = (model_hash, window)
                score = _DIRECT_CACHE.get(key)
                if score is not None:
                    _DIRECT_CACHE.move_to_end(key)
            if score is None:
                missing.append(window)
            else:
                scores[window] = score
        if missing:
            predictions = dict(_predict_windows(path, stamp, missing))
            scores.update(predictions)
            for window, score in predictions.items():
                _DIRECT_CACHE[model_hash, window] = score
                if len(_DIRECT_CACHE) > _DIRECT_CACHE_CAPACITY:
                    _DIRECT_CACHE.popitem(last=False)
    return scores


def describe_n_terminus(
    sequence: str,
    *,
    model_path: str | Path | None = None,
    lookup_path: str | Path | None = None,
    prepared_scores: Mapping[str, float] | None = None,
) -> dict[str, str | float | None]:
    """Score a normalized construct before and, for MSG, after Met excision.

    ``prepared_scores`` supplies the complete window scores from a prepared
    batch and bypasses model and lookup reads.
    """
    windows = nfive_windows(sequence)
    scores = prepared_scores
    if scores is None:
        scores = _score_windows(windows.values(), model_path, lookup_path)

    metrics = {}
    for suffix in _SCORE_SUFFIXES:
        window = windows.get(suffix)
        metrics[f"n_degron_sequence{suffix}"] = window
        try:
            metrics[f"n_degron_stability_index{suffix}"] = None if window is None else scores[window]
        except KeyError as exc:
            raise ValueError(f"missing prepared N-FIVE score for {window!r}") from exc
    return metrics


def build_msg_lookup(output: str | Path, *, model_path: str | Path | None = None) -> Path:
    """Build a non-pickle NPZ of MSGXX and SGXXX windows from a trusted model."""
    output = Path(output).expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"N-FIVE lookup output already exists: {output}")
    path, stamp, model_hash = _model_source(model_path)
    tables = {}
    for prefix in _LOOKUP_PREFIXES:
        windows = [prefix + "".join(suffix) for suffix in product(_ALPHABET, repeat=WINDOW_LENGTH - len(prefix))]
        scores = np.fromiter((score for _, score in _predict_windows(path, stamp, windows)), dtype=np.float32)
        tables[prefix.lower()] = scores.reshape(_lookup_shape(prefix))

    with TemporaryDirectory(dir=output.parent) as directory:
        staged = Path(directory) / "lookup.npz"
        np.savez(staged, version=_LOOKUP_VERSION, alphabet=_ALPHABET, model_sha256=model_hash, **tables)
        # Publish complete contents without overwriting an existing path.
        output.hardlink_to(staged)
    return output


def _main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Generate an N-FIVE fixed-MSG lookup from a trusted model.")
    parser.add_argument("--build-msg-lookup", required=True, type=Path, metavar="PATH", help="New output NPZ path")
    parser.add_argument("--nfive-model", type=Path, help="Trusted model path. Defaults to the pinned checkout model.")
    args = parser.parse_args()
    try:
        print(build_msg_lookup(args.build_msg_lookup, model_path=args.nfive_model))
    except (OSError, ValueError) as exc:
        parser.exit(2, f"{parser.prog}: error: {exc}\n")


if __name__ == "__main__":
    _main()
