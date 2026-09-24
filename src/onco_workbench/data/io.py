"""Read and write expression matrices, sample metadata, ground truth, and data manifests.

Files are UTF-8 CSVs with LF line endings. They begin with ``#`` comment lines that
carry the DEMONSTRATION / SYNTHETIC label, and loaders read them back with
``comment="#"``. Writers use a fixed float format so that output is byte-identical
across operating systems.
"""

from __future__ import annotations

import csv
import dataclasses
import hashlib
import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from onco_workbench.config import SyntheticConfig
from onco_workbench.data.synthetic import SAMPLE_ID_COLUMN, SyntheticDataset
from onco_workbench.disclaimers import DATA_LABEL, SHORT_DISCLAIMER, as_comment_block

DATA_LICENSE: dict[str, str] = {
    "spdx": "CC0-1.0",
    "name": "Creative Commons Zero v1.0 Universal (public-domain dedication)",
    "url": "https://creativecommons.org/publicdomain/zero/1.0/",
    "scope": (
        "Applies only to the synthetic data files listed in this manifest. The source "
        "code is licensed separately under the MIT License (see LICENSE); third-party "
        "dependencies keep their own licenses."
    ),
}


class DataLoadError(ValueError):
    """Raised when an input file is missing, empty, or cannot be parsed."""


def _require_file(path: Path, description: str) -> None:
    if not path.is_file():
        raise DataLoadError(
            f"{description} file not found: {path}. "
            "Generate the synthetic demo data with `obw generate-data`, or pass a valid path."
        )


def _read_header(path: Path) -> list[str]:
    """Return the first non-comment CSV row. Duplicate names are preserved."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = csv.reader(line for line in handle if not line.startswith("#"))
        return next(rows, [])


def load_expression(path: str | Path) -> pd.DataFrame:
    """Load an expression matrix with samples as rows and genes as columns.

    The first column holds sample identifiers. Column names are taken verbatim from
    the header, so duplicated gene IDs are preserved for validation to report
    instead of being silently renamed by pandas.

    Args:
        path: CSV file path.

    Returns:
        DataFrame indexed by ``sample_id`` with one column per gene.

    Raises:
        DataLoadError: If the file is missing, empty, or has no gene columns.
    """
    path = Path(path)
    _require_file(path, "Expression matrix")
    header = _read_header(path)
    if len(header) < 2:
        raise DataLoadError(
            f"Expression matrix {path} must have a sample ID column followed by at least "
            "one gene column."
        )
    try:
        frame = pd.read_csv(path, comment="#", index_col=0, converters={0: str})
    except (pd.errors.ParserError, pd.errors.EmptyDataError, UnicodeDecodeError) as exc:
        raise DataLoadError(f"Could not parse expression matrix {path}: {exc}") from exc
    if frame.shape[1] != len(header) - 1:
        raise DataLoadError(
            f"Expression matrix {path}: header has {len(header) - 1} gene columns but "
            f"{frame.shape[1]} were parsed."
        )
    frame.columns = pd.Index(header[1:], dtype=object)
    frame.index.name = SAMPLE_ID_COLUMN
    return frame


def load_metadata(path: str | Path) -> pd.DataFrame:
    """Load sample metadata as strings. Only empty cells are treated as missing.

    Args:
        path: CSV file path.

    Returns:
        DataFrame with one row per sample.

    Raises:
        DataLoadError: If the file is missing or cannot be parsed.
    """
    path = Path(path)
    _require_file(path, "Sample metadata")
    try:
        return pd.read_csv(path, comment="#", dtype=str, keep_default_na=False, na_values=[""])
    except (pd.errors.ParserError, pd.errors.EmptyDataError, UnicodeDecodeError) as exc:
        raise DataLoadError(f"Could not parse sample metadata {path}: {exc}") from exc


def load_ground_truth(path: str | Path) -> pd.DataFrame:
    """Load the synthetic ground-truth table (planted signal genes).

    Args:
        path: CSV file path.

    Returns:
        DataFrame with ``gene_id``, ``is_signal``, ``direction``, and ``true_shift``.

    Raises:
        DataLoadError: If the file is missing or cannot be parsed.
    """
    path = Path(path)
    _require_file(path, "Ground-truth")
    try:
        return pd.read_csv(path, comment="#", dtype={"gene_id": str, "direction": str})
    except (pd.errors.ParserError, pd.errors.EmptyDataError, UnicodeDecodeError) as exc:
        raise DataLoadError(f"Could not parse ground-truth table {path}: {exc}") from exc


def file_sha256(path: str | Path) -> str:
    """Return the hexadecimal SHA-256 digest of a file's bytes.

    Args:
        path: File path.

    Returns:
        64-character lowercase hex digest.
    """
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_labeled_csv(
    frame: pd.DataFrame,
    path: str | Path,
    *,
    index: bool,
    float_decimals: int | None = None,
    float_format: str | None = None,
    extra_comments: Sequence[str] = (),
) -> Path:
    """Write a CSV that starts with the synthetic-data comment header.

    Args:
        frame: Table to write.
        path: Destination file; parent directories are created.
        index: Whether to write the DataFrame index as the first column.
        float_decimals: Fixed number of decimals for floats (e.g. ``3`` -> ``%.3f``).
        float_format: Explicit printf-style float format; used when ``float_decimals``
            is not given. Defaults to ``"%.6g"``, which keeps very small p-values readable.
        extra_comments: Additional ``#`` comment lines written after the standard header.

    Returns:
        The written path.

    Raises:
        ValueError: If both ``float_decimals`` and ``float_format`` are given.
    """
    if float_decimals is not None and float_format is not None:
        raise ValueError("Pass either float_decimals or float_format, not both.")
    fmt = f"%.{float_decimals}f" if float_decimals is not None else (float_format or "%.6g")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(as_comment_block("# "))
        for line in extra_comments:
            handle.write(f"# {line}\n")
        frame.to_csv(handle, index=index, float_format=fmt, lineterminator="\n", na_rep="")
    return path


def write_dataset(
    dataset: SyntheticDataset,
    config: SyntheticConfig,
    *,
    expression_path: str | Path,
    metadata_path: str | Path,
    ground_truth_path: str | Path,
    manifest_path: str | Path,
) -> dict[str, Path]:
    """Write a synthetic dataset and its provenance manifest.

    The manifest records the seed, generator parameters, CC0-1.0 data dedication,
    and SHA-256 checksums. It contains no timestamps or environment details, so
    regenerating with the same seed reproduces it byte for byte.

    Args:
        dataset: Generated dataset.
        config: Generator parameters that produced ``dataset``.
        expression_path: Destination of the expression matrix.
        metadata_path: Destination of the sample metadata.
        ground_truth_path: Destination of the ground-truth table.
        manifest_path: Destination of the JSON manifest.

    Returns:
        Mapping of ``expression``, ``metadata``, ``ground_truth``, and ``manifest``
        to the written paths.
    """
    decimals = config.float_decimals
    written = {
        "expression": write_labeled_csv(
            dataset.expression, expression_path, index=True, float_decimals=decimals
        ),
        "metadata": write_labeled_csv(
            dataset.metadata, metadata_path, index=False, float_decimals=decimals
        ),
        "ground_truth": write_labeled_csv(
            dataset.ground_truth, ground_truth_path, index=False, float_decimals=decimals
        ),
    }
    shapes = {
        "expression": dataset.expression.shape,
        "metadata": dataset.metadata.shape,
        "ground_truth": dataset.ground_truth.shape,
    }
    manifest: dict[str, Any] = {
        "data_label": DATA_LABEL,
        "disclaimer": SHORT_DISCLAIMER,
        "description": (
            "Fully synthetic demonstration data simulated by "
            "onco_workbench.data.synthetic.generate_synthetic_dataset. Not patient data. "
            "Gene identifiers do not correspond to real genes; group labels are not diagnoses."
        ),
        "license": DATA_LICENSE,
        "seed": dataset.seed,
        "random_generator": "numpy.random.default_rng (PCG64)",
        "synthetic_config": dataclasses.asdict(config),
        "files": {
            written[key].name: {
                "role": key,
                "rows": shapes[key][0],
                "columns": shapes[key][1],
                "bytes": written[key].stat().st_size,
                "sha256": file_sha256(written[key]),
            }
            for key in ("expression", "metadata", "ground_truth")
        },
    }
    manifest_file = Path(manifest_path)
    manifest_file.parent.mkdir(parents=True, exist_ok=True)
    with manifest_file.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    written["manifest"] = manifest_file
    return written


def load_manifest(path: str | Path) -> dict[str, Any]:
    """Load a data manifest written by :func:`write_dataset`.

    Args:
        path: JSON file path.

    Returns:
        The parsed manifest.

    Raises:
        DataLoadError: If the file is missing or is not valid JSON.
    """
    path = Path(path)
    _require_file(path, "Data manifest")
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise DataLoadError(f"Data manifest {path} is not valid JSON: {exc}") from exc
