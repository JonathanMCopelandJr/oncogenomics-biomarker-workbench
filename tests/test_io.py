from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pandas as pd
import pytest

from onco_workbench.config import SyntheticConfig
from onco_workbench.data.io import (
    DataLoadError,
    file_sha256,
    load_expression,
    load_ground_truth,
    load_manifest,
    load_metadata,
    write_dataset,
)
from onco_workbench.data.synthetic import SyntheticDataset
from onco_workbench.disclaimers import DATA_LABEL


def _write(tmp_path: Path, dataset: SyntheticDataset, config: SyntheticConfig) -> dict[str, Path]:
    return write_dataset(
        dataset,
        config,
        expression_path=tmp_path / "expr.csv",
        metadata_path=tmp_path / "meta.csv",
        ground_truth_path=tmp_path / "truth.csv",
        manifest_path=tmp_path / "manifest.json",
    )


def test_round_trip(
    tmp_path: Path, small_dataset: SyntheticDataset, small_config: SyntheticConfig
) -> None:
    paths = _write(tmp_path, small_dataset, small_config)
    expression = load_expression(paths["expression"])
    pd.testing.assert_frame_equal(expression, small_dataset.expression, check_column_type=False)
    metadata = load_metadata(paths["metadata"])
    pd.testing.assert_frame_equal(metadata, small_dataset.metadata, check_dtype=False)
    truth = load_ground_truth(paths["ground_truth"])
    assert truth["is_signal"].sum() == small_dataset.ground_truth["is_signal"].sum()
    pd.testing.assert_series_equal(truth["true_shift"], small_dataset.ground_truth["true_shift"])


def test_files_carry_synthetic_label_and_lf_endings(
    tmp_path: Path, small_dataset: SyntheticDataset, small_config: SyntheticConfig
) -> None:
    paths = _write(tmp_path, small_dataset, small_config)
    for key in ("expression", "metadata", "ground_truth"):
        raw = paths[key].read_bytes()
        assert raw.startswith(f"# {DATA_LABEL}".encode())
        assert b"\r\n" not in raw


def test_manifest_records_license_seed_and_checksums(
    tmp_path: Path, small_dataset: SyntheticDataset, small_config: SyntheticConfig
) -> None:
    paths = _write(tmp_path, small_dataset, small_config)
    manifest = load_manifest(paths["manifest"])
    assert manifest["data_label"] == DATA_LABEL
    assert manifest["license"]["spdx"] == "CC0-1.0"
    assert "MIT" in manifest["license"]["scope"]
    assert manifest["seed"] == 123
    for key in ("expression", "metadata", "ground_truth"):
        entry = manifest["files"][paths[key].name]
        assert entry["sha256"] == file_sha256(paths[key])
        assert entry["bytes"] == paths[key].stat().st_size


def test_missing_values_round_trip_as_nan(tmp_path: Path) -> None:
    path = tmp_path / "expr.csv"
    path.write_text("sample_id,G1,G2\nS1,1.0,\nS2,,2.0\n", encoding="utf-8")
    expression = load_expression(path)
    assert expression.isna().sum().sum() == 2


def test_duplicate_gene_columns_are_preserved(tmp_path: Path) -> None:
    path = tmp_path / "expr.csv"
    path.write_text("# comment\nsample_id,G1,G1,G2\nS1,1,2,3\n", encoding="utf-8")
    expression = load_expression(path)
    assert list(expression.columns) == ["G1", "G1", "G2"]


def test_sample_ids_keep_leading_zeros(tmp_path: Path) -> None:
    path = tmp_path / "expr.csv"
    path.write_text("sample_id,G1\n001,1.0\n002,2.0\n", encoding="utf-8")
    assert list(load_expression(path).index) == ["001", "002"]


def test_metadata_treats_only_blank_as_missing(tmp_path: Path) -> None:
    path = tmp_path / "meta.csv"
    path.write_text("sample_id,group\nS1,NA\nS2,\n", encoding="utf-8")
    metadata = load_metadata(path)
    assert metadata.loc[0, "group"] == "NA"
    assert pd.isna(metadata.loc[1, "group"])


@pytest.mark.parametrize(
    "loader", [load_expression, load_metadata, load_ground_truth, load_manifest]
)
def test_missing_file_has_helpful_message(tmp_path: Path, loader: Callable[[Path], object]) -> None:
    with pytest.raises(DataLoadError, match="obw generate-data"):
        loader(tmp_path / "absent.csv")


def test_expression_without_gene_columns_rejected(tmp_path: Path) -> None:
    path = tmp_path / "expr.csv"
    path.write_text("sample_id\nS1\n", encoding="utf-8")
    with pytest.raises(DataLoadError, match="at least one gene column"):
        load_expression(path)


def test_invalid_manifest_json(tmp_path: Path) -> None:
    path = tmp_path / "manifest.json"
    path.write_text("{not json", encoding="utf-8")
    with pytest.raises(DataLoadError, match="not valid JSON"):
        load_manifest(path)
