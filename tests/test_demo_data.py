"""Tests for the committed synthetic demo files in data/synthetic/.

These guarantee that the committed files are small, labeled, CC0-documented,
valid, and exactly reproducible from the configured seed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from onco_workbench.config import load_config
from onco_workbench.data.io import (
    file_sha256,
    load_expression,
    load_manifest,
    load_metadata,
    write_dataset,
)
from onco_workbench.data.synthetic import generate_synthetic_dataset
from onco_workbench.data.validation import validate_dataset
from onco_workbench.disclaimers import DATA_LABEL

MAX_TOTAL_BYTES = 300 * 1024


@pytest.fixture(scope="module")
def committed_paths() -> dict[str, Path]:
    paths = load_config().paths
    result = {
        "expression": paths.expression_file,
        "metadata": paths.metadata_file,
        "ground_truth": paths.ground_truth_file,
        "manifest": paths.manifest_file,
    }
    missing = [str(p) for p in result.values() if not p.is_file()]
    if missing:
        pytest.fail(f"Committed demo files missing (run `obw generate-data`): {missing}")
    return result


def test_demo_files_are_small(committed_paths: dict[str, Path]) -> None:
    total = sum(path.stat().st_size for path in committed_paths.values())
    assert total <= MAX_TOTAL_BYTES, f"demo data is {total:,} bytes"


def test_demo_files_are_labeled_synthetic(committed_paths: dict[str, Path]) -> None:
    for key in ("expression", "metadata", "ground_truth"):
        first_line = committed_paths[key].read_text(encoding="utf-8").splitlines()[0]
        assert DATA_LABEL in first_line


def test_manifest_matches_committed_files(committed_paths: dict[str, Path]) -> None:
    manifest = load_manifest(committed_paths["manifest"])
    config = load_config()
    assert manifest["seed"] == config.seed
    assert manifest["license"]["spdx"] == "CC0-1.0"
    for key in ("expression", "metadata", "ground_truth"):
        path = committed_paths[key]
        assert manifest["files"][path.name]["sha256"] == file_sha256(path)


def test_regeneration_from_configured_seed_is_byte_identical(
    committed_paths: dict[str, Path], tmp_path: Path
) -> None:
    config = load_config()
    dataset = generate_synthetic_dataset(config.synthetic, config.seed)
    regenerated = write_dataset(
        dataset,
        config.synthetic,
        expression_path=tmp_path / committed_paths["expression"].name,
        metadata_path=tmp_path / committed_paths["metadata"].name,
        ground_truth_path=tmp_path / committed_paths["ground_truth"].name,
        manifest_path=tmp_path / committed_paths["manifest"].name,
    )
    for key, path in committed_paths.items():
        assert file_sha256(regenerated[key]) == file_sha256(path), (
            f"{path.name} differs from a fresh regeneration; run `obw generate-data` "
            "and review the change."
        )


def test_committed_demo_data_pass_validation(committed_paths: dict[str, Path]) -> None:
    config = load_config()
    report = validate_dataset(
        load_expression(committed_paths["expression"]),
        load_metadata(committed_paths["metadata"]),
        allowed_groups=config.validation.allowed_groups,
        min_samples_per_group=config.validation.min_samples_per_group,
        max_missing_fraction_per_gene=config.validation.max_missing_fraction_per_gene,
        max_missing_fraction_per_sample=config.validation.max_missing_fraction_per_sample,
    )
    assert report.is_valid, report.summary()
