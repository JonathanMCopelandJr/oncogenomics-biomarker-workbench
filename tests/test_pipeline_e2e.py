"""End-to-end tests of the analysis pipeline on SYNTHETIC data only."""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
import pytest

from onco_workbench.config import SyntheticConfig, load_config
from onco_workbench.data.io import file_sha256, load_manifest, write_dataset
from onco_workbench.data.synthetic import SyntheticDataset
from onco_workbench.data.validation import DataValidationError
from onco_workbench.disclaimers import DATA_LABEL
from onco_workbench.pipeline import (
    FIGURES_DIR,
    MANIFEST_NAME,
    REPORT_NAME,
    RESULT_COLUMN_ORDER,
    TABLES_DIR,
    PipelineRun,
    run_pipeline,
)

FIXED_TIME = datetime(2026, 1, 1, tzinfo=UTC)


@pytest.fixture
def small_files(
    tmp_path: Path, small_dataset: SyntheticDataset, small_config: SyntheticConfig
) -> dict[str, Path]:
    return write_dataset(
        small_dataset,
        small_config,
        expression_path=tmp_path / "in" / "expr.csv",
        metadata_path=tmp_path / "in" / "meta.csv",
        ground_truth_path=tmp_path / "in" / "truth.csv",
        manifest_path=tmp_path / "in" / "manifest.json",
    )


def _run(files: dict[str, Path], out: Path, **kwargs: object) -> PipelineRun:
    return run_pipeline(
        load_config(),
        expression_path=files["expression"],
        metadata_path=files["metadata"],
        output_dir=out,
        generated_at=FIXED_TIME,
        **kwargs,  # type: ignore[arg-type]
    )


def test_full_run_writes_labeled_outputs(small_files: dict[str, Path], tmp_path: Path) -> None:
    out = tmp_path / "out"
    run = _run(small_files, out, ground_truth_path=small_files["ground_truth"])

    expected_tables = {
        "qc_sample_summary",
        "qc_gene_summary",
        "pca_scores",
        "group_comparison_results",
        "top_ranked_genes",
    }
    assert {p.stem for p in (out / TABLES_DIR).iterdir()} == expected_tables
    assert len(list((out / FIGURES_DIR).glob("*.png"))) == 7
    for csv in (out / TABLES_DIR).glob("*.csv"):
        assert csv.read_text(encoding="utf-8").startswith(f"# {DATA_LABEL}")

    report = (out / REPORT_NAME).read_text(encoding="utf-8")
    assert DATA_LABEL in report and "## 5. Workflow check" in report
    for figure in (out / FIGURES_DIR).iterdir():
        assert f"{FIGURES_DIR}/{figure.name}" in report

    results = pd.read_csv(out / TABLES_DIR / "group_comparison_results.csv", comment="#")
    assert tuple(results.columns) == RESULT_COLUMN_ORDER
    assert len(results) == 40
    assert run.comparison is not None and run.comparison.recovery is not None
    assert run.comparison.recovery.n_planted == 8


def test_manifest_matches_written_files(small_files: dict[str, Path], tmp_path: Path) -> None:
    out = tmp_path / "out"
    run = _run(small_files, out)
    manifest = load_manifest(out / MANIFEST_NAME)
    assert manifest["data_label"] == DATA_LABEL
    assert manifest["generated_at"] == "2026-01-01T00:00:00+00:00"
    assert manifest["inputs"]["expression"]["sha256"] == file_sha256(small_files["expression"])
    for role, entry in manifest["outputs"].items():
        assert entry["sha256"] == file_sha256(run.files[role])
    assert manifest["parameters"]["comparison"]["group_b"] == "Group_B"


def test_custom_inputs_skip_ground_truth_unless_given(
    small_files: dict[str, Path], tmp_path: Path
) -> None:
    run = _run(small_files, tmp_path / "out")
    assert run.comparison is not None and run.comparison.recovery is None
    report = (tmp_path / "out" / REPORT_NAME).read_text(encoding="utf-8")
    assert "## 5." not in report


def test_qc_only_run(small_files: dict[str, Path], tmp_path: Path) -> None:
    out = tmp_path / "qc"
    run = _run(small_files, out, include_comparison=False)
    assert run.comparison is None
    assert not (out / TABLES_DIR / "group_comparison_results.csv").exists()
    assert not (out / FIGURES_DIR / "volcano.png").exists()
    assert len(list((out / FIGURES_DIR).glob("*.png"))) == 4


def test_invalid_inputs_raise_before_writing(small_files: dict[str, Path], tmp_path: Path) -> None:
    metadata = pd.read_csv(small_files["metadata"], comment="#")
    metadata.loc[0, "group"] = "Group_Z"
    metadata.to_csv(small_files["metadata"], index=False)
    out = tmp_path / "out"
    with pytest.raises(DataValidationError, match="invalid_group"):
        _run(small_files, out)
    assert not out.exists()


def test_unknown_comparison_group_raises(small_files: dict[str, Path], tmp_path: Path) -> None:
    config = load_config().with_comparison(group_b="Group_Q")
    with pytest.raises(ValueError, match="Group_Q"):
        run_pipeline(
            config,
            expression_path=small_files["expression"],
            metadata_path=small_files["metadata"],
            output_dir=tmp_path / "out",
        )


@pytest.mark.e2e
def test_demo_data_end_to_end(tmp_path: Path) -> None:
    """Full run on the committed synthetic demo data, including the ground-truth check.

    With 40 samples per group and planted shifts of 0.8-2.0 against noise SD 0.5,
    every planted gene has a large standardized effect, so near-complete recovery is
    expected. BH at 5% bounds the expected share of spurious flags.
    """
    run = run_pipeline(load_config(), output_dir=tmp_path / "demo", generated_at=FIXED_TIME)
    recovery = run.comparison.recovery  # type: ignore[union-attr]
    assert recovery is not None
    assert recovery.n_planted == 40
    assert recovery.true_positives >= 38
    assert recovery.false_discovery_proportion <= 0.10
    assert recovery.direction_matches == recovery.true_positives
    assert run.qc.summary.n_samples == 80 and run.qc.summary.n_genes == 500
