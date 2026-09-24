from __future__ import annotations

import json
from pathlib import Path

import pytest

from onco_workbench.analysis.qc import compute_qc_summary
from onco_workbench.config import load_config
from onco_workbench.data.io import file_sha256
from onco_workbench.data.synthetic import SyntheticDataset
from onco_workbench.data.validation import validate_dataset
from onco_workbench.disclaimers import DATA_LABEL, FULL_DISCLAIMER
from onco_workbench.pipeline import run_comparison
from onco_workbench.reporting.manifest import (
    build_run_manifest,
    display_path,
    git_state,
    write_run_manifest,
)
from onco_workbench.reporting.markdown_report import (
    ComparisonSection,
    ReportContext,
    md_table,
    render_markdown_report,
)


def _context(dataset: SyntheticDataset, *, with_comparison: bool) -> ReportContext:
    config = load_config()
    comparison = None
    if with_comparison:
        result = run_comparison(dataset.expression, dataset.metadata, config, dataset.ground_truth)
        comparison = ComparisonSection(
            config.comparison, config.ranking, result.results, result.recovery
        )
    return ReportContext(
        generated_at="2026-01-01T00:00:00+00:00",
        config_path="config/default.yaml",
        seed=42,
        inputs={"expression": "expr.csv", "metadata": "meta.csv"},
        validation=validate_dataset(dataset.expression, dataset.metadata),
        qc=compute_qc_summary(dataset.expression, dataset.metadata),
        pca_variance=[0.2, 0.1],
        normalization_steps=["step one"],
        figures={"volcano": "figures/volcano.png", "pca": "figures/qc_pca.png"},
        tables={"group_comparison_results": "tables/group_comparison_results.csv"},
        comparison=comparison,
    )


def test_md_table_escapes_pipes() -> None:
    table = md_table(["a", "b"], [("x|y", 1)])
    assert table.splitlines()[0] == "| a | b |"
    assert "x\\|y" in table


def test_report_is_labeled_and_complete(small_dataset: SyntheticDataset) -> None:
    text = render_markdown_report(_context(small_dataset, with_comparison=True))
    assert text.startswith(f"# Analysis report - {DATA_LABEL}")
    assert FULL_DISCLAIMER in text
    assert text.rstrip().endswith("not for clinical use.*")
    for heading in (
        "## 1. Run information",
        "## 2. Data quality control",
        "## 3. Normalization",
        "## 4. Group comparison: Group_B vs Group_A",
        "## 5. Workflow check against the synthetic ground truth",
        "## Interpretation limits",
    ):
        assert heading in text
    assert "![Volcano-style plot](figures/volcano.png)" in text
    assert "No output of this workbench identifies or validates a biomarker." in text


def test_qc_only_report_has_no_comparison(small_dataset: SyntheticDataset) -> None:
    text = render_markdown_report(_context(small_dataset, with_comparison=False))
    assert "## 4." not in text and "## 5." not in text
    assert "## 2. Data quality control" in text


def test_display_path(tmp_path: Path) -> None:
    inside = tmp_path / "a" / "b.csv"
    assert display_path(inside, tmp_path) == "a/b.csv"
    outside = tmp_path.parent / "elsewhere" / "secret_dir" / "c.csv"
    assert display_path(outside, tmp_path) == "c.csv"


def test_git_state_outside_repository(tmp_path: Path) -> None:
    state = git_state(tmp_path)
    assert set(state) == {"commit", "dirty"}
    if state["commit"] is None:
        assert state["dirty"] is None


def test_manifest_checksums_and_labels(tmp_path: Path, repo_root: Path) -> None:
    config_file = repo_root / "config" / "default.yaml"
    output = tmp_path / "out.csv"
    output.write_text("x\n1\n", encoding="utf-8")
    manifest = build_run_manifest(
        root=repo_root,
        config_path=config_file,
        seed=42,
        parameters={"k": 1},
        inputs={"config_copy": config_file},
        outputs={"table": output},
        generated_at="2026-01-01T00:00:00+00:00",
    )
    assert manifest["data_label"] == DATA_LABEL
    assert manifest["config"]["path"] == "config/default.yaml"
    assert manifest["outputs"]["table"]["sha256"] == file_sha256(output)
    assert manifest["outputs"]["table"]["path"] == "out.csv"  # outside root: name only
    assert "numpy" in manifest["packages"]
    path = write_run_manifest(manifest, tmp_path / "m" / "run_manifest.json")
    assert json.loads(path.read_text(encoding="utf-8"))["seed"] == 42
    assert b"\r\n" not in path.read_bytes()


@pytest.mark.parametrize("with_comparison", [True, False])
def test_report_has_no_crlf(small_dataset: SyntheticDataset, with_comparison: bool) -> None:
    assert "\r" not in render_markdown_report(
        _context(small_dataset, with_comparison=with_comparison)
    )
