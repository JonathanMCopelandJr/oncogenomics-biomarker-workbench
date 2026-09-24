from __future__ import annotations

from pathlib import Path

import pytest

from onco_workbench import cli
from onco_workbench.data.io import file_sha256
from onco_workbench.disclaimers import DATA_LABEL


def test_generate_data_to_output_dir(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    assert cli.main(["generate-data", "--output-dir", str(tmp_path)]) == cli.EXIT_OK
    out = capsys.readouterr().out
    assert DATA_LABEL in out
    assert "80 samples x 500 genes" in out
    for name in ("demo_expression.csv", "demo_metadata.csv", "demo_ground_truth.csv"):
        assert (tmp_path / name).is_file()
    assert (tmp_path / "demo_manifest.json").is_file()


def test_generate_data_is_deterministic(tmp_path: Path) -> None:
    first, second = tmp_path / "a", tmp_path / "b"
    cli.main(["generate-data", "--output-dir", str(first)])
    cli.main(["generate-data", "--output-dir", str(second)])
    for path in first.iterdir():
        assert file_sha256(path) == file_sha256(second / path.name)


def test_generate_data_seed_override_changes_output(tmp_path: Path) -> None:
    cli.main(["generate-data", "--output-dir", str(tmp_path / "a")])
    cli.main(["generate-data", "--output-dir", str(tmp_path / "b"), "--seed", "7"])
    name = "demo_expression.csv"
    assert file_sha256(tmp_path / "a" / name) != file_sha256(tmp_path / "b" / name)


def test_validate_generated_files(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    cli.main(["generate-data", "--output-dir", str(tmp_path)])
    code = cli.main(
        [
            "validate",
            "--expression",
            str(tmp_path / "demo_expression.csv"),
            "--metadata",
            str(tmp_path / "demo_metadata.csv"),
        ]
    )
    assert code == cli.EXIT_OK
    assert "Validation PASSED" in capsys.readouterr().out


def test_validate_reports_failures(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    expression = tmp_path / "expr.csv"
    metadata = tmp_path / "meta.csv"
    expression.write_text("sample_id,G1,G1\nS1,1,2\nS2,3,x\n", encoding="utf-8")
    metadata.write_text("sample_id,group\nS1,Group_A\nS2,Group_Z\n", encoding="utf-8")
    code = cli.main(["validate", "--expression", str(expression), "--metadata", str(metadata)])
    assert code == cli.EXIT_FAILURE
    out = capsys.readouterr().out
    assert "Validation FAILED" in out
    for expected in ("duplicate_gene_id", "non_numeric_values", "invalid_group"):
        assert expected in out


def test_validate_missing_file(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    code = cli.main(["validate", "--expression", str(tmp_path / "none.csv")])
    assert code == cli.EXIT_FAILURE
    assert "not found" in capsys.readouterr().err


def test_bad_config_path(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    code = cli.main(["generate-data", "--config", str(tmp_path / "none.yaml")])
    assert code == cli.EXIT_FAILURE
    assert "Configuration file not found" in capsys.readouterr().err


@pytest.mark.e2e
def test_run_analysis_on_demo_data(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    assert cli.main(["run-analysis", "--output-dir", str(tmp_path)]) == cli.EXIT_OK
    out = capsys.readouterr().out
    assert DATA_LABEL in out
    assert "Comparison Group_B vs Group_A: 500 genes tested" in out
    assert "Synthetic ground-truth check" in out
    assert (tmp_path / "report.md").is_file()
    assert (tmp_path / "run_manifest.json").is_file()
    assert (tmp_path / "tables" / "group_comparison_results.csv").is_file()


def test_qc_command(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    assert cli.main(["qc", "--output-dir", str(tmp_path)]) == cli.EXIT_OK
    out = capsys.readouterr().out
    assert "QC: 80 samples x 500 genes" in out
    assert "Comparison" not in out
    assert not (tmp_path / "tables" / "group_comparison_results.csv").exists()


def test_run_analysis_unknown_group(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    code = cli.main(["run-analysis", "--output-dir", str(tmp_path), "--group-b", "Group_Q"])
    assert code == cli.EXIT_FAILURE
    assert "Group_Q" in capsys.readouterr().err


def test_run_analysis_invalid_threshold(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    code = cli.main(["run-analysis", "--output-dir", str(tmp_path), "--fdr-threshold", "2"])
    assert code == cli.EXIT_FAILURE
    assert "fdr_threshold" in capsys.readouterr().err


def test_run_analysis_invalid_inputs(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    expression = tmp_path / "expr.csv"
    metadata = tmp_path / "meta.csv"
    expression.write_text("sample_id,G1\nS1,1\nS2,2\n", encoding="utf-8")
    metadata.write_text("sample_id,group\nS1,Group_A\nS2,Group_B\n", encoding="utf-8")
    code = cli.main(
        [
            "run-analysis",
            "--expression",
            str(expression),
            "--metadata",
            str(metadata),
            "--output-dir",
            str(tmp_path / "out"),
        ]
    )
    assert code == cli.EXIT_FAILURE
    assert "group_too_small" in capsys.readouterr().err
    assert not (tmp_path / "out").exists()
