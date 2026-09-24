"""Phase 1 smoke tests: package metadata, disclaimer text, CLI skeleton, repo layout."""

from __future__ import annotations

from pathlib import Path

import pytest

import onco_workbench
from onco_workbench import cli, disclaimers


def test_version_is_resolved_from_installed_metadata() -> None:
    assert onco_workbench.__version__ == "0.1.0"


def test_disclaimer_states_nonclinical_boundary() -> None:
    text = disclaimers.FULL_DISCLAIMER.lower()
    assert "not" in text and "medical device" in text
    assert "synthetic" in text
    assert "not for clinical use" in disclaimers.SHORT_DISCLAIMER.lower()


def test_comment_block_prefixes_every_line() -> None:
    block = disclaimers.as_comment_block("# ")
    lines = block.splitlines()
    assert len(lines) == 2
    assert all(line.startswith("# ") for line in lines)
    assert disclaimers.DATA_LABEL in block
    assert block.endswith("\n")


def test_cli_version(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        cli.main(["--version"])
    assert exc.value.code == 0
    assert "0.1.0" in capsys.readouterr().out


def test_cli_without_command_prints_help(capsys: pytest.CaptureFixture[str]) -> None:
    assert cli.main([]) == 0
    assert "research and education only" in capsys.readouterr().out.lower()


def test_cli_disclaimer(capsys: pytest.CaptureFixture[str]) -> None:
    assert cli.main(["disclaimer"]) == 0
    out = capsys.readouterr().out
    assert disclaimers.DATA_LABEL in out
    assert disclaimers.FULL_DISCLAIMER in out


@pytest.mark.parametrize(
    "relative_path",
    [
        "config/default.yaml",
        "data/raw",
        "data/synthetic",
        "data/processed",
        "outputs",
        "docs/plans/initial_implementation_plan.md",
        "CLAUDE.md",
        "LICENSE",
    ],
)
def test_expected_layout_exists(repo_root: Path, relative_path: str) -> None:
    assert (repo_root / relative_path).exists(), f"missing {relative_path}"
