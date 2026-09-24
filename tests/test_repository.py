"""Repository hygiene checks: required files, CI and citation metadata, disclaimers."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
import yaml

REQUIRED_FILES = [
    "README.md",
    "LICENSE",
    "CITATION.cff",
    "CONTRIBUTING.md",
    "CODE_OF_CONDUCT.md",
    "SECURITY.md",
    "CHANGELOG.md",
    "CLAUDE.md",
    "docs/methodology.md",
    "docs/data_dictionary.md",
    "docs/architecture.md",
    "docs/reproducibility.md",
    "docs/limitations_and_ethics.md",
    "docs/portfolio_talking_points.md",
    "docs/plans/initial_implementation_plan.md",
    "docs/release_checklist.md",
    ".github/workflows/ci.yml",
    ".github/pull_request_template.md",
    ".github/ISSUE_TEMPLATE/bug_report.md",
    ".github/ISSUE_TEMPLATE/feature_request.md",
    "data/synthetic/README.md",
    "data/raw/README.md",
]


@pytest.mark.parametrize("relative_path", REQUIRED_FILES)
def test_required_file_exists_and_is_not_empty(repo_root: Path, relative_path: str) -> None:
    path = repo_root / relative_path
    assert path.is_file(), f"missing {relative_path}"
    assert path.stat().st_size > 0, f"empty {relative_path}"


def test_citation_metadata(repo_root: Path) -> None:
    citation = yaml.safe_load((repo_root / "CITATION.cff").read_text(encoding="utf-8"))
    for key in ("cff-version", "message", "title", "authors", "license", "version"):
        assert key in citation
    assert citation["license"] == "MIT"
    assert citation["authors"] == [{"family-names": "Copeland", "given-names": "Jonathan"}]
    assert "orcid" not in str(citation).lower()  # only added if a verified iD is supplied
    assert "synthetic" in citation["message"] and "clinical" in citation["message"]


def test_citation_version_matches_package(repo_root: Path) -> None:
    import onco_workbench

    citation = yaml.safe_load((repo_root / "CITATION.cff").read_text(encoding="utf-8"))
    assert str(citation["version"]) == onco_workbench.__version__


def test_ci_workflow(repo_root: Path) -> None:
    text = (repo_root / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    workflow = yaml.safe_load(text)
    triggers = workflow[True] if True in workflow else workflow["on"]  # YAML 1.1 "on" -> True
    assert {"push", "pull_request"} <= set(triggers)
    assert workflow["permissions"] == {"contents": "read"}
    matrix = workflow["jobs"]["test"]["strategy"]["matrix"]["include"]
    assert {"os": "ubuntu-latest", "python-version": "3.11"} in matrix
    assert {m["os"] for m in matrix} == {"ubuntu-latest", "windows-latest", "macos-latest"}
    assert "secrets." not in text, "CI must not need secrets"
    # Only executable content matters; comments may mention what CI does *not* do.
    code = "\n".join(line.split("#", 1)[0] for line in text.lower().splitlines())
    for forbidden in ("deploy", "publish", "gh release", "twine", "upload"):
        assert forbidden not in code


@pytest.mark.parametrize(
    "page", sorted(p.name for p in (Path(__file__).parents[1] / "app" / "pages").glob("*.py"))
)
def test_every_dashboard_page_has_banner_and_footer(repo_root: Path, page: str) -> None:
    source = (repo_root / "app" / "pages" / page).read_text(encoding="utf-8")
    assert "page_setup(" in source and "footer()" in source


@pytest.mark.parametrize(
    "document",
    ["README.md", "docs/limitations_and_ethics.md", "data/synthetic/README.md", "SECURITY.md"],
)
def test_documents_state_nonclinical_boundary(repo_root: Path, document: str) -> None:
    raw = (repo_root / document).read_text(encoding="utf-8").lower()
    # Ignore Markdown emphasis, block-quote markers, and line wrapping.
    text = " ".join(raw.replace("*", "").replace(">", " ").split())
    assert "synthetic" in text
    assert (
        "not for clinical use" in text or "not patient data" in text or "real patient data" in text
    )


def test_data_dictionary_documents_every_validation_code(repo_root: Path) -> None:
    source = (repo_root / "src" / "onco_workbench" / "data" / "validation.py").read_text(
        encoding="utf-8"
    )
    # Codes follow a severity argument in report.add(...); drop the Literal type's members.
    codes = set(re.findall(r'"(?:error|warning)",\s*"([a-z_]+)"', source)) - {"error", "warning"}
    codes |= {
        f"{kind}_{suffix}"
        for kind in ("gene", "sample")
        for suffix in ("all_missing", "high_missingness")
    }
    dictionary = (repo_root / "docs" / "data_dictionary.md").read_text(encoding="utf-8")
    assert len(codes) == 23
    missing = sorted(code for code in codes if f"`{code}`" not in dictionary)
    assert not missing, f"undocumented validation codes: {missing}"


def test_gitignore_protects_data_and_secrets(repo_root: Path) -> None:
    rules = (repo_root / ".gitignore").read_text(encoding="utf-8").splitlines()
    for rule in (".venv/", ".env", ".streamlit/secrets.toml", "outputs/*", "data/raw/*"):
        assert rule in rules
