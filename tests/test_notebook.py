"""Checks for the walkthrough notebook (no Jupyter dependency needed).

Structure tests read the ``.ipynb`` JSON directly. The ``e2e`` test executes every code
cell in order in one namespace (IPython magics stripped), which catches API drift
between the notebook and the package without requiring a Jupyter kernel in CI.
"""

from __future__ import annotations

import ast
import json
from pathlib import Path
from typing import Any

import pytest

NOTEBOOK = (
    Path(__file__).resolve().parents[1] / "notebooks" / "oncogenomics_workbench_walkthrough.ipynb"
)
REQUIRED_DISCLAIMER = (
    "This notebook uses fully synthetic transcriptomics-style data for educational and "
    "research demonstration purposes only. It is not intended for diagnosis, prognosis, "
    "treatment selection, clinical decision-making, cancer prediction, biological discovery, "
    "or biomarker validation."
)
ALLOWED_IMPORT_ROOTS = {"onco_workbench", "pandas", "numpy"}
REQUIRED_SECTIONS = (
    "## 1. Project overview",
    "## 2. Synthetic data",
    "## 3. Dataset validation summary",
    "## 4. Compact quality-control summary",
    "## 5. PCA visualization",
    "## 6. Group comparison",
    "## 7. Ranked simulated candidate genes",
    "## 8. Volcano-style plot",
    "## 9. Top-ranked genes across samples",
    "## 10. Optional ML demonstration",
    "## 11. Reproducibility, limitations",
)


@pytest.fixture(scope="module")
def notebook() -> dict[str, Any]:
    return json.loads(NOTEBOOK.read_text(encoding="utf-8"))


def _source(cell: dict[str, Any]) -> str:
    return "".join(cell["source"])


def _code_without_magics(cell: dict[str, Any]) -> str:
    return "\n".join(
        line for line in _source(cell).splitlines() if not line.lstrip().startswith("%")
    )


def _code_cells(notebook: dict[str, Any]) -> list[dict[str, Any]]:
    return [c for c in notebook["cells"] if c["cell_type"] == "code"]


def _markdown(notebook: dict[str, Any]) -> str:
    return "\n".join(_source(c) for c in notebook["cells"] if c["cell_type"] == "markdown")


def test_notebook_format(notebook: dict[str, Any]) -> None:
    assert notebook["nbformat"] == 4
    assert notebook["metadata"]["kernelspec"]["name"] == "python3"
    ids = [cell["id"] for cell in notebook["cells"]]
    assert len(ids) == len(set(ids))


def test_first_cell_is_the_required_disclaimer(notebook: dict[str, Any]) -> None:
    first = notebook["cells"][0]
    assert first["cell_type"] == "markdown"
    assert REQUIRED_DISCLAIMER in _source(first)
    assert "NOT FOR CLINICAL USE" in _source(first)


def test_required_sections_present_in_order(notebook: dict[str, Any]) -> None:
    text = _markdown(notebook)
    positions = [text.find(section) for section in REQUIRED_SECTIONS]
    assert all(p >= 0 for p in positions), dict(zip(REQUIRED_SECTIONS, positions, strict=True))
    assert positions == sorted(positions)


def test_results_are_framed_as_implanted_synthetic_signals(notebook: dict[str, Any]) -> None:
    text = " ".join(_markdown(notebook).split())  # ignore Markdown line wrapping
    assert "intentionally implanted synthetic" in text
    assert "sanity check, not a biological benchmark" in text
    assert "NOT a clinical prediction model" in text
    assert "(../README.md)" in text


def test_outputs_are_stripped(notebook: dict[str, Any]) -> None:
    for cell in _code_cells(notebook):
        assert cell["outputs"] == []
        assert cell["execution_count"] is None


def test_no_core_logic_is_duplicated_in_cells(notebook: dict[str, Any]) -> None:
    for cell in _code_cells(notebook):
        tree = ast.parse(_code_without_magics(cell))
        for node in ast.walk(tree):
            assert not isinstance(node, ast.FunctionDef | ast.ClassDef | ast.Lambda), (
                "notebook cells must call package functions, not define logic"
            )
            if isinstance(node, ast.Import):
                roots = {alias.name.split(".")[0] for alias in node.names}
            elif isinstance(node, ast.ImportFrom):
                roots = {(node.module or "").split(".")[0]}
            else:
                continue
            assert roots <= ALLOWED_IMPORT_ROOTS, f"unexpected import: {roots}"


def test_no_file_writes_or_network_in_cells(notebook: dict[str, Any]) -> None:
    code = "\n".join(_source(c) for c in _code_cells(notebook))
    for forbidden in ("to_csv(", "savefig(", "write_text(", "open(", "urlopen", "requests", "http"):
        assert forbidden not in code, forbidden


@pytest.mark.e2e
def test_notebook_code_cells_execute_in_order(notebook: dict[str, Any]) -> None:
    namespace: dict[str, Any] = {"__name__": "__notebook__"}
    for number, cell in enumerate(_code_cells(notebook), start=1):
        code = compile(_code_without_magics(cell), f"<notebook code cell {number}>", "exec")
        exec(code, namespace)  # executes the repository's own notebook cells
    assert namespace["same_ids"] and namespace["same_values"]
    assert len(namespace["results"]) == 500
    assert namespace["comparison"].recovery.n_planted == 40
    assert namespace["ml"].labels.positive == "Group_B"
