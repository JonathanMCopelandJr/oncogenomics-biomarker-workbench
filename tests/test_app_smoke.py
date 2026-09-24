"""Smoke tests for every Streamlit page, using Streamlit's built-in AppTest (no browser).

Each page must run without exceptions and show the research-only / synthetic-data
banner and the footer disclaimer.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from onco_workbench.disclaimers import DATA_LABEL

APP_DIR = Path(__file__).resolve().parents[1] / "app"
PAGES = [
    "Home.py",
    "pages/1_Data_Explorer.py",
    "pages/2_Quality_Control.py",
    "pages/3_Group_Comparison.py",
    "pages/4_Model_Demonstration.py",
    "pages/5_Methods_and_Limitations.py",
]
TIMEOUT = 60


def _run(page: str) -> AppTest:
    app = AppTest.from_file(str(APP_DIR / page), default_timeout=TIMEOUT)
    app.run()
    return app


def _all_text(app: AppTest) -> str:
    parts = [e.value for e in app.warning] + [e.value for e in app.error]
    parts += [e.value for e in app.markdown] + [e.value for e in app.caption]
    parts += [e.value for e in app.info]
    return "\n".join(str(p) for p in parts)


@pytest.mark.parametrize("page", PAGES)
def test_page_runs_and_shows_disclaimers(page: str) -> None:
    app = _run(page)
    assert not app.exception, [e.value for e in app.exception]
    banner = [w.value for w in app.warning]
    assert any("RESEARCH AND EDUCATION ONLY - NOT FOR CLINICAL USE" in w for w in banner)
    assert any(DATA_LABEL in w for w in banner)
    assert any("not for clinical use" in c.value.lower() for c in app.caption)  # footer


def test_group_comparison_defaults_and_download() -> None:
    app = _run("pages/3_Group_Comparison.py")
    assert [s.value for s in app.selectbox][:2] == ["Group_A", "Group_B"]
    metrics = {m.label: m.value for m in app.metric}
    assert metrics["Genes tested"] == "500"
    assert "not a biomarker" in _all_text(app)


def test_group_comparison_same_group_shows_error() -> None:
    app = _run("pages/3_Group_Comparison.py")
    app.selectbox[1].set_value("Group_A").run()
    assert not app.exception
    assert any("two different groups" in e.value for e in app.error)


def test_group_comparison_stricter_threshold_flags_fewer() -> None:
    app = _run("pages/3_Group_Comparison.py")
    before = int({m.label: m.value for m in app.metric}["Meeting both thresholds"])
    app.number_input[0].set_value(3.0).run()
    after = int({m.label: m.value for m in app.metric}["Meeting both thresholds"])
    assert not app.exception
    assert after <= before


def test_group_comparison_filter_only_flagged() -> None:
    app = _run("pages/3_Group_Comparison.py")
    app.checkbox[0].check().run()
    assert not app.exception
    assert "genes shown" in _all_text(app)


def test_quality_control_batch_coloring() -> None:
    app = _run("pages/2_Quality_Control.py")
    app.radio[0].set_value("batch").run()
    assert not app.exception


def test_data_explorer_search_no_match() -> None:
    app = _run("pages/1_Data_Explorer.py")
    app.text_input[0].input("does-not-exist").run()
    assert not app.exception
    assert any("No synthetic gene IDs match" in i.value for i in app.info)


def test_model_page_makes_no_metric_claims() -> None:
    app = _run("pages/4_Model_Demonstration.py")
    text = _all_text(app)
    assert "NOT a clinical prediction model" in text
    assert "Not yet implemented" in text
    assert len(app.metric) == 0
