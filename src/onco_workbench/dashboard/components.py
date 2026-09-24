"""Streamlit building blocks shared by every dashboard page.

Every page calls :func:`page_setup` first. It renders the research-only /
synthetic-data banner above the page title and a matching notice in the sidebar.
Every page ends with :func:`footer`, which repeats the full disclaimer.

Data access goes through cached wrappers around
:mod:`onco_workbench.dashboard.data`, which in turn call the same package functions
as the CLI pipeline.
"""

from __future__ import annotations

import streamlit as st

from onco_workbench.config import ComparisonConfig, ConfigError
from onco_workbench.dashboard.data import DemoData, compare, load_demo_data
from onco_workbench.data.io import DataLoadError
from onco_workbench.data.validation import DataValidationError
from onco_workbench.disclaimers import DATA_LABEL, FULL_DISCLAIMER, SHORT_DISCLAIMER
from onco_workbench.pipeline import ComparisonResult, QCResult, run_qc

APP_NAME = "Oncogenomics Biomarker Workbench"


def disclaimer_banner() -> None:
    """Render the prominent research-only / synthetic-data warning."""
    st.warning(
        f"**{DATA_LABEL} - RESEARCH AND EDUCATION ONLY - NOT FOR CLINICAL USE.**  \n"
        "This dashboard is not a medical device, diagnostic tool, or clinical "
        "decision-support system. Every value shown is simulated. Synthetic gene IDs do "
        "not correspond to real genes, and nothing here identifies or validates a biomarker.",
        icon=":material/warning:",
    )


def page_setup(title: str) -> None:
    """Configure the page, then show the banner, sidebar notice, and title.

    Args:
        title: Page title.
    """
    st.set_page_config(
        page_title=f"{title} - {APP_NAME} (synthetic demo)",
        page_icon=":material/science:",
        layout="wide",
    )
    disclaimer_banner()
    st.sidebar.markdown(f"**{DATA_LABEL}**")
    st.sidebar.caption(SHORT_DISCLAIMER)
    st.title(title)


def footer() -> None:
    """Repeat the full disclaimer at the bottom of the page."""
    st.divider()
    st.caption(f"**{DATA_LABEL}.** {FULL_DISCLAIMER}")


@st.cache_resource(show_spinner="Loading synthetic demo data...")
def get_demo_data() -> DemoData:
    """Load and validate the synthetic demo data once per server process."""
    return load_demo_data()


def require_demo_data() -> DemoData:
    """Return the demo data, or show an actionable error and stop the page.

    Returns:
        The loaded :class:`DemoData`.
    """
    try:
        return get_demo_data()
    except (ConfigError, DataLoadError, DataValidationError) as exc:
        st.error(
            "The synthetic demo data could not be loaded. Regenerate it with "
            "`obw generate-data`, then reload this page.\n\n"
            f"Details: {exc}",
            icon=":material/error:",
        )
        footer()
        st.stop()
        raise  # unreachable: st.stop() ends the script run


@st.cache_data(show_spinner="Computing QC summaries...")
def get_qc() -> QCResult:
    """QC summaries, PCA, and correlation for the demo data (cached)."""
    data = get_demo_data()
    return run_qc(data.expression, data.metadata, data.config)


@st.cache_data(show_spinner="Comparing groups...")
def get_comparison(
    group_a: str, group_b: str, fdr_threshold: float, effect_size_threshold: float
) -> tuple[ComparisonConfig, ComparisonResult]:
    """Group comparison for the chosen groups and thresholds (cached per combination)."""
    return compare(
        get_demo_data(),
        group_a=group_a,
        group_b=group_b,
        fdr_threshold=fdr_threshold,
        effect_size_threshold=effect_size_threshold,
    )
