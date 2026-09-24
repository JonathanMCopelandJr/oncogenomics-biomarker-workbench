"""Dashboard home page: overview, disclaimer, quick start.

Run from the repository root with ``obw dashboard`` (or ``streamlit run app/Home.py``).
"""

import streamlit as st

from onco_workbench.dashboard.components import APP_NAME, footer, page_setup, require_demo_data
from onco_workbench.disclaimers import FULL_DISCLAIMER

page_setup(APP_NAME)

st.error(f"**Important:** {FULL_DISCLAIMER}", icon=":material/gpp_maybe:")

st.markdown(
    """
This workbench shows, step by step, how a transcriptomics-style dataset can be
**checked, summarized, and compared between two groups of samples**. It is built for
learning and for demonstrating transparent, reproducible analysis practice.

The data are **simulated with a known, planted answer**: a small set of synthetic genes
was deliberately made to differ between `Group_A` and `Group_B`. A correct workflow
should recover them. This makes every step checkable. It also means the results say
nothing about real biology.
"""
)

data = require_demo_data()
truth = data.ground_truth
cols = st.columns(4)
cols[0].metric("Samples", f"{data.expression.shape[0]}")
cols[1].metric("Synthetic genes", f"{data.expression.shape[1]}")
cols[2].metric("Groups", ", ".join(data.groups))
cols[3].metric(
    "Planted signal genes", "n/a" if truth is None else f"{int(truth['is_signal'].sum())}"
)

st.subheader("Pages")
st.markdown(
    """
- **Data Explorer**: dataset shape, sample metadata, gene search, summary statistics.
- **Quality Control**: validation results, PCA, expression distributions, missingness.
- **Group Comparison**: choose two groups and display thresholds, then explore the
  volcano-style plot, filter and sort the results table, and download a labeled CSV.
- **Model Demonstration**: planned educational example; not yet implemented.
- **Methods and Limitations**: plain-language methods, limitations, data provenance,
  and reproducibility details.
"""
)

st.subheader("Quick start (command line)")
st.code(
    """python -m pip install -r requirements-dev.txt -e .   # install (inside a virtual env)
obw validate                                        # check the synthetic demo data
obw run-analysis                                    # tables, figures, report -> outputs/
obw dashboard                                       # this dashboard at http://localhost:8501""",
    language="bash",
)

footer()
