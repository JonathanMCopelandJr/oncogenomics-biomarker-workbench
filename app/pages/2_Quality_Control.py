"""Quality Control: validation, PCA, expression distributions, missingness."""

import pandas as pd
import streamlit as st

from onco_workbench.dashboard.components import footer, get_qc, page_setup, require_demo_data
from onco_workbench.viz.interactive import expression_box_figure, missingness_figure, pca_figure

page_setup("Quality Control")
data = require_demo_data()
qc = get_qc()

report = data.validation
if report.warnings:
    st.warning(report.summary(), icon=":material/info:")
else:
    st.success(report.summary().splitlines()[0], icon=":material/check_circle:")

st.subheader("Dataset summary")
st.dataframe(
    pd.DataFrame(qc.summary.as_rows(), columns=["Measure", "Value"]),
    hide_index=True,
    width="stretch",
)

st.subheader("PCA")
options = ["group", "batch"] if "batch" in qc.pca.scores.columns else ["group"]
color_by = st.radio("Color samples by", options, horizontal=True)
st.plotly_chart(pca_figure(qc.pca, color_by), width="stretch")
st.caption(
    "For visual QC only: missing values are replaced with gene means and genes are "
    "scaled before PCA. In this synthetic dataset a simulated batch offset is visible. "
    "Batches are balanced within groups, so batch is not confounded with group."
)

st.subheader("Expression distributions")
st.plotly_chart(expression_box_figure(data.expression, data.aligned), width="stretch")

st.subheader("Missingness")
summary = qc.summary
cols = st.columns(4)
cols[0].metric("Missing values", f"{summary.n_missing:,}")
cols[1].metric("Missing fraction", f"{summary.missing_fraction:.2%}")
cols[2].metric("Genes with any missing", summary.genes_with_missing)
cols[3].metric("Samples with any missing", summary.samples_with_missing)
st.plotly_chart(missingness_figure(qc.sample_table), width="stretch")
st.markdown("**Genes with the most missing values**")
st.dataframe(
    qc.gene_table.sort_values(["n_missing", "gene_id"], ascending=[False, True]).head(10),
    hide_index=True,
    width="stretch",
)

with st.expander("Per-sample QC table"):
    st.dataframe(qc.sample_table, hide_index=True, width="stretch")

footer()
