"""Data Explorer: dataset shape, metadata, gene search, summary statistics."""

import streamlit as st

from onco_workbench.analysis.qc import gene_qc_table
from onco_workbench.dashboard.components import footer, page_setup, require_demo_data
from onco_workbench.dashboard.data import gene_group_summary, search_genes, summary_statistics
from onco_workbench.viz.interactive import gene_by_group_figure

page_setup("Data Explorer")
data = require_demo_data()
expression, metadata = data.expression, data.metadata

st.caption(
    "All values are simulated, continuous, log2-like numbers. Gene IDs such as "
    "`SYN_G0001` are synthetic and do not correspond to real genes."
)

cols = st.columns(4)
cols[0].metric("Samples (rows)", expression.shape[0])
cols[1].metric("Synthetic genes (columns)", expression.shape[1])
cols[2].metric("Groups", len(data.groups))
cols[3].metric("Batches", metadata["batch"].nunique() if "batch" in metadata else 0)

st.subheader("Sample metadata")
st.dataframe(metadata, hide_index=True, width="stretch", height=280)

st.subheader("Gene search")
query = st.text_input("Search synthetic gene IDs", placeholder="for example SYN_G004")
matches = search_genes(list(expression.columns), query)
if not matches:
    st.info("No synthetic gene IDs match that search.")
else:
    st.caption(f"{len(matches)} match(es) shown (limit 200).")
    table = gene_qc_table(expression)
    st.dataframe(
        table[table["gene_id"].isin(matches)],
        hide_index=True,
        width="stretch",
        height=240,
        column_config={
            "mean": st.column_config.NumberColumn(format="%.3f"),
            "sd": st.column_config.NumberColumn(format="%.3f"),
            "missing_fraction": st.column_config.NumberColumn(format="%.3f"),
        },
    )
    gene = st.selectbox("Inspect one gene by group", matches)
    left, right = st.columns([2, 1])
    left.plotly_chart(gene_by_group_figure(expression, data.aligned, gene), width="stretch")
    right.dataframe(gene_group_summary(data, gene), hide_index=True, width="stretch")

st.subheader("Summary statistics")
left, right = st.columns(2)
left.markdown("**All observed expression values**")
left.dataframe(summary_statistics(expression), width="stretch")
right.markdown("**Samples per group and batch**")
if "batch" in metadata:
    right.dataframe(
        metadata.groupby(["group", "batch"]).size().unstack(fill_value=0), width="stretch"
    )
else:
    right.dataframe(metadata["group"].value_counts().rename("samples"), width="stretch")

footer()
