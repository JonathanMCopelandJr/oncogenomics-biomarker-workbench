"""Group Comparison: choose groups and thresholds; volcano plot, results table, CSV download."""

import pandas as pd
import streamlit as st

from onco_workbench.config import ConfigError
from onco_workbench.dashboard.components import (
    footer,
    get_comparison,
    page_setup,
    require_demo_data,
)
from onco_workbench.dashboard.data import filter_results, results_csv_bytes
from onco_workbench.viz.interactive import top_genes_figure, volcano_figure

page_setup("Group Comparison")
data = require_demo_data()
defaults = data.config.comparison
groups = data.groups

st.markdown(
    "Per-gene Welch t-test, Cohen's d, and Benjamini-Hochberg adjusted p-values, computed "
    "on the synthetic data. Every difference is **group B minus group A**."
)

cols = st.columns(4)
group_a = cols[0].selectbox(
    "Group A (reference)",
    groups,
    index=groups.index(defaults.group_a) if defaults.group_a in groups else 0,
)
group_b = cols[1].selectbox(
    "Group B (comparison)",
    groups,
    index=groups.index(defaults.group_b) if defaults.group_b in groups else min(1, len(groups) - 1),
)
fdr_threshold = cols[2].slider(
    "BH-adjusted p-value threshold",
    0.001,
    0.25,
    float(defaults.fdr_threshold),
    step=0.001,
    format="%.3f",
)
effect_threshold = cols[3].number_input(
    "|Cohen's d| threshold", 0.0, 5.0, float(defaults.effect_size_threshold), step=0.1
)

if group_a == group_b:
    st.error("Choose two different groups to compare.", icon=":material/error:")
    footer()
    st.stop()

try:
    comparison, result = get_comparison(group_a, group_b, fdr_threshold, effect_threshold)
except (ConfigError, ValueError) as exc:
    st.error(f"The comparison could not be run: {exc}", icon=":material/error:")
    footer()
    st.stop()

results = result.results
flagged = results["meets_thresholds"].astype(bool)
cols = st.columns(4)
cols[0].metric("Genes tested", int(results["p_value"].notna().sum()))
cols[1].metric("Meeting both thresholds", int(flagged.sum()))
cols[2].metric(
    f"Higher in {group_b}", int((flagged & (results["direction"] == "higher_in_b")).sum())
)
cols[3].metric(f"Lower in {group_b}", int((flagged & (results["direction"] == "lower_in_b")).sum()))
st.info(
    "Thresholds are display settings you choose. Meeting them only means that a simulated "
    "difference was detected by this workflow. It is not a biomarker, and it is not a "
    "biological or clinical finding.",
    icon=":material/info:",
)

st.plotly_chart(
    volcano_figure(
        results,
        group_a=group_a,
        group_b=group_b,
        fdr_threshold=fdr_threshold,
        effect_size_threshold=effect_threshold,
    ),
    width="stretch",
)
st.plotly_chart(
    top_genes_figure(results, group_a=group_a, group_b=group_b, top_n=data.config.ranking.top_n),
    width="stretch",
)
st.caption(
    f"Bars to the right (red) are higher in {group_b}; bars to the left (blue) are lower. "
    "Rank uses |d| x -log10(BH adjusted p), a documented sorting heuristic."
)

st.subheader("Results table")
left, right = st.columns([3, 1])
query = left.text_input("Filter by gene ID", placeholder="for example SYN_G01")
only_flagged = right.checkbox("Only genes meeting both thresholds", value=False)
view = filter_results(results, query=query, only_meeting_thresholds=only_flagged)
st.caption(f"{len(view)} of {len(results)} genes shown. Click a column header to sort.")
scientific = st.column_config.NumberColumn(format="%.2e")
fixed = st.column_config.NumberColumn(format="%.3f")
st.dataframe(
    view,
    hide_index=True,
    width="stretch",
    height=380,
    column_config={
        "p_value": scientific,
        "p_adj": scientific,
        "mean_a": fixed,
        "mean_b": fixed,
        "mean_diff": fixed,
        "sd_a": fixed,
        "sd_b": fixed,
        "cohens_d": fixed,
        "t_statistic": fixed,
        "df": st.column_config.NumberColumn(format="%.1f"),
        "ranking_score": st.column_config.NumberColumn(format="%.2f"),
    },
)
st.download_button(
    "Download shown rows (CSV, labeled synthetic demonstration)",
    data=results_csv_bytes(view, comparison),
    file_name=f"synthetic_demo_{group_b}_vs_{group_a}_results.csv",
    mime="text/csv",
    icon=":material/download:",
)

with st.expander("Workflow check against the synthetic ground truth"):
    if result.recovery is None:
        st.write(
            "The planted differences were simulated as Group_B relative to Group_A, so this "
            "check is only shown for that comparison order."
        )
    else:
        st.write(
            "The demo data contain a known set of planted differences. This table counts "
            "how many the workflow recovered at the current thresholds. It verifies the "
            "code only and says nothing about real data. High recovery is expected here "
            "because the planted shifts are large relative to the simulated noise."
        )
        st.dataframe(
            pd.DataFrame(result.recovery.as_rows(), columns=["Check", "Result"]),
            hide_index=True,
            width="stretch",
        )

footer()
