"""Methods and Limitations: plain-language workflow, limitations, provenance, reproducibility."""

import pandas as pd
import streamlit as st

from onco_workbench.dashboard.components import footer, page_setup, require_demo_data
from onco_workbench.dashboard.data import environment_info

page_setup("Methods and Limitations")
data = require_demo_data()
c, r = data.config.comparison, data.config.ranking

st.subheader("How the workflow works")
st.markdown(
    f"""
1. **Simulate data.** A generator creates {data.expression.shape[0]} samples and
   {data.expression.shape[1]} synthetic genes from a fixed random seed. A few genes get a
   deliberate shift in one group. A small share of values is removed to exercise the
   missing-data handling.
2. **Validate.** Checks for duplicate or blank identifiers, non-numeric values, missing
   data, unexpected group labels, samples without metadata, group sizes, and whether
   batch is confounded with group.
3. **Quality control.** Counts, missingness, per-sample distributions, PCA, and
   sample-to-sample correlation. The PCA and heatmaps use gene-scaled values for display only.
4. **Normalize.** Each sample's median is aligned to a common level, which removes
   per-sample offsets while keeping the original scale.
5. **Compare two groups.** For every gene: group means, the difference (B minus A),
   Cohen's d, and a Welch t-test. The p-values are adjusted for testing many genes at once
   with the Benjamini-Hochberg procedure.
6. **Flag and rank.** Genes meet the display thresholds when the adjusted p-value is at
   most {c.fdr_threshold:g} and |d| is at least {c.effect_size_threshold:g} (by default;
   adjustable on the Group Comparison page). The ranking score |d| x -log10(max(adjusted
   p, {r.p_floor:g})) only sorts the table.
7. **Check against the planted answer.** On the demo data, flagged genes are compared
   with the genes that were deliberately shifted. This checks the code, not biology.
8. **Optional ML demonstration** (separate from steps 1-7). A logistic-regression
   classifier is trained on a stratified training split of the synthetic groups and
   evaluated on held-out samples. All preprocessing happens inside a pipeline fitted on
   training data only. It is compared with a permuted-label baseline. It is an
   educational example, not a clinical prediction model.

Full formulas are in `docs/methodology.md`, and field definitions are in
`docs/data_dictionary.md`.
"""
)

st.subheader("Limitations")
st.markdown(
    """
- **All data are simulated.** Gene IDs do not correspond to real genes, and `Group_A` /
  `Group_B` are not diagnoses or clinical categories.
- **The statistics are a teaching simplification.** A Welch t-test on continuous values
  is not appropriate for real RNA-seq counts, which need count-aware models,
  normalization, and design-aware handling of covariates.
- **Recovering planted signal only shows the code works.** It is not evidence of
  biological relevance, diagnostic value, or clinical validity.
- **Thresholds are display choices.** Meeting them is not a scientific conclusion.
- **False discovery control has assumptions.** Benjamini-Hochberg control depends on the
  dependence between tests. Real genes are correlated in ways these data do not simulate.
- **Nothing in this workbench identifies or validates a biomarker**, and nothing in it
  may be used to inform any decision about any person's health.
"""
)

st.subheader("Data provenance")
manifest = data.data_manifest
if manifest is None:
    st.warning("The demo data manifest (`demo_manifest.json`) was not found.")
else:
    st.markdown(
        f"""
- {manifest["description"]}
- Generated with seed **{manifest["seed"]}** using `{manifest["random_generator"]}`.
- Data license: **{manifest["license"]["spdx"]}** ({manifest["license"]["name"]}).
  {manifest["license"]["scope"]}
- No external or patient data are used anywhere in this project.
"""
    )

st.subheader("Reproducibility")
st.markdown(
    "Re-running with the same configuration and seed reproduces the demo data byte for byte, "
    "and the analysis tables exactly. `obw run-analysis` also writes `run_manifest.json`, "
    "which records versions, parameters, the git commit, and file checksums."
)
st.dataframe(
    pd.DataFrame(environment_info(data), columns=["Item", "Value"]),
    hide_index=True,
    width="stretch",
)

footer()
