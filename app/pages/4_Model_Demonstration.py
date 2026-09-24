"""Model Demonstration: optional educational ML example on SYNTHETIC data only.

Not a clinical prediction model. All logic lives in ``onco_workbench.ml.classifier_demo``.
"""

import pandas as pd
import streamlit as st

from onco_workbench.dashboard.components import footer, get_ml_demo, page_setup, require_demo_data
from onco_workbench.ml.classifier_demo import (
    HEADLINE_WARNING,
    LIMITATIONS,
    InsufficientDataError,
    output_statements,
)
from onco_workbench.viz.interactive import confusion_matrix_figure

page_setup("Model Demonstration")

st.error(f"**{HEADLINE_WARNING}**", icon=":material/gpp_maybe:")
st.warning(
    "Small samples and simulated data make every evaluation metric unstable and optimistic. "
    "Metrics from this demonstration must never be quoted as predictive performance.",
    icon=":material/warning:",
)

data = require_demo_data()
if not data.config.ml_demo.enabled:
    st.info("The ML demonstration is disabled in the configuration (`ml_demo.enabled: false`).")
    footer()
    st.stop()

try:
    result = get_ml_demo()
except InsufficientDataError as exc:
    st.error(f"The ML demonstration cannot run on these data: {exc}", icon=":material/error:")
    footer()
    st.stop()


def fmt(value: float | None) -> str:
    """Format a metric, showing undefined values explicitly."""
    return "undefined" if value is None else f"{value:.2f}"


labels = result.labels
st.info(
    "\n".join(f"- **{statement}**" for statement in output_statements(labels.positive)),
    icon=":material/info:",
)
st.subheader("Setup")
st.markdown(
    f"""
- **Classes (synthetic, arbitrary):** positive class **{labels.positive}** (`comparison.group_b`),
  negative class **{labels.negative}** (`comparison.group_a`). These are not diagnoses or
  clinical categories.
- **Features:** all {result.n_features} raw synthetic gene values. No feature selection.
- **Split:** stratified, seed {result.seed}. Training {result.train_class_counts},
  test {result.test_class_counts}.
- **Pipeline:** mean imputation, then standardization, then logistic regression
  (C = {result.settings.logistic_regression_c:g}). **Fitted on the training split only.**
  Test samples only pass through the already fitted pipeline.
- This preprocessing is intentionally separate from the Group Comparison page's per-sample
  median centering, so each workflow's transformations stay independent.
"""
)

st.subheader("Held-out synthetic test set")
st.caption(
    f"Positive class: {labels.positive}. Precision, recall, and F1 refer to {labels.positive}."
)
m = result.metrics
cols = st.columns(5)
cols[0].metric("Accuracy", fmt(m.accuracy))
cols[1].metric(f"Precision ({labels.positive})", fmt(m.precision))
cols[2].metric(f"Recall ({labels.positive})", fmt(m.recall))
cols[3].metric(f"F1 ({labels.positive})", fmt(m.f1))
cols[4].metric("ROC-AUC", fmt(m.roc_auc))
for note in m.notes:
    st.caption(note)
left, right = st.columns([3, 2])
left.plotly_chart(
    confusion_matrix_figure(m.confusion, positive_label=labels.positive), width="stretch"
)
right.markdown("**Confusion matrix** (rows = true group, columns = predicted group)")
right.dataframe(m.confusion, width="stretch")

st.subheader("Cross-validation (training split only)")
st.write(result.cv.reason)
if result.cv.ran:
    st.dataframe(
        pd.DataFrame(
            {
                "fold": range(1, result.cv.n_folds + 1),
                "validation samples": result.cv.fold_sizes,
                "accuracy": result.cv.accuracies,
                "ROC-AUC": result.cv.roc_aucs,
            }
        ),
        hide_index=True,
        width="stretch",
    )

st.subheader("Correct labels vs permuted-label baseline")
perm = result.permutation
if perm.ran:
    cols = st.columns(2)
    cols[0].metric("Test accuracy, correct synthetic labels", fmt(m.accuracy))
    cols[1].metric(
        f"Mean test accuracy, shuffled labels ({perm.n_permutations} runs)",
        fmt(perm.accuracy_mean),
    )
st.info(result.comparison_text(), icon=":material/info:")

st.subheader("Limitations")
st.warning("\n".join(f"- {text}" for text in LIMITATIONS), icon=":material/warning:")

footer()
