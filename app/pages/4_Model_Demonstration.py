"""Model Demonstration: placeholder until the educational ML example is approved and built."""

import streamlit as st

from onco_workbench.dashboard.components import footer, page_setup

page_setup("Model Demonstration")

st.error(
    "**Synthetic-data educational example only. This is NOT a clinical prediction model.** "
    "No model on this page predicts anything about any person. Any classifier trained on "
    "these simulated data will look accurate **by construction**, because the group "
    "differences were planted. Such numbers say nothing about performance on real data.",
    icon=":material/gpp_maybe:",
)
st.warning(
    "Small samples and simulated data make any evaluation metric unstable and optimistic. "
    "Metrics from this demonstration must never be quoted as predictive performance.",
    icon=":material/warning:",
)

st.info(
    "**Not yet implemented.** The educational classifier is planned but has not been built, "
    "and its scheduling is pending the project owner's approval. No model results are shown "
    "here. None have been computed, and none are simulated or invented.",
    icon=":material/construction:",
)

st.subheader("Planned design")
st.markdown(
    """
- A stratified train/test split of the synthetic samples.
- A scikit-learn `Pipeline` of `StandardScaler` and `LogisticRegression`, fitted on the
  **training data only**, so no information leaks from the test set.
- No feature selection on the full dataset (that would be double-dipping).
- Accuracy, precision, recall, F1, ROC-AUC (when defined), and a confusion matrix, all on
  held-out data.
- Cross-validation only when there are enough samples per class. A label-permutation
  check shows what chance-level performance looks like.
- Automated tests that check the no-leakage property.
"""
)

footer()
