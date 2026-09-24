"""Support code for the Streamlit dashboard in ``app/``.

- :mod:`onco_workbench.dashboard.data` is pure, Streamlit-free logic, fully unit tested.
- :mod:`onco_workbench.dashboard.components` holds the Streamlit widgets and cached loaders.

All analysis goes through the same package functions as the CLI pipeline. The dashboard
reads only the configured SYNTHETIC demo files. It has no upload, external data
ingestion, authentication, or network access.
"""
