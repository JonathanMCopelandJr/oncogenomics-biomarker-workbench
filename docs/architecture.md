# Architecture

**Status:** Reflects Phases 1–4. Expanded in Phase 5.

## Principles

- **One source of logic.** All computation lives in the `onco_workbench` package. The CLI
  (`obw`), the notebooks, and the Streamlit app only call package functions.
- **Configuration over constants.** Seeds, thresholds, and paths come from
  `config/default.yaml`.
- **Pure functions where possible.** Analysis functions take DataFrames and return
  DataFrames or figures. File I/O is kept in `data/io.py` and `reporting/`.
- **One disclaimer source.** `onco_workbench.disclaimers` supplies all nonclinical wording.

## Package layout

| Module | Responsibility | Phase |
|---|---|---|
| `cli.py` | `obw` command and cross-platform task runner (`disclaimer`, `generate-data`, `validate`, `qc`, `run-analysis`, `ml-demo`, `dashboard`) | done |
| `disclaimers.py` | Nonclinical disclaimer text and data label | 1 (done) |
| `config.py` | Load and validate YAML into typed, frozen dataclasses; resolve paths | 2–3 (done) |
| `data/synthetic.py` | Seeded synthetic data generator | 2 (done) |
| `data/io.py` | Labeled CSV and manifest reading/writing, SHA-256 checksums | 2 (done) |
| `data/validation.py` | Schema and data-quality checks; `align_metadata` helper | 2 (done) |
| `analysis/qc.py` | QC summaries, per-sample and per-gene tables, PCA, sample correlation | 3 (done) |
| `analysis/normalization.py` | Median centering; plot-only imputation and z-scoring | 3 (done) |
| `analysis/differential.py` | Per-gene means, Cohen's d, vectorized Welch t-test | 3 (done) |
| `analysis/multiple_testing.py` | Benjamini-Hochberg (NaN-aware) | 3 (done) |
| `analysis/ranking.py` | Threshold flags, direction, transparent ranking score | 3 (done) |
| `analysis/recovery.py` | Synthetic ground-truth workflow check | 3 (done) |
| `viz/static.py` | matplotlib figures (object-oriented API, no global state), validated palette | 3 (done) |
| `reporting/markdown_report.py` | Markdown summary report | 3 (done) |
| `reporting/manifest.py` | Run manifest: versions, git state, parameters, checksums | 3 (done) |
| `pipeline.py` | Step functions (`load_and_validate`, `run_qc`, `run_comparison`) and `run_pipeline` | 3 (done) |
| `viz/interactive.py` | plotly figures for the dashboard (same palette, hover tooltips, labeled) | 4 (done) |
| `dashboard/data.py` | Streamlit-free dashboard helpers: load demo data, compare, filter, labeled CSV download | 4 (done) |
| `dashboard/components.py` | Streamlit banner, page setup, footer, cached loaders | 4 (done) |
| `ml/classifier_demo.py` | Optional educational classifier (synthetic only, NOT clinical). Stratified seeded split; imputer, scaler, and logistic-regression `Pipeline` fitted on training data only; gated CV; metrics; permuted-label baseline; labeled outputs. Independent of `pipeline.py` and `run-analysis` | done (optional) |
| `notebooks/oncogenomics_workbench_walkthrough.ipynb` (outside the package) | Walkthrough that only calls package functions. Committed with outputs stripped. Checked and executed by `tests/test_notebook.py` | done |

## Dashboard

```
app/Home.py, app/pages/*.py        UI only: layout, widgets, text
        │  import
        ▼
onco_workbench.dashboard.components  banner, page_setup(), footer(), @st.cache_* loaders
        │  call
        ▼
onco_workbench.dashboard.data        pure helpers (unit tested without Streamlit)
        │  call
        ▼
onco_workbench.pipeline / analysis / viz.interactive   the same code the CLI uses
```

- The demo data are loaded once per server process (`st.cache_resource`). Comparisons are
  cached per group/threshold combination (`st.cache_data`).
- The dashboard only reads the configured synthetic demo files. It never writes files.
  Its one download (the results CSV) is built in memory with the same labeled format
  as the pipeline exports.
- `.streamlit/config.toml` binds the server to `localhost` and turns off usage
  statistics.

## Data flow

See the Mermaid diagram in the [README](../README.md#architecture).

`run_pipeline` works in this order:

1. Load and validate the inputs. It stops before writing anything if validation fails.
2. Run QC.
3. Median-center the samples, compare the groups, apply BH correction, flag, and rank.
4. Optionally run the ground-truth check.
5. Write the tables, figures, report, and manifest (in that order, so the manifest can
   checksum everything else).

The step functions are public so that the dashboard can call them on in-memory data
without writing files.
