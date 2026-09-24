# Architecture

**Status:** Planned structure (Phase 1). Expanded in Phase 5.

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
| `cli.py` | `obw` command and cross-platform task runner | 1 (skeleton) |
| `disclaimers.py` | Nonclinical disclaimer text and data label | 1 |
| `config.py` | Load and validate YAML configuration, resolve paths | 2 |
| `data/synthetic.py` | Seeded synthetic data generator | 2 |
| `data/io.py` | Read and write the expression, metadata, and ground-truth files | 2 |
| `data/validation.py` | Schema and data-quality checks with explicit errors | 2 |
| `analysis/qc.py` | QC summaries and PCA | 3 |
| `analysis/normalization.py` | Median centering, z-scoring | 3 |
| `analysis/differential.py` | Per-gene group statistics | 3 |
| `analysis/multiple_testing.py` | Benjamini-Hochberg wrapper | 3 |
| `analysis/ranking.py` | Transparent ranking score | 3 |
| `ml/classifier_demo.py` | Leakage-safe educational classifier | 3 |
| `viz/static.py`, `viz/interactive.py` | matplotlib PNGs, plotly figures | 3 |
| `reporting/` | Markdown report, run manifest | 3 |
| `pipeline.py` | End-to-end orchestration | 3 |

## Data flow

See the Mermaid diagram in the [README](../README.md#architecture).
