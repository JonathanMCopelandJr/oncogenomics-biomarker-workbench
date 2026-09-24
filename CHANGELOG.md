# Changelog

All notable changes to this project are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-09-23

Initial public release. A reproducible, **synthetic-data, research-and-education-only**
workbench for exploring transcriptomics-style analysis workflows. Its scope is
nonclinical: it is **not** a medical device or a clinical decision-support system. The
release includes documentation, a walkthrough notebook, a local Streamlit dashboard and
command-line workbench, a test suite, and cross-platform CI (Ubuntu, Windows, and macOS).
Code is released under the MIT License. The synthetic demo data files are dedicated under
CC0-1.0.

### Added

- **Foundation:** src-layout package `onco_workbench`, `obw` CLI, pinned dependencies,
  ruff and pytest configuration, `.gitignore`, `CLAUDE.md`, and baseline documentation.
  A single source for the nonclinical disclaimer text.
- **Synthetic data and validation:**
  - A seeded generator of synthetic transcriptomics-style data with planted signal genes
    and neutral IDs (`SYN_G0001`, ...).
  - Typed, validated YAML configuration.
  - Labeled CSV input/output and a data manifest with SHA-256 checksums.
  - Validation that reports every problem at once (23 issue codes: 18 errors, 5 warnings).
  - Committed demo data (80 samples x 500 genes, about 267 KB, CC0-1.0) with a
    byte-identical regeneration test.
- **Analysis pipeline:**
  - QC summaries, PCA, and sample correlation.
  - Median-centering normalization.
  - Per-gene means, Cohen's d, and a Welch t-test (cross-checked against SciPy).
  - Benjamini-Hochberg adjustment, threshold flags, and a transparent ranking score.
  - A synthetic ground-truth workflow check.
  - Labeled CSV tables and seven matplotlib figures with a CVD-validated palette.
  - A Markdown report and a run manifest.
  - `obw qc` and `obw run-analysis`.
- **Dashboard:** a local Streamlit app (`obw dashboard`) with interactive plotly
  charts, a filterable results table, and a labeled CSV download. It shows the
  research-only / synthetic-data banner on every page and has AppTest smoke tests.
- **Quality and project files:** GitHub Actions CI (lint, plus tests on Ubuntu with
  Python 3.11/3.12, Windows, and macOS), issue and pull-request templates,
  `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `SECURITY.md`, `CITATION.cff`,
  a release-readiness checklist, and repository hygiene tests.
- **Optional educational ML demonstration** (synthetic data only, NOT a clinical
  prediction model):
  - `onco_workbench.ml.classifier_demo` and a new `obw ml-demo` command writing labeled
    outputs to `outputs/ml_demo/`.
  - A stratified, seeded train/test split, and an imputer → scaler → logistic-regression
    `Pipeline` fitted on training data only.
  - Cross-validation that runs only when training class sizes meet
    `ml_demo.min_train_per_class_for_cv`.
  - Test-set metrics with explicit "undefined" handling, and a confusion matrix.
  - A deterministic permuted-label sanity-check baseline.
  - A typed `ml_demo` configuration section (new key `min_samples_per_class`).
  - The dashboard's Model Demonstration page now shows the results with the limitations.
  - Tests for determinism, no preprocessing leakage, metric behavior, and insufficient
    class counts.
  - `obw run-analysis` and the group-comparison analysis are unchanged.
- **Walkthrough notebook** (`notebooks/oncogenomics_workbench_walkthrough.ipynb`):
  - A 5–10 minute tour on synthetic data only, opening with the required disclaimer.
  - Every cell calls the `onco_workbench` package; no logic is duplicated.
  - Committed with outputs stripped.
  - `tests/test_notebook.py` checks the structure, disclaimer, section order, stripped
    outputs, and absence of duplicated logic, and executes every cell without Jupyter.
  - It adds no new dependencies; it uses the existing optional `notebooks` extra.
- **Phase 6A documentation: related work and a future roadmap (documentation only).**
  - New `docs/related_work.md` and a README "Related Work and Future Extensions" section.
  - Records that the public list `hussius/deeplearning-biology` was consulted only as a
    discovery index. Nothing was copied, reproduced, vendored, depended on, or
    benchmarked against.
  - Describes a possible future synthetic-data autoencoder as **not implemented**.
  - A matching guardrail was added to `CLAUDE.md`.
  - No code, data, models, dependencies, tests, configuration, or application behavior
    changed.
