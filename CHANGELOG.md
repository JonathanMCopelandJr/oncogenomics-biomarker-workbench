# Changelog

All notable changes to this project are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

No version has been released yet. The package metadata declares `0.1.0` as the first
planned version.

## [Unreleased]

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
- **Dashboard:** a six-page local Streamlit app (`obw dashboard`) with interactive plotly
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

### Not yet implemented

- The walkthrough notebook (planned as a separate, later deliverable).
