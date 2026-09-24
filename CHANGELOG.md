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

### Not yet implemented

- The educational machine-learning demonstration and the walkthrough notebook. The
  dashboard's Model Demonstration page states this explicitly.
