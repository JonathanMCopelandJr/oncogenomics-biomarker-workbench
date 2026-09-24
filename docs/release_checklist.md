# Release-readiness checklist

> Research and education software on **synthetic** data. It is not for clinical use.

Status key: **Done** = verified in the Phase 5 audit (2026-09-23, Windows 11, Python 3.12.10).
**Owner** = needs an action or decision by the repository owner. No push, release, tag, or
GitHub-settings change has been made.

## Code quality

| Check | Status | Evidence |
|---|---|---|
| `ruff format --check src tests app` | Done | 56 files already formatted (after the ML demonstration was added) |
| `ruff check src tests app` | Done | All checks passed |
| Strict one-off audit: `ruff --select ANN,S,PERF,RUF` | Done | 10 findings. 1 fixed (`RUF046`). 9 reviewed and accepted: `ANN401` on the YAML parsers, which take untrusted `Any` values and check them at runtime; `S603`/`S607` on the fixed-argument, no-shell `streamlit` and `git` subprocess calls |
| No `TODO`/`FIXME`, no stray `print` outside the CLI | Done | grep: none |
| Public functions typed and documented | Done | ruff `D` (Google style) enforced; `ANN` audit shows no missing annotations |
| Full test suite | Done | 300 passed, including 26 ML-demo tests. Line/branch coverage 96% (after the ML demonstration was added) |
| `pip check` | Done | No broken requirements |

## Reproducibility

| Check | Status | Evidence |
|---|---|---|
| Fresh environment from pinned requirements (a new venv; files copied as a fresh clone would be) | Done | `pip install -r requirements-dev.txt -e .` succeeded with exactly the pinned versions. `pip check` clean. (The venv had to be in a short path because of Windows `WinError 206`; the README now documents this) |
| Lint and tests in the fresh environment | Done (re-run after the ML demo and notebook) | Clean copy and a separate venv, both in short paths outside the development repository; Python 3.12.10; README Windows steps plus `-r requirements-notebooks.txt`. `pip check` clean; 57 files formatted; lint clean; **309 passed** |
| Full workflow in the fresh environment | Done | `obw generate-data` (byte-identical to committed), `obw run-analysis`, and `obw ml-demo` all exit 0. The notebook executed with `python -m nbconvert --execute`: 12/12 cells, 0 errors, 0 stderr, no kernel registration needed. The dashboard launched: `/_stcore/health` returned 200 `ok`, `/` returned 200, and it stopped cleanly |
| Demo data regenerate byte-identically from the seed | Done | All 4 files (3 CSVs and the manifest) have SHA-256 identical to the committed files |
| Analysis tables identical across environments | Done | All 5 CSV tables from `obw run-analysis` are byte-identical between the development and fresh environments. Both runs: 40/40 planted genes recovered, 0 unplanted genes flagged (expected by construction) |
| Analysis tables identical across repeated runs | Done | Phase 3: all 5 CSV tables byte-identical |
| Python 3.11 (NumPy 2.4.6 / SciPy 1.17.1 pins) | Owner | Only Python 3.12 is installed locally. It is covered by the CI matrix once the repository is pushed |
| Linux and macOS | Owner | Covered by the CI matrix once pushed |

## Repository setup

| Check | Status | Evidence |
|---|---|---|
| Required docs and community files present | Done | `tests/test_repository.py` |
| CI workflow valid; read-only permissions; no secrets; no deploy or publish steps | Done | `tests/test_repository.py` parses and checks `ci.yml` |
| `CITATION.cff` valid; author Jonathan Copeland; no ORCID | Done | `tests/test_repository.py` |
| `.gitignore` covers venvs, caches, `outputs/`, `data/raw/`, `data/processed/`, `.env`, `.streamlit/secrets.toml` | Done | `git check-ignore` on 9 sample paths: all ignored. Demo data, `.streamlit/config.toml`, and `.gitkeep` files are still tracked |
| No secrets, emails, or absolute user paths in tracked files | Done | `git grep` scans: none |
| No CRLF line endings in the index | Done | `git ls-files --eol`: none |
| Largest tracked file | Done | `demo_expression.csv`, 249,430 bytes. Demo data total 267,334 bytes |
| Package builds | Done | `pip wheel .` produced `oncogenomics_biomarker_workbench-0.1.0-py3-none-any.whl` (33 entries: package code and LICENSE only) |
| Licenses | Done | Code MIT. Synthetic data CC0-1.0 (data files only). Direct dependencies are under permissive licenses (BSD, MIT, Apache-2.0, PSF-style) |

## Scientific and ethical boundary

| Check | Status | Evidence |
|---|---|---|
| Disclaimer on README, docs, CSV exports, figures, report, manifest, and every dashboard page | Done | Tests in `test_package`, `test_io`, `test_viz_static`, `test_interactive`, `test_reporting`, `test_app_smoke`, `test_repository` |
| No real gene symbols; neutral synthetic IDs only | Done | `test_synthetic.py` |
| No clinical, biomarker, or performance claims. ML outputs and page carry the "NOT a clinical prediction model" warning, the positive-class statement, and the limitations | Done | `test_ml_demo.py::test_outputs_are_labeled_and_contain_no_model_weights`, `test_app_smoke.py::test_model_page_shows_labeled_results_and_limitations`, `test_cli.py::test_ml_demo_command` |
| ML demonstration is leakage-safe (preprocessing fitted on the training split only; test-set changes do not affect the fit or CV) | Done | `test_ml_demo.py` leakage and invariance tests |
| No external data, downloads, credentials, authentication, or deployment | Done | Code review. Dashboard bound to `localhost` |

## Owner actions before a public release

1. **Push and confirm CI passes** on all four matrix jobs, including Python 3.11.
2. **ML demonstration and walkthrough notebook:** both are implemented. The notebook
   executes cleanly with `python -m nbconvert --execute` and in `tests/test_notebook.py`.
   The fresh-environment reproduction has been re-run to include both (see Reproducibility).
3. **Add dashboard screenshots** to `docs/assets/screenshots/` and link them from the
   README. They weren't captured in the audit because no browser was available.
4. **Add the repository URL** to `CITATION.cff` (`repository-code`) and, if wanted, a CI
   status badge to the README. The URL isn't known locally.
5. **GitHub settings** (not changed by this work): consider enabling private
   vulnerability reporting, which `SECURITY.md` refers to, and branch protection on `main`.
6. **Code of Conduct contact:** confirm that "contact the maintainer through their
   GitHub profile" is acceptable, or add a dedicated contact.
7. **To release:** move the CHANGELOG `[Unreleased]` entries under a dated `[0.1.0]`
   heading, add `date-released` to `CITATION.cff`, then tag `v0.1.0`.
