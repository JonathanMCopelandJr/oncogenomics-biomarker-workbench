# Initial implementation plan

Approved in Phase 0 (2026-09-23). Work proceeds in phases. Each phase ends with a
summary, verification commands, open decisions, and a wait for explicit approval.

## Decisions (Phase 0 defaults, accepted)

| Topic | Decision |
|---|---|
| Code license | MIT |
| Demo-data license | CC0-1.0 (to be confirmed in Phase 2) |
| Dependency tool | pip + `pyproject.toml` + pinned `requirements*.txt` |
| Task runner | `obw` CLI (cross-platform) + thin `Makefile` for macOS/Linux |
| Import package / CLI | `onco_workbench` / `obw` |
| Gene identifiers | Neutral synthetic IDs (`SYN_G0001`), never real symbols |

**Pending user approval:** rename `master` to `main`; local commit per phase; commit
the small demo CSVs (Phase 2); citation author metadata and ORCID (Phase 5).

## Phases

1. **Foundation.** Skeleton, packaging, `.gitignore`, config structure, `CLAUDE.md`,
   baseline docs, README with placeholders, a CLI skeleton, smoke tests.
2. **Synthetic data and validation.** Generator, config loader, loaders, schema
   validation, data dictionary, small demo files, tests.
3. **Analysis pipeline.** QC, normalization, group comparison, BH correction, ranking,
   figures, CSV and Markdown exports, run manifest, CLI subcommands, ML demo with
   leakage tests, walkthrough notebook, end-to-end run.
4. **Streamlit dashboard.** Six pages on top of the package functions, disclaimers on
   every page, an AppTest smoke test.
5. **Quality, CI, and polish.** GitHub Actions (Ubuntu 3.11/3.12, Windows and macOS on
   3.12), community files, `CITATION.cff`, `CHANGELOG.md`, docs completion, a
   release-readiness checklist, and a final audit. No push or publish.

## Dependencies

- **Runtime:** numpy, pandas, scipy, statsmodels, scikit-learn, plotly, matplotlib,
  streamlit, pyyaml.
- **Dev:** pytest, pytest-cov, ruff.
- **Optional `[notebooks]`:** jupyterlab, ipykernel, nbstripout.
- **Deliberately excluded:** seaborn, pandera/pydantic, click/typer, kaleido. Docker is
  deferred until after Phase 5.

## Key risks and mitigations

See `docs/limitations_and_ethics.md`. In summary: the nonclinical boundary is enforced
through a single disclaimer source, neutral IDs, and `CLAUDE.md` rules. Double-dipping
in the ML demo is prevented by a pipeline and a test. Batch confounding is handled by a
balanced generator plus a validation warning. Dependency drift is handled by pins.
Repository hygiene is handled by `.gitignore` rules for data, outputs, and secrets.
