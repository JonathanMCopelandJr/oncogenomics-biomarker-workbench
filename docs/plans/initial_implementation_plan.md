# Initial implementation plan

Approved in Phase 0 (2026-09-23). Work proceeds in phases. Each phase ends with a
summary, verification commands, open decisions, and a wait for explicit approval.

## Decisions (Phase 0 defaults, accepted)

| Topic | Decision |
|---|---|
| Code license | MIT |
| Demo-data license | CC0-1.0, for the synthetic data files only (confirmed) |
| Default branch | `main` (renamed from `master` before the first commit) |
| Commits | One local commit per approved phase; never pushed |
| Demo CSVs | Committed (≤ ~300 KB), seed-reproducible, checked by checksum test |
| Citation | Author "Jonathan Copeland"; ORCID omitted unless a verified iD is supplied |
| Dependency tool | pip + `pyproject.toml` + pinned `requirements*.txt` |
| Task runner | `obw` CLI (cross-platform) + thin `Makefile` for macOS/Linux |
| Import package / CLI | `onco_workbench` / `obw` |
| Gene identifiers | Neutral synthetic IDs (`SYN_G0001`), never real symbols |


## Phases

1. **Foundation.** Skeleton, packaging, `.gitignore`, config structure, `CLAUDE.md`,
   baseline docs, README with placeholders, a CLI skeleton, smoke tests.
2. **Synthetic data and validation.** Generator, config loader, loaders, schema
   validation, data dictionary, small demo files, tests.
3. **Analysis pipeline.** QC, normalization, group comparison, BH correction, ranking,
   figures, CSV and Markdown exports, run manifest, CLI subcommands, end-to-end run.
   *As delivered:* the Phase 3 approval scoped the work to "QC and group-comparison",
   so the ML demo (with leakage tests) and the walkthrough notebook were **not** built
   in Phase 3. Their scheduling is pending the user's decision. Plotly builders
   (`viz/interactive.py`) move to Phase 4, where the dashboard uses them.
4. **Streamlit dashboard.** Six pages on top of the package functions, disclaimers on
   every page, an AppTest smoke test. *As delivered:* it includes the plotly builders
   (`viz/interactive.py`). The Model Demonstration page shows the required warnings and
   an explicit "not yet implemented" notice, because the ML demo is still unscheduled.
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
