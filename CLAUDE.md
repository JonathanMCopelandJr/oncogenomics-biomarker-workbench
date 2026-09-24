# CLAUDE.md

Guidance for Claude Code (and other AI assistants) working in this repository.

## What this project is

`oncogenomics-biomarker-workbench` is a **research and education portfolio project**
that runs a transparent validation, QC, group-comparison, and visualization workflow
on **SYNTHETIC** transcriptomics-style data. It is **not** a medical device,
diagnostic tool, clinical decision-support system, treatment recommender, or
validated biomarker-discovery system.

## Non-negotiable rules

1. **Never fabricate** scientific claims, citations, DOIs, real datasets, results,
   benchmark outcomes, or performance metrics. If a number is reported, it must come
   from code that was actually run in this repository.
2. **Preserve the research-only / nonclinical boundary.** No diagnostic, prognostic,
   survival, drug-efficacy, or patient-level claims or recommendations anywhere
   (code, UI, docs, commit messages). Reuse the text in
   `src/onco_workbench/disclaimers.py` instead of rewording it.
3. **Ask for clarification** before adding any medical claim, real-gene biological
   interpretation, or restricted / controlled-access / patient data. Do not proceed
   on assumption.
4. **Never add credentials**, API keys, tokens, `.env` files with secrets, or private data.
5. **Data in Git:** keep raw data out of Git unless it is small, openly licensed,
   documented (source + license in a README), and explicitly approved. Never commit
   large generated data or anything under `outputs/` or `data/processed/`.
6. **Synthetic gene IDs stay neutral** (`SYN_G0001`, ...). Do not map them to real
   gene symbols.
7. Label demo outputs **DEMONSTRATION / SYNTHETIC** (see `disclaimers.DATA_LABEL`).

## Working practices

- Prefer **small, reviewable changes**. Keep business logic in `src/onco_workbench/`;
  Streamlit pages (`app/`) and notebooks only call package functions.
- Every new major feature ships with **tests and documentation** updates.
- When behavior changes, update the relevant files in `docs/`, the README, and
  (from Phase 5) `CHANGELOG.md`.
- Public functions need type hints and Google-style docstrings (ruff `D` rules enforce this).
- Use `pathlib` and paths from `config/default.yaml`; no hard-coded absolute paths.
- Do not push, create releases, change repository visibility, or deploy.
- Do not run destructive shell commands.

## Before proposing completion

Run, and report the real outcome of:

```bash
python -m ruff format --check src tests
python -m ruff check src tests
python -m pytest
```

If anything fails, say so and show the output. Do not claim success without running them.

## Commands

| Task | Command |
|---|---|
| Install (in an activated venv) | `python -m pip install -r requirements-dev.txt -e .` |
| Format | `python -m ruff format src tests` |
| Lint | `python -m ruff check src tests` |
| Test | `python -m pytest` |
| CLI help | `obw --help` |
| Regenerate committed demo data | `obw generate-data` (then review the diff; `tests/test_demo_data.py` checks checksums) |
| Validate data | `obw validate [--expression X.csv --metadata Y.csv]` |

If a change alters generator output, regenerate the demo files in the same change and
say so explicitly. Never edit files in `data/synthetic/` by hand.

## Project phases

Work is delivered in approved phases (see `docs/plans/initial_implementation_plan.md`).
Stop at the end of each phase, summarize changes and verification commands, and wait
for explicit approval before starting the next phase.
