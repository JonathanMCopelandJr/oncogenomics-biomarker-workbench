# Contributing

Thank you for your interest. This is a **research and education** project that works
on **synthetic** data. Contributions must keep it that way.

## Ground rules

1. **No clinical claims.** Do not add diagnostic, prognostic, survival, treatment,
   drug-efficacy, or patient-level claims or recommendations. Nothing in this project
   identifies or validates a biomarker.
2. **No real patient data, restricted data, or credentials.** Never commit protected
   health information, controlled-access data, API keys, tokens, or `.env` files.
3. **No fabricated science.** Do not add invented citations, datasets, results, or
   performance numbers. Any number in the docs must come from code in this repository.
4. **Public data only with approval.** Small, openly licensed public data may be proposed
   in an issue first, with its source and license documented.
5. Keep the **DEMONSTRATION / SYNTHETIC** labels and the disclaimer on every output and
   dashboard page. Reuse `src/onco_workbench/disclaimers.py` rather than rewording it.

Please also follow the [Code of Conduct](CODE_OF_CONDUCT.md).

## Development setup

```bash
python -m venv .venv
# macOS/Linux: source .venv/bin/activate    Windows: .\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements-dev.txt -e .
```

## Workflow

1. Open an issue describing the change (bug or feature template).
2. Create a branch and make a **small, reviewable** change.
3. Put logic in `src/onco_workbench/`. The CLI (`cli.py`), notebooks, and the Streamlit
   pages in `app/` only call package functions.
4. Add or update tests. Public functions need type hints and Google-style docstrings.
5. Update the docs and `CHANGELOG.md` when behavior changes.
6. Run the checks:

   ```bash
   python -m ruff format src tests app
   python -m ruff check src tests app
   python -m pytest
   ```

7. Open a pull request and complete its checklist. CI runs the same checks on Linux,
   Windows, and macOS.

Commit messages follow the conventional style: `feat: ...`, `fix: ...`, `docs: ...`,
`test: ...`, `chore: ...`.

## Changing the synthetic data generator

If a change alters generator output, regenerate the committed demo files with
`obw generate-data` in the same pull request and explain why.
`tests/test_demo_data.py` fails if the committed files and a fresh regeneration differ.
