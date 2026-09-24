# Reproducibility

**Status:** Phase 1 baseline. Extended as features are added.

## Environment

- Python 3.11 or newer.
- **Dependency tool: pip** (decision recorded in Phase 0). pip needs no extra installation
  on any OS and is familiar to most reviewers. `pyproject.toml` is standards-based, so
  tools such as `uv` can also install from it.
- `pyproject.toml` declares compatible lower bounds. `requirements.txt` (runtime) and
  `requirements-dev.txt` (runtime + dev tools) pin the **direct** dependencies to the
  exact versions used for local verification. Transitive dependencies are not pinned,
  to avoid lock files that break on other operating systems.

```bash
python -m pip install -r requirements-dev.txt -e .
```

## Determinism

- The random seed lives in `config/default.yaml` (`seed`). The synthetic generator
  (Phase 2) uses `numpy.random.default_rng(seed)` exclusively. It does not use global
  random state.
- ML steps (Phase 3) receive `random_state` derived from the same seed.

## Cross-platform

- `pathlib` for every path. No hard-coded separators or absolute paths.
- `.gitattributes` normalizes line endings to LF.
- The `obw` CLI is the task runner on every OS. The `Makefile` is an optional shortcut
  for macOS and Linux.

## Run provenance (Phase 3)

Each analysis run will write `outputs/run_manifest.json` with the config hash, seed,
package versions, git commit (if available), and timestamp.
