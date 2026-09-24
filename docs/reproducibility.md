# Reproducibility

**Status:** Covers Phases 1–3 (environment, synthetic data, analysis pipeline).

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

- The random seed lives in `config/default.yaml` (`seed: 42`). The synthetic generator
  uses one `numpy.random.default_rng(seed)` generator in a fixed order and never uses
  global random state.
- CSVs are written with a fixed float format (`%.3f`), UTF-8 encoding, and LF line
  endings. The data manifest contains no timestamps or environment details.
- The analysis pipeline involves no randomness. PCA uses a full SVD, and heatmap
  sample selection is a deterministic, group-balanced choice by sample ID. Result
  ordering breaks ties by `gene_id`. Re-running on the same inputs gives the same
  tables. Only the timestamp and git state in the run manifest change.
- The planned ML demonstration will receive `random_state` derived from the same seed.

## Committed demo data and checksums

The small synthetic demo files in `data/synthetic/` (about 267 KB in total) are committed
so that the project works right after cloning. `demo_manifest.json` records the SHA-256
checksum of each CSV. `tests/test_demo_data.py` regenerates the dataset from the
configured seed and asserts that the new files match the committed ones byte for byte.
To regenerate the files deliberately, run `obw generate-data`.

**Caveat:** NumPy does not guarantee identical random streams across NumPy versions.
The committed files were generated with the pinned NumPy for Python 3.12 (2.5.3). If a
different NumPy version ever changes the stream, the regeneration test fails loudly
instead of silently shipping different data.

## Cross-platform

- `pathlib` for every path. No hard-coded separators or absolute paths.
- `.gitattributes` normalizes line endings to LF.
- The `obw` CLI is the task runner on every OS. The `Makefile` is an optional shortcut
  for macOS and Linux.

## Run provenance

Every `obw qc` and `obw run-analysis` run writes `run_manifest.json` next to its
outputs. The manifest records:

- the UTC timestamp, workbench version, Python version, and OS family
- versions of the key scientific packages
- the git commit, and whether the working tree had uncommitted changes (`dirty`)
- the configuration file path and SHA-256
- the seed and all analysis parameters
- SHA-256 checksums of every input and output file

Paths are stored relative to the project root, or as bare file names outside it, so
manifests do not reveal local directory names.

To reproduce a run: check out the recorded commit, install the pinned requirements,
and run the same command with the same configuration. Then compare the output
checksums with the manifest. PNG checksums may differ across matplotlib versions or
platforms because of font rendering. The CSV tables are the reproducible record.
