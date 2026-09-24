"""Run manifests: a machine-readable provenance record for each analysis run.

A manifest records what was run (configuration, parameters, seed), on what (input
file checksums), with what (Python and package versions, git commit), and what it
produced (output file checksums). Paths are stored relative to the project root when
possible, so manifests do not expose local directory names.
"""

from __future__ import annotations

import json
import platform
import subprocess
from collections.abc import Mapping
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from onco_workbench import __version__
from onco_workbench.data.io import file_sha256
from onco_workbench.disclaimers import DATA_LABEL, SHORT_DISCLAIMER

TRACKED_PACKAGES: tuple[str, ...] = (
    "numpy",
    "pandas",
    "scipy",
    "statsmodels",
    "scikit-learn",
    "matplotlib",
    "PyYAML",
)


def display_path(path: str | Path, root: str | Path) -> str:
    """Return ``path`` relative to ``root`` (POSIX style) or, failing that, its file name.

    Args:
        path: Path to display.
        root: Project root.

    Returns:
        A string that does not reveal directories outside the project.
    """
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(Path(root).resolve()).as_posix()
    except ValueError:
        return resolved.name


def package_versions(names: tuple[str, ...] = TRACKED_PACKAGES) -> dict[str, str | None]:
    """Return installed versions of selected packages (None if not installed).

    Args:
        names: Distribution names.

    Returns:
        Mapping of name to version string or None.
    """
    versions: dict[str, str | None] = {}
    for name in names:
        try:
            versions[name] = version(name)
        except PackageNotFoundError:
            versions[name] = None
    return versions


def git_state(root: str | Path) -> dict[str, Any]:
    """Return the current git commit and whether the tree has uncommitted changes.

    Never raises. Returns ``{"commit": None, "dirty": None}`` if git is unavailable
    or ``root`` is not a repository.

    Args:
        root: Directory inside the repository.

    Returns:
        Mapping with ``commit`` (str or None) and ``dirty`` (bool or None).
    """

    def run(*args: str) -> str | None:
        try:
            completed = subprocess.run(
                ["git", *args], cwd=root, capture_output=True, text=True, timeout=10, check=False
            )
        except (OSError, subprocess.SubprocessError):
            return None
        return completed.stdout.strip() if completed.returncode == 0 else None

    commit = run("rev-parse", "HEAD")
    status = run("status", "--porcelain") if commit else None
    return {"commit": commit, "dirty": None if status is None else bool(status)}


def build_run_manifest(
    *,
    root: str | Path,
    config_path: str | Path,
    seed: int,
    parameters: Mapping[str, Any],
    inputs: Mapping[str, Path],
    outputs: Mapping[str, Path],
    generated_at: str,
) -> dict[str, Any]:
    """Assemble a JSON-serializable manifest for an analysis run.

    Args:
        root: Project root, used to make paths relative.
        config_path: Configuration file used.
        seed: Configured random seed.
        parameters: Analysis parameters (thresholds, groups, normalization switches).
        inputs: Role -> input file path.
        outputs: Role -> output file path.
        generated_at: ISO-8601 UTC timestamp.

    Returns:
        The manifest dictionary.
    """

    def describe(files: Mapping[str, Path]) -> dict[str, dict[str, Any]]:
        return {
            role: {
                "path": display_path(path, root),
                "bytes": Path(path).stat().st_size,
                "sha256": file_sha256(path),
            }
            for role, path in sorted(files.items())
        }

    return {
        "data_label": DATA_LABEL,
        "disclaimer": SHORT_DISCLAIMER,
        "generated_at": generated_at,
        "workbench_version": __version__,
        "python": platform.python_version(),
        "platform": platform.system(),
        "packages": package_versions(),
        "git": git_state(root),
        "config": {
            "path": display_path(config_path, root),
            "sha256": file_sha256(config_path),
        },
        "seed": seed,
        "parameters": dict(parameters),
        "inputs": describe(inputs),
        "outputs": describe(outputs),
    }


def write_run_manifest(manifest: Mapping[str, Any], path: str | Path) -> Path:
    """Write a manifest as indented, key-sorted JSON with LF line endings.

    Args:
        manifest: Manifest dictionary.
        path: Destination file.

    Returns:
        The written path.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")
    return path
