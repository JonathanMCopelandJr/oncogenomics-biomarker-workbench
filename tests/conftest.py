"""Shared pytest fixtures.

Fixtures for small, in-memory synthetic datasets are added in Phase 2 so that
unit tests stay fast and never depend on generated files on disk.
"""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def repo_root() -> Path:
    """Absolute path to the repository root."""
    return Path(__file__).resolve().parents[1]
