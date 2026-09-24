"""Shared pytest fixtures.

Unit tests use small in-memory synthetic datasets so they stay fast and never
depend on generated files on disk (except tests of the committed demo files).
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pandas as pd
import pytest

from onco_workbench.config import SyntheticConfig
from onco_workbench.data.synthetic import SyntheticDataset, generate_synthetic_dataset


@pytest.fixture(scope="session")
def repo_root() -> Path:
    """Absolute path to the repository root."""
    return Path(__file__).resolve().parents[1]


def make_synthetic_config(**overrides: object) -> SyntheticConfig:
    """Small, valid generator configuration; keyword arguments override fields."""
    params: dict[str, object] = {
        "n_samples": 20,
        "n_genes": 40,
        "group_labels": ("Group_A", "Group_B"),
        "group_proportions": (0.5, 0.5),
        "n_batches": 2,
        "batch_shift_sd": 0.3,
        "baseline_mean": 8.0,
        "baseline_sd": 1.5,
        "noise_sd": 0.5,
        "n_signal_genes_up": 4,
        "n_signal_genes_down": 4,
        "signal_effect_range": (1.5, 2.5),
        "missing_fraction": 0.01,
        "gene_id_prefix": "SYN_G",
        "sample_id_prefix": "SYN_S",
        "float_decimals": 3,
    }
    params.update(overrides)
    return SyntheticConfig(**params)  # type: ignore[arg-type]


@pytest.fixture
def config_factory() -> Callable[..., SyntheticConfig]:
    """Return :func:`make_synthetic_config` so tests can build variants."""
    return make_synthetic_config


@pytest.fixture
def small_config() -> SyntheticConfig:
    return make_synthetic_config()


@pytest.fixture
def small_dataset(small_config: SyntheticConfig) -> SyntheticDataset:
    return generate_synthetic_dataset(small_config, seed=123)


@pytest.fixture
def tiny_expression() -> pd.DataFrame:
    """Hand-written 6 x 3 valid expression matrix for validation tests."""
    return pd.DataFrame(
        {
            "SYN_G0001": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            "SYN_G0002": [2.0, 2.5, 3.0, 3.5, 4.0, 4.5],
            "SYN_G0003": [9.0, 8.0, 7.0, 6.0, 5.0, 4.0],
        },
        index=pd.Index([f"S{i}" for i in range(1, 7)], name="sample_id"),
    )


@pytest.fixture
def tiny_metadata() -> pd.DataFrame:
    """Metadata matching :func:`tiny_expression`: 3 samples per group, balanced batches."""
    return pd.DataFrame(
        {
            "sample_id": [f"S{i}" for i in range(1, 7)],
            "group": ["Group_A"] * 3 + ["Group_B"] * 3,
            "batch": ["Batch_1", "Batch_2", "Batch_1", "Batch_2", "Batch_1", "Batch_2"],
        }
    )
