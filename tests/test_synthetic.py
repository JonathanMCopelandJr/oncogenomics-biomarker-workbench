from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd
import pytest

from onco_workbench.config import SyntheticConfig, load_config
from onco_workbench.data.synthetic import SyntheticDataset, generate_synthetic_dataset, make_ids
from onco_workbench.data.validation import validate_dataset

Factory = Callable[..., SyntheticConfig]


def test_shapes_and_identifiers(small_dataset: SyntheticDataset) -> None:
    expr = small_dataset.expression
    assert expr.shape == (20, 40)
    assert expr.index.name == "sample_id"
    assert list(expr.index) == list(small_dataset.metadata["sample_id"])
    assert expr.index.is_unique and expr.columns.is_unique
    assert list(small_dataset.ground_truth["gene_id"]) == list(expr.columns)


def test_gene_ids_are_neutral_synthetic_identifiers(small_dataset: SyntheticDataset) -> None:
    assert small_dataset.expression.columns.str.fullmatch(r"SYN_G\d{4}").all()
    assert small_dataset.metadata["sample_id"].str.fullmatch(r"SYN_S\d{3}").all()


def test_make_ids_widens_for_large_counts() -> None:
    assert make_ids("X", 3, min_width=4) == ["X0001", "X0002", "X0003"]
    assert make_ids("X", 12345, min_width=4)[-1] == "X12345"


def test_same_seed_is_identical(small_config: SyntheticConfig) -> None:
    first = generate_synthetic_dataset(small_config, seed=5)
    second = generate_synthetic_dataset(small_config, seed=5)
    pd.testing.assert_frame_equal(first.expression, second.expression)
    pd.testing.assert_frame_equal(first.metadata, second.metadata)
    pd.testing.assert_frame_equal(first.ground_truth, second.ground_truth)


def test_different_seed_differs(small_config: SyntheticConfig) -> None:
    first = generate_synthetic_dataset(small_config, seed=5)
    second = generate_synthetic_dataset(small_config, seed=6)
    assert not first.expression.equals(second.expression)


def test_group_sizes_and_labels(small_dataset: SyntheticDataset) -> None:
    counts = small_dataset.metadata["group"].value_counts()
    assert counts.to_dict() == {"Group_A": 10, "Group_B": 10}


def test_unequal_group_proportions(config_factory: Factory) -> None:
    config = config_factory(group_proportions=(0.25, 0.75))
    dataset = generate_synthetic_dataset(config, seed=1)
    assert dataset.metadata["group"].value_counts().to_dict() == {"Group_A": 5, "Group_B": 15}


def test_batches_are_balanced_within_each_group(small_dataset: SyntheticDataset) -> None:
    table = pd.crosstab(small_dataset.metadata["group"], small_dataset.metadata["batch"])
    assert (table.to_numpy() == 5).all()


def test_ground_truth_signal_counts(small_dataset: SyntheticDataset) -> None:
    truth = small_dataset.ground_truth
    assert truth["is_signal"].sum() == 8
    assert (truth["direction"] == "up").sum() == 4
    assert (truth["direction"] == "down").sum() == 4
    signal_shift = truth.loc[truth["is_signal"], "true_shift"].abs()
    assert signal_shift.between(1.5, 2.5).all()
    assert (truth.loc[~truth["is_signal"], "true_shift"] == 0).all()
    assert (truth.loc[~truth["is_signal"], "direction"] == "none").all()


def test_planted_shifts_are_visible_in_group_means(small_dataset: SyntheticDataset) -> None:
    """Generator sanity check: the planted shift direction appears in the raw group means."""
    expr, meta, truth = small_dataset.expression, small_dataset.metadata, small_dataset.ground_truth
    in_b = (meta["group"] == "Group_B").to_numpy()
    observed = expr[in_b].mean() - expr[~in_b].mean()
    signal = truth.set_index("gene_id")
    for gene_id, row in signal[signal["is_signal"]].iterrows():
        assert np.sign(observed[gene_id]) == np.sign(row["true_shift"])


def test_values_are_rounded(small_dataset: SyntheticDataset) -> None:
    values = small_dataset.expression.to_numpy()
    finite = values[~np.isnan(values)]
    assert np.allclose(finite, np.round(finite, 3))


def test_missing_fraction_close_to_config() -> None:
    config = load_config()
    dataset = generate_synthetic_dataset(config.synthetic, config.seed)
    observed = float(dataset.expression.isna().to_numpy().mean())
    assert 0.002 < observed < 0.008  # configured 0.005 over 40,000 cells


def test_zero_missing_fraction_gives_complete_matrix(config_factory: Factory) -> None:
    dataset = generate_synthetic_dataset(config_factory(missing_fraction=0.0), seed=3)
    assert not dataset.expression.isna().any().any()


def test_no_signal_genes_is_supported(config_factory: Factory) -> None:
    config = config_factory(n_signal_genes_up=0, n_signal_genes_down=0)
    dataset = generate_synthetic_dataset(config, seed=3)
    assert not dataset.ground_truth["is_signal"].any()


def test_generated_default_dataset_passes_validation() -> None:
    config = load_config()
    dataset = generate_synthetic_dataset(config.synthetic, config.seed)
    report = validate_dataset(
        dataset.expression,
        dataset.metadata,
        allowed_groups=config.validation.allowed_groups,
        min_samples_per_group=config.validation.min_samples_per_group,
    )
    assert report.is_valid, report.summary()
    assert not report.warnings, report.summary()


def test_negative_seed_rejected(small_config: SyntheticConfig) -> None:
    with pytest.raises(ValueError, match="non-negative"):
        generate_synthetic_dataset(small_config, seed=-1)
