from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from onco_workbench.analysis.qc import (
    compute_qc_summary,
    gene_qc_table,
    run_pca,
    sample_correlation,
    sample_qc_table,
    select_heatmap_samples,
)
from onco_workbench.data.synthetic import SyntheticDataset


def test_qc_summary_counts(tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame) -> None:
    tiny_expression.iloc[0, 0] = np.nan
    tiny_expression.iloc[1, 0] = np.nan
    summary = compute_qc_summary(tiny_expression, tiny_metadata)
    assert (summary.n_samples, summary.n_genes, summary.n_values) == (6, 3, 18)
    assert summary.n_missing == 2
    assert summary.missing_fraction == pytest.approx(2 / 18)
    assert summary.genes_with_missing == 1
    assert summary.samples_with_missing == 2
    assert summary.group_counts == {"Group_A": 3, "Group_B": 3}
    assert summary.batch_counts == {"Batch_1": 3, "Batch_2": 3}
    labels = [label for label, _ in summary.as_rows()]
    assert "Samples per batch" in labels


def test_qc_summary_without_batch(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    summary = compute_qc_summary(tiny_expression, tiny_metadata.drop(columns="batch"))
    assert summary.batch_counts is None


def test_sample_qc_table(tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame) -> None:
    tiny_expression.iloc[0, 1] = np.nan
    table = sample_qc_table(tiny_expression, tiny_metadata)
    assert list(table["sample_id"]) == list(tiny_expression.index)
    assert table.loc[0, "n_missing"] == 1
    assert table.loc[0, "total_expression"] == pytest.approx(1.0 + 9.0)
    assert table.loc[1, "median"] == pytest.approx(2.5)
    assert {"group", "batch", "mean", "sd"} <= set(table.columns)


def test_sample_qc_table_follows_expression_order(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    shuffled = tiny_metadata.iloc[::-1].reset_index(drop=True)
    table = sample_qc_table(tiny_expression, shuffled)
    assert list(table["group"]) == ["Group_A"] * 3 + ["Group_B"] * 3


def test_gene_qc_table(tiny_expression: pd.DataFrame) -> None:
    tiny_expression.iloc[0:3, 2] = np.nan
    table = gene_qc_table(tiny_expression)
    assert list(table["gene_id"]) == list(tiny_expression.columns)
    row = table.set_index("gene_id").loc["SYN_G0003"]
    assert row["n_missing"] == 3
    assert row["missing_fraction"] == pytest.approx(0.5)
    assert row["mean"] == pytest.approx(5.0)


def test_pca_shapes_and_variance(small_dataset: SyntheticDataset) -> None:
    result = run_pca(small_dataset.expression, small_dataset.metadata, n_components=3)
    assert list(result.scores.columns) == ["sample_id", "group", "batch", "PC1", "PC2", "PC3"]
    assert len(result.scores) == 20
    ratios = result.explained_variance_ratio
    assert np.all(np.diff(ratios) <= 1e-12)
    assert 0 < ratios.sum() <= 1
    assert result.n_genes_used == 40
    assert result.standardized


def test_pca_is_deterministic(small_dataset: SyntheticDataset) -> None:
    first = run_pca(small_dataset.expression, small_dataset.metadata)
    second = run_pca(small_dataset.expression, small_dataset.metadata)
    pd.testing.assert_frame_equal(first.scores, second.scores)


def test_pca_components_capped_by_data(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    result = run_pca(tiny_expression, tiny_metadata, n_components=10)
    assert result.explained_variance_ratio.size <= 3


def test_pca_skips_constant_genes(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    tiny_expression["SYN_G0004"] = 1.0
    assert run_pca(tiny_expression, tiny_metadata).n_genes_used == 3


def test_select_heatmap_samples_balanced(small_dataset: SyntheticDataset) -> None:
    samples = select_heatmap_samples(small_dataset.expression, small_dataset.metadata, 8)
    groups = small_dataset.metadata.set_index("sample_id").loc[samples, "group"]
    assert groups.value_counts().to_dict() == {"Group_A": 4, "Group_B": 4}
    assert list(groups) == ["Group_A"] * 4 + ["Group_B"] * 4
    assert samples == select_heatmap_samples(small_dataset.expression, small_dataset.metadata, 8)


def test_select_heatmap_samples_returns_all_when_small(small_dataset: SyntheticDataset) -> None:
    samples = select_heatmap_samples(small_dataset.expression, small_dataset.metadata, 100)
    assert sorted(samples) == sorted(small_dataset.expression.index)


def test_sample_correlation_properties(small_dataset: SyntheticDataset) -> None:
    samples = list(small_dataset.expression.index[:6])
    corr = sample_correlation(small_dataset.expression, samples)
    assert list(corr.index) == samples and list(corr.columns) == samples
    assert np.allclose(np.diag(corr), 1.0)
    assert np.allclose(corr.to_numpy(), corr.to_numpy().T)
    assert (corr.abs() <= 1 + 1e-12).all().all()
