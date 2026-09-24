from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from onco_workbench.analysis.differential import RESULT_COLUMNS, compare_groups, welch_t_test


def _two_group_data(a: np.ndarray, b: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build expression/metadata frames; a and b are (samples x genes) arrays."""
    values = np.vstack([a, b])
    ids = [f"S{i}" for i in range(values.shape[0])]
    expression = pd.DataFrame(
        values,
        index=pd.Index(ids, name="sample_id"),
        columns=[f"G{j}" for j in range(values.shape[1])],
    )
    metadata = pd.DataFrame(
        {"sample_id": ids, "group": ["Group_A"] * len(a) + ["Group_B"] * len(b)}
    )
    return expression, metadata


def test_hand_computed_example() -> None:
    # A = [1, 2, 3], B = [4, 5, 6]: means 2 and 5, variances 1 and 1.
    expression, metadata = _two_group_data(
        np.array([[1.0], [2.0], [3.0]]), np.array([[4.0], [5.0], [6.0]])
    )
    row = compare_groups(expression, metadata, "Group_A", "Group_B").iloc[0]
    assert row["mean_a"] == pytest.approx(2.0)
    assert row["mean_b"] == pytest.approx(5.0)
    assert row["mean_diff"] == pytest.approx(3.0)
    assert row["sd_a"] == pytest.approx(1.0)
    assert row["cohens_d"] == pytest.approx(3.0)  # 3 / pooled SD of 1
    assert row["t_statistic"] == pytest.approx(3.0 / np.sqrt(2 / 3))
    assert row["df"] == pytest.approx(4.0)


def test_matches_scipy_welch_on_complete_data() -> None:
    rng = np.random.default_rng(0)
    a = rng.normal(0, 1, size=(8, 25))
    b = rng.normal(0.5, 2, size=(11, 25))
    expression, metadata = _two_group_data(a, b)
    results = compare_groups(expression, metadata, "Group_A", "Group_B")
    reference = stats.ttest_ind(b, a, equal_var=False, axis=0)
    assert np.allclose(results["t_statistic"], reference.statistic)
    assert np.allclose(results["p_value"], reference.pvalue)


def test_matches_scipy_with_missing_values() -> None:
    rng = np.random.default_rng(1)
    a = rng.normal(0, 1, size=(10, 6))
    b = rng.normal(1, 1, size=(10, 6))
    a[0, 0] = a[3, 2] = np.nan
    b[5, 0] = b[1, 4] = b[2, 4] = np.nan
    expression, metadata = _two_group_data(a, b)
    results = compare_groups(expression, metadata, "Group_A", "Group_B")
    for j in range(6):
        col_a, col_b = a[:, j], b[:, j]
        ref = stats.ttest_ind(col_b[~np.isnan(col_b)], col_a[~np.isnan(col_a)], equal_var=False)
        assert results.loc[j, "p_value"] == pytest.approx(ref.pvalue)
        assert results.loc[j, "n_a"] == (~np.isnan(col_a)).sum()
        assert results.loc[j, "n_b"] == (~np.isnan(col_b)).sum()


def test_direction_is_b_minus_a() -> None:
    expression, metadata = _two_group_data(
        np.array([[10.0], [11.0], [12.0]]), np.array([[1.0], [2.0], [3.0]])
    )
    row = compare_groups(expression, metadata, "Group_A", "Group_B").iloc[0]
    assert row["mean_diff"] < 0 and row["cohens_d"] < 0 and row["t_statistic"] < 0
    swapped = compare_groups(expression, metadata, "Group_B", "Group_A").iloc[0]
    assert swapped["mean_diff"] == pytest.approx(-row["mean_diff"])
    assert swapped["p_value"] == pytest.approx(row["p_value"])


def test_genes_with_too_few_observations_are_not_tested() -> None:
    a = np.array([[1.0, 1.0], [2.0, np.nan], [3.0, np.nan]])
    b = np.array([[4.0, 5.0], [5.0, 6.0], [6.0, 7.0]])
    expression, metadata = _two_group_data(a, b)
    results = compare_groups(
        expression, metadata, "Group_A", "Group_B", min_non_missing_per_group=3
    )
    assert results.loc[1, "n_a"] == 1
    for column in ("cohens_d", "t_statistic", "df", "p_value"):
        assert np.isnan(results.loc[1, column])
    assert not np.isnan(results.loc[0, "p_value"])


def test_zero_variance_gives_nan_statistics() -> None:
    expression, metadata = _two_group_data(np.full((3, 1), 2.0), np.full((3, 1), 2.0))
    row = compare_groups(expression, metadata, "Group_A", "Group_B").iloc[0]
    assert row["mean_diff"] == 0
    assert np.isnan(row["p_value"]) and np.isnan(row["cohens_d"])


def test_result_columns_and_order(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    results = compare_groups(tiny_expression, tiny_metadata, "Group_A", "Group_B")
    assert tuple(results.columns) == RESULT_COLUMNS
    assert list(results["gene_id"]) == list(tiny_expression.columns)


def test_invalid_groups_raise(tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame) -> None:
    with pytest.raises(ValueError, match="must differ"):
        compare_groups(tiny_expression, tiny_metadata, "Group_A", "Group_A")
    with pytest.raises(ValueError, match="Group_C.*Available groups: Group_A, Group_B"):
        compare_groups(tiny_expression, tiny_metadata, "Group_A", "Group_C")


def test_missing_metadata_row_raises(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    with pytest.raises(ValueError, match="no metadata row"):
        compare_groups(tiny_expression, tiny_metadata.iloc[1:], "Group_A", "Group_B")


def test_welch_undefined_cases_are_nan() -> None:
    zeros = np.zeros(1)
    t, df, p = welch_t_test(zeros, zeros, zeros, zeros, np.full(1, 3), np.full(1, 3))
    assert np.isnan(t).all() and np.isnan(df).all() and np.isnan(p).all()
