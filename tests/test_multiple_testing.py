from __future__ import annotations

import numpy as np
import pytest

from onco_workbench.analysis.multiple_testing import benjamini_hochberg


def test_hand_computed_reference() -> None:
    # Sorted p: 0.005, 0.01, 0.03, 0.04 (m = 4). Raw p * m / rank: 0.02, 0.02, 0.04, 0.04.
    # Cumulative minimum from the largest rank down leaves them unchanged.
    adjusted = benjamini_hochberg([0.01, 0.04, 0.03, 0.005])
    assert np.allclose(adjusted, [0.02, 0.04, 0.04, 0.02])


def test_step_up_monotonicity_example() -> None:
    # m = 3. p * m / rank = 0.01*3/1, 0.0125*3/2, 0.02*3/3 = 0.03, 0.01875, 0.02.
    # Step-up (cumulative minimum from the largest rank down): 0.01875, 0.01875, 0.02.
    # The smallest p inherits the smaller value from rank 2, so adjusted p stays monotone.
    assert np.allclose(benjamini_hochberg([0.01, 0.0125, 0.02]), [0.01875, 0.01875, 0.02])


def test_adjusted_values_are_capped_at_one() -> None:
    assert np.allclose(benjamini_hochberg([0.9, 0.95, 1.0]), [1.0, 1.0, 1.0])


def test_properties_on_random_pvalues() -> None:
    rng = np.random.default_rng(3)
    p = rng.uniform(size=200)
    adjusted = benjamini_hochberg(p)
    assert np.all(adjusted >= p - 1e-15)
    assert np.all(adjusted <= 1.0)
    order = np.argsort(p)
    assert np.all(np.diff(adjusted[order]) >= -1e-15)  # ranking by p is preserved


def test_nan_values_are_excluded_from_m() -> None:
    adjusted = benjamini_hochberg([0.01, np.nan, 0.04])
    assert np.isnan(adjusted[1])
    # With m = 2 (not 3): 0.01 * 2 / 1 = 0.02, 0.04 * 2 / 2 = 0.04.
    assert np.allclose(adjusted[[0, 2]], [0.02, 0.04])


def test_all_nan_and_empty_inputs() -> None:
    assert np.isnan(benjamini_hochberg([np.nan, np.nan])).all()
    assert benjamini_hochberg([]).size == 0


def test_single_pvalue_is_unchanged() -> None:
    assert benjamini_hochberg([0.03])[0] == pytest.approx(0.03)


@pytest.mark.parametrize("bad", [[-0.1, 0.5], [0.2, 1.5]])
def test_out_of_range_raises(bad: list[float]) -> None:
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        benjamini_hochberg(bad)


def test_two_dimensional_input_raises() -> None:
    with pytest.raises(ValueError, match="one-dimensional"):
        benjamini_hochberg([[0.1, 0.2]])
