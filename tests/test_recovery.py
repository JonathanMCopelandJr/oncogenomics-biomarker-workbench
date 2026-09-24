from __future__ import annotations

import math

import pandas as pd
import pytest

from onco_workbench.analysis.recovery import compare_with_ground_truth


def _results() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_id": ["G1", "G2", "G3", "G4", "G5"],
            "meets_thresholds": [True, True, False, True, False],
            "direction": ["higher_in_b", "higher_in_b", "lower_in_b", "lower_in_b", "lower_in_b"],
        }
    )


def _truth() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_id": ["G1", "G2", "G3", "G4", "G5"],
            "is_signal": [True, True, True, False, False],
            "direction": ["up", "down", "down", "none", "none"],
        }
    )


def test_counts() -> None:
    summary = compare_with_ground_truth(_results(), _truth())
    assert summary.n_planted == 3
    assert summary.n_flagged == 3
    assert summary.true_positives == 2
    assert summary.false_negatives == 1
    assert summary.false_positives == 1
    assert summary.direction_matches == 1  # G2 was planted "down" but observed higher
    assert summary.sensitivity == pytest.approx(2 / 3)
    assert summary.false_discovery_proportion == pytest.approx(1 / 3)


def test_rates_are_nan_when_undefined() -> None:
    results = _results().assign(meets_thresholds=False)
    truth = _truth().assign(is_signal=False)
    summary = compare_with_ground_truth(results, truth)
    assert math.isnan(summary.sensitivity)
    assert math.isnan(summary.false_discovery_proportion)
    assert "n/a" in dict(summary.as_rows())["Planted genes recovered"]


def test_mismatched_genes_raise() -> None:
    with pytest.raises(ValueError, match="do not match"):
        compare_with_ground_truth(_results(), _truth().iloc[:3])
