from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from onco_workbench.analysis.ranking import (
    DIRECTION_HIGHER,
    DIRECTION_LOWER,
    DIRECTION_NOT_TESTED,
    add_ranking,
    flag_results,
    ranking_score,
)


def _results() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_id": ["G1", "G2", "G3", "G4", "G5"],
            "mean_diff": [1.0, -2.0, 0.1, 0.5, 0.3],
            "cohens_d": [2.0, -1.0, 0.1, 0.5, np.nan],
            "p_value": [0.001, 0.002, 0.5, 0.01, np.nan],
            "p_adj": [0.01, 0.01, 0.5, 0.05, np.nan],
        }
    )


def test_score_formula() -> None:
    assert ranking_score([2.0], [0.01], 1e-300)[0] == pytest.approx(4.0)
    assert ranking_score([-2.0], [0.01], 1e-300)[0] == pytest.approx(4.0)


def test_score_floor_prevents_infinity() -> None:
    assert ranking_score([1.0], [0.0], 1e-10)[0] == pytest.approx(10.0)


def test_score_nan_propagates() -> None:
    assert np.isnan(ranking_score([np.nan, 1.0], [0.1, np.nan], 1e-300)).all()


def test_score_invalid_floor() -> None:
    with pytest.raises(ValueError, match="p_floor"):
        ranking_score([1.0], [0.1], 0.0)


def test_flags_and_directions() -> None:
    flagged = flag_results(_results(), fdr_threshold=0.05, effect_size_threshold=0.5)
    assert flagged["meets_thresholds"].tolist() == [True, True, False, True, False]
    assert flagged["direction"].tolist() == [
        DIRECTION_HIGHER,
        DIRECTION_LOWER,
        DIRECTION_HIGHER,
        DIRECTION_HIGHER,
        DIRECTION_NOT_TESTED,
    ]


def test_flag_boundaries_are_inclusive() -> None:
    flagged = flag_results(_results(), fdr_threshold=0.05, effect_size_threshold=0.5)
    assert bool(flagged.loc[3, "meets_thresholds"])  # p_adj == 0.05 and |d| == 0.5


def test_flag_does_not_modify_input() -> None:
    results = _results()
    flag_results(results, fdr_threshold=0.05, effect_size_threshold=0.5)
    assert "meets_thresholds" not in results.columns


def test_ranking_order_ties_and_untested() -> None:
    ranked = add_ranking(_results(), p_floor=1e-300)
    # G1 score 2*2=4, G2 1*2=2, G4 0.5*1.30=0.65, G3 0.1*0.30=0.03, G5 untested.
    assert ranked["gene_id"].tolist() == ["G1", "G2", "G4", "G3", "G5"]
    assert ranked["rank"].tolist()[:4] == [1, 2, 3, 4]
    assert pd.isna(ranked["rank"].iloc[4])


def test_ranking_ties_broken_by_gene_id() -> None:
    results = pd.DataFrame({"gene_id": ["B", "A"], "cohens_d": [1.0, 1.0], "p_adj": [0.01, 0.01]})
    assert add_ranking(results, p_floor=1e-300)["gene_id"].tolist() == ["A", "B"]
