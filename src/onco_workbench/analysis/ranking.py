"""Transparent flagging and ranking of group-comparison results.

The ranking score is a documented **sorting heuristic**, not a statistical test:

    ranking_score = |cohens_d| * -log10(max(p_adj, p_floor))

It rewards genes that have both a large standardized difference and strong evidence
after multiple-testing correction. ``p_floor`` prevents infinite scores when an
adjusted p-value underflows to 0.

``meets_thresholds`` marks genes with ``p_adj <= fdr_threshold`` **and**
``|cohens_d| >= effect_size_threshold``. Both thresholds are user-chosen display
settings, not scientific conclusions.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike

DIRECTION_HIGHER = "higher_in_b"
DIRECTION_LOWER = "lower_in_b"
DIRECTION_NONE = "no_difference"
DIRECTION_NOT_TESTED = "not_tested"


def ranking_score(cohens_d: ArrayLike, p_adj: ArrayLike, p_floor: float) -> np.ndarray:
    """Compute ``|d| * -log10(max(p_adj, p_floor))``.

    Args:
        cohens_d: Effect sizes.
        p_adj: Adjusted p-values.
        p_floor: Smallest p-value used in the logarithm, in (0, 1).

    Returns:
        Scores. NaN where either input is NaN.

    Raises:
        ValueError: If ``p_floor`` is not in (0, 1).
    """
    if not 0.0 < p_floor < 1.0:
        raise ValueError(f"p_floor must be in (0, 1), got {p_floor}.")
    d = np.asarray(cohens_d, dtype=float)
    p = np.asarray(p_adj, dtype=float)
    with np.errstate(invalid="ignore"):
        return np.abs(d) * -np.log10(np.maximum(p, p_floor))


def flag_results(
    results: pd.DataFrame, *, fdr_threshold: float, effect_size_threshold: float
) -> pd.DataFrame:
    """Add ``direction`` and ``meets_thresholds`` columns.

    Args:
        results: Output of the group comparison with ``mean_diff``, ``cohens_d``,
            ``p_value``, and ``p_adj``.
        fdr_threshold: Maximum adjusted p-value.
        effect_size_threshold: Minimum absolute Cohen's d.

    Returns:
        A copy of ``results`` with the new columns.
    """
    out = results.copy()
    tested = out["p_value"].notna().to_numpy()
    diff = out["mean_diff"].to_numpy(dtype=float)
    direction = np.where(
        diff > 0, DIRECTION_HIGHER, np.where(diff < 0, DIRECTION_LOWER, DIRECTION_NONE)
    )
    out["direction"] = np.where(tested, direction, DIRECTION_NOT_TESTED)
    p_adj = out["p_adj"].to_numpy(dtype=float)
    d = out["cohens_d"].to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        meets = (p_adj <= fdr_threshold) & (np.abs(d) >= effect_size_threshold)
    out["meets_thresholds"] = meets & tested
    return out


def add_ranking(results: pd.DataFrame, *, p_floor: float) -> pd.DataFrame:
    """Add ``ranking_score`` and ``rank`` columns, sorted best first.

    Ties are broken by ``gene_id`` so the order is deterministic. Untested genes get
    a missing score and rank and are placed last.

    Args:
        results: Group-comparison results with ``cohens_d`` and ``p_adj``.
        p_floor: Floor applied to adjusted p-values inside the logarithm.

    Returns:
        A new, re-sorted DataFrame with a fresh index.
    """
    out = results.copy()
    out["ranking_score"] = ranking_score(out["cohens_d"], out["p_adj"], p_floor)
    out = out.sort_values(
        ["ranking_score", "gene_id"], ascending=[False, True], na_position="last", kind="stable"
    ).reset_index(drop=True)
    scored = out["ranking_score"].notna()
    rank = pd.Series(pd.NA, index=out.index, dtype="Int64")
    rank[scored] = np.arange(1, int(scored.sum()) + 1)
    out["rank"] = rank
    return out
