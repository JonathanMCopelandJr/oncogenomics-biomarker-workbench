"""Multiple-testing correction.

Benjamini-Hochberg (BH) controls the expected false discovery rate when the tests
are independent or positively dependent. With ``m`` tested genes and p-values sorted
ascending, the adjusted value at rank ``i`` is ``min over j >= i of (p_(j) * m / j)``,
capped at 1. Genes that were not tested (NaN p-value) are excluded from ``m`` and keep
a NaN adjusted value.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from statsmodels.stats.multitest import multipletests


def benjamini_hochberg(p_values: ArrayLike) -> np.ndarray:
    """Return Benjamini-Hochberg adjusted p-values.

    Args:
        p_values: One-dimensional raw p-values in ``[0, 1]``. NaN marks untested genes.

    Returns:
        Adjusted p-values with the same shape as the input. NaN stays NaN.

    Raises:
        ValueError: If the input is not one-dimensional or contains values outside [0, 1].
    """
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1:
        raise ValueError(f"p_values must be one-dimensional, got shape {p.shape}.")
    tested = ~np.isnan(p)
    if np.any((p[tested] < 0) | (p[tested] > 1)):
        raise ValueError("p-values must lie within [0, 1].")
    adjusted = np.full(p.shape, np.nan)
    if tested.any():
        adjusted[tested] = multipletests(p[tested], method="fdr_bh")[1]
    return adjusted
