"""Per-gene two-group comparison: group means, mean difference, Cohen's d, and a Welch t-test.

Educational, transparent statistics for SYNTHETIC continuous values. Every formula is
written out below so it can be checked by hand. Direction convention: every
difference is **group B minus group A**.

- ``mean_diff = mean_b - mean_a``
- ``cohens_d = mean_diff / pooled_sd`` with
  ``pooled_sd = sqrt(((n_a - 1) * var_a + (n_b - 1) * var_b) / (n_a + n_b - 2))``
- Welch's t-test (unequal variances): ``se = sqrt(var_a / n_a + var_b / n_b)``,
  ``t = mean_diff / se``, and Welch-Satterthwaite degrees of freedom
  ``df = se**4 / ((var_a / n_a)**2 / (n_a - 1) + (var_b / n_b)**2 / (n_b - 1))``.
  The two-sided ``p = 2 * t_sf(|t|, df)``.

Missing values are ignored gene by gene. A gene is tested only if both groups have at
least ``min_non_missing_per_group`` observed values; otherwise its statistics are NaN.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from onco_workbench.data.synthetic import GROUP_COLUMN
from onco_workbench.data.validation import align_metadata

RESULT_COLUMNS: tuple[str, ...] = (
    "gene_id",
    "n_a",
    "n_b",
    "mean_a",
    "mean_b",
    "mean_diff",
    "sd_a",
    "sd_b",
    "cohens_d",
    "t_statistic",
    "df",
    "p_value",
)


def _group_moments(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return per-column observed count, mean, and sample variance (ddof=1)."""
    observed = ~np.isnan(values)
    n = observed.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(observed, values, 0.0).sum(axis=0) / n
        squared = np.where(observed, (values - mean) ** 2, 0.0).sum(axis=0)
        var = squared / (n - 1)
    return n, mean, var


def welch_t_test(
    mean_a: np.ndarray,
    mean_b: np.ndarray,
    var_a: np.ndarray,
    var_b: np.ndarray,
    n_a: np.ndarray,
    n_b: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized Welch two-sample t-test of B versus A.

    Args:
        mean_a: Group A means.
        mean_b: Group B means.
        var_a: Group A sample variances (ddof=1).
        var_b: Group B sample variances (ddof=1).
        n_a: Group A observed counts.
        n_b: Group B observed counts.

    Returns:
        ``(t_statistic, df, p_value)``. Entries are NaN where the test is undefined
        (for example zero variance in both groups).
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        term_a = var_a / n_a
        term_b = var_b / n_b
        se = np.sqrt(term_a + term_b)
        t = (mean_b - mean_a) / se
        df = (term_a + term_b) ** 2 / (term_a**2 / (n_a - 1) + term_b**2 / (n_b - 1))
        p = 2.0 * stats.t.sf(np.abs(t), df)
    undefined = ~np.isfinite(t) | ~np.isfinite(df) | (se == 0)
    t = np.where(undefined, np.nan, t)
    df = np.where(undefined, np.nan, df)
    p = np.where(undefined, np.nan, p)
    return t, df, p


def cohens_d(
    mean_a: np.ndarray,
    mean_b: np.ndarray,
    var_a: np.ndarray,
    var_b: np.ndarray,
    n_a: np.ndarray,
    n_b: np.ndarray,
) -> np.ndarray:
    """Cohen's d of B versus A using the pooled standard deviation.

    Args:
        mean_a: Group A means.
        mean_b: Group B means.
        var_a: Group A sample variances (ddof=1).
        var_b: Group B sample variances (ddof=1).
        n_a: Group A observed counts.
        n_b: Group B observed counts.

    Returns:
        Effect sizes. NaN where the pooled SD is zero or undefined.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        pooled = np.sqrt(((n_a - 1) * var_a + (n_b - 1) * var_b) / (n_a + n_b - 2))
        d = (mean_b - mean_a) / pooled
    return np.where(np.isfinite(d) & (pooled > 0), d, np.nan)


def compare_groups(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    group_a: str,
    group_b: str,
    *,
    min_non_missing_per_group: int = 3,
) -> pd.DataFrame:
    """Compute per-gene statistics comparing ``group_b`` against ``group_a``.

    Args:
        expression: Samples as rows, genes as columns (numeric).
        metadata: Sample metadata with ``sample_id`` and ``group``.
        group_a: Reference group label.
        group_b: Comparison group label.
        min_non_missing_per_group: Minimum observed values per group for a gene to be tested.

    Returns:
        One row per gene with the columns in :data:`RESULT_COLUMNS`, in gene order.

    Raises:
        ValueError: If the groups are identical, or either group has no samples.
    """
    if group_a == group_b:
        raise ValueError(f"group_a and group_b must differ, both are {group_a!r}.")
    if min_non_missing_per_group < 2:
        raise ValueError("min_non_missing_per_group must be at least 2.")
    labels = align_metadata(expression, metadata)[GROUP_COLUMN]
    available = sorted(labels.dropna().astype(str).unique())
    for group in (group_a, group_b):
        if group not in available:
            raise ValueError(
                f"Group {group!r} has no samples. Available groups: {', '.join(available)}."
            )

    in_a = (labels == group_a).to_numpy()
    in_b = (labels == group_b).to_numpy()
    values = expression.to_numpy(dtype=float)
    n_a, mean_a, var_a = _group_moments(values[in_a])
    n_b, mean_b, var_b = _group_moments(values[in_b])

    tested = (n_a >= min_non_missing_per_group) & (n_b >= min_non_missing_per_group)
    t, df, p = welch_t_test(mean_a, mean_b, var_a, var_b, n_a, n_b)
    d = cohens_d(mean_a, mean_b, var_a, var_b, n_a, n_b)

    def only_tested(array: np.ndarray) -> np.ndarray:
        return np.where(tested, array, np.nan)

    return pd.DataFrame(
        {
            "gene_id": expression.columns.to_list(),
            "n_a": n_a,
            "n_b": n_b,
            "mean_a": mean_a,
            "mean_b": mean_b,
            "mean_diff": mean_b - mean_a,
            "sd_a": np.sqrt(var_a),
            "sd_b": np.sqrt(var_b),
            "cohens_d": only_tested(d),
            "t_statistic": only_tested(t),
            "df": only_tested(df),
            "p_value": only_tested(p),
        },
        columns=list(RESULT_COLUMNS),
    )
