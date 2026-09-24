"""Transparent normalization helpers for SYNTHETIC continuous, log2-like values.

These steps are deliberately simple so every transformation can be read and checked.
They are **not** appropriate for real RNA-seq counts, which need count-aware methods
(library-size normalization and variance modelling) that are out of scope for this
educational project.

- :func:`median_center_samples` is used before the group comparison.
- :func:`prepare_for_plots` (mean imputation plus per-gene centering or z-scoring) is
  used **only** for PCA and heatmaps, never for statistical testing.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def median_center_samples(expression: pd.DataFrame) -> pd.DataFrame:
    """Remove per-sample offsets by aligning every sample's median.

    Each value becomes ``x - median(sample) + median(all sample medians)``. Adding back
    the median of medians keeps values on their original, interpretable scale.
    Missing values stay missing and are ignored when computing medians.

    Args:
        expression: Samples as rows, genes as columns.

    Returns:
        A new DataFrame with the same shape, index, and columns.
    """
    sample_medians = expression.median(axis=1, skipna=True)
    reference = float(sample_medians.median())
    return expression.sub(sample_medians, axis=0) + reference


def impute_gene_means(expression: pd.DataFrame) -> pd.DataFrame:
    """Replace missing values with the mean of their gene (plots only).

    Args:
        expression: Samples as rows, genes as columns.

    Returns:
        A new DataFrame without missing values, except in genes that have no observed
        values at all.
    """
    return expression.fillna(expression.mean(axis=0, skipna=True))


def zscore_genes(expression: pd.DataFrame) -> pd.DataFrame:
    """Standardize every gene to mean 0 and sample standard deviation 1.

    Genes with zero variance become 0 (they carry no information). Missing values
    stay missing.

    Args:
        expression: Samples as rows, genes as columns.

    Returns:
        A new DataFrame of z-scores.
    """
    means = expression.mean(axis=0, skipna=True)
    sds = expression.std(axis=0, skipna=True, ddof=1)
    constant = ~(sds > 0)
    z = expression.sub(means, axis=1).div(sds.where(~constant), axis=1)
    if constant.any():
        observed = expression.loc[:, constant].notna()
        z.loc[:, constant] = np.where(observed, 0.0, np.nan)
    return z


def prepare_for_plots(expression: pd.DataFrame, *, zscore: bool = True) -> pd.DataFrame:
    """Prepare values for PCA or heatmaps: impute gene means, then z-score or center.

    Args:
        expression: Samples as rows, genes as columns.
        zscore: If True, standardize each gene. If False, only subtract each gene's mean.

    Returns:
        A new DataFrame with no missing values in genes that have observations.
    """
    imputed = impute_gene_means(expression)
    if zscore:
        return zscore_genes(imputed)
    return imputed.sub(imputed.mean(axis=0), axis=1)
