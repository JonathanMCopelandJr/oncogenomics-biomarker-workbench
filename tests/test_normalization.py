from __future__ import annotations

import numpy as np
import pandas as pd

from onco_workbench.analysis.normalization import (
    impute_gene_means,
    median_center_samples,
    prepare_for_plots,
    zscore_genes,
)


def _frame() -> pd.DataFrame:
    return pd.DataFrame(
        {"G1": [1.0, 2.0, 3.0], "G2": [4.0, np.nan, 6.0], "G3": [7.0, 8.0, 12.0]},
        index=["S1", "S2", "S3"],
    )


def test_median_centering_aligns_sample_medians() -> None:
    frame = _frame()
    centered = median_center_samples(frame)
    medians = centered.median(axis=1)
    assert np.allclose(medians, medians.iloc[0])
    # Median of original sample medians is kept as the common level.
    assert np.isclose(medians.iloc[0], frame.median(axis=1).median())


def test_median_centering_preserves_missing_and_shape() -> None:
    frame = _frame()
    centered = median_center_samples(frame)
    assert centered.shape == frame.shape
    assert centered.isna().equals(frame.isna())


def test_median_centering_preserves_within_sample_differences() -> None:
    frame = _frame()
    centered = median_center_samples(frame)
    assert np.allclose(
        (centered["G3"] - centered["G1"]).to_numpy(), (frame["G3"] - frame["G1"]).to_numpy()
    )


def test_impute_gene_means() -> None:
    imputed = impute_gene_means(_frame())
    assert imputed.loc["S2", "G2"] == 5.0
    assert not imputed.isna().any().any()


def test_zscore_genes_mean_zero_sd_one() -> None:
    z = zscore_genes(_frame())
    assert np.allclose(z.mean(), 0.0)
    assert np.allclose(z.std(ddof=1), 1.0)
    assert np.isnan(z.loc["S2", "G2"])


def test_zscore_constant_gene_becomes_zero() -> None:
    frame = pd.DataFrame({"G1": [5.0, 5.0, np.nan], "G2": [1.0, 2.0, 3.0]})
    z = zscore_genes(frame)
    assert z["G1"].iloc[:2].tolist() == [0.0, 0.0]
    assert np.isnan(z["G1"].iloc[2])


def test_prepare_for_plots_centers_without_scaling() -> None:
    prepared = prepare_for_plots(_frame(), zscore=False)
    assert np.allclose(prepared.mean(), 0.0)
    assert not prepared.isna().any().any()
    assert not np.allclose(prepared.std(ddof=1), 1.0)


def test_inputs_are_not_modified() -> None:
    frame = _frame()
    before = frame.copy()
    median_center_samples(frame)
    zscore_genes(frame)
    prepare_for_plots(frame)
    pd.testing.assert_frame_equal(frame, before)
