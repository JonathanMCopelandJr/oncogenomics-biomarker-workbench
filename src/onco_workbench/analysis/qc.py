"""Quality-control summaries for SYNTHETIC input data.

Covers counts, missingness, per-sample and per-gene statistics, PCA, and
sample-to-sample correlation. These are data-quality checks, not findings.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from onco_workbench.analysis.normalization import prepare_for_plots
from onco_workbench.data.synthetic import BATCH_COLUMN, GROUP_COLUMN, SAMPLE_ID_COLUMN
from onco_workbench.data.validation import align_metadata


@dataclass(frozen=True)
class QCSummary:
    """Dataset-level QC counts."""

    n_samples: int
    n_genes: int
    n_values: int
    n_missing: int
    missing_fraction: float
    genes_with_missing: int
    samples_with_missing: int
    group_counts: dict[str, int]
    batch_counts: dict[str, int] | None

    def as_rows(self) -> list[tuple[str, str]]:
        """Return ``(label, value)`` pairs for tables in reports and dashboards."""
        rows = [
            ("Samples", f"{self.n_samples:,}"),
            ("Genes", f"{self.n_genes:,}"),
            (
                "Missing values",
                f"{self.n_missing:,} of {self.n_values:,} ({self.missing_fraction:.2%})",
            ),
            ("Genes with any missing value", f"{self.genes_with_missing:,}"),
            ("Samples with any missing value", f"{self.samples_with_missing:,}"),
            ("Samples per group", ", ".join(f"{k}: {v}" for k, v in self.group_counts.items())),
        ]
        if self.batch_counts is not None:
            rows.append(
                ("Samples per batch", ", ".join(f"{k}: {v}" for k, v in self.batch_counts.items()))
            )
        return rows


@dataclass(frozen=True)
class PCAResult:
    """Principal component scores and explained variance.

    Attributes:
        scores: One row per sample with ``sample_id``, ``group``, optional ``batch``,
            and ``PC1``, ``PC2``, ... columns.
        explained_variance_ratio: Fraction of variance explained by each component.
        n_genes_used: Genes with non-zero variance that entered the PCA.
        standardized: Whether genes were z-scored (True) or only centered (False).
    """

    scores: pd.DataFrame
    explained_variance_ratio: np.ndarray
    n_genes_used: int
    standardized: bool


def _counts(series: pd.Series) -> dict[str, int]:
    return {str(k): int(v) for k, v in series.value_counts().sort_index().items()}


def compute_qc_summary(expression: pd.DataFrame, metadata: pd.DataFrame) -> QCSummary:
    """Count samples, genes, missing values, and samples per group and batch.

    Args:
        expression: Samples as rows, genes as columns.
        metadata: Sample metadata with ``sample_id`` and ``group`` columns.

    Returns:
        A :class:`QCSummary`.
    """
    aligned = align_metadata(expression, metadata)
    missing = expression.isna()
    n_values = int(expression.size)
    n_missing = int(missing.to_numpy().sum())
    return QCSummary(
        n_samples=int(expression.shape[0]),
        n_genes=int(expression.shape[1]),
        n_values=n_values,
        n_missing=n_missing,
        missing_fraction=n_missing / n_values if n_values else 0.0,
        genes_with_missing=int(missing.any(axis=0).sum()),
        samples_with_missing=int(missing.any(axis=1).sum()),
        group_counts=_counts(aligned[GROUP_COLUMN]),
        batch_counts=_counts(aligned[BATCH_COLUMN]) if BATCH_COLUMN in aligned.columns else None,
    )


def sample_qc_table(expression: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Per-sample summaries: missing count, total (sum), mean, median, and SD.

    ``total_expression`` sums the observed log2-like values of a sample. It is a simple
    sample-level summary for spotting outliers, not a library size.

    Args:
        expression: Samples as rows, genes as columns.
        metadata: Sample metadata.

    Returns:
        One row per sample, in expression-row order.
    """
    aligned = align_metadata(expression, metadata)
    table = pd.DataFrame(
        {
            SAMPLE_ID_COLUMN: expression.index.to_list(),
            GROUP_COLUMN: aligned[GROUP_COLUMN].to_numpy(),
        }
    )
    if BATCH_COLUMN in aligned.columns:
        table[BATCH_COLUMN] = aligned[BATCH_COLUMN].to_numpy()
    table["n_missing"] = expression.isna().sum(axis=1).to_numpy()
    table["total_expression"] = expression.sum(axis=1, skipna=True).to_numpy()
    table["mean"] = expression.mean(axis=1, skipna=True).to_numpy()
    table["median"] = expression.median(axis=1, skipna=True).to_numpy()
    table["sd"] = expression.std(axis=1, skipna=True, ddof=1).to_numpy()
    return table


def gene_qc_table(expression: pd.DataFrame) -> pd.DataFrame:
    """Per-gene mean, SD, and missingness.

    Args:
        expression: Samples as rows, genes as columns.

    Returns:
        One row per gene, in column order.
    """
    n_missing = expression.isna().sum(axis=0)
    return pd.DataFrame(
        {
            "gene_id": expression.columns.to_list(),
            "mean": expression.mean(axis=0, skipna=True).to_numpy(),
            "sd": expression.std(axis=0, skipna=True, ddof=1).to_numpy(),
            "n_missing": n_missing.to_numpy(),
            "missing_fraction": (n_missing / expression.shape[0]).to_numpy(),
        }
    )


def run_pca(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    n_components: int = 2,
    standardize: bool = True,
) -> PCAResult:
    """Run PCA on samples after mean imputation and per-gene scaling (plots only).

    Args:
        expression: Samples as rows, genes as columns.
        metadata: Sample metadata (used to attach group and batch to the scores).
        n_components: Number of components (capped by the data dimensions).
        standardize: Z-score genes (True) or only center them (False).

    Returns:
        A :class:`PCAResult`.

    Raises:
        ValueError: If fewer than two samples or no variable genes are available.
    """
    aligned = align_metadata(expression, metadata)
    prepared = prepare_for_plots(expression, zscore=standardize)
    variable = prepared.columns[(prepared.std(axis=0, ddof=1) > 0).to_numpy()]
    prepared = prepared.loc[:, variable].dropna(axis=1, how="any")
    if prepared.shape[0] < 2 or prepared.shape[1] < 1:
        raise ValueError("PCA needs at least two samples and one variable gene.")
    k = min(n_components, prepared.shape[0], prepared.shape[1])
    pca = PCA(n_components=k, svd_solver="full")
    values = pca.fit_transform(prepared.to_numpy())

    scores = pd.DataFrame(
        {
            SAMPLE_ID_COLUMN: expression.index.to_list(),
            GROUP_COLUMN: aligned[GROUP_COLUMN].to_numpy(),
        }
    )
    if BATCH_COLUMN in aligned.columns:
        scores[BATCH_COLUMN] = aligned[BATCH_COLUMN].to_numpy()
    for i in range(k):
        scores[f"PC{i + 1}"] = values[:, i]
    return PCAResult(scores, pca.explained_variance_ratio_, int(prepared.shape[1]), standardize)


def select_heatmap_samples(
    expression: pd.DataFrame, metadata: pd.DataFrame, max_samples: int
) -> list[str]:
    """Choose a deterministic, group-balanced subset of samples for heatmaps.

    Samples are ordered by group, then by sample ID. If there are more than
    ``max_samples``, each group contributes the same number of samples (its first ones
    by ID), so no randomness is involved.

    Args:
        expression: Samples as rows.
        metadata: Sample metadata.
        max_samples: Upper bound on the number of samples returned.

    Returns:
        Sample IDs in display order.
    """
    aligned = align_metadata(expression, metadata)
    ordered = aligned.reset_index().sort_values([GROUP_COLUMN, SAMPLE_ID_COLUMN], kind="stable")
    if len(ordered) <= max_samples:
        return ordered[SAMPLE_ID_COLUMN].tolist()
    quota = max(1, math.floor(max_samples / ordered[GROUP_COLUMN].nunique()))
    return ordered.groupby(GROUP_COLUMN, sort=True).head(quota)[SAMPLE_ID_COLUMN].tolist()


def sample_correlation(
    expression: pd.DataFrame, sample_ids: list[str], *, standardize: bool = True
) -> pd.DataFrame:
    """Pearson correlation between samples, computed on per-gene centered/scaled values.

    Centering each gene first means the correlation reflects sample-to-sample
    similarity in deviations from the gene mean, instead of being dominated by
    baseline expression levels.

    Args:
        expression: Samples as rows, genes as columns.
        sample_ids: Samples to include, in display order.
        standardize: Z-score genes (True) or only center them (False).

    Returns:
        A square DataFrame indexed and labeled by ``sample_ids``.
    """
    prepared = prepare_for_plots(expression, zscore=standardize).loc[sample_ids]
    prepared = prepared.dropna(axis=1, how="any")
    matrix = np.corrcoef(prepared.to_numpy())
    return pd.DataFrame(matrix, index=sample_ids, columns=sample_ids)
