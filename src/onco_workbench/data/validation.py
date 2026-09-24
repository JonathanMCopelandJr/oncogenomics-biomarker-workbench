"""Schema and data-quality validation for expression matrices and sample metadata.

Validation collects every problem into a :class:`ValidationReport` instead of
stopping at the first one, so users see everything that needs fixing at once.
Problems that make analysis impossible or misleading are **errors**. Problems a
user should look at, but that do not block analysis, are **warnings**.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import pandas as pd
from pandas.api.types import is_bool_dtype, is_numeric_dtype

from onco_workbench.data.synthetic import BATCH_COLUMN, GROUP_COLUMN, SAMPLE_ID_COLUMN

Severity = Literal["error", "warning"]

_PREVIEW_LIMIT = 5


@dataclass(frozen=True)
class ValidationIssue:
    """One validation finding.

    Attributes:
        severity: ``"error"`` (blocks analysis) or ``"warning"`` (review recommended).
        code: Stable machine-readable identifier, e.g. ``"duplicate_gene_id"``.
        message: Human-readable explanation, including offending identifiers.
    """

    severity: Severity
    code: str
    message: str


@dataclass
class ValidationReport:
    """All findings from validating one dataset."""

    n_samples: int
    n_genes: int
    issues: list[ValidationIssue] = field(default_factory=list)

    def add(self, severity: Severity, code: str, message: str) -> None:
        """Record a finding.

        Args:
            severity: ``"error"`` or ``"warning"``.
            code: Machine-readable identifier.
            message: Human-readable explanation.
        """
        self.issues.append(ValidationIssue(severity, code, message))

    @property
    def errors(self) -> list[ValidationIssue]:
        """Findings that block analysis."""
        return [issue for issue in self.issues if issue.severity == "error"]

    @property
    def warnings(self) -> list[ValidationIssue]:
        """Findings that should be reviewed but do not block analysis."""
        return [issue for issue in self.issues if issue.severity == "warning"]

    @property
    def is_valid(self) -> bool:
        """True when there are no errors."""
        return not self.errors

    def codes(self) -> set[str]:
        """Return the set of issue codes, which is convenient for tests and filtering."""
        return {issue.code for issue in self.issues}

    def summary(self) -> str:
        """Return a multi-line, human-readable summary of the report."""
        status = "PASSED" if self.is_valid else "FAILED"
        lines = [
            f"Validation {status}: {self.n_samples} samples x {self.n_genes} genes, "
            f"{len(self.errors)} error(s), {len(self.warnings)} warning(s)."
        ]
        lines += [f"  [{i.severity.upper()}] {i.code}: {i.message}" for i in self.issues]
        return "\n".join(lines)

    def raise_if_invalid(self) -> None:
        """Raise :class:`DataValidationError` if the report contains errors."""
        if not self.is_valid:
            raise DataValidationError(self)


class DataValidationError(ValueError):
    """Raised when data fail validation. The full report is available as ``.report``."""

    def __init__(self, report: ValidationReport) -> None:
        super().__init__(report.summary())
        self.report = report


def _preview(items: Iterable[object]) -> str:
    values = [str(item) for item in items]
    shown = ", ".join(values[:_PREVIEW_LIMIT])
    extra = len(values) - _PREVIEW_LIMIT
    return f"{shown} (+{extra} more)" if extra > 0 else shown


def _is_blank(value: object) -> bool:
    return value is None or (isinstance(value, float) and np.isnan(value)) or not str(value).strip()


def _check_expression(
    expression: pd.DataFrame,
    report: ValidationReport,
    max_missing_gene: float,
    max_missing_sample: float,
) -> None:
    if expression.shape[0] == 0 or expression.shape[1] == 0:
        report.add("error", "empty_expression", "Expression matrix has no samples or no genes.")
        return

    blank_samples = [pos for pos, v in enumerate(expression.index) if _is_blank(v)]
    if blank_samples:
        report.add(
            "error",
            "blank_sample_id",
            f"{len(blank_samples)} expression row(s) have a blank sample ID "
            f"(row positions: {_preview(blank_samples)}).",
        )
    duplicated_samples = expression.index[expression.index.duplicated()].unique()
    if len(duplicated_samples):
        report.add(
            "error",
            "duplicate_sample_id",
            f"Duplicated sample IDs in expression matrix: {_preview(duplicated_samples)}.",
        )

    blank_genes = [pos for pos, v in enumerate(expression.columns) if _is_blank(v)]
    if blank_genes:
        report.add(
            "error",
            "blank_gene_id",
            f"{len(blank_genes)} gene column(s) have a blank ID "
            f"(column positions: {_preview(blank_genes)}).",
        )
    duplicated_genes = expression.columns[expression.columns.duplicated()].unique()
    if len(duplicated_genes):
        report.add(
            "error",
            "duplicate_gene_id",
            f"Duplicated gene IDs in expression matrix: {_preview(duplicated_genes)}.",
        )

    dtypes = list(expression.dtypes)
    numeric_positions = [
        pos for pos, dt in enumerate(dtypes) if is_numeric_dtype(dt) and not is_bool_dtype(dt)
    ]
    non_numeric = [
        expression.columns[pos] for pos in range(len(dtypes)) if pos not in numeric_positions
    ]
    if non_numeric:
        report.add(
            "error",
            "non_numeric_values",
            f"{len(non_numeric)} gene column(s) contain non-numeric values: "
            f"{_preview(non_numeric)}. Expression values must be numbers (blank = missing).",
        )
    if numeric_positions:
        values = expression.iloc[:, numeric_positions].to_numpy(dtype=float)
        n_infinite = int(np.isinf(values).sum())
        if n_infinite:
            report.add(
                "error",
                "non_finite_values",
                f"Expression matrix contains {n_infinite} infinite value(s).",
            )

    missing = expression.isna().to_numpy()
    gene_fraction = missing.mean(axis=0)
    sample_fraction = missing.mean(axis=1)
    _missingness_findings(report, "gene", expression.columns, gene_fraction, max_missing_gene)
    _missingness_findings(report, "sample", expression.index, sample_fraction, max_missing_sample)


def _missingness_findings(
    report: ValidationReport,
    kind: str,
    labels: pd.Index,
    fractions: np.ndarray,
    threshold: float,
) -> None:
    all_missing = labels[fractions >= 1.0]
    if len(all_missing):
        report.add(
            "error",
            f"{kind}_all_missing",
            f"{len(all_missing)} {kind}(s) have no observed values: {_preview(all_missing)}.",
        )
    high = labels[(fractions > threshold) & (fractions < 1.0)]
    if len(high):
        report.add(
            "warning",
            f"{kind}_high_missingness",
            f"{len(high)} {kind}(s) exceed the missingness threshold of {threshold:.0%}: "
            f"{_preview(high)}.",
        )


def _check_metadata(
    metadata: pd.DataFrame,
    report: ValidationReport,
    allowed_groups: Sequence[str] | None,
) -> bool:
    """Check metadata structure and values. Returns False if alignment checks cannot run."""
    missing_columns = [c for c in (SAMPLE_ID_COLUMN, GROUP_COLUMN) if c not in metadata.columns]
    if missing_columns:
        report.add(
            "error",
            "metadata_missing_column",
            f"Sample metadata is missing required column(s): {', '.join(missing_columns)}. "
            f"Found columns: {_preview(metadata.columns)}.",
        )
        return False
    if metadata.empty:
        report.add("error", "empty_metadata", "Sample metadata has no rows.")
        return False

    ids = metadata[SAMPLE_ID_COLUMN]
    blank_ids = [pos for pos, v in enumerate(ids) if _is_blank(v)]
    if blank_ids:
        report.add(
            "error",
            "metadata_blank_sample_id",
            f"{len(blank_ids)} metadata row(s) have a blank sample_id "
            f"(row positions: {_preview(blank_ids)}).",
        )
    duplicated = ids[ids.duplicated() & ~ids.isna()].unique()
    if len(duplicated):
        report.add(
            "error",
            "metadata_duplicate_sample_id",
            f"Duplicated sample IDs in metadata: {_preview(duplicated)}.",
        )

    groups = metadata[GROUP_COLUMN]
    missing_group = metadata.loc[groups.map(_is_blank), SAMPLE_ID_COLUMN]
    if len(missing_group):
        report.add(
            "error",
            "missing_group",
            f"{len(missing_group)} sample(s) have no group label: {_preview(missing_group)}.",
        )
    if allowed_groups is not None:
        present = groups[~groups.map(_is_blank)]
        invalid = sorted(set(present) - set(allowed_groups))
        if invalid:
            report.add(
                "error",
                "invalid_group",
                f"Unexpected group label(s): {_preview(invalid)}. "
                f"Allowed labels: {', '.join(allowed_groups)}.",
            )

    if BATCH_COLUMN in metadata.columns:
        missing_batch = metadata.loc[metadata[BATCH_COLUMN].map(_is_blank), SAMPLE_ID_COLUMN]
        if len(missing_batch):
            report.add(
                "warning",
                "missing_batch",
                f"{len(missing_batch)} sample(s) have no batch label: {_preview(missing_batch)}.",
            )
    return True


def _check_alignment_and_groups(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    report: ValidationReport,
    allowed_groups: Sequence[str] | None,
    min_samples_per_group: int,
) -> None:
    expression_ids = {str(v) for v in expression.index if not _is_blank(v)}
    metadata_ids = {str(v) for v in metadata[SAMPLE_ID_COLUMN] if not _is_blank(v)}

    only_expression = sorted(expression_ids - metadata_ids)
    if only_expression:
        report.add(
            "error",
            "sample_not_in_metadata",
            f"{len(only_expression)} expression sample(s) have no metadata row: "
            f"{_preview(only_expression)}.",
        )
    only_metadata = sorted(metadata_ids - expression_ids)
    if only_metadata:
        report.add(
            "warning",
            "sample_not_in_expression",
            f"{len(only_metadata)} metadata sample(s) are absent from the expression matrix "
            f"and will be ignored: {_preview(only_metadata)}.",
        )

    shared = metadata[metadata[SAMPLE_ID_COLUMN].astype(str).isin(expression_ids)]
    shared = shared.drop_duplicates(subset=SAMPLE_ID_COLUMN)
    shared = shared[~shared[GROUP_COLUMN].map(_is_blank)]
    counts = shared[GROUP_COLUMN].value_counts()

    groups_to_check = list(allowed_groups) if allowed_groups is not None else list(counts.index)
    if allowed_groups is None and len(counts) < 2:
        report.add(
            "error",
            "too_few_groups",
            f"At least two groups are required for comparison; found {len(counts)}.",
        )
    for group in groups_to_check:
        n = int(counts.get(group, 0))
        if n < min_samples_per_group:
            report.add(
                "error",
                "group_too_small",
                f"Group '{group}' has {n} sample(s); at least {min_samples_per_group} required.",
            )

    if BATCH_COLUMN in shared.columns:
        with_batch = shared[~shared[BATCH_COLUMN].map(_is_blank)]
        if with_batch[BATCH_COLUMN].nunique() > 1 and with_batch[GROUP_COLUMN].nunique() > 1:
            groups_per_batch = with_batch.groupby(BATCH_COLUMN)[GROUP_COLUMN].nunique()
            if (groups_per_batch == 1).all():
                report.add(
                    "warning",
                    "batch_confounded_with_group",
                    "Every batch contains samples from only one group, so batch effects "
                    "cannot be separated from group differences.",
                )


def validate_dataset(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    allowed_groups: Sequence[str] | None = None,
    min_samples_per_group: int = 3,
    max_missing_fraction_per_gene: float = 0.2,
    max_missing_fraction_per_sample: float = 0.2,
) -> ValidationReport:
    """Validate an expression matrix and its sample metadata.

    Checks: non-empty matrix; blank or duplicated sample and gene IDs; numeric,
    finite expression values; genes or samples that are entirely or heavily
    missing; required metadata columns; blank or duplicated metadata IDs;
    missing or unexpected group labels; missing batch labels; sample-ID agreement
    between the matrix and the metadata; minimum group sizes; and batch/group
    confounding.

    Args:
        expression: Samples as rows (index = sample IDs), genes as columns.
        metadata: Table with ``sample_id``, ``group``, and optionally ``batch`` columns.
        allowed_groups: Permitted group labels, or ``None`` to accept any labels
            (at least two groups are then required).
        min_samples_per_group: Minimum number of samples each group must have.
        max_missing_fraction_per_gene: Warn for genes above this missing fraction.
        max_missing_fraction_per_sample: Warn for samples above this missing fraction.

    Returns:
        A :class:`ValidationReport`. Call :meth:`ValidationReport.raise_if_invalid`
        to turn errors into an exception.
    """
    report = ValidationReport(n_samples=expression.shape[0], n_genes=expression.shape[1])
    _check_expression(
        expression, report, max_missing_fraction_per_gene, max_missing_fraction_per_sample
    )
    if _check_metadata(metadata, report, allowed_groups) and expression.shape[0] > 0:
        _check_alignment_and_groups(
            expression, metadata, report, allowed_groups, min_samples_per_group
        )
    return report


def align_metadata(expression: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Return metadata rows for the expression samples, in expression-row order.

    Args:
        expression: Samples as rows (index = sample IDs).
        metadata: Table with a ``sample_id`` column.

    Returns:
        Metadata indexed by ``sample_id``, with one row per expression sample and the
        same order as the expression matrix.

    Raises:
        ValueError: If metadata lacks ``sample_id`` or any expression sample has no row.
    """
    if SAMPLE_ID_COLUMN not in metadata.columns:
        raise ValueError(f"Sample metadata must contain a '{SAMPLE_ID_COLUMN}' column.")
    indexed = metadata.drop_duplicates(subset=SAMPLE_ID_COLUMN).set_index(SAMPLE_ID_COLUMN)
    missing = [sample for sample in expression.index if sample not in indexed.index]
    if missing:
        raise ValueError(
            f"{len(missing)} expression sample(s) have no metadata row: {_preview(missing)}. "
            "Run `obw validate` for a full report."
        )
    return indexed.loc[list(expression.index)]
