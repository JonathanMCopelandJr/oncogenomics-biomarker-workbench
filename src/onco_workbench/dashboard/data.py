"""Streamlit-free data helpers for the dashboard.

Everything here is deterministic and unit tested. Streamlit pages call these functions
(through cached wrappers in :mod:`onco_workbench.dashboard.components`) instead of
re-implementing any analysis.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import pandas as pd

from onco_workbench.config import ComparisonConfig, WorkbenchConfig, load_config
from onco_workbench.data.io import labeled_csv_text, load_ground_truth, load_manifest
from onco_workbench.data.synthetic import GROUP_COLUMN
from onco_workbench.data.validation import ValidationReport, align_metadata
from onco_workbench.pipeline import (
    ComparisonResult,
    load_and_validate,
    results_header_lines,
    run_comparison,
)
from onco_workbench.reporting.manifest import display_path, package_versions


@dataclass(frozen=True)
class DemoData:
    """The validated synthetic demo dataset plus its provenance.

    Attributes:
        config: Configuration used to locate and validate the data.
        expression: Samples as rows, synthetic genes as columns.
        metadata: Raw sample metadata table.
        aligned: Metadata indexed by ``sample_id`` in expression-row order.
        validation: Validation report (always valid; invalid data raise on load).
        ground_truth: Planted-signal table, or None if the file is absent.
        data_manifest: Contents of ``demo_manifest.json``, or None if absent.
    """

    config: WorkbenchConfig
    expression: pd.DataFrame
    metadata: pd.DataFrame
    aligned: pd.DataFrame
    validation: ValidationReport
    ground_truth: pd.DataFrame | None
    data_manifest: dict[str, Any] | None

    @property
    def groups(self) -> list[str]:
        """Sorted group labels present in the data."""
        return sorted(self.aligned[GROUP_COLUMN].astype(str).unique())


def load_demo_data(config_path: str | Path | None = None) -> DemoData:
    """Load and validate the configured synthetic demo files.

    Args:
        config_path: Optional configuration file (default: ``config/default.yaml``).

    Returns:
        A :class:`DemoData` bundle.

    Raises:
        ConfigError: If the configuration is invalid.
        DataLoadError: If a demo file is missing or unreadable.
        DataValidationError: If the demo data fail validation.
    """
    config = load_config(config_path)
    paths = config.paths
    expression, metadata, validation = load_and_validate(
        config, paths.expression_file, paths.metadata_file
    )
    truth = (
        load_ground_truth(paths.ground_truth_file) if paths.ground_truth_file.is_file() else None
    )
    manifest = load_manifest(paths.manifest_file) if paths.manifest_file.is_file() else None
    return DemoData(
        config=config,
        expression=expression,
        metadata=metadata,
        aligned=align_metadata(expression, metadata),
        validation=validation,
        ground_truth=truth,
        data_manifest=manifest,
    )


def compare(
    data: DemoData,
    *,
    group_a: str,
    group_b: str,
    fdr_threshold: float,
    effect_size_threshold: float,
) -> tuple[ComparisonConfig, ComparisonResult]:
    """Run the group comparison with user-chosen groups and display thresholds.

    The synthetic ground-truth check only runs when the groups are in the order the
    signal was planted (second synthetic label versus the first), because planted
    directions are defined relative to that order.

    Args:
        data: Loaded demo data.
        group_a: Reference group.
        group_b: Comparison group.
        fdr_threshold: Adjusted p-value display threshold.
        effect_size_threshold: |Cohen's d| display threshold.

    Returns:
        ``(comparison settings used, comparison result)``.

    Raises:
        ConfigError: If the groups are identical or a threshold is out of range.
        ValueError: If a group is not present in the data.
    """
    config = data.config.with_comparison(
        group_a=group_a,
        group_b=group_b,
        fdr_threshold=fdr_threshold,
        effect_size_threshold=effect_size_threshold,
    )
    planted_order = (group_a, group_b) == tuple(data.config.synthetic.group_labels)
    truth = data.ground_truth if planted_order else None
    return config.comparison, run_comparison(data.expression, data.metadata, config, truth)


def search_genes(gene_ids: list[str], query: str, limit: int = 200) -> list[str]:
    """Return gene IDs containing ``query`` (case-insensitive), in original order.

    Args:
        gene_ids: Candidate identifiers.
        query: Substring to match; empty matches everything.
        limit: Maximum number of matches returned.

    Returns:
        Matching identifiers.
    """
    needle = query.strip().lower()
    matches = [g for g in gene_ids if needle in g.lower()]
    return matches[:limit]


def filter_results(
    results: pd.DataFrame, *, query: str = "", only_meeting_thresholds: bool = False
) -> pd.DataFrame:
    """Filter comparison results by gene-ID substring and/or threshold status.

    Args:
        results: Ranked comparison results.
        query: Case-insensitive gene-ID substring.
        only_meeting_thresholds: Keep only rows with ``meets_thresholds`` True.

    Returns:
        A filtered copy, preserving rank order.
    """
    mask = pd.Series(True, index=results.index)
    if query.strip():
        mask &= results["gene_id"].str.lower().str.contains(query.strip().lower(), regex=False)
    if only_meeting_thresholds:
        mask &= results["meets_thresholds"].astype(bool)
    return results[mask].copy()


def results_csv_bytes(results: pd.DataFrame, comparison: ComparisonConfig) -> bytes:
    """Render results as a labeled CSV download, identical in format to pipeline exports.

    Args:
        results: Results to export (possibly filtered).
        comparison: Settings that produced them (written into the header comments).

    Returns:
        UTF-8 bytes beginning with the DEMONSTRATION / SYNTHETIC comment lines.
    """
    text = labeled_csv_text(results, index=False, extra_comments=results_header_lines(comparison))
    return text.encode("utf-8")


def summary_statistics(expression: pd.DataFrame) -> pd.DataFrame:
    """Describe the distribution of all observed expression values.

    Args:
        expression: Samples as rows, genes as columns.

    Returns:
        One-column table (count, mean, std, min, quartiles, max).
    """
    values = pd.Series(expression.to_numpy().ravel()).dropna()
    return values.describe().to_frame(name="all observed values")


def gene_group_summary(data: DemoData, gene_id: str) -> pd.DataFrame:
    """Per-group count, mean, SD, and missing count for one gene.

    Args:
        data: Loaded demo data.
        gene_id: Gene to summarize.

    Returns:
        One row per group.

    Raises:
        KeyError: If the gene is not present.
    """
    if gene_id not in data.expression.columns:
        raise KeyError(f"Unknown gene ID: {gene_id}")
    values = data.expression[gene_id]
    groups = data.aligned[GROUP_COLUMN].astype(str)
    frame = pd.DataFrame({"group": groups.to_numpy(), "value": values.to_numpy()})
    summary = frame.groupby("group")["value"].agg(
        n_observed="count", mean="mean", sd="std", n_missing=lambda s: int(s.isna().sum())
    )
    return summary.reset_index()


def environment_info(data: DemoData) -> list[tuple[str, str]]:
    """Reproducibility details for display: config, seed, data checksums, versions.

    Args:
        data: Loaded demo data.

    Returns:
        ``(label, value)`` pairs.
    """
    config = data.config
    rows = [
        ("Configuration", display_path(config.source, config.root)),
        ("Seed", str(config.seed)),
    ]
    if data.data_manifest is not None:
        rows.append(("Demo data license", data.data_manifest["license"]["spdx"]))
        for name, entry in sorted(data.data_manifest["files"].items()):
            rows.append((f"SHA-256 {name}", entry["sha256"]))
    versions = dict(package_versions())
    try:
        versions["streamlit"] = version("streamlit")
    except PackageNotFoundError:  # pragma: no cover - streamlit is a runtime dependency
        versions["streamlit"] = None
    rows += [(f"{name} version", str(value)) for name, value in versions.items()]
    return rows
