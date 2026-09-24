"""End-to-end analysis pipeline: load, validate, QC, normalize, compare, rank, export.

Each step is also available as a separate function, so the CLI, notebooks, and the
dashboard reuse exactly the same logic. Every output (tables, figures, report,
manifest) is labeled DEMONSTRATION / SYNTHETIC.

Output layout (under the configured ``outputs_dir`` or ``output_dir``)::

    tables/   qc_sample_summary.csv, qc_gene_summary.csv, pca_scores.csv,
              group_comparison_results.csv, top_ranked_genes.csv
    figures/  qc_pca, qc_expression_distributions, qc_sample_totals,
              qc_sample_correlation, volcano, top_genes, top_genes_heatmap
    report.md
    run_manifest.json
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from onco_workbench.analysis.differential import compare_groups
from onco_workbench.analysis.multiple_testing import benjamini_hochberg
from onco_workbench.analysis.normalization import median_center_samples
from onco_workbench.analysis.qc import (
    PCAResult,
    QCSummary,
    compute_qc_summary,
    gene_qc_table,
    run_pca,
    sample_correlation,
    sample_qc_table,
    select_heatmap_samples,
)
from onco_workbench.analysis.ranking import add_ranking, flag_results
from onco_workbench.analysis.recovery import RecoverySummary, compare_with_ground_truth
from onco_workbench.config import WorkbenchConfig
from onco_workbench.data.io import (
    load_expression,
    load_ground_truth,
    load_metadata,
    write_labeled_csv,
)
from onco_workbench.data.validation import (
    DataValidationError,
    ValidationReport,
    align_metadata,
    validate_dataset,
)
from onco_workbench.reporting.manifest import build_run_manifest, display_path, write_run_manifest
from onco_workbench.reporting.markdown_report import (
    ComparisonSection,
    ReportContext,
    render_markdown_report,
)
from onco_workbench.viz import static

TABLES_DIR = "tables"
FIGURES_DIR = "figures"
REPORT_NAME = "report.md"
MANIFEST_NAME = "run_manifest.json"

RESULT_COLUMN_ORDER: tuple[str, ...] = (
    "rank",
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
    "p_adj",
    "ranking_score",
    "direction",
    "meets_thresholds",
)


@dataclass(frozen=True)
class QCResult:
    """Outputs of the QC step."""

    summary: QCSummary
    sample_table: pd.DataFrame
    gene_table: pd.DataFrame
    pca: PCAResult
    correlation: pd.DataFrame


@dataclass(frozen=True)
class ComparisonResult:
    """Outputs of the group-comparison step."""

    results: pd.DataFrame
    normalized: pd.DataFrame
    recovery: RecoverySummary | None


@dataclass
class PipelineRun:
    """Everything produced by :func:`run_pipeline`."""

    output_dir: Path
    validation: ValidationReport
    qc: QCResult
    comparison: ComparisonResult | None
    files: dict[str, Path] = field(default_factory=dict)


def load_and_validate(
    config: WorkbenchConfig, expression_path: Path, metadata_path: Path
) -> tuple[pd.DataFrame, pd.DataFrame, ValidationReport]:
    """Load both input files and validate them with the configured thresholds.

    Args:
        config: Workbench configuration.
        expression_path: Expression matrix CSV.
        metadata_path: Sample metadata CSV.

    Returns:
        ``(expression, metadata, report)``.

    Raises:
        DataValidationError: If validation finds errors.
    """
    expression = load_expression(expression_path)
    metadata = load_metadata(metadata_path)
    v = config.validation
    report = validate_dataset(
        expression,
        metadata,
        allowed_groups=v.allowed_groups,
        min_samples_per_group=v.min_samples_per_group,
        max_missing_fraction_per_gene=v.max_missing_fraction_per_gene,
        max_missing_fraction_per_sample=v.max_missing_fraction_per_sample,
    )
    if not report.is_valid:
        raise DataValidationError(report)
    return expression, metadata, report


def run_qc(expression: pd.DataFrame, metadata: pd.DataFrame, config: WorkbenchConfig) -> QCResult:
    """Compute QC summaries, PCA, and a group-balanced sample correlation matrix.

    Args:
        expression: Validated expression matrix.
        metadata: Validated metadata.
        config: Workbench configuration.

    Returns:
        A :class:`QCResult`.
    """
    standardize = config.normalization.zscore_genes_for_plots
    samples = select_heatmap_samples(expression, metadata, config.qc.heatmap_max_samples)
    return QCResult(
        summary=compute_qc_summary(expression, metadata),
        sample_table=sample_qc_table(expression, metadata),
        gene_table=gene_qc_table(expression),
        pca=run_pca(
            expression, metadata, n_components=config.qc.pca_components, standardize=standardize
        ),
        correlation=sample_correlation(expression, samples, standardize=standardize),
    )


def run_comparison(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    config: WorkbenchConfig,
    ground_truth: pd.DataFrame | None = None,
) -> ComparisonResult:
    """Normalize, compare the configured groups, adjust p-values, flag, and rank.

    Args:
        expression: Validated expression matrix.
        metadata: Validated metadata.
        config: Workbench configuration (``comparison``, ``ranking``, ``normalization``).
        ground_truth: Optional synthetic ground truth for the workflow check.

    Returns:
        A :class:`ComparisonResult` with results sorted by rank.
    """
    c = config.comparison
    normalized = (
        median_center_samples(expression)
        if config.normalization.median_center_samples
        else expression
    )
    results = compare_groups(
        normalized,
        metadata,
        c.group_a,
        c.group_b,
        min_non_missing_per_group=c.min_non_missing_per_group,
    )
    results["p_adj"] = benjamini_hochberg(results["p_value"])
    results = flag_results(
        results, fdr_threshold=c.fdr_threshold, effect_size_threshold=c.effect_size_threshold
    )
    results = add_ranking(results, p_floor=config.ranking.p_floor)
    results = results[list(RESULT_COLUMN_ORDER)]
    recovery = (
        compare_with_ground_truth(results, ground_truth) if ground_truth is not None else None
    )
    return ComparisonResult(results=results, normalized=normalized, recovery=recovery)


def _normalization_steps(config: WorkbenchConfig, include_comparison: bool) -> list[str]:
    steps = []
    if include_comparison:
        steps.append(
            "Group comparison: per-sample median centering "
            "(x - sample median + median of sample medians)."
            if config.normalization.median_center_samples
            else "Group comparison: no normalization (values used as provided)."
        )
    scaling = "z-scored" if config.normalization.zscore_genes_for_plots else "centered"
    steps.append(f"PCA and heatmaps only: missing values replaced by gene means, genes {scaling}.")
    steps.append(
        "These steps suit simulated continuous values only. They are not appropriate for "
        "real RNA-seq counts."
    )
    return steps


def run_pipeline(
    config: WorkbenchConfig,
    *,
    expression_path: str | Path | None = None,
    metadata_path: str | Path | None = None,
    ground_truth_path: str | Path | None = None,
    output_dir: str | Path | None = None,
    include_comparison: bool = True,
    generated_at: datetime | None = None,
) -> PipelineRun:
    """Run the analysis and write labeled tables, figures, a report, and a manifest.

    When no input paths are given, the configured synthetic demo files are used, and
    the synthetic ground truth (if present) is used for the workflow check. When a
    custom expression file is given, the ground truth is only used if passed explicitly.

    Args:
        config: Workbench configuration.
        expression_path: Expression CSV (default: configured demo file).
        metadata_path: Metadata CSV (default: configured demo file).
        ground_truth_path: Optional ground-truth CSV for the workflow check.
        output_dir: Destination directory (default: configured ``outputs_dir``).
        include_comparison: If False, only QC outputs are produced.
        generated_at: Timestamp to record (default: now, UTC).

    Returns:
        A :class:`PipelineRun` describing everything produced.

    Raises:
        DataValidationError: If the inputs fail validation.
        ValueError: If the configured groups are not present in the data.
    """
    paths = config.paths
    using_demo = expression_path is None
    expression_file = Path(expression_path) if expression_path else paths.expression_file
    metadata_file = Path(metadata_path) if metadata_path else paths.metadata_file
    truth_file: Path | None = Path(ground_truth_path) if ground_truth_path else None
    if truth_file is None and using_demo and paths.ground_truth_file.is_file():
        truth_file = paths.ground_truth_file
    out = Path(output_dir) if output_dir else paths.outputs_dir
    timestamp = (generated_at or datetime.now(UTC)).isoformat(timespec="seconds")

    expression, metadata, validation = load_and_validate(config, expression_file, metadata_file)
    aligned = align_metadata(expression, metadata)
    qc = run_qc(expression, metadata, config)
    comparison = None
    if include_comparison:
        truth = load_ground_truth(truth_file) if truth_file else None
        comparison = run_comparison(expression, metadata, config, truth)

    run = PipelineRun(out, validation, qc, comparison)
    tables_dir, figures_dir = out / TABLES_DIR, out / FIGURES_DIR
    c = config.comparison

    tables = {
        "qc_sample_summary": (qc.sample_table, []),
        "qc_gene_summary": (qc.gene_table, []),
        "pca_scores": (qc.pca.scores, []),
    }
    if comparison is not None:
        header = [
            f"Comparison: {c.group_b} (b) vs {c.group_a} (a); "
            "mean_diff and cohens_d are b minus a.",
            f"meets_thresholds: p_adj <= {c.fdr_threshold:g} and |cohens_d| >= "
            f"{c.effect_size_threshold:g} (display thresholds, not conclusions).",
        ]
        tables["group_comparison_results"] = (comparison.results, header)
        top = comparison.results[comparison.results["rank"].notna()].nsmallest(
            config.ranking.top_n, "rank"
        )
        tables["top_ranked_genes"] = (top, header)
    for name, (frame, comments) in tables.items():
        run.files[f"table:{name}"] = write_labeled_csv(
            frame, tables_dir / f"{name}.csv", index=False, extra_comments=comments
        )

    fmt, dpi = config.figures.format, config.figures.dpi
    standardize = config.normalization.zscore_genes_for_plots
    figures = {
        "pca": ("qc_pca", static.plot_pca(qc.pca)),
        "distributions": (
            "qc_expression_distributions",
            static.plot_expression_distributions(expression, aligned),
        ),
        "sample_totals": ("qc_sample_totals", static.plot_sample_totals(qc.sample_table)),
        "correlation": (
            "qc_sample_correlation",
            static.plot_correlation_heatmap(qc.correlation, aligned),
        ),
    }
    if comparison is not None:
        results = comparison.results
        top_ids = results[results["rank"].notna()].nsmallest(config.ranking.top_n, "rank")
        figures["volcano"] = (
            "volcano",
            static.plot_volcano(
                results,
                group_a=c.group_a,
                group_b=c.group_b,
                fdr_threshold=c.fdr_threshold,
                effect_size_threshold=c.effect_size_threshold,
            ),
        )
        figures["top_genes"] = (
            "top_genes",
            static.plot_top_genes(
                results, group_a=c.group_a, group_b=c.group_b, top_n=config.ranking.top_n
            ),
        )
        figures["top_heatmap"] = (
            "top_genes_heatmap",
            static.plot_top_gene_heatmap(
                comparison.normalized, aligned, top_ids["gene_id"].tolist(), standardize=standardize
            ),
        )
    figure_links: dict[str, str] = {}
    for key, (stem, fig) in figures.items():
        path = static.save_figure(fig, figures_dir / f"{stem}.{fmt}", dpi=dpi)
        run.files[f"figure:{key}"] = path
        figure_links[key] = f"{FIGURES_DIR}/{path.name}"

    inputs = {"expression": expression_file, "metadata": metadata_file}
    if comparison is not None and truth_file is not None:
        inputs["ground_truth"] = truth_file
    report = ReportContext(
        generated_at=timestamp,
        config_path=display_path(config.source, config.root),
        seed=config.seed,
        inputs={role: display_path(p, config.root) for role, p in inputs.items()},
        validation=validation,
        qc=qc.summary,
        pca_variance=list(qc.pca.explained_variance_ratio),
        normalization_steps=_normalization_steps(config, include_comparison),
        figures=figure_links,
        tables={
            key.split(":", 1)[1]: f"{TABLES_DIR}/{path.name}"
            for key, path in run.files.items()
            if key.startswith("table:")
        },
        comparison=(
            ComparisonSection(c, config.ranking, comparison.results, comparison.recovery)
            if comparison is not None
            else None
        ),
    )
    report_path = out / REPORT_NAME
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(render_markdown_report(report), encoding="utf-8", newline="\n")
    run.files["report"] = report_path

    parameters = {
        "include_comparison": include_comparison,
        "normalization": asdict(config.normalization),
        "comparison": asdict(c) if include_comparison else None,
        "ranking": asdict(config.ranking) if include_comparison else None,
        "qc": asdict(config.qc),
        "validation": asdict(config.validation),
    }
    manifest = build_run_manifest(
        root=config.root,
        config_path=config.source,
        seed=config.seed,
        parameters=parameters,
        inputs=inputs,
        outputs=run.files,
        generated_at=timestamp,
    )
    run.files["manifest"] = write_run_manifest(manifest, out / MANIFEST_NAME)
    return run
