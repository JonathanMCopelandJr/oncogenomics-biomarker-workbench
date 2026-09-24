"""Command-line interface for the workbench (``obw``).

This is the project's cross-platform task runner: the Makefile and README
commands delegate here. Commands: ``disclaimer``, ``generate-data``, ``validate``,
``qc``, ``run-analysis``, and ``dashboard``.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from onco_workbench import __version__
from onco_workbench.config import PROJECT_ROOT, ConfigError, WorkbenchConfig, load_config
from onco_workbench.data.io import DataLoadError, load_expression, load_metadata, write_dataset
from onco_workbench.data.synthetic import generate_synthetic_dataset
from onco_workbench.data.validation import DataValidationError, ValidationReport, validate_dataset
from onco_workbench.disclaimers import DATA_LABEL, FULL_DISCLAIMER, SHORT_DISCLAIMER

EXIT_OK = 0
EXIT_FAILURE = 1

DASHBOARD_ENTRY = Path("app") / "Home.py"


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level argument parser.

    Returns:
        The configured :class:`argparse.ArgumentParser`.
    """
    parser = argparse.ArgumentParser(
        prog="obw",
        description=(
            "oncogenomics-biomarker-workbench - research and education only, not for clinical use."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="<command>")

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--config",
        type=Path,
        default=None,
        help="YAML configuration file (default: config/default.yaml).",
    )

    subparsers.add_parser("disclaimer", help="Print the full nonclinical disclaimer.")

    generate = subparsers.add_parser(
        "generate-data",
        parents=[common],
        help="Generate the small synthetic demo dataset.",
    )
    generate.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Write files here instead of the configured synthetic data directory.",
    )
    generate.add_argument(
        "--seed", type=int, default=None, help="Override the configured random seed."
    )

    validate = subparsers.add_parser(
        "validate",
        parents=[common],
        help="Validate an expression matrix and its sample metadata.",
    )
    validate.add_argument(
        "--expression", type=Path, default=None, help="Expression matrix CSV (samples x genes)."
    )
    validate.add_argument("--metadata", type=Path, default=None, help="Sample metadata CSV.")

    inputs = argparse.ArgumentParser(add_help=False)
    inputs.add_argument(
        "--expression", type=Path, default=None, help="Expression CSV (default: demo data)."
    )
    inputs.add_argument(
        "--metadata", type=Path, default=None, help="Metadata CSV (default: demo data)."
    )
    inputs.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for tables, figures, report, and manifest (default: outputs/).",
    )

    subparsers.add_parser(
        "qc",
        parents=[common, inputs],
        help="Validate inputs and write QC tables, figures, and a QC report.",
    )

    analysis = subparsers.add_parser(
        "run-analysis",
        parents=[common, inputs],
        help="Run QC and the two-group comparison; write tables, figures, and a report.",
    )
    analysis.add_argument(
        "--ground-truth",
        type=Path,
        default=None,
        help="Synthetic ground-truth CSV for the workflow check (default: demo file when "
        "using demo data).",
    )
    analysis.add_argument("--group-a", default=None, help="Reference group label.")
    analysis.add_argument("--group-b", default=None, help="Comparison group label.")
    analysis.add_argument(
        "--fdr-threshold", type=float, default=None, help="Adjusted p-value display threshold."
    )
    analysis.add_argument(
        "--effect-size-threshold",
        type=float,
        default=None,
        help="Absolute Cohen's d display threshold.",
    )

    dashboard = subparsers.add_parser(
        "dashboard",
        help="Launch the local Streamlit dashboard (synthetic demo data only).",
    )
    dashboard.add_argument(
        "--port", type=int, default=8501, help="Local port to serve on (default: 8501)."
    )
    dashboard.add_argument(
        "--headless",
        action="store_true",
        help="Do not open a browser window automatically.",
    )

    return parser


def _validate_with_config(
    config: WorkbenchConfig, expression_path: Path, metadata_path: Path
) -> ValidationReport:
    expression = load_expression(expression_path)
    metadata = load_metadata(metadata_path)
    return _validate_frames(config, expression, metadata)


def _validate_frames(
    config: WorkbenchConfig, expression: pd.DataFrame, metadata: pd.DataFrame
) -> ValidationReport:
    v = config.validation
    return validate_dataset(
        expression,
        metadata,
        allowed_groups=v.allowed_groups,
        min_samples_per_group=v.min_samples_per_group,
        max_missing_fraction_per_gene=v.max_missing_fraction_per_gene,
        max_missing_fraction_per_sample=v.max_missing_fraction_per_sample,
    )


def _cmd_generate_data(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    if args.seed is not None:
        config = config.with_seed(args.seed)
    paths = config.paths
    targets = {
        "expression_path": paths.expression_file,
        "metadata_path": paths.metadata_file,
        "ground_truth_path": paths.ground_truth_file,
        "manifest_path": paths.manifest_file,
    }
    if args.output_dir is not None:
        targets = {key: args.output_dir / path.name for key, path in targets.items()}

    dataset = generate_synthetic_dataset(config.synthetic, config.seed)
    report = _validate_frames(config, dataset.expression, dataset.metadata)
    if not report.is_valid:
        print(report.summary(), file=sys.stderr)
        print("Generated data failed validation; no files were written.", file=sys.stderr)
        return EXIT_FAILURE

    written = write_dataset(dataset, config.synthetic, **targets)
    n_samples, n_genes = dataset.expression.shape
    n_signal = int(dataset.ground_truth["is_signal"].sum())
    print(DATA_LABEL)
    print(
        f"Generated synthetic demo data (seed={config.seed}): {n_samples} samples x "
        f"{n_genes} genes, {n_signal} planted signal genes."
    )
    for path in written.values():
        print(f"  wrote {path} ({path.stat().st_size:,} bytes)")
    print(report.summary())
    print(SHORT_DISCLAIMER)
    return EXIT_OK


def _cmd_validate(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    expression_path = args.expression or config.paths.expression_file
    metadata_path = args.metadata or config.paths.metadata_file
    report = _validate_with_config(config, expression_path, metadata_path)
    print(f"Expression: {expression_path}")
    print(f"Metadata:   {metadata_path}")
    print(report.summary())
    return EXIT_OK if report.is_valid else EXIT_FAILURE


def _run_pipeline_command(args: argparse.Namespace, *, include_comparison: bool) -> int:
    # Imported lazily so `obw --version` and `obw validate` stay fast.
    from onco_workbench.pipeline import run_pipeline

    config = load_config(args.config)
    overrides = {
        key: value
        for key, value in {
            "group_a": getattr(args, "group_a", None),
            "group_b": getattr(args, "group_b", None),
            "fdr_threshold": getattr(args, "fdr_threshold", None),
            "effect_size_threshold": getattr(args, "effect_size_threshold", None),
        }.items()
        if value is not None
    }
    if overrides:
        config = config.with_comparison(**overrides)
    try:
        run = run_pipeline(
            config,
            expression_path=args.expression,
            metadata_path=args.metadata,
            ground_truth_path=getattr(args, "ground_truth", None),
            output_dir=args.output_dir,
            include_comparison=include_comparison,
        )
    except DataValidationError as exc:
        print(exc.report.summary(), file=sys.stderr)
        print("Inputs failed validation; no analysis outputs were written.", file=sys.stderr)
        return EXIT_FAILURE
    except ValueError as exc:
        if isinstance(exc, ConfigError | DataLoadError):
            raise
        print(f"obw {args.command}: {exc}", file=sys.stderr)
        return EXIT_FAILURE

    print(DATA_LABEL)
    qc = run.qc.summary
    print(f"{run.validation.summary().splitlines()[0]}")
    print(f"QC: {qc.n_samples} samples x {qc.n_genes} genes, {qc.missing_fraction:.2%} missing.")
    if run.comparison is not None:
        c = config.comparison
        results = run.comparison.results
        print(
            f"Comparison {c.group_b} vs {c.group_a}: {int(results['p_value'].notna().sum())} "
            f"genes tested, {int(results['meets_thresholds'].sum())} meet the display "
            f"thresholds (adj. p <= {c.fdr_threshold:g}, |d| >= {c.effect_size_threshold:g})."
        )
        if run.comparison.recovery is not None:
            r = run.comparison.recovery
            print(
                f"Synthetic ground-truth check: {r.true_positives}/{r.n_planted} planted genes "
                f"recovered, {r.false_positives} flagged genes were not planted."
            )
    print(f"Outputs written to {run.output_dir}:")
    for path in run.files.values():
        print(f"  {path.relative_to(run.output_dir).as_posix()}")
    print(SHORT_DISCLAIMER)
    return EXIT_OK


def _cmd_qc(args: argparse.Namespace) -> int:
    return _run_pipeline_command(args, include_comparison=False)


def _cmd_run_analysis(args: argparse.Namespace) -> int:
    return _run_pipeline_command(args, include_comparison=True)


def dashboard_command(port: int, *, headless: bool, root: Path = PROJECT_ROOT) -> list[str]:
    """Build the ``streamlit run`` command for the local dashboard.

    The server binds to ``localhost`` only and never sends usage statistics.

    Args:
        port: Local port.
        headless: If True, do not open a browser.
        root: Project root containing ``app/Home.py``.

    Returns:
        The command as a list of arguments.
    """
    command = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(root / DASHBOARD_ENTRY),
        "--server.address",
        "localhost",
        "--server.port",
        str(port),
        "--browser.gatherUsageStats",
        "false",
    ]
    if headless:
        command += ["--server.headless", "true"]
    return command


def _cmd_dashboard(args: argparse.Namespace) -> int:
    entry = PROJECT_ROOT / DASHBOARD_ENTRY
    if not entry.is_file():
        print(f"obw dashboard: app entry point not found: {entry}", file=sys.stderr)
        return EXIT_FAILURE
    print(DATA_LABEL)
    print(f"Starting the local dashboard at http://localhost:{args.port}  (Ctrl+C to stop)")
    print(SHORT_DISCLAIMER)
    try:
        return subprocess.call(
            dashboard_command(args.port, headless=args.headless), cwd=PROJECT_ROOT
        )
    except KeyboardInterrupt:
        return EXIT_OK


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI.

    Args:
        argv: Arguments excluding the program name. Defaults to ``sys.argv[1:]``.

    Returns:
        Process exit code: 0 on success, 1 on failure.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return EXIT_OK

    if args.command == "disclaimer":
        print(DATA_LABEL)
        print(FULL_DISCLAIMER)
        return EXIT_OK

    handlers = {
        "generate-data": _cmd_generate_data,
        "validate": _cmd_validate,
        "qc": _cmd_qc,
        "run-analysis": _cmd_run_analysis,
        "dashboard": _cmd_dashboard,
    }
    try:
        return handlers[args.command](args)
    except (ConfigError, DataLoadError) as exc:
        print(f"obw {args.command}: {exc}", file=sys.stderr)
        return EXIT_FAILURE
