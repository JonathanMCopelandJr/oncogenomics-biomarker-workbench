"""Command-line interface for the workbench (``obw``).

This is the project's cross-platform task runner: the Makefile and README
commands delegate here. Implemented: ``disclaimer``, ``generate-data``, and
``validate``. The remaining subcommands are registered so the interface is
visible, but they exit with a clear message until later phases implement them.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from onco_workbench import __version__
from onco_workbench.config import ConfigError, WorkbenchConfig, load_config
from onco_workbench.data.io import DataLoadError, load_expression, load_metadata, write_dataset
from onco_workbench.data.synthetic import generate_synthetic_dataset
from onco_workbench.data.validation import ValidationReport, validate_dataset
from onco_workbench.disclaimers import DATA_LABEL, FULL_DISCLAIMER, SHORT_DISCLAIMER

# Subcommand name -> (help text, phase in which it will be implemented).
_PLANNED_COMMANDS: dict[str, tuple[str, int]] = {
    "run-analysis": ("Run QC, group comparison, figures, and reports.", 3),
    "dashboard": ("Launch the Streamlit dashboard.", 4),
}

EXIT_OK = 0
EXIT_FAILURE = 1
NOT_IMPLEMENTED_EXIT_CODE = 2


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

    for name, (help_text, phase) in _PLANNED_COMMANDS.items():
        subparsers.add_parser(name, help=f"{help_text} [not yet implemented: Phase {phase}]")

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


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI.

    Args:
        argv: Arguments excluding the program name. Defaults to ``sys.argv[1:]``.

    Returns:
        Process exit code: 0 on success, 1 on failure, 2 for not-yet-implemented commands.
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

    handlers = {"generate-data": _cmd_generate_data, "validate": _cmd_validate}
    if args.command in handlers:
        try:
            return handlers[args.command](args)
        except (ConfigError, DataLoadError) as exc:
            print(f"obw {args.command}: {exc}", file=sys.stderr)
            return EXIT_FAILURE

    _, phase = _PLANNED_COMMANDS[args.command]
    print(
        f"obw {args.command}: not implemented yet (planned for Phase {phase}).",
        file=sys.stderr,
    )
    return NOT_IMPLEMENTED_EXIT_CODE
