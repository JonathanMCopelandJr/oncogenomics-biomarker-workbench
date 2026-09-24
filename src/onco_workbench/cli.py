"""Command-line interface for the workbench (``obw``).

This is the project's cross-platform task runner: the Makefile and README
commands delegate here. In Phase 1 only ``disclaimer`` is functional; the
pipeline subcommands are registered so the interface is visible, but they exit
with a clear message until they are implemented in later phases.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from onco_workbench import __version__
from onco_workbench.disclaimers import DATA_LABEL, FULL_DISCLAIMER

# Subcommand name -> (help text, phase in which it will be implemented).
_PLANNED_COMMANDS: dict[str, tuple[str, int]] = {
    "generate-data": ("Generate the small synthetic demo dataset.", 2),
    "validate": ("Validate an expression matrix and its sample metadata.", 2),
    "run-analysis": ("Run QC, group comparison, figures, and reports.", 3),
    "dashboard": ("Launch the Streamlit dashboard.", 4),
}

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

    subparsers.add_parser("disclaimer", help="Print the full nonclinical disclaimer.")
    for name, (help_text, phase) in _PLANNED_COMMANDS.items():
        subparsers.add_parser(name, help=f"{help_text} [not yet implemented: Phase {phase}]")

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI.

    Args:
        argv: Arguments excluding the program name. Defaults to ``sys.argv[1:]``.

    Returns:
        Process exit code (0 on success).
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 0

    if args.command == "disclaimer":
        print(DATA_LABEL)
        print(FULL_DISCLAIMER)
        return 0

    _, phase = _PLANNED_COMMANDS[args.command]
    print(
        f"obw {args.command}: not implemented yet (planned for Phase {phase}).",
        file=sys.stderr,
    )
    return NOT_IMPLEMENTED_EXIT_CODE
