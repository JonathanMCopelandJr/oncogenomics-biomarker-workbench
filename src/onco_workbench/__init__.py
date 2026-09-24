"""oncogenomics-biomarker-workbench: a research and education workbench.

Explores SYNTHETIC transcriptomics-style data through a transparent validation,
QC, group-comparison, and visualization workflow. Not for clinical use; see
:mod:`onco_workbench.disclaimers`.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__: str = version("oncogenomics-biomarker-workbench")
except PackageNotFoundError:  # pragma: no cover - only when running from an uninstalled tree
    __version__ = "0.0.0+unknown"

__all__ = ["__version__"]
