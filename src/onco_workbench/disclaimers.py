"""Single source of truth for the project's nonclinical disclaimer text.

Every user-facing surface (README, dashboard, reports, CSV exports, run manifests)
should import from this module rather than restating the wording, so the boundary
cannot drift between surfaces.
"""

from __future__ import annotations

DATA_LABEL: str = "DEMONSTRATION / SYNTHETIC DATA"

SHORT_DISCLAIMER: str = (
    "RESEARCH AND EDUCATION ONLY - NOT FOR CLINICAL USE. "
    "Outputs are generated from SYNTHETIC data and do not describe real patients."
)

FULL_DISCLAIMER: str = (
    "This software is a research and education portfolio project only. It is NOT a "
    "medical device, diagnostic tool, clinical decision-support system, treatment "
    "recommender, or validated biomarker-discovery system. It must not be used to "
    "inform any decision about any person's health. Demonstration outputs are "
    "generated from SYNTHETIC data created by this repository; synthetic gene "
    "identifiers (e.g. SYN_G0001) do not correspond to real genes, and no result "
    "should be interpreted as evidence about real biology or as a validated biomarker."
)


def as_comment_block(prefix: str = "# ") -> str:
    """Return the data label and short disclaimer as prefixed comment lines.

    Intended for headers of exported text files (e.g. CSVs read back with
    ``pandas.read_csv(..., comment="#")``).

    Args:
        prefix: String prepended to each line, such as ``"# "``.

    Returns:
        A newline-terminated block of comment lines.
    """
    lines = [DATA_LABEL, SHORT_DISCLAIMER]
    return "".join(f"{prefix}{line}\n" for line in lines)
