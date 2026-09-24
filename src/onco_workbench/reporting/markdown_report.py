"""Render the Markdown summary report for an analysis run.

The report describes what was computed on SYNTHETIC demonstration data, in neutral
language. It never describes genes as biomarkers and never makes clinical or
biological claims. The disclaimer appears at the top and again at the bottom.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import pandas as pd

from onco_workbench.analysis.qc import QCSummary
from onco_workbench.analysis.ranking import DIRECTION_HIGHER, DIRECTION_LOWER
from onco_workbench.analysis.recovery import RecoverySummary
from onco_workbench.config import ComparisonConfig, RankingConfig
from onco_workbench.data.validation import ValidationReport
from onco_workbench.disclaimers import DATA_LABEL, FULL_DISCLAIMER


@dataclass(frozen=True)
class ComparisonSection:
    """Inputs for the group-comparison part of the report."""

    comparison: ComparisonConfig
    ranking: RankingConfig
    results: pd.DataFrame
    recovery: RecoverySummary | None


@dataclass(frozen=True)
class ReportContext:
    """Everything the report needs. Paths are relative to the report file."""

    generated_at: str
    config_path: str
    seed: int
    inputs: Mapping[str, str]
    validation: ValidationReport
    qc: QCSummary
    pca_variance: Sequence[float]
    normalization_steps: Sequence[str]
    figures: Mapping[str, str]
    tables: Mapping[str, str]
    comparison: ComparisonSection | None = None
    notes: Sequence[str] = field(default_factory=tuple)


def md_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    """Render a GitHub-flavored Markdown table.

    Args:
        headers: Column headers.
        rows: Row values, converted with ``str``. Pipe characters are escaped.

    Returns:
        The table as a string.
    """

    def cell(value: object) -> str:
        return str(value).replace("|", "\\|")

    lines = [
        "| " + " | ".join(cell(h) for h in headers) + " |",
        "|" + "|".join("---" for _ in headers) + "|",
    ]
    lines += ["| " + " | ".join(cell(v) for v in row) + " |" for row in rows]
    return "\n".join(lines)


def _fmt_p(value: float) -> str:
    return (
        "n/a"
        if value is None or (isinstance(value, float) and math.isnan(value))
        else f"{value:.2e}"
    )


def _figure(figures: Mapping[str, str], key: str, alt: str) -> str:
    return f"![{alt}]({figures[key]})" if key in figures else ""


def _qc_section(ctx: ReportContext) -> list[str]:
    variance = ", ".join(f"PC{i + 1} {v:.1%}" for i, v in enumerate(ctx.pca_variance))
    return [
        "## 2. Data quality control",
        "",
        md_table(["Measure", "Value"], ctx.qc.as_rows()),
        "",
        f"Validation: **{'passed' if ctx.validation.is_valid else 'failed'}** with "
        f"{len(ctx.validation.errors)} error(s) and {len(ctx.validation.warnings)} warning(s).",
        *[f"- `{i.code}`: {i.message}" for i in ctx.validation.issues],
        "",
        f"PCA explained variance: {variance}. PCA and heatmaps use mean-imputed, per-gene "
        "scaled values. These are for visual QC only and are not used for testing.",
        "",
        _figure(ctx.figures, "pca", "PCA of samples colored by group and batch"),
        "",
        _figure(ctx.figures, "distributions", "Expression distribution per sample"),
        "",
        _figure(ctx.figures, "sample_totals", "Per-sample total of observed values"),
        "",
        _figure(ctx.figures, "correlation", "Sample-to-sample correlation heatmap"),
        "",
    ]


def _comparison_section(section: ComparisonSection, figures: Mapping[str, str]) -> list[str]:
    c, r, results = section.comparison, section.ranking, section.results
    tested = results["p_value"].notna()
    flagged = results["meets_thresholds"].astype(bool)
    n_higher = int((flagged & (results["direction"] == DIRECTION_HIGHER)).sum())
    n_lower = int((flagged & (results["direction"] == DIRECTION_LOWER)).sum())
    top = results[results["rank"].notna()].nsmallest(r.top_n, "rank")
    rows = [
        (
            int(row["rank"]),
            row["gene_id"],
            f"{row['mean_diff']:+.3f}",
            f"{row['cohens_d']:+.2f}",
            _fmt_p(row["p_value"]),
            _fmt_p(row["p_adj"]),
            f"{row['ranking_score']:.2f}",
            "yes" if row["meets_thresholds"] else "no",
        )
        for _, row in top.iterrows()
    ]
    lines = [
        f"## 4. Group comparison: {c.group_b} vs {c.group_a}",
        "",
        "Per-gene statistics (every difference is "
        f"**{c.group_b} minus {c.group_a}**): group means, mean difference, Cohen's d "
        "(pooled SD), and Welch's two-sample t-test. Raw p-values are adjusted across all "
        "tested genes with the Benjamini-Hochberg procedure.",
        "",
        md_table(
            ["Setting", "Value"],
            [
                ("Minimum observed values per group", c.min_non_missing_per_group),
                ("Adjusted p-value threshold", f"{c.fdr_threshold:g}"),
                ("|Cohen's d| threshold", f"{c.effect_size_threshold:g}"),
                ("Ranking score", f"abs(d) x -log10(max(p_adj, {r.p_floor:g}))"),
            ],
        ),
        "",
        f"- Genes tested: {int(tested.sum())} of {len(results)}",
        f"- Genes meeting both display thresholds: **{int(flagged.sum())}** "
        f"({n_higher} higher in {c.group_b}, {n_lower} lower in {c.group_b})",
        "",
        "The thresholds are display settings chosen by the user. Meeting them means only that "
        "a simulated difference was detected by this workflow.",
        "",
        _figure(figures, "volcano", "Volcano-style plot"),
        "",
        f"### Top {len(top)} genes by ranking score",
        "",
        md_table(
            [
                "Rank",
                "Gene",
                "Mean diff",
                "Cohen's d",
                "p",
                "BH adj. p",
                "Score",
                "Meets thresholds",
            ],
            rows,
        ),
        "",
        _figure(figures, "top_genes", "Top-ranked genes by effect size"),
        "",
        _figure(figures, "top_heatmap", "Heatmap of top-ranked genes"),
        "",
    ]
    if section.recovery is not None:
        lines += [
            "## 5. Workflow check against the synthetic ground truth",
            "",
            "The demo data were simulated with a known set of planted differences, so we can "
            "check whether the workflow recovers them. This verifies the code only. It says "
            "nothing about performance on real data. High recovery is expected whenever the "
            "planted shifts are large relative to the simulated noise, as in the default "
            "configuration.",
            "",
            md_table(["Check", "Result"], section.recovery.as_rows()),
            "",
        ]
    return lines


def render_markdown_report(ctx: ReportContext) -> str:
    """Render the full report.

    Args:
        ctx: Report inputs.

    Returns:
        Markdown text ending with a newline.
    """
    lines = [
        f"# Analysis report - {DATA_LABEL}",
        "",
        "> **RESEARCH AND EDUCATION ONLY - NOT FOR CLINICAL USE.**",
        f"> {FULL_DISCLAIMER}",
        "",
        "## 1. Run information",
        "",
        md_table(
            ["Item", "Value"],
            [
                ("Generated (UTC)", ctx.generated_at),
                ("Configuration", f"`{ctx.config_path}`"),
                ("Seed", ctx.seed),
                *[(f"Input: {role}", f"`{path}`") for role, path in ctx.inputs.items()],
            ],
        ),
        "",
        *_qc_section(ctx),
        "## 3. Normalization",
        "",
        *[f"- {step}" for step in ctx.normalization_steps],
        "",
    ]
    if ctx.comparison is not None:
        lines += _comparison_section(ctx.comparison, ctx.figures)
    lines += [
        "## Output files",
        "",
        md_table(["Table", "File"], [(k, f"`{v}`") for k, v in ctx.tables.items()]),
        "",
        "Every file is labeled DEMONSTRATION / SYNTHETIC. `run_manifest.json` records "
        "checksums, parameters, and software versions.",
        "",
        "## Interpretation limits",
        "",
        "- All values are simulated. Gene identifiers do not correspond to real genes, and "
        "the group labels are not diagnoses.",
        "- The Welch t-test on continuous values is a teaching simplification. Real "
        "transcriptomics data need count-aware models and design-aware analysis.",
        "- On simulated data, recovering the planted differences shows that the code works "
        "as designed. It is not evidence of biological relevance or clinical validity.",
        "- No output of this workbench identifies or validates a biomarker.",
        *[f"- {note}" for note in ctx.notes],
        "",
        "---",
        "",
        f"*{DATA_LABEL}. Research and education only - not for clinical use.*",
        "",
    ]
    return "\n".join(line for line in lines if line is not None)
