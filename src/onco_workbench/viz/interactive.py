"""Interactive plotly figures for the dashboard.

These use the same CVD-validated palette and conventions as :mod:`onco_workbench.viz.static`:

- categorical hues in fixed order, with marker shapes as a second identity cue
- blue/red for lower/higher in group B, and neutral gray for "not meeting thresholds"
- hover tooltips on every mark
- a DEMONSTRATION / SYNTHETIC annotation on every figure

Figures leave the background transparent, so Streamlit's light and dark themes both apply.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from onco_workbench.analysis.qc import PCAResult
from onco_workbench.analysis.ranking import DIRECTION_HIGHER, DIRECTION_LOWER
from onco_workbench.data.synthetic import BATCH_COLUMN, GROUP_COLUMN
from onco_workbench.disclaimers import DATA_LABEL
from onco_workbench.viz.static import (
    BATCH_PALETTE,
    COLOR_HIGHER,
    COLOR_LOWER,
    GROUP_PALETTE,
    NEUTRAL,
    category_styles,
)

# Plotly marker names equivalent to the matplotlib markers used in static figures.
_PLOTLY_MARKERS = {
    "o": "circle",
    "^": "triangle-up",
    "s": "square",
    "D": "diamond",
    "v": "triangle-down",
    "P": "cross",
    "X": "x",
    "*": "star",
}


def _label(fig: go.Figure, title: str, *, height: int = 460) -> go.Figure:
    fig.update_layout(
        title={"text": title, "x": 0.0, "xanchor": "left"},
        height=height,
        margin={"l": 60, "r": 20, "t": 60, "b": 70},
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.0, "x": 1.0, "xanchor": "right"},
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        hoverlabel={"namelength": -1},
    )
    fig.add_annotation(
        text=f"{DATA_LABEL} - research and education only, not for clinical use",
        xref="paper",
        yref="paper",
        x=0.0,
        y=-0.2,
        xanchor="left",
        showarrow=False,
        font={"size": 10},
        opacity=0.75,
    )
    return fig


def pca_figure(pca: PCAResult, color_by: str = GROUP_COLUMN) -> go.Figure:
    """PC1 vs PC2 scatter colored (and shaped) by group or batch.

    Args:
        pca: Result of :func:`onco_workbench.analysis.qc.run_pca`.
        color_by: ``"group"`` or ``"batch"``.

    Returns:
        The figure.

    Raises:
        ValueError: If ``color_by`` is not an available column.
    """
    scores = pca.scores
    if color_by not in (GROUP_COLUMN, BATCH_COLUMN) or color_by not in scores.columns:
        raise ValueError(f"Cannot color PCA by {color_by!r}.")
    palette = GROUP_PALETTE if color_by == GROUP_COLUMN else BATCH_PALETTE
    styles = category_styles(scores[color_by].dropna(), palette)
    ratio = pca.explained_variance_ratio
    fig = go.Figure()
    extra = [c for c in (GROUP_COLUMN, BATCH_COLUMN) if c in scores.columns]
    for label, (color, marker) in styles.items():
        subset = scores[scores[color_by] == label]
        fig.add_trace(
            go.Scatter(
                x=subset["PC1"],
                y=subset["PC2"],
                mode="markers",
                name=f"{label} (n={len(subset)})",
                marker={
                    "color": color,
                    "symbol": _PLOTLY_MARKERS[marker],
                    "size": 10,
                    "line": {"width": 1, "color": "white"},
                },
                customdata=subset[["sample_id", *extra]].to_numpy(),
                hovertemplate="<b>%{customdata[0]}</b><br>"
                + "".join(f"{c}: %{{customdata[{i + 1}]}}<br>" for i, c in enumerate(extra))
                + "PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>",
            )
        )
    fig.update_xaxes(title=f"PC1 ({ratio[0]:.1%} of variance)", zeroline=False)
    fig.update_yaxes(
        title=f"PC2 ({ratio[1]:.1%} of variance)" if len(ratio) > 1 else "PC2", zeroline=False
    )
    return _label(fig, f"PCA of samples, colored by {color_by}")


def expression_box_figure(expression: pd.DataFrame, aligned: pd.DataFrame) -> go.Figure:
    """One box per sample (samples ordered by group), colored by group.

    Args:
        expression: Samples as rows, genes as columns.
        aligned: Metadata indexed by sample ID in expression-row order.

    Returns:
        The figure.
    """
    groups = aligned[GROUP_COLUMN].astype(str)
    order = groups.sort_values(kind="stable").index.tolist()
    styles = category_styles(groups, GROUP_PALETTE)
    fig = go.Figure()
    for label, (color, _) in styles.items():
        samples = [s for s in order if groups[s] == label]
        block = expression.loc[samples]
        # Long form built directly from the array (row-major: each sample's genes in turn).
        long = pd.DataFrame(
            {
                "sample_id": np.repeat(block.index.to_numpy(), block.shape[1]),
                "value": block.to_numpy(dtype=float).ravel(),
            }
        ).dropna()
        fig.add_trace(
            go.Box(
                x=long["sample_id"],
                y=long["value"],
                name=label,
                marker_color=color,
                boxpoints=False,
                line={"width": 1},
            )
        )
    fig.update_xaxes(
        title=f"Samples (n={len(order)}), ordered by group",
        showticklabels=False,
        categoryorder="array",
        categoryarray=order,
    )
    fig.update_yaxes(title="Expression value (log2-like units)")
    return _label(fig, "Expression distribution per sample", height=420)


def missingness_figure(sample_table: pd.DataFrame) -> go.Figure:
    """Bar chart of missing values per sample, colored by group.

    Args:
        sample_table: Output of :func:`onco_workbench.analysis.qc.sample_qc_table`.

    Returns:
        The figure.
    """
    table = sample_table.sort_values([GROUP_COLUMN, "sample_id"], kind="stable")
    styles = category_styles(table[GROUP_COLUMN], GROUP_PALETTE)
    fig = go.Figure()
    for label, (color, _) in styles.items():
        subset = table[table[GROUP_COLUMN] == label]
        fig.add_trace(
            go.Bar(
                x=subset["sample_id"],
                y=subset["n_missing"],
                name=label,
                marker_color=color,
                hovertemplate="<b>%{x}</b><br>missing values: %{y}<extra></extra>",
            )
        )
    fig.update_layout(bargap=0.15)
    fig.update_xaxes(
        title="Samples, ordered by group",
        showticklabels=False,
        categoryorder="array",
        categoryarray=table["sample_id"].tolist(),
    )
    fig.update_yaxes(title="Missing values per sample", rangemode="tozero")
    return _label(fig, "Missing values per sample", height=360)


def gene_by_group_figure(
    expression: pd.DataFrame, aligned: pd.DataFrame, gene_id: str
) -> go.Figure:
    """Box plus individual points for one gene, by group.

    Args:
        expression: Samples as rows, genes as columns.
        aligned: Metadata indexed by sample ID in expression-row order.
        gene_id: Gene to show.

    Returns:
        The figure.
    """
    groups = aligned[GROUP_COLUMN].astype(str)
    styles = category_styles(groups, GROUP_PALETTE)
    fig = go.Figure()
    for label, (color, _) in styles.items():
        values = expression.loc[groups == label, gene_id]
        fig.add_trace(
            go.Box(
                y=values,
                name=label,
                marker_color=color,
                boxpoints="all",
                jitter=0.4,
                pointpos=0,
                text=values.index.tolist(),
                hovertemplate="<b>%{text}</b><br>value: %{y:.3f}<extra></extra>",
            )
        )
    fig.update_layout(showlegend=False)
    fig.update_yaxes(title="Expression value (log2-like units)")
    return _label(fig, f"{gene_id} by group", height=380)


def volcano_figure(
    results: pd.DataFrame,
    *,
    group_a: str,
    group_b: str,
    fdr_threshold: float,
    effect_size_threshold: float,
) -> go.Figure:
    """Interactive volcano-style plot: Cohen's d vs -log10(BH-adjusted p).

    Args:
        results: Flagged, ranked comparison results.
        group_a: Reference group.
        group_b: Comparison group.
        fdr_threshold: Adjusted p-value threshold (horizontal line).
        effect_size_threshold: |d| threshold (vertical lines).

    Returns:
        The figure.
    """
    tested = results[results["p_adj"].notna()].copy()
    tested["neg_log_p"] = -np.log10(np.maximum(tested["p_adj"].to_numpy(dtype=float), 1e-300))
    flagged = tested["meets_thresholds"].astype(bool)
    categories = [
        (~flagged, NEUTRAL, "Not meeting thresholds"),
        (flagged & (tested["direction"] == DIRECTION_LOWER), COLOR_LOWER, f"Lower in {group_b}"),
        (flagged & (tested["direction"] == DIRECTION_HIGHER), COLOR_HIGHER, f"Higher in {group_b}"),
    ]
    fig = go.Figure()
    for mask, color, label in categories:
        subset = tested[mask]
        fig.add_trace(
            go.Scatter(
                x=subset["cohens_d"],
                y=subset["neg_log_p"],
                mode="markers",
                name=f"{label} ({len(subset)})",
                marker={"color": color, "size": 8, "line": {"width": 0.5, "color": "white"}},
                customdata=subset[["gene_id", "mean_diff", "p_value", "p_adj", "rank"]]
                .astype(object)
                .to_numpy(),
                hovertemplate="<b>%{customdata[0]}</b><br>"
                "rank: %{customdata[4]}<br>"
                "Cohen's d: %{x:.2f}<br>"
                "mean difference: %{customdata[1]:.3f}<br>"
                "p: %{customdata[2]:.2e}<br>"
                "BH adj. p: %{customdata[3]:.2e}<extra></extra>",
            )
        )
    line = {"dash": "dash", "width": 1, "color": "gray"}
    fig.add_hline(y=-np.log10(fdr_threshold), line=line)
    fig.add_vline(x=-effect_size_threshold, line=line)
    fig.add_vline(x=effect_size_threshold, line=line)
    fig.update_xaxes(title=f"Cohen's d ({group_b} minus {group_a})", zeroline=False)
    fig.update_yaxes(title="-log10(BH-adjusted p-value)", zeroline=False)
    return _label(
        fig,
        f"{group_b} vs {group_a} "
        f"(lines: adj. p = {fdr_threshold:g}, |d| = {effect_size_threshold:g})",
        height=520,
    )


def top_genes_figure(results: pd.DataFrame, *, group_a: str, group_b: str, top_n: int) -> go.Figure:
    """Horizontal bars of Cohen's d for the top-ranked genes (rank 1 at top).

    Args:
        results: Ranked comparison results.
        group_a: Reference group.
        group_b: Comparison group.
        top_n: Number of genes to show.

    Returns:
        The figure.
    """
    top = results[results["rank"].notna()].nsmallest(top_n, "rank")
    colors = [
        COLOR_HIGHER if d == DIRECTION_HIGHER else COLOR_LOWER if d == DIRECTION_LOWER else NEUTRAL
        for d in top["direction"]
    ]
    fig = go.Figure(
        go.Bar(
            x=top["cohens_d"],
            y=top["gene_id"],
            orientation="h",
            marker_color=colors,
            customdata=top[["rank", "p_adj"]].astype(object).to_numpy(),
            hovertemplate="<b>%{y}</b><br>rank: %{customdata[0]}<br>Cohen's d: %{x:.2f}<br>"
            "BH adj. p: %{customdata[1]:.2e}<extra></extra>",
            showlegend=False,
        )
    )
    fig.update_yaxes(autorange="reversed", title=None)
    fig.update_xaxes(title=f"Cohen's d ({group_b} minus {group_a})", zeroline=True)
    return _label(
        fig, f"Top {len(top)} genes by ranking score", height=max(320, 22 * len(top) + 140)
    )
