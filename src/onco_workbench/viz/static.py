"""Static matplotlib figures for QC and group-comparison results.

Figures are built with the object-oriented API (``matplotlib.figure.Figure``), not
``pyplot``, so they need no display backend and hold no global state. Every figure
carries the DEMONSTRATION / SYNTHETIC label as a footer.

Colors come from a CVD-validated categorical palette, used in fixed order and never
cycled. Groups also get distinct marker shapes, so identity never depends on color
alone. Directions use a blue/red diverging pair, with neutral gray for "not meeting
thresholds".
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from onco_workbench.analysis.normalization import prepare_for_plots
from onco_workbench.analysis.qc import PCAResult
from onco_workbench.analysis.ranking import DIRECTION_HIGHER, DIRECTION_LOWER
from onco_workbench.data.synthetic import BATCH_COLUMN, GROUP_COLUMN
from onco_workbench.disclaimers import DATA_LABEL

SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRIDLINE = "#e1e0d9"
NEUTRAL = "#b8b6ae"
DIVERGING_MID = "#f0efec"
COLOR_LOWER = "#2a78d6"  # blue pole: lower in group B
COLOR_HIGHER = "#e34948"  # red pole: higher in group B

# Categorical slots in fixed order (validated palette). Groups use slots 1-2 onward.
GROUP_PALETTE: tuple[str, ...] = (
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#e87ba4",
    "#008300",
    "#4a3aa7",
    "#e34948",
)
# Batches are a different entity from groups, so they use a separately validated pair.
BATCH_PALETTE: tuple[str, ...] = ("#1baf7a", "#4a3aa7", "#eda100", "#e87ba4", "#008300")
MARKERS: tuple[str, ...] = ("o", "^", "s", "D", "v", "P", "X", "*")

DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "obw_diverging", [COLOR_LOWER, DIVERGING_MID, COLOR_HIGHER]
)
FOOTER = f"{DATA_LABEL} - research and education only, not for clinical use"


def category_styles(labels: Sequence[str], palette: Sequence[str]) -> dict[str, tuple[str, str]]:
    """Map sorted unique labels to ``(color, marker)`` in fixed palette order.

    Labels beyond the palette length fall back to neutral gray with their own marker,
    so hues are never generated or cycled.

    Args:
        labels: Category labels (duplicates allowed).
        palette: Ordered hex colors.

    Returns:
        Mapping of label to ``(color, marker)``.
    """
    unique = sorted({str(label) for label in labels})
    return {
        label: (
            palette[i] if i < len(palette) else NEUTRAL,
            MARKERS[i % len(MARKERS)],
        )
        for i, label in enumerate(unique)
    }


def _new_figure(width: float, height: float) -> Figure:
    return Figure(figsize=(width, height), facecolor=SURFACE, layout="constrained")


def _style_axes(ax: Axes, *, grid: bool = True) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(TEXT_SECONDARY)
        ax.spines[side].set_linewidth(0.6)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=8, width=0.6)
    ax.xaxis.label.set_color(TEXT_PRIMARY)
    ax.yaxis.label.set_color(TEXT_PRIMARY)
    ax.title.set_color(TEXT_PRIMARY)
    if grid:
        ax.grid(color=GRIDLINE, linewidth=0.6)
        ax.set_axisbelow(True)


def _finish(fig: Figure, title: str) -> Figure:
    fig.suptitle(title, color=TEXT_PRIMARY, fontsize=11, x=0.01, ha="left")
    fig.supxlabel(FOOTER, color=TEXT_SECONDARY, fontsize=7, x=0.01, ha="left")
    return fig


def _legend(
    ax: Axes,
    handles: list,
    *,
    title: str | None = None,
    loc: str = "best",
    outside: bool = False,
) -> None:
    placement = {"loc": "upper left", "bbox_to_anchor": (1.01, 1.0)} if outside else {"loc": loc}
    legend = ax.legend(
        handles=handles, title=title, fontsize=8, title_fontsize=8, frameon=False, **placement
    )
    for text in legend.get_texts():
        text.set_color(TEXT_PRIMARY)


def _ordered_samples(metadata: pd.DataFrame, samples: Sequence[str]) -> list[str]:
    frame = metadata.loc[list(samples)].reset_index()
    id_column = frame.columns[0]
    return frame.sort_values([GROUP_COLUMN, id_column], kind="stable")[id_column].tolist()


def _group_strip(ax: Axes, groups: Sequence[str], styles: dict[str, tuple[str, str]]) -> None:
    order = list(styles)
    codes = np.array([[order.index(str(g)) for g in groups]])
    ax.imshow(
        codes,
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap([styles[label][0] for label in order]),
        vmin=0,
        vmax=max(len(order) - 1, 1),
    )
    ax.set_yticks([0], labels=["group"], fontsize=7, color=TEXT_SECONDARY)
    # The strip shares its x-axis with the heatmap; hide its ticks without removing
    # the heatmap's tick labels.
    ax.tick_params(axis="x", length=0, labelbottom=False, labeltop=False)
    ax.tick_params(axis="y", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_pca(pca: PCAResult) -> Figure:
    """Scatter of PC1 vs PC2, colored by group and, when present, by batch.

    Args:
        pca: Result of :func:`onco_workbench.analysis.qc.run_pca`.

    Returns:
        The figure.
    """
    scores = pca.scores
    has_batch = BATCH_COLUMN in scores.columns and scores[BATCH_COLUMN].notna().any()
    fig = _new_figure(10.0 if has_batch else 5.5, 4.6)
    axes = fig.subplots(1, 2 if has_batch else 1, squeeze=False)[0]
    ratio = pca.explained_variance_ratio
    panels = [(GROUP_COLUMN, GROUP_PALETTE, "Colored by group")]
    if has_batch:
        panels.append((BATCH_COLUMN, BATCH_PALETTE, "Colored by batch"))
    for ax, (column, palette, subtitle) in zip(axes, panels, strict=True):
        styles = category_styles(scores[column].dropna(), palette)
        handles = []
        for label, (color, marker) in styles.items():
            subset = scores[scores[column] == label]
            ax.scatter(
                subset["PC1"],
                subset["PC2"],
                s=30,
                c=color,
                marker=marker,
                edgecolors=SURFACE,
                linewidths=0.6,
            )
            handles.append(
                Line2D(
                    [],
                    [],
                    marker=marker,
                    linestyle="",
                    color=color,
                    markersize=7,
                    label=f"{label} (n={len(subset)})",
                )
            )
        _style_axes(ax)
        ax.set_xlabel(f"PC1 ({ratio[0]:.1%} of variance)")
        ax.set_ylabel(f"PC2 ({ratio[1]:.1%} of variance)" if len(ratio) > 1 else "PC2")
        ax.set_title(subtitle, fontsize=9, loc="left")
        _legend(ax, handles)
    scaling = "z-scored" if pca.standardized else "centered"
    return _finish(fig, f"PCA of samples ({pca.n_genes_used} genes, {scaling}, mean-imputed)")


def plot_expression_distributions(expression: pd.DataFrame, metadata: pd.DataFrame) -> Figure:
    """Box plot of expression values per sample, ordered and colored by group.

    Args:
        expression: Samples as rows, genes as columns.
        metadata: Metadata indexed by sample ID (aligned to ``expression``).

    Returns:
        The figure.
    """
    order = _ordered_samples(metadata, expression.index)
    groups = metadata.loc[order, GROUP_COLUMN].astype(str)
    styles = category_styles(groups, GROUP_PALETTE)
    data = [expression.loc[s].dropna().to_numpy() for s in order]
    fig = _new_figure(10.0, 4.2)
    ax = fig.subplots()
    boxes = ax.boxplot(
        data,
        patch_artist=True,
        showfliers=False,
        widths=0.7,
        medianprops={"color": TEXT_PRIMARY, "linewidth": 1.0},
        whiskerprops={"color": TEXT_SECONDARY, "linewidth": 0.6},
        capprops={"color": TEXT_SECONDARY, "linewidth": 0.6},
        boxprops={"linewidth": 0.4, "edgecolor": SURFACE},
    )
    for patch, group in zip(boxes["boxes"], groups, strict=True):
        patch.set_facecolor(styles[group][0])
    _style_axes(ax)
    ax.grid(axis="x", visible=False)
    ax.set_xticks([])
    ax.set_xlabel(f"Samples (n={len(order)}), ordered by group")
    ax.set_ylabel("Expression value (log2-like units)")
    _legend(ax, [Patch(color=c, label=g) for g, (c, _) in styles.items()], outside=True)
    return _finish(fig, "Expression distribution per sample (boxes: IQR; whiskers: 1.5 x IQR)")


def plot_sample_totals(sample_table: pd.DataFrame) -> Figure:
    """Bar chart of the per-sample sum of observed values, colored by group.

    Args:
        sample_table: Output of :func:`onco_workbench.analysis.qc.sample_qc_table`.

    Returns:
        The figure.
    """
    table = sample_table.sort_values([GROUP_COLUMN, "sample_id"], kind="stable")
    styles = category_styles(table[GROUP_COLUMN], GROUP_PALETTE)
    fig = _new_figure(10.0, 3.8)
    ax = fig.subplots()
    ax.bar(
        np.arange(len(table)),
        table["total_expression"],
        width=0.8,
        color=[styles[str(g)][0] for g in table[GROUP_COLUMN]],
    )
    _style_axes(ax)
    ax.grid(axis="x", visible=False)
    ax.set_xticks([])
    ax.set_xlabel(f"Samples (n={len(table)}), ordered by group")
    ax.set_ylabel("Sum of observed values")
    _legend(ax, [Patch(color=c, label=g) for g, (c, _) in styles.items()], outside=True)
    return _finish(fig, "Per-sample total of observed expression values")


def plot_correlation_heatmap(correlation: pd.DataFrame, metadata: pd.DataFrame) -> Figure:
    """Heatmap of sample-to-sample Pearson correlation with a group color strip.

    Args:
        correlation: Square matrix from :func:`onco_workbench.analysis.qc.sample_correlation`.
        metadata: Metadata indexed by sample ID.

    Returns:
        The figure.
    """
    samples = list(correlation.index)
    groups = metadata.loc[samples, GROUP_COLUMN].astype(str).tolist()
    styles = category_styles(groups, GROUP_PALETTE)
    fig = _new_figure(7.0, 6.6)
    strip_ax, ax = fig.subplots(2, 1, height_ratios=[0.35, 10], sharex=True)
    _group_strip(strip_ax, groups, styles)
    image = ax.imshow(
        correlation.to_numpy(),
        cmap=DIVERGING_CMAP,
        vmin=-1,
        vmax=1,
        aspect="auto",
        interpolation="nearest",
    )
    _style_axes(ax, grid=False)
    if len(samples) <= 40:
        ax.set_xticks(range(len(samples)), labels=samples, rotation=90, fontsize=6)
        ax.set_yticks(range(len(samples)), labels=samples, fontsize=6)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    colorbar = fig.colorbar(image, ax=[strip_ax, ax], shrink=0.8)
    colorbar.set_label("Pearson r (genes centered and scaled)", fontsize=8, color=TEXT_PRIMARY)
    colorbar.ax.tick_params(labelsize=7, colors=TEXT_SECONDARY)
    _legend(
        strip_ax,
        [Patch(color=c, label=g) for g, (c, _) in styles.items()],
        loc="lower left",
    )
    strip_ax.get_legend().set_bbox_to_anchor((0.0, 1.2))
    return _finish(fig, f"Sample correlation ({len(samples)} samples, balanced by group)")


def plot_volcano(
    results: pd.DataFrame,
    *,
    group_a: str,
    group_b: str,
    fdr_threshold: float,
    effect_size_threshold: float,
    n_labels: int = 5,
) -> Figure:
    """Volcano-style plot: Cohen's d against -log10(BH-adjusted p).

    Args:
        results: Flagged, ranked comparison results.
        group_a: Reference group label.
        group_b: Comparison group label.
        fdr_threshold: Adjusted p-value threshold (horizontal line).
        effect_size_threshold: |d| threshold (vertical lines).
        n_labels: Number of top-ranked genes labeled directly.

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
    fig = _new_figure(6.8, 5.2)
    ax = fig.subplots()
    handles = []
    for mask, color, label in categories:
        subset = tested[mask]
        ax.scatter(
            subset["cohens_d"],
            subset["neg_log_p"],
            s=16,
            c=color,
            edgecolors=SURFACE,
            linewidths=0.4,
        )
        handles.append(
            Line2D(
                [],
                [],
                marker="o",
                linestyle="",
                color=color,
                markersize=6,
                label=f"{label} ({len(subset)})",
            )
        )
    line = {"color": TEXT_SECONDARY, "linewidth": 0.8, "linestyle": (0, (4, 3))}
    ax.axhline(-np.log10(fdr_threshold), **line)
    ax.axvline(-effect_size_threshold, **line)
    ax.axvline(effect_size_threshold, **line)
    # Direct labels for the top-ranked genes, placed inward and staggered per side so
    # they never collide with each other or run off the plot edge.
    labeled = tested.nsmallest(n_labels, "rank").sort_values("neg_log_p", ascending=False)
    placed = {True: 0, False: 0}
    for _, row in labeled.iterrows():
        right = bool(row["cohens_d"] > 0)
        step = placed[right]
        placed[right] += 1
        ax.annotate(
            row["gene_id"],
            (row["cohens_d"], row["neg_log_p"]),
            xytext=(-14 if right else 14, -6 - 13 * step),
            textcoords="offset points",
            ha="right" if right else "left",
            va="center",
            fontsize=7,
            color=TEXT_PRIMARY,
            arrowprops={"arrowstyle": "-", "color": TEXT_SECONDARY, "linewidth": 0.5},
        )
    _style_axes(ax)
    ax.set_xlabel(f"Cohen's d ({group_b} minus {group_a})")
    ax.set_ylabel("-log10(BH-adjusted p-value)")
    _legend(ax, handles, loc="upper center")
    return _finish(
        fig,
        f"Group comparison: {group_b} vs {group_a} "
        f"(lines: adj. p = {fdr_threshold:g}, |d| = {effect_size_threshold:g})",
    )


def plot_top_genes(results: pd.DataFrame, *, group_a: str, group_b: str, top_n: int) -> Figure:
    """Horizontal bar chart of Cohen's d for the top-ranked genes.

    Args:
        results: Ranked comparison results.
        group_a: Reference group label.
        group_b: Comparison group label.
        top_n: Number of genes to show.

    Returns:
        The figure.
    """
    top = results[results["rank"].notna()].nsmallest(top_n, "rank")
    colors = [
        COLOR_HIGHER if d == DIRECTION_HIGHER else COLOR_LOWER if d == DIRECTION_LOWER else NEUTRAL
        for d in top["direction"]
    ]
    fig = _new_figure(6.4, max(3.0, 0.24 * len(top) + 1.4))
    ax = fig.subplots()
    ax.barh(np.arange(len(top)), top["cohens_d"], height=0.7, color=colors)
    ax.axvline(0, color=TEXT_SECONDARY, linewidth=0.8)
    ax.set_yticks(np.arange(len(top)), labels=top["gene_id"], fontsize=7)
    ax.invert_yaxis()
    _style_axes(ax)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel(f"Cohen's d ({group_b} minus {group_a})")
    _legend(
        ax,
        [
            Patch(color=COLOR_HIGHER, label=f"Higher in {group_b}"),
            Patch(color=COLOR_LOWER, label=f"Lower in {group_b}"),
        ],
        loc="lower right",
    )
    return _finish(fig, f"Top {len(top)} genes by ranking score (rank 1 at top)")


def plot_top_gene_heatmap(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    gene_ids: Sequence[str],
    *,
    standardize: bool = True,
) -> Figure:
    """Heatmap of the top-ranked genes across samples ordered by group.

    Values are mean-imputed and z-scored per gene (or centered), then clipped to +/-3.

    Args:
        expression: Samples as rows, genes as columns.
        metadata: Metadata indexed by sample ID.
        gene_ids: Genes to show, top rank first.
        standardize: Z-score genes (True) or only center them (False).

    Returns:
        The figure.
    """
    order = _ordered_samples(metadata, expression.index)
    groups = metadata.loc[order, GROUP_COLUMN].astype(str).tolist()
    styles = category_styles(groups, GROUP_PALETTE)
    values = prepare_for_plots(expression, zscore=standardize).loc[order, list(gene_ids)]
    fig = _new_figure(9.0, max(3.5, 0.2 * len(gene_ids) + 1.8))
    strip_ax, ax = fig.subplots(2, 1, height_ratios=[0.35, 10], sharex=True)
    _group_strip(strip_ax, groups, styles)
    image = ax.imshow(
        np.clip(values.to_numpy().T, -3, 3),
        cmap=DIVERGING_CMAP,
        vmin=-3,
        vmax=3,
        aspect="auto",
        interpolation="nearest",
    )
    _style_axes(ax, grid=False)
    ax.set_yticks(range(len(gene_ids)), labels=list(gene_ids), fontsize=7)
    ax.set_xticks([])
    ax.set_xlabel(f"Samples (n={len(order)}), ordered by group")
    colorbar = fig.colorbar(image, ax=[strip_ax, ax], shrink=0.8)
    label = "Gene z-score (clipped to +/-3)" if standardize else "Centered value (clipped)"
    colorbar.set_label(label, fontsize=8, color=TEXT_PRIMARY)
    colorbar.ax.tick_params(labelsize=7, colors=TEXT_SECONDARY)
    _legend(strip_ax, [Patch(color=c, label=g) for g, (c, _) in styles.items()], loc="lower left")
    strip_ax.get_legend().set_bbox_to_anchor((0.0, 1.2))
    return _finish(fig, f"Top {len(gene_ids)} ranked genes across samples")


def save_figure(fig: Figure, path: str | Path, *, dpi: int = 150) -> Path:
    """Save a figure, creating parent directories. The format follows the file suffix.

    Args:
        fig: Figure to save.
        path: Destination path (e.g. ``outputs/figures/volcano.png``).
        dpi: Resolution for raster formats.

    Returns:
        The written path.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, facecolor=SURFACE)
    return path
