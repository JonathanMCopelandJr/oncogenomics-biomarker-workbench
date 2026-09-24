from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go
import pytest

from onco_workbench.analysis.qc import run_pca, sample_qc_table
from onco_workbench.config import load_config
from onco_workbench.data.synthetic import SyntheticDataset
from onco_workbench.data.validation import align_metadata
from onco_workbench.disclaimers import DATA_LABEL
from onco_workbench.pipeline import run_comparison
from onco_workbench.viz import interactive


def _is_labeled(fig: go.Figure) -> bool:
    return any(DATA_LABEL in (a.text or "") for a in fig.layout.annotations)


@pytest.fixture
def results(small_dataset: SyntheticDataset) -> pd.DataFrame:
    return run_comparison(small_dataset.expression, small_dataset.metadata, load_config()).results


def test_all_figures_are_labeled(small_dataset: SyntheticDataset, results: pd.DataFrame) -> None:
    expr, meta = small_dataset.expression, small_dataset.metadata
    aligned = align_metadata(expr, meta)
    pca = run_pca(expr, meta)
    figures = [
        interactive.pca_figure(pca, "group"),
        interactive.pca_figure(pca, "batch"),
        interactive.expression_box_figure(expr, aligned),
        interactive.missingness_figure(sample_qc_table(expr, meta)),
        interactive.gene_by_group_figure(expr, aligned, expr.columns[0]),
        interactive.volcano_figure(
            results,
            group_a="Group_A",
            group_b="Group_B",
            fdr_threshold=0.05,
            effect_size_threshold=0.5,
        ),
        interactive.top_genes_figure(results, group_a="Group_A", group_b="Group_B", top_n=5),
    ]
    for fig in figures:
        assert isinstance(fig, go.Figure)
        assert _is_labeled(fig)


def test_pca_uses_shapes_and_counts(small_dataset: SyntheticDataset) -> None:
    fig = interactive.pca_figure(run_pca(small_dataset.expression, small_dataset.metadata))
    assert [t.name for t in fig.data] == ["Group_A (n=10)", "Group_B (n=10)"]
    assert [t.marker.symbol for t in fig.data] == ["circle", "triangle-up"]
    assert all("<extra></extra>" in t.hovertemplate for t in fig.data)


def test_pca_invalid_color_column(small_dataset: SyntheticDataset) -> None:
    pca = run_pca(small_dataset.expression, small_dataset.metadata.drop(columns="batch"))
    with pytest.raises(ValueError, match="batch"):
        interactive.pca_figure(pca, "batch")


def test_volcano_trace_counts_match_flags(results: pd.DataFrame) -> None:
    fig = interactive.volcano_figure(
        results, group_a="Group_A", group_b="Group_B", fdr_threshold=0.05, effect_size_threshold=0.5
    )
    total = sum(len(t.x) for t in fig.data)
    assert total == int(results["p_adj"].notna().sum())
    flagged = int(results["meets_thresholds"].sum())
    assert len(fig.data[1].x) + len(fig.data[2].x) == flagged


def test_expression_box_covers_all_observed_values(small_dataset: SyntheticDataset) -> None:
    expr = small_dataset.expression
    fig = interactive.expression_box_figure(expr, align_metadata(expr, small_dataset.metadata))
    assert sum(len(t.y) for t in fig.data) == int(expr.notna().to_numpy().sum())


def test_top_genes_order(results: pd.DataFrame) -> None:
    fig = interactive.top_genes_figure(results, group_a="Group_A", group_b="Group_B", top_n=3)
    assert list(fig.data[0].y) == results.nsmallest(3, "rank")["gene_id"].tolist()


def test_confusion_matrix_figure() -> None:
    confusion = pd.DataFrame(
        [[9, 1], [2, 8]],
        index=["true: Group_A", "true: Group_B"],
        columns=["predicted: Group_A", "predicted: Group_B"],
    )
    fig = interactive.confusion_matrix_figure(confusion, positive_label="Group_B")
    assert _is_labeled(fig)
    assert list(fig.data[0].x) == ["Group_A", "Group_B"]
    assert fig.data[0].z.tolist() == [[9, 1], [2, 8]]
    assert "positive class: Group_B" in fig.layout.title.text
    cell_text = sorted(a.text for a in fig.layout.annotations if a.text.startswith("<b>"))
    assert cell_text == ["<b>1</b>", "<b>2</b>", "<b>8</b>", "<b>9</b>"]
