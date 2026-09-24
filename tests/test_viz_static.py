from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from matplotlib.figure import Figure
from matplotlib.text import Text

from onco_workbench.analysis.qc import run_pca, sample_correlation, sample_qc_table
from onco_workbench.config import load_config
from onco_workbench.data.synthetic import SyntheticDataset
from onco_workbench.data.validation import align_metadata
from onco_workbench.disclaimers import DATA_LABEL
from onco_workbench.pipeline import run_comparison
from onco_workbench.viz import static


def _texts(fig: Figure) -> str:
    return " ".join(t.get_text() for t in fig.findobj(Text))


@pytest.fixture
def results(small_dataset: SyntheticDataset) -> pd.DataFrame:
    return run_comparison(small_dataset.expression, small_dataset.metadata, load_config()).results


def test_category_styles_fixed_order_and_fallback() -> None:
    styles = static.category_styles(["b", "a", "b"], ["#111111", "#222222"])
    assert styles == {"a": ("#111111", "o"), "b": ("#222222", "^")}
    many = static.category_styles(["x", "y", "z"], ["#111111"])
    assert many["z"][0] == static.NEUTRAL


def test_every_figure_is_labeled(small_dataset: SyntheticDataset, results: pd.DataFrame) -> None:
    expr, meta = small_dataset.expression, small_dataset.metadata
    aligned = align_metadata(expr, meta)
    figures = [
        static.plot_pca(run_pca(expr, meta)),
        static.plot_expression_distributions(expr, aligned),
        static.plot_sample_totals(sample_qc_table(expr, meta)),
        static.plot_correlation_heatmap(sample_correlation(expr, list(expr.index)), aligned),
        static.plot_volcano(
            results,
            group_a="Group_A",
            group_b="Group_B",
            fdr_threshold=0.05,
            effect_size_threshold=0.5,
        ),
        static.plot_top_genes(results, group_a="Group_A", group_b="Group_B", top_n=10),
        static.plot_top_gene_heatmap(expr, aligned, results["gene_id"].head(5).tolist()),
    ]
    for fig in figures:
        assert isinstance(fig, Figure)
        assert DATA_LABEL in _texts(fig)


def test_pca_without_batch_has_one_panel(small_dataset: SyntheticDataset) -> None:
    pca = run_pca(small_dataset.expression, small_dataset.metadata.drop(columns="batch"))
    assert len(static.plot_pca(pca).axes) == 1


def test_volcano_legend_counts_match_flags(results: pd.DataFrame) -> None:
    fig = static.plot_volcano(
        results, group_a="Group_A", group_b="Group_B", fdr_threshold=0.05, effect_size_threshold=0.5
    )
    text = _texts(fig)
    n_not = int((~results["meets_thresholds"]).sum())
    assert f"Not meeting thresholds ({n_not})" in text


def test_save_figure(tmp_path: Path, small_dataset: SyntheticDataset) -> None:
    fig = static.plot_pca(run_pca(small_dataset.expression, small_dataset.metadata))
    path = static.save_figure(fig, tmp_path / "nested" / "pca.png", dpi=60)
    assert path.is_file() and path.stat().st_size > 1000
    assert path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
