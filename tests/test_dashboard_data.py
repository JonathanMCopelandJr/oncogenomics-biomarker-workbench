from __future__ import annotations

import io

import pandas as pd
import pytest

from onco_workbench.config import ConfigError, load_config
from onco_workbench.dashboard.data import (
    DemoData,
    compare,
    environment_info,
    filter_results,
    gene_group_summary,
    load_demo_data,
    results_csv_bytes,
    search_genes,
    summary_statistics,
)
from onco_workbench.disclaimers import DATA_LABEL
from onco_workbench.pipeline import run_comparison


@pytest.fixture(scope="module")
def demo() -> DemoData:
    return load_demo_data()


def test_load_demo_data(demo: DemoData) -> None:
    assert demo.expression.shape == (80, 500)
    assert demo.validation.is_valid
    assert demo.groups == ["Group_A", "Group_B"]
    assert list(demo.aligned.index) == list(demo.expression.index)
    assert demo.ground_truth is not None
    assert demo.data_manifest is not None and demo.data_manifest["license"]["spdx"] == "CC0-1.0"


def test_compare_matches_pipeline(demo: DemoData) -> None:
    settings, result = compare(
        demo, group_a="Group_A", group_b="Group_B", fdr_threshold=0.05, effect_size_threshold=0.5
    )
    reference = run_comparison(demo.expression, demo.metadata, load_config(), demo.ground_truth)
    pd.testing.assert_frame_equal(result.results, reference.results)
    assert settings.group_b == "Group_B"
    assert result.recovery == reference.recovery


def test_compare_thresholds_change_flags_not_statistics(demo: DemoData) -> None:
    _, loose = compare(
        demo, group_a="Group_A", group_b="Group_B", fdr_threshold=0.2, effect_size_threshold=0.0
    )
    _, strict = compare(
        demo, group_a="Group_A", group_b="Group_B", fdr_threshold=0.001, effect_size_threshold=2.0
    )
    assert strict.results["meets_thresholds"].sum() <= loose.results["meets_thresholds"].sum()
    by_gene = strict.results.set_index("gene_id")["p_adj"]
    assert by_gene.equals(loose.results.set_index("gene_id")["p_adj"].loc[by_gene.index])


def test_swapped_groups_skip_ground_truth_check(demo: DemoData) -> None:
    _, result = compare(
        demo, group_a="Group_B", group_b="Group_A", fdr_threshold=0.05, effect_size_threshold=0.5
    )
    assert result.recovery is None


def test_compare_rejects_identical_groups(demo: DemoData) -> None:
    with pytest.raises(ConfigError, match="must differ"):
        compare(
            demo,
            group_a="Group_A",
            group_b="Group_A",
            fdr_threshold=0.05,
            effect_size_threshold=0.5,
        )


def test_search_genes() -> None:
    genes = ["SYN_G0001", "SYN_G0002", "SYN_G0010"]
    assert search_genes(genes, "") == genes
    assert search_genes(genes, " syn_g000 ") == ["SYN_G0001", "SYN_G0002"]
    assert search_genes(genes, "zzz") == []
    assert search_genes(genes, "", limit=1) == ["SYN_G0001"]


def test_filter_results() -> None:
    results = pd.DataFrame(
        {
            "gene_id": ["SYN_G0001", "SYN_G0002", "SYN_G0013"],
            "meets_thresholds": [True, False, True],
        }
    )
    assert filter_results(results)["gene_id"].tolist() == results["gene_id"].tolist()
    assert filter_results(results, query="g000")["gene_id"].tolist() == ["SYN_G0001", "SYN_G0002"]
    assert filter_results(results, query="13")["gene_id"].tolist() == ["SYN_G0013"]
    assert filter_results(results, only_meeting_thresholds=True)["gene_id"].tolist() == [
        "SYN_G0001",
        "SYN_G0013",
    ]
    assert filter_results(results, query="0002", only_meeting_thresholds=True).empty


def test_results_csv_bytes_is_labeled_and_parseable(demo: DemoData) -> None:
    settings, result = compare(
        demo, group_a="Group_A", group_b="Group_B", fdr_threshold=0.05, effect_size_threshold=0.5
    )
    payload = results_csv_bytes(result.results.head(5), settings)
    text = payload.decode("utf-8")
    assert text.startswith(f"# {DATA_LABEL}")
    assert "# Comparison: Group_B (b) vs Group_A (a)" in text
    assert "\r\n" not in text
    parsed = pd.read_csv(io.StringIO(text), comment="#")
    assert parsed["gene_id"].tolist() == result.results["gene_id"].head(5).tolist()


def test_summary_statistics(demo: DemoData) -> None:
    stats = summary_statistics(demo.expression)
    assert stats.loc["count"].iloc[0] == demo.expression.notna().to_numpy().sum()


def test_gene_group_summary(demo: DemoData) -> None:
    summary = gene_group_summary(demo, "SYN_G0001")
    assert summary["group"].tolist() == ["Group_A", "Group_B"]
    assert (summary["n_observed"] + summary["n_missing"]).tolist() == [40, 40]
    with pytest.raises(KeyError, match="Unknown gene"):
        gene_group_summary(demo, "NOT_A_GENE")


def test_environment_info(demo: DemoData) -> None:
    info = dict(environment_info(demo))
    assert info["Seed"] == "42"
    assert info["Configuration"] == "config/default.yaml"
    assert info["Demo data license"] == "CC0-1.0"
    assert "streamlit version" in info
