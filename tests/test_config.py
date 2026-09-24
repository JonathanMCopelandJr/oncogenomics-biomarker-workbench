from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
import yaml

from onco_workbench.config import DEFAULT_CONFIG_PATH, ConfigError, load_config

Mutator = Callable[[dict[str, Any]], object]


def _write_config(tmp_path: Path, mutate: Mutator) -> Path:
    raw = yaml.safe_load(DEFAULT_CONFIG_PATH.read_text(encoding="utf-8"))
    mutate(raw)
    config_dir = tmp_path / "config"
    config_dir.mkdir()
    path = config_dir / "custom.yaml"
    path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    return path


def test_default_config_loads_with_expected_values(repo_root: Path) -> None:
    config = load_config()
    assert config.seed == 42
    assert config.root == repo_root
    assert config.synthetic.n_samples == 80
    assert config.synthetic.n_genes == 500
    assert config.synthetic.group_labels == ("Group_A", "Group_B")
    assert config.synthetic.group_sizes == (40, 40)
    assert config.validation.allowed_groups == ("Group_A", "Group_B")
    assert config.comparison.group_a == "Group_A"
    assert config.comparison.fdr_threshold == 0.05
    assert config.ranking.top_n == 20
    assert config.figures.format == "png"
    assert config.ml_demo.test_size == 0.25
    assert config.ml_demo.min_samples_per_class == 10
    assert config.ml_demo.n_label_permutations == 20
    assert "ml_demo" not in config.sections and "comparison" not in config.sections


def test_with_comparison_revalidates() -> None:
    config = load_config()
    changed = config.with_comparison(fdr_threshold=0.1, group_b="Group_A", group_a="Group_B")
    assert changed.comparison.fdr_threshold == 0.1
    assert config.comparison.fdr_threshold == 0.05
    with pytest.raises(ConfigError, match="must differ"):
        config.with_comparison(group_b="Group_A")
    with pytest.raises(ConfigError, match="fdr_threshold"):
        config.with_comparison(fdr_threshold=0.0)


def test_paths_are_absolute_and_under_root(repo_root: Path) -> None:
    config = load_config()
    assert config.paths.expression_file == repo_root / "data" / "synthetic" / "demo_expression.csv"
    for value in vars(config.paths).values():
        assert value.is_absolute()


def test_relative_paths_resolve_against_config_parent(tmp_path: Path) -> None:
    path = _write_config(tmp_path, lambda raw: None)
    config = load_config(path)
    assert config.root == tmp_path.resolve()
    assert config.paths.outputs_dir == tmp_path.resolve() / "outputs"


def test_with_seed_returns_new_config() -> None:
    config = load_config()
    changed = config.with_seed(7)
    assert changed.seed == 7 and config.seed == 42
    with pytest.raises(ConfigError, match="seed"):
        config.with_seed(-1)


def test_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="not found"):
        load_config(tmp_path / "nope.yaml")


def test_invalid_yaml_raises(tmp_path: Path) -> None:
    path = tmp_path / "bad.yaml"
    path.write_text("seed: [unclosed", encoding="utf-8")
    with pytest.raises(ConfigError, match="not valid YAML"):
        load_config(path)


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda raw: raw.pop("seed"), "seed"),
        (lambda raw: raw.update(seed=-3), "seed must be >= 0"),
        (lambda raw: raw.update(seed=1.5), "seed must be an integer"),
        (lambda raw: raw.pop("synthetic"), "section 'synthetic'"),
        (lambda raw: raw["synthetic"].update(n_sampels=10), "unknown key"),
        (lambda raw: raw["synthetic"].pop("noise_sd"), "missing required key"),
        (lambda raw: raw["synthetic"].update(n_genes="many"), "n_genes must be an integer"),
        (lambda raw: raw["synthetic"].update(n_signal_genes_up=600), "exceeds n_genes"),
        (lambda raw: raw["synthetic"].update(group_proportions=[0.7, 0.7]), "sum to 1"),
        (lambda raw: raw["synthetic"].update(group_labels=["A", "A"]), "two distinct"),
        (lambda raw: raw["synthetic"].update(missing_fraction=0.9), "missing_fraction"),
        (lambda raw: raw["synthetic"].update(gene_id_prefix="G,"), "must not contain"),
        (lambda raw: raw["validation"].update(min_samples_per_group=1), "min_samples_per_group"),
        (lambda raw: raw["validation"].update(allowed_groups=["Only"]), "allowed_groups"),
        (lambda raw: raw["paths"].update(outputs_dir=""), "outputs_dir"),
        (lambda raw: raw["comparison"].update(fdr_method="bonferroni"), "fdr_method"),
        (lambda raw: raw["comparison"].update(effect_size_threshold=-1), "effect_size"),
        (lambda raw: raw["comparison"].update(group_b="Group_A"), "must differ"),
        (lambda raw: raw["normalization"].update(median_center_samples="yes"), "true or false"),
        (lambda raw: raw["ranking"].update(p_floor=0), "p_floor"),
        (lambda raw: raw["ranking"].update(top_n=0), "top_n"),
        (lambda raw: raw["qc"].update(pca_components=1), "pca_components"),
        (lambda raw: raw["figures"].update(format="gif"), "format"),
        (lambda raw: raw["figures"].update(dpi=5), "dpi"),
        (lambda raw: raw.pop("comparison"), "section 'comparison'"),
        (lambda raw: raw.pop("ml_demo"), "section 'ml_demo'"),
        (lambda raw: raw["ml_demo"].update(test_size=1.0), "test_size"),
        (lambda raw: raw["ml_demo"].update(cv_folds=1), "cv_folds"),
        (lambda raw: raw["ml_demo"].update(min_train_per_class_for_cv=3), ">= cv_folds"),
        (lambda raw: raw["ml_demo"].update(logistic_regression_c=0), "logistic_regression_c"),
        (lambda raw: raw["ml_demo"].update(n_label_permutations=-1), "n_label_permutations"),
        (lambda raw: raw["ml_demo"].update(min_samples_per_class=2), "min_samples_per_class"),
        (lambda raw: raw["ml_demo"].update(enabled="yes"), "true or false"),
        (lambda raw: raw["ml_demo"].update(max_iter=10), "max_iter"),
        (lambda raw: raw["ml_demo"].update(seed=1), "unknown key"),
    ],
)
def test_invalid_values_raise_clear_errors(tmp_path: Path, mutate: Mutator, message: str) -> None:
    path = _write_config(tmp_path, mutate)
    with pytest.raises(ConfigError, match=message):
        load_config(path)


def test_synthetic_config_validates_on_direct_construction(
    config_factory: Callable[..., object],
) -> None:
    with pytest.raises(ConfigError, match="n_batches"):
        config_factory(n_batches=50)
    with pytest.raises(ConfigError, match="noise_sd"):
        config_factory(noise_sd=0.0)
    with pytest.raises(ConfigError, match="signal_effect_range"):
        config_factory(signal_effect_range=(2.0, 1.0))
