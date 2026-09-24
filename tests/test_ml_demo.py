"""Tests for the educational ML demonstration (synthetic data only)."""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, roc_auc_score

from onco_workbench.config import MLDemoConfig, WorkbenchConfig, load_config
from onco_workbench.dashboard.data import DemoData, load_demo_data
from onco_workbench.data.synthetic import SyntheticDataset
from onco_workbench.disclaimers import DATA_LABEL
from onco_workbench.ml import classifier_demo as ml
from onco_workbench.ml.classifier_demo import (
    ClassLabels,
    InsufficientDataError,
    MLDemoResult,
    check_class_counts,
    compute_metrics,
    cross_validate_training,
    output_statements,
    permutation_baseline,
    prepare_features_and_labels,
    run_ml_demo,
    split_samples,
    write_ml_demo_outputs,
)

REQUIRED_STATEMENTS = (
    "Positive synthetic class: Group_B",
    "Dataset: fully synthetic transcriptomics-style demonstration data",
    "Purpose: educational machine-learning workflow demonstration",
    "Not for clinical, diagnostic, prognostic, or treatment use",
)

LABELS = ClassLabels(negative="Group_A", positive="Group_B")


@pytest.fixture(scope="module")
def demo() -> DemoData:
    return load_demo_data()


@pytest.fixture(scope="module")
def demo_result(demo: DemoData) -> MLDemoResult:
    return run_ml_demo(demo.expression, demo.metadata, demo.config)


def _with_ml(config: WorkbenchConfig, **changes: object) -> WorkbenchConfig:
    return dataclasses.replace(config, ml_demo=dataclasses.replace(config.ml_demo, **changes))


# --- deterministic, stratified splitting --------------------------------------------------


def test_split_is_deterministic_for_a_seed(demo: DemoData) -> None:
    _, y = prepare_features_and_labels(demo.expression, demo.metadata, LABELS)
    first = split_samples(y, test_size=0.25, seed=42)
    assert first == split_samples(y, test_size=0.25, seed=42)
    assert first != split_samples(y, test_size=0.25, seed=7)


def test_split_is_stratified_disjoint_and_complete(demo_result: MLDemoResult) -> None:
    r = demo_result
    assert set(r.train_ids).isdisjoint(r.test_ids)
    assert len(r.train_ids) + len(r.test_ids) == r.n_samples == 80
    assert r.train_class_counts == {"Group_A": 30, "Group_B": 30}
    assert r.test_class_counts == {"Group_A": 10, "Group_B": 10}


def test_split_uses_configured_seed(demo: DemoData, demo_result: MLDemoResult) -> None:
    _, y = prepare_features_and_labels(demo.expression, demo.metadata, LABELS)
    assert demo_result.seed == demo.config.seed
    train, test = split_samples(y, test_size=demo.config.ml_demo.test_size, seed=demo.config.seed)
    assert (train, test) == (demo_result.train_ids, demo_result.test_ids)


def test_labels_encode_group_b_as_positive(demo: DemoData) -> None:
    _, y = prepare_features_and_labels(demo.expression, demo.metadata, LABELS)
    groups = demo.aligned["group"]
    assert (y[groups == "Group_B"] == 1).all() and (y[groups == "Group_A"] == 0).all()


# --- no preprocessing leakage --------------------------------------------------------------


def test_imputer_and_scaler_learn_from_training_rows_only(
    demo: DemoData, demo_result: MLDemoResult
) -> None:
    x_train = demo.expression.loc[demo_result.train_ids]
    imputer = demo_result.pipeline.named_steps["impute"]
    scaler = demo_result.pipeline.named_steps["scale"]
    assert np.allclose(imputer.statistics_, x_train.mean().to_numpy())
    assert not np.allclose(imputer.statistics_, demo.expression.mean().to_numpy())
    imputed_train = x_train.fillna(x_train.mean())
    assert np.allclose(scaler.mean_, imputed_train.mean().to_numpy())
    assert np.allclose(scaler.scale_, imputed_train.std(ddof=0).to_numpy())


def test_changing_test_rows_does_not_change_the_fitted_pipeline_or_cv(
    demo: DemoData, demo_result: MLDemoResult
) -> None:
    tampered = demo.expression.copy()
    tampered.loc[demo_result.test_ids] = tampered.loc[demo_result.test_ids] + 1000.0
    tampered.loc[demo_result.test_ids[:3], tampered.columns[:50]] = np.nan
    other = run_ml_demo(tampered, demo.metadata, demo.config)
    assert other.train_ids == demo_result.train_ids
    for step in ("impute", "scale", "model"):
        a, b = demo_result.pipeline.named_steps[step], other.pipeline.named_steps[step]
        for attribute in ("statistics_", "mean_", "scale_", "coef_", "intercept_"):
            if hasattr(a, attribute):
                assert np.array_equal(getattr(a, attribute), getattr(b, attribute)), attribute
    assert other.cv.accuracies == demo_result.cv.accuracies
    assert other.cv.roc_aucs == demo_result.cv.roc_aucs


def test_all_features_used_without_selection(demo_result: MLDemoResult) -> None:
    assert demo_result.n_features == 500
    assert demo_result.pipeline.named_steps["model"].n_features_in_ == 500


def test_cv_folds_are_built_from_training_rows_only(
    monkeypatch: pytest.MonkeyPatch, demo: DemoData, demo_result: MLDemoResult
) -> None:
    seen: list[pd.Index] = []
    original = ml.StratifiedKFold

    class Recording(original):  # type: ignore[misc, valid-type]
        def split(self, x, y=None, groups=None):
            seen.append(x.index)
            return super().split(x, y, groups)

    monkeypatch.setattr(ml, "StratifiedKFold", Recording)
    run_ml_demo(demo.expression, demo.metadata, demo.config)
    assert seen, "cross-validation did not run"
    assert all(set(index) == set(demo_result.train_ids) for index in seen)
    assert sum(demo_result.cv.fold_sizes) == len(demo_result.train_ids)


# --- metric behavior -----------------------------------------------------------------------


def test_metrics_match_hand_counts_and_sklearn() -> None:
    y_true = np.array([1, 1, 1, 0, 0, 0, 0, 1])
    y_pred = np.array([1, 0, 1, 0, 1, 0, 0, 1])
    y_score = np.array([0.9, 0.4, 0.8, 0.1, 0.6, 0.2, 0.3, 0.7])
    m = compute_metrics(y_true, y_pred, y_score, LABELS)
    assert (m.accuracy, m.precision, m.recall, m.f1) == (0.75, 0.75, 0.75, 0.75)
    assert m.accuracy == accuracy_score(y_true, y_pred)
    assert m.precision == precision_score(y_true, y_pred)
    assert m.recall == recall_score(y_true, y_pred)
    assert m.f1 == pytest.approx(f1_score(y_true, y_pred))
    assert m.roc_auc == pytest.approx(roc_auc_score(y_true, y_score))
    assert m.notes == ()


def test_confusion_matrix_layout_and_labels() -> None:
    m = compute_metrics(
        np.array([0, 0, 1, 1, 1]), np.array([0, 1, 1, 1, 0]), np.zeros(5) + 0.5, LABELS
    )
    assert list(m.confusion.index) == ["true: Group_A", "true: Group_B"]
    assert list(m.confusion.columns) == ["predicted: Group_A", "predicted: Group_B"]
    assert m.confusion.to_numpy().tolist() == [[1, 1], [1, 2]]  # [[TN, FP], [FN, TP]]


def test_roc_auc_undefined_for_single_class_test_set() -> None:
    m = compute_metrics(np.array([1, 1, 1]), np.array([1, 0, 1]), np.array([0.9, 0.2, 0.8]), LABELS)
    assert m.roc_auc is None
    assert any("ROC-AUC undefined" in note for note in m.notes)


def test_precision_and_f1_undefined_without_positive_predictions() -> None:
    m = compute_metrics(np.array([1, 0, 1, 0]), np.zeros(4, dtype=int), np.full(4, 0.2), LABELS)
    assert m.precision is None and m.f1 is None
    assert m.recall == 0.0
    assert any("Precision undefined" in note for note in m.notes)


def test_recall_undefined_without_positive_samples() -> None:
    m = compute_metrics(np.zeros(4, dtype=int), np.array([0, 1, 0, 0]), np.full(4, 0.4), LABELS)
    assert m.recall is None and m.f1 is None and m.roc_auc is None
    assert m.precision == 0.0


def test_f1_is_zero_when_precision_and_recall_are_zero() -> None:
    m = compute_metrics(np.array([1, 0]), np.array([0, 1]), np.array([0.2, 0.8]), LABELS)
    assert (m.precision, m.recall, m.f1) == (0.0, 0.0, 0.0)


# --- insufficient data and cross-validation eligibility ------------------------------------


def test_single_class_raises(tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame) -> None:
    y = pd.Series([0] * 6, index=tiny_expression.index)
    with pytest.raises(InsufficientDataError, match="'Group_B' has none"):
        check_class_counts(y, LABELS, min_samples_per_class=4, test_size=0.25)


def test_class_below_minimum_raises_with_count(small_dataset: SyntheticDataset) -> None:
    config = load_config()
    metadata = small_dataset.metadata.copy()
    b_samples = metadata.index[metadata["group"] == "Group_B"][:5]
    keep = ~metadata.index.isin(b_samples)
    expression = small_dataset.expression.loc[metadata.loc[keep, "sample_id"]]
    with pytest.raises(InsufficientDataError, match=r"'Group_B' has 5 samples.*at least 10"):
        run_ml_demo(expression, metadata[keep], config)


def test_split_too_small_for_test_size_raises() -> None:
    y = pd.Series([0] * 5 + [1] * 5, index=[f"S{i}" for i in range(10)])
    with pytest.raises(InsufficientDataError, match="too few for a stratified split"):
        check_class_counts(y, LABELS, min_samples_per_class=4, test_size=0.1)


def test_cv_skipped_when_training_classes_are_small(small_dataset: SyntheticDataset) -> None:
    result = run_ml_demo(small_dataset.expression, small_dataset.metadata, load_config())
    assert not result.cv.ran
    assert "at least 25 per class" in result.cv.reason
    assert result.cv.accuracy_mean is None


def test_cv_runs_with_documented_settings(demo_result: MLDemoResult) -> None:
    cv = demo_result.cv
    assert cv.ran and cv.n_folds == 5 and len(cv.accuracies) == 5
    assert cv.fold_sizes == (12, 12, 12, 12, 12)


def test_cross_validate_training_threshold_is_configurable(demo: DemoData) -> None:
    x, y = prepare_features_and_labels(demo.expression, demo.metadata, LABELS)
    strict = dataclasses.replace(demo.config.ml_demo, min_train_per_class_for_cv=50)
    assert not cross_validate_training(x, y, strict, seed=1).ran


# --- permuted-label baseline ---------------------------------------------------------------


def test_permutation_baseline_is_deterministic(demo: DemoData, demo_result: MLDemoResult) -> None:
    x, y = prepare_features_and_labels(demo.expression, demo.metadata, LABELS)
    train, test = demo_result.train_ids, demo_result.test_ids
    settings: MLDemoConfig = dataclasses.replace(demo.config.ml_demo, n_label_permutations=5)
    args = (x.loc[train], y.loc[train], x.loc[test], y.loc[test], settings)
    first = permutation_baseline(*args, seed=42)
    assert first == permutation_baseline(*args, seed=42)
    assert first.n_permutations == 5
    assert first.accuracies != permutation_baseline(*args, seed=3).accuracies


def test_permutation_baseline_can_be_disabled(demo: DemoData) -> None:
    result = run_ml_demo(
        demo.expression, demo.metadata, _with_ml(demo.config, n_label_permutations=0)
    )
    assert not result.permutation.ran
    assert "Permuted-label baseline not run" in result.comparison_text()


# --- outputs -------------------------------------------------------------------------------


def test_outputs_are_labeled_and_contain_no_model_weights(
    tmp_path: Path, demo_result: MLDemoResult
) -> None:
    files = write_ml_demo_outputs(demo_result, tmp_path / "ml_demo")
    assert {p.name for p in files.values()} == {
        "ml_demo_results.json",
        "confusion_matrix.csv",
        "test_predictions.csv",
        "summary.md",
    }
    payload = json.loads(files["results_json"].read_text(encoding="utf-8"))
    assert payload["data_label"] == DATA_LABEL
    assert payload["classes"]["positive_class_1"] == "Group_B"
    assert payload["classes"]["negative_class_0"] == "Group_A"
    assert "NOT a clinical prediction model" in payload["headline_warning"]
    assert payload["permuted_label_baseline"]["purpose"].startswith("sanity-check")
    text = json.dumps(payload).lower()
    assert "coef" not in text and "intercept" not in text
    for key in ("confusion_matrix", "test_predictions"):
        content = files[key].read_text(encoding="utf-8")
        assert content.startswith(f"# {DATA_LABEL}")
        assert "positive class (1) = Group_B" in content
    summary = files["summary"].read_text(encoding="utf-8")
    assert "NOT a clinical prediction model" in summary
    assert "Positive class: **Group_B**" in summary
    assert "sanity check, not a biological benchmark" in summary


def test_output_statements_follow_configured_positive_class() -> None:
    statements = output_statements("Group_X")
    assert statements[0] == "Positive synthetic class: Group_X (comparison.group_b)"
    assert statements[1:] == REQUIRED_STATEMENTS[1:]


def test_every_output_file_carries_the_required_statements(
    tmp_path: Path, demo_result: MLDemoResult
) -> None:
    files = write_ml_demo_outputs(demo_result, tmp_path)
    for path in files.values():
        text = path.read_text(encoding="utf-8")
        for statement in REQUIRED_STATEMENTS:
            assert statement in text, f"{path.name} is missing: {statement}"


@pytest.mark.e2e
def test_demo_run_beats_permuted_baseline_by_construction(demo_result: MLDemoResult) -> None:
    """On the synthetic demo data the planted groups separate easily, so the correctly
    labeled model should clearly beat the permuted-label baseline. This checks the code,
    not any real-world performance."""
    perm = demo_result.permutation
    assert perm.ran and perm.n_permutations == 20
    assert demo_result.metrics.accuracy > perm.accuracy_mean
    assert demo_result.metrics.roc_auc is not None
