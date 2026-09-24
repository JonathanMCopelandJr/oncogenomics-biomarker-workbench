from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from onco_workbench.data.validation import (
    DataValidationError,
    ValidationReport,
    validate_dataset,
)

GROUPS = ("Group_A", "Group_B")


def _validate(
    expression: pd.DataFrame, metadata: pd.DataFrame, **kwargs: object
) -> ValidationReport:
    kwargs.setdefault("allowed_groups", GROUPS)
    return validate_dataset(expression, metadata, **kwargs)


def test_valid_data_passes(tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame) -> None:
    report = _validate(tiny_expression, tiny_metadata)
    assert report.is_valid
    assert report.issues == []
    assert "PASSED" in report.summary()
    report.raise_if_invalid()  # does not raise


def test_empty_expression(tiny_metadata: pd.DataFrame) -> None:
    report = _validate(pd.DataFrame(), tiny_metadata)
    assert "empty_expression" in report.codes()
    assert not report.is_valid


def test_duplicate_sample_ids(tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame) -> None:
    tiny_expression.index = pd.Index(["S1", "S1", "S3", "S4", "S5", "S6"], name="sample_id")
    report = _validate(tiny_expression, tiny_metadata)
    assert "duplicate_sample_id" in report.codes()
    assert "S1" in report.summary()


def test_blank_sample_id(tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame) -> None:
    tiny_expression.index = pd.Index(["", "S2", "S3", "S4", "S5", "S6"], name="sample_id")
    assert "blank_sample_id" in _validate(tiny_expression, tiny_metadata).codes()


def test_duplicate_gene_ids(tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame) -> None:
    tiny_expression.columns = ["SYN_G0001", "SYN_G0001", "SYN_G0003"]
    report = _validate(tiny_expression, tiny_metadata)
    assert "duplicate_gene_id" in report.codes()


def test_blank_gene_id(tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame) -> None:
    tiny_expression.columns = ["SYN_G0001", " ", "SYN_G0003"]
    assert "blank_gene_id" in _validate(tiny_expression, tiny_metadata).codes()


def test_non_numeric_values(tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame) -> None:
    tiny_expression["SYN_G0002"] = ["1.0", "abc", "3", "4", "5", "6"]
    report = _validate(tiny_expression, tiny_metadata)
    assert "non_numeric_values" in report.codes()
    assert "SYN_G0002" in report.summary()


def test_boolean_column_is_not_numeric_expression(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    tiny_expression["SYN_G0002"] = [True, False, True, False, True, False]
    assert "non_numeric_values" in _validate(tiny_expression, tiny_metadata).codes()


def test_infinite_values(tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame) -> None:
    tiny_expression.iloc[0, 0] = np.inf
    assert "non_finite_values" in _validate(tiny_expression, tiny_metadata).codes()


def test_all_missing_gene_is_error(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    tiny_expression["SYN_G0002"] = np.nan
    report = _validate(tiny_expression, tiny_metadata)
    assert "gene_all_missing" in report.codes()
    assert not report.is_valid


def test_high_missingness_gene_is_warning(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    tiny_expression.iloc[0:2, 1] = np.nan  # 2 of 6 = 33% > 20%
    report = _validate(tiny_expression, tiny_metadata)
    assert "gene_high_missingness" in report.codes()
    assert report.is_valid  # warnings do not block


def test_all_missing_sample_is_error(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    tiny_expression.iloc[0, :] = np.nan
    assert "sample_all_missing" in _validate(tiny_expression, tiny_metadata).codes()


def test_high_missingness_sample_is_warning(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    tiny_expression.iloc[0, 0] = np.nan  # 1 of 3 genes = 33%
    report = _validate(tiny_expression, tiny_metadata)
    assert "sample_high_missingness" in report.codes()


def test_metadata_missing_required_column(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    report = _validate(tiny_expression, tiny_metadata.drop(columns="group"))
    assert "metadata_missing_column" in report.codes()
    assert "group" in report.summary()


def test_metadata_duplicate_and_blank_ids(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    tiny_metadata.loc[1, "sample_id"] = "S1"
    tiny_metadata.loc[2, "sample_id"] = np.nan
    codes = _validate(tiny_expression, tiny_metadata).codes()
    assert {"metadata_duplicate_sample_id", "metadata_blank_sample_id"} <= codes


def test_missing_group_label(tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame) -> None:
    tiny_metadata.loc[0, "group"] = np.nan
    assert "missing_group" in _validate(tiny_expression, tiny_metadata).codes()


def test_unexpected_group_label(tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame) -> None:
    tiny_metadata.loc[0, "group"] = "Group_C"
    report = _validate(tiny_expression, tiny_metadata)
    assert "invalid_group" in report.codes()
    assert "Group_C" in report.summary()


def test_any_labels_accepted_when_allowed_groups_is_none(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    tiny_metadata["group"] = ["X"] * 3 + ["Y"] * 3
    assert _validate(tiny_expression, tiny_metadata, allowed_groups=None).is_valid


def test_single_group_rejected_when_allowed_groups_is_none(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    tiny_metadata["group"] = "X"
    report = _validate(tiny_expression, tiny_metadata, allowed_groups=None)
    assert "too_few_groups" in report.codes()


def test_group_too_small(tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame) -> None:
    report = _validate(tiny_expression, tiny_metadata, min_samples_per_group=4)
    assert "group_too_small" in report.codes()
    assert "at least 4 required" in report.summary()


def test_expression_sample_without_metadata_is_error(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    report = _validate(tiny_expression, tiny_metadata.iloc[1:])
    assert "sample_not_in_metadata" in report.codes()
    assert not report.is_valid


def test_extra_metadata_sample_is_warning(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    extra = pd.DataFrame({"sample_id": ["S99"], "group": ["Group_A"], "batch": ["Batch_1"]})
    report = _validate(tiny_expression, pd.concat([tiny_metadata, extra], ignore_index=True))
    assert "sample_not_in_expression" in report.codes()
    assert report.is_valid


def test_missing_batch_is_warning(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    tiny_metadata.loc[0, "batch"] = np.nan
    report = _validate(tiny_expression, tiny_metadata)
    assert "missing_batch" in report.codes()
    assert report.is_valid


def test_batch_confounded_with_group_warns(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    tiny_metadata["batch"] = ["Batch_1"] * 3 + ["Batch_2"] * 3
    assert "batch_confounded_with_group" in _validate(tiny_expression, tiny_metadata).codes()


def test_batch_column_is_optional(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    assert _validate(tiny_expression, tiny_metadata.drop(columns="batch")).is_valid


def test_multiple_problems_are_all_reported(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    tiny_expression.columns = ["SYN_G0001", "SYN_G0001", "SYN_G0003"]
    tiny_metadata.loc[0, "group"] = "Group_C"
    report = _validate(tiny_expression, tiny_metadata)
    assert {"duplicate_gene_id", "invalid_group", "group_too_small"} <= report.codes()


def test_raise_if_invalid_carries_report(
    tiny_expression: pd.DataFrame, tiny_metadata: pd.DataFrame
) -> None:
    tiny_expression.iloc[0, 0] = np.inf
    report = _validate(tiny_expression, tiny_metadata)
    with pytest.raises(DataValidationError) as exc:
        report.raise_if_invalid()
    assert exc.value.report is report
    assert "non_finite_values" in str(exc.value)


def test_long_id_lists_are_truncated(tiny_metadata: pd.DataFrame) -> None:
    expression = pd.DataFrame(
        np.full((6, 12), np.nan),
        index=pd.Index([f"S{i}" for i in range(1, 7)], name="sample_id"),
        columns=[f"G{i:02d}" for i in range(12)],
    )
    expression["G00"] = 1.0  # keep samples partially observed
    summary = _validate(expression, tiny_metadata).summary()
    assert "11 gene(s) have no observed values" in summary
    assert "G01, G02, G03, G04, G05 (+6 more)" in summary
