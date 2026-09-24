"""Optional educational ML demonstration on SYNTHETIC data. NOT a clinical prediction model.

A logistic-regression classifier learns to tell the two arbitrary synthetic groups
apart. The point is to show leakage-safe evaluation practice, not to predict anything
about any person:

- **Classes:** the positive class (1) is ``comparison.group_b`` and the negative class
  (0) is ``comparison.group_a``. Every output states this mapping.
- **Features:** all raw synthetic gene values. There is no feature selection, and
  nothing from the group comparison or ranking is used, so the result cannot be
  inflated by selecting on the full data (double-dipping).
- **Split:** stratified train/test split seeded by the configured ``seed``.
- **Pipeline:** ``SimpleImputer(mean) -> StandardScaler -> LogisticRegression``, fitted
  on the training split only. Test data pass through the already fitted pipeline.
- **Cross-validation:** runs on the training split only, and only when every training
  class has at least ``ml_demo.min_train_per_class_for_cv`` samples. Each fold refits the
  whole pipeline.
- **Permuted-label baseline:** refits on shuffled training labels (seeded), to show the
  scores expected when labels carry no information. It is a sanity check, not a
  biological benchmark.

This preprocessing is intentionally separate from the group comparison's per-sample
median centering. The model's own pipeline must learn every transformation from
training data alone, and keeping it separate means changing one workflow never
silently changes the other.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from onco_workbench.config import MLDemoConfig, WorkbenchConfig
from onco_workbench.data.io import labeled_csv_text
from onco_workbench.data.synthetic import GROUP_COLUMN, SAMPLE_ID_COLUMN
from onco_workbench.data.validation import align_metadata
from onco_workbench.disclaimers import DATA_LABEL, FULL_DISCLAIMER, SHORT_DISCLAIMER

HEADLINE_WARNING = (
    "Synthetic-data educational example only. This is NOT a clinical prediction model. "
    "No model here predicts anything about any person. The classifier distinguishes two "
    "arbitrary synthetic groups that were simulated with deliberately planted differences. "
    "It does not classify cancer, predict any disease or outcome, or validate any biomarker."
)

LIMITATIONS: tuple[str, ...] = (
    "High scores are expected by construction. The synthetic groups differ in planted "
    "genes with large simulated effects, so near-perfect separation says only that the "
    "code works. It says nothing about real data.",
    "The test set is tiny. With 20 test samples, each misclassified sample changes accuracy "
    "by 5 percentage points. Metrics from sets this small are unstable and must not be "
    "quoted as performance.",
    "The data are simulated: independent genes, uniform noise, continuous values, balanced "
    "batches. Real transcriptomics data are far more complex, and a model like this would "
    "need independent external validation that this project does not and cannot provide.",
    "This is not a diagnostic, prognostic, or treatment tool. The metrics are not evidence "
    "of clinical performance and must not inform any decision about any person's health.",
    "The permuted-label baseline shows the scores a model gets when labels carry no "
    "information. It is a sanity check, not a biological benchmark. Compare the "
    "correctly labeled result with it, not with an absolute standard.",
)

OUTPUT_SUBDIR = "ml_demo"

DATASET_STATEMENT = "Dataset: fully synthetic transcriptomics-style demonstration data"
PURPOSE_STATEMENT = "Purpose: educational machine-learning workflow demonstration"
USE_STATEMENT = "Not for clinical, diagnostic, prognostic, or treatment use"


def output_statements(positive_label: str) -> tuple[str, ...]:
    """The four statements every ML output must carry.

    Args:
        positive_label: The positive-class group label (``comparison.group_b``).

    Returns:
        Positive class, dataset, purpose, and permitted-use statements, in that order.
    """
    return (
        f"Positive synthetic class: {positive_label} (comparison.group_b)",
        DATASET_STATEMENT,
        PURPOSE_STATEMENT,
        USE_STATEMENT,
    )


class InsufficientDataError(ValueError):
    """Raised when the synthetic data have too few samples per class for the demo."""


@dataclass(frozen=True)
class ClassLabels:
    """Mapping between synthetic group labels and model classes (0/1)."""

    negative: str
    positive: str

    def name(self, code: int) -> str:
        """Return the group label for class code 0 or 1."""
        return self.positive if int(code) == 1 else self.negative

    def describe(self) -> str:
        """Human-readable statement of the class mapping."""
        return (
            f"positive class (1) = {self.positive} (comparison.group_b); "
            f"negative class (0) = {self.negative} (comparison.group_a)"
        )


@dataclass(frozen=True)
class TestMetrics:
    """Held-out test-set metrics. ``None`` means undefined; ``notes`` explain why."""

    n_test: int
    accuracy: float
    precision: float | None
    recall: float | None
    f1: float | None
    roc_auc: float | None
    confusion: pd.DataFrame
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class CVResult:
    """Cross-validation on the training split, or the reason it was skipped."""

    ran: bool
    reason: str
    n_folds: int = 0
    fold_sizes: tuple[int, ...] = ()
    accuracies: tuple[float, ...] = ()
    roc_aucs: tuple[float, ...] = ()

    @property
    def accuracy_mean(self) -> float | None:
        """Mean fold accuracy, or None if CV did not run."""
        return float(np.mean(self.accuracies)) if self.accuracies else None

    @property
    def accuracy_sd(self) -> float | None:
        """Sample SD of fold accuracy, or None if CV did not run."""
        return float(np.std(self.accuracies, ddof=1)) if len(self.accuracies) > 1 else None

    @property
    def roc_auc_mean(self) -> float | None:
        """Mean fold ROC-AUC, or None if CV did not run."""
        return float(np.mean(self.roc_aucs)) if self.roc_aucs else None


@dataclass(frozen=True)
class PermutationBaseline:
    """Scores after refitting on shuffled training labels (seeded)."""

    ran: bool
    reason: str
    accuracies: tuple[float, ...] = ()
    roc_aucs: tuple[float, ...] = ()

    @property
    def n_permutations(self) -> int:
        """Number of permutations run."""
        return len(self.accuracies)

    @property
    def accuracy_mean(self) -> float | None:
        """Mean permuted-label test accuracy, or None if not run."""
        return float(np.mean(self.accuracies)) if self.accuracies else None

    @property
    def accuracy_sd(self) -> float | None:
        """Sample SD of permuted-label test accuracy, or None."""
        return float(np.std(self.accuracies, ddof=1)) if len(self.accuracies) > 1 else None

    @property
    def roc_auc_mean(self) -> float | None:
        """Mean permuted-label test ROC-AUC, or None."""
        return float(np.mean(self.roc_aucs)) if self.roc_aucs else None


@dataclass
class MLDemoResult:
    """Everything produced by :func:`run_ml_demo`.

    The fitted pipeline is kept in memory for inspection and tests. It is never written
    to disk.
    """

    labels: ClassLabels
    seed: int
    settings: MLDemoConfig
    n_samples: int
    n_features: int
    train_ids: list[str]
    test_ids: list[str]
    train_class_counts: dict[str, int]
    test_class_counts: dict[str, int]
    metrics: TestMetrics
    cv: CVResult
    permutation: PermutationBaseline
    predictions: pd.DataFrame
    pipeline: Pipeline = field(repr=False)

    def comparison_text(self) -> str:
        """Plain-language comparison of the real run with the permuted-label baseline."""
        real = f"Correctly labeled run: test accuracy {self.metrics.accuracy:.2f}"
        if self.metrics.roc_auc is not None:
            real += f", ROC-AUC {self.metrics.roc_auc:.2f}"
        perm = self.permutation
        if not perm.ran:
            return f"{real}. Permuted-label baseline not run: {perm.reason}"
        sd = f" (SD {perm.accuracy_sd:.2f})" if perm.accuracy_sd is not None else ""
        return (
            f"{real}. Permuted-label baseline ({perm.n_permutations} runs with shuffled "
            f"training labels): mean test accuracy {perm.accuracy_mean:.2f}{sd}. "
            "A large gap shows the model is using the planted synthetic group structure. "
            "It is a sanity check, not a biological benchmark."
        )

    def to_dict(self) -> dict[str, Any]:
        """JSON-serializable summary (no model weights)."""
        m, cv, perm = self.metrics, self.cv, self.permutation
        return {
            "data_label": DATA_LABEL,
            "disclaimer": SHORT_DISCLAIMER,
            "statements": list(output_statements(self.labels.positive)),
            "headline_warning": HEADLINE_WARNING,
            "limitations": list(LIMITATIONS),
            "classes": {
                "positive_class_1": self.labels.positive,
                "negative_class_0": self.labels.negative,
                "description": self.labels.describe(),
            },
            "setup": {
                "seed": self.seed,
                "n_samples": self.n_samples,
                "n_features": self.n_features,
                "features": "all raw synthetic gene values (no feature selection)",
                "pipeline": "SimpleImputer(mean) -> StandardScaler -> LogisticRegression",
                "fitted_on": "training split only",
                "test_size": self.settings.test_size,
                "logistic_regression_c": self.settings.logistic_regression_c,
                "train_class_counts": self.train_class_counts,
                "test_class_counts": self.test_class_counts,
            },
            "test_metrics": {
                "n_test": m.n_test,
                "accuracy": m.accuracy,
                "precision": m.precision,
                "recall": m.recall,
                "f1": m.f1,
                "roc_auc": m.roc_auc,
                "notes": list(m.notes),
                "confusion_matrix": {
                    "rows_true": list(m.confusion.index),
                    "columns_predicted": list(m.confusion.columns),
                    "counts": m.confusion.to_numpy().tolist(),
                },
            },
            "cross_validation": {
                "ran": cv.ran,
                "reason": cv.reason,
                "n_folds": cv.n_folds,
                "fold_sizes": list(cv.fold_sizes),
                "accuracies": list(cv.accuracies),
                "roc_aucs": list(cv.roc_aucs),
                "accuracy_mean": cv.accuracy_mean,
                "accuracy_sd": cv.accuracy_sd,
                "roc_auc_mean": cv.roc_auc_mean,
            },
            "permuted_label_baseline": {
                "purpose": "sanity-check baseline, not a biological benchmark",
                "ran": perm.ran,
                "reason": perm.reason,
                "n_permutations": perm.n_permutations,
                "accuracies": list(perm.accuracies),
                "roc_aucs": list(perm.roc_aucs),
                "accuracy_mean": perm.accuracy_mean,
                "accuracy_sd": perm.accuracy_sd,
                "roc_auc_mean": perm.roc_auc_mean,
            },
            "comparison": self.comparison_text(),
        }


def prepare_features_and_labels(
    expression: pd.DataFrame, metadata: pd.DataFrame, labels: ClassLabels
) -> tuple[pd.DataFrame, pd.Series]:
    """Select samples from the two groups and encode labels (group B = 1).

    Args:
        expression: Samples as rows, synthetic genes as columns (raw values).
        metadata: Sample metadata with ``sample_id`` and ``group``.
        labels: Class mapping.

    Returns:
        ``(x, y)`` indexed by sample ID, in expression-row order.
    """
    groups = align_metadata(expression, metadata)[GROUP_COLUMN].astype(str)
    keep = groups.isin([labels.negative, labels.positive])
    x = expression.loc[keep.to_numpy()]
    y = (groups[keep] == labels.positive).astype(int).rename("label")
    y.index.name = SAMPLE_ID_COLUMN
    return x, y


def check_class_counts(
    y: pd.Series, labels: ClassLabels, *, min_samples_per_class: int, test_size: float
) -> None:
    """Refuse to run when classes are missing or too small for a stratified split.

    Args:
        y: Encoded labels (0/1).
        labels: Class mapping (for messages).
        min_samples_per_class: Minimum samples required in each class.
        test_size: Planned test fraction.

    Raises:
        InsufficientDataError: With a message naming the class and its count.
    """
    counts = {code: int((y == code).sum()) for code in (0, 1)}
    for code, count in counts.items():
        name = labels.name(code)
        if count == 0:
            raise InsufficientDataError(
                f"The ML demonstration needs samples from both groups, but '{name}' has none."
            )
        if count < min_samples_per_class:
            raise InsufficientDataError(
                f"Group '{name}' has {count} samples; the ML demonstration requires at least "
                f"{min_samples_per_class} per class (ml_demo.min_samples_per_class)."
            )
        n_test = math.floor(count * test_size)
        if n_test < 1 or count - math.ceil(count * test_size) < 2:
            raise InsufficientDataError(
                f"Group '{name}' has {count} samples, too few for a stratified split with "
                f"test_size={test_size} (each class needs >= 1 test and >= 2 training samples)."
            )


def split_samples(y: pd.Series, *, test_size: float, seed: int) -> tuple[list[str], list[str]]:
    """Stratified, seeded train/test split of sample IDs.

    Args:
        y: Encoded labels indexed by sample ID.
        test_size: Test fraction.
        seed: Random seed (from configuration).

    Returns:
        ``(train_ids, test_ids)``, each in the original sample order.
    """
    train, test = train_test_split(
        y.index.to_numpy(), test_size=test_size, stratify=y.to_numpy(), random_state=seed
    )
    order = {sample: i for i, sample in enumerate(y.index)}
    return sorted(train, key=order.__getitem__), sorted(test, key=order.__getitem__)


def build_pipeline(settings: MLDemoConfig, seed: int) -> Pipeline:
    """Imputer -> scaler -> logistic regression. Fit it on training data only.

    Args:
        settings: ML demo settings.
        seed: Random seed for the model.

    Returns:
        An unfitted :class:`~sklearn.pipeline.Pipeline`.
    """
    return Pipeline(
        [
            ("impute", SimpleImputer(strategy="mean", keep_empty_features=True)),
            ("scale", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    C=settings.logistic_regression_c,
                    max_iter=settings.max_iter,
                    random_state=seed,
                ),
            ),
        ]
    )


def compute_metrics(
    y_true: np.ndarray, y_pred: np.ndarray, y_score: np.ndarray, labels: ClassLabels
) -> TestMetrics:
    """Compute test-set metrics from the confusion counts (positive class = 1).

    Undefined metrics are ``None`` (never silently 0), with an explanatory note:
    precision when nothing is predicted positive, recall when no true positives exist,
    F1 when either is undefined, and ROC-AUC when the test set has only one class.

    Args:
        y_true: True labels (0/1).
        y_pred: Predicted labels (0/1).
        y_score: Predicted probability of the positive class.
        labels: Class mapping (for confusion-matrix labels).

    Returns:
        A :class:`TestMetrics`.
    """
    y_true = np.asarray(y_true, dtype=int)
    y_pred = np.asarray(y_pred, dtype=int)
    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())
    notes: list[str] = []

    precision = tp / (tp + fp) if tp + fp else None
    if precision is None:
        notes.append(f"Precision undefined: no test sample was predicted as {labels.positive}.")
    recall = tp / (tp + fn) if tp + fn else None
    if recall is None:
        notes.append(f"Recall undefined: the test set contains no {labels.positive} samples.")
    if precision is None or recall is None:
        f1 = None
        notes.append("F1 undefined because precision or recall is undefined.")
    else:
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    if len(np.unique(y_true)) < 2:
        roc_auc = None
        notes.append("ROC-AUC undefined: the test set contains only one class.")
    else:
        roc_auc = float(roc_auc_score(y_true, y_score))

    confusion = pd.DataFrame(
        [[tn, fp], [fn, tp]],
        index=[f"true: {labels.negative}", f"true: {labels.positive}"],
        columns=[f"predicted: {labels.negative}", f"predicted: {labels.positive}"],
    )
    return TestMetrics(
        n_test=len(y_true),
        accuracy=(tp + tn) / len(y_true),
        precision=precision,
        recall=recall,
        f1=f1,
        roc_auc=roc_auc,
        confusion=confusion,
        notes=tuple(notes),
    )


def cross_validate_training(
    x_train: pd.DataFrame, y_train: pd.Series, settings: MLDemoConfig, seed: int
) -> CVResult:
    """Stratified K-fold CV on the training split only, if class sizes allow.

    Every fold clones and refits the whole pipeline on that fold's training part, so
    imputation and scaling never see the fold's validation rows. The held-out test set
    is not an input to this function.

    Args:
        x_train: Training features.
        y_train: Training labels (0/1).
        settings: ML demo settings.
        seed: Random seed for fold shuffling.

    Returns:
        A :class:`CVResult` (with ``ran=False`` and a reason if skipped).
    """
    counts = [int((y_train == code).sum()) for code in (0, 1)]
    smallest = min(counts)
    required = settings.min_train_per_class_for_cv
    if smallest < required:
        return CVResult(
            ran=False,
            reason=(
                f"Skipped: the smallest training class has {smallest} samples; "
                f"cross-validation requires at least {required} per class "
                "(ml_demo.min_train_per_class_for_cv)."
            ),
        )
    folds = StratifiedKFold(n_splits=settings.cv_folds, shuffle=True, random_state=seed)
    accuracies, roc_aucs, sizes = [], [], []
    template = build_pipeline(settings, seed)
    for fit_idx, val_idx in folds.split(x_train, y_train):
        model = clone(template).fit(x_train.iloc[fit_idx], y_train.iloc[fit_idx])
        y_val = y_train.iloc[val_idx].to_numpy()
        accuracies.append(float((model.predict(x_train.iloc[val_idx]) == y_val).mean()))
        roc_aucs.append(
            float(roc_auc_score(y_val, model.predict_proba(x_train.iloc[val_idx])[:, 1]))
        )
        sizes.append(len(val_idx))
    return CVResult(
        ran=True,
        reason=(
            f"Ran {settings.cv_folds}-fold stratified CV on the training split "
            f"(every training class has >= {required} samples)."
        ),
        n_folds=settings.cv_folds,
        fold_sizes=tuple(sizes),
        accuracies=tuple(accuracies),
        roc_aucs=tuple(roc_aucs),
    )


def permutation_baseline(
    x_train: pd.DataFrame,
    y_train: pd.Series,
    x_test: pd.DataFrame,
    y_test: pd.Series,
    settings: MLDemoConfig,
    seed: int,
) -> PermutationBaseline:
    """Refit on shuffled training labels and score on the unchanged test set.

    Shuffling uses ``numpy.random.default_rng(seed)``, so results are deterministic.

    Args:
        x_train: Training features.
        y_train: True training labels (shuffled here).
        x_test: Test features.
        y_test: True test labels.
        settings: ML demo settings (``n_label_permutations``).
        seed: Random seed.

    Returns:
        A :class:`PermutationBaseline`.
    """
    n = settings.n_label_permutations
    if n == 0:
        return PermutationBaseline(ran=False, reason="Disabled (ml_demo.n_label_permutations = 0).")
    rng = np.random.default_rng(seed)
    template = build_pipeline(settings, seed)
    y_true = y_test.to_numpy()
    accuracies, roc_aucs = [], []
    for _ in range(n):
        shuffled = rng.permutation(y_train.to_numpy())
        model = clone(template).fit(x_train, shuffled)
        accuracies.append(float((model.predict(x_test) == y_true).mean()))
        roc_aucs.append(float(roc_auc_score(y_true, model.predict_proba(x_test)[:, 1])))
    return PermutationBaseline(
        ran=True,
        reason=f"{n} refits with randomly shuffled training labels (seed {seed}).",
        accuracies=tuple(accuracies),
        roc_aucs=tuple(roc_aucs),
    )


def run_ml_demo(
    expression: pd.DataFrame, metadata: pd.DataFrame, config: WorkbenchConfig
) -> MLDemoResult:
    """Run the full educational demonstration on synthetic data.

    Args:
        expression: Validated synthetic expression matrix (raw values).
        metadata: Validated metadata.
        config: Workbench configuration (``seed``, ``comparison``, ``ml_demo``).

    Returns:
        An :class:`MLDemoResult`.

    Raises:
        InsufficientDataError: If either class is missing or too small.
    """
    settings = config.ml_demo
    labels = ClassLabels(negative=config.comparison.group_a, positive=config.comparison.group_b)
    x, y = prepare_features_and_labels(expression, metadata, labels)
    check_class_counts(
        y,
        labels,
        min_samples_per_class=settings.min_samples_per_class,
        test_size=settings.test_size,
    )
    train_ids, test_ids = split_samples(y, test_size=settings.test_size, seed=config.seed)
    x_train, y_train = x.loc[train_ids], y.loc[train_ids]
    x_test, y_test = x.loc[test_ids], y.loc[test_ids]

    pipeline = build_pipeline(settings, config.seed).fit(x_train, y_train)
    y_pred = pipeline.predict(x_test)
    y_score = pipeline.predict_proba(x_test)[:, 1]
    metrics = compute_metrics(y_test.to_numpy(), y_pred, y_score, labels)

    predictions = pd.DataFrame(
        {
            SAMPLE_ID_COLUMN: test_ids,
            "true_group": [labels.name(v) for v in y_test],
            "predicted_group": [labels.name(v) for v in y_pred],
            f"score_{labels.positive}": y_score,
            "correct": y_pred == y_test.to_numpy(),
        }
    )

    def counts(values: pd.Series) -> dict[str, int]:
        return {labels.name(code): int((values == code).sum()) for code in (0, 1)}

    return MLDemoResult(
        labels=labels,
        seed=config.seed,
        settings=settings,
        n_samples=len(y),
        n_features=x.shape[1],
        train_ids=train_ids,
        test_ids=test_ids,
        train_class_counts=counts(y_train),
        test_class_counts=counts(y_test),
        metrics=metrics,
        cv=cross_validate_training(x_train, y_train, settings, config.seed),
        permutation=permutation_baseline(x_train, y_train, x_test, y_test, settings, config.seed),
        predictions=predictions,
        pipeline=pipeline,
    )


def _fmt(value: float | None) -> str:
    return "undefined" if value is None else f"{value:.3f}"


def render_ml_summary(result: MLDemoResult) -> str:
    """Render a Markdown summary with the warnings, class mapping, and results.

    Args:
        result: Demo result.

    Returns:
        Markdown text ending with a newline.
    """
    m, cv = result.metrics, result.cv
    lines = [
        f"# ML demonstration - {DATA_LABEL}",
        "",
        "> **RESEARCH AND EDUCATION ONLY - NOT FOR CLINICAL USE.**",
        f"> **{HEADLINE_WARNING}**",
        "",
        *[f"- **{statement}**" for statement in output_statements(result.labels.positive)],
        "",
        "## Setup",
        "",
        f"- Classes: {result.labels.describe()}.",
        f"- {result.n_samples} synthetic samples, {result.n_features} synthetic genes "
        "(all used, no feature selection).",
        f"- Stratified split (seed {result.seed}): train {result.train_class_counts}, "
        f"test {result.test_class_counts}.",
        "- Pipeline: SimpleImputer(mean) -> StandardScaler -> LogisticRegression "
        f"(C = {result.settings.logistic_regression_c:g}), fitted on the training split only.",
        "",
        "## Held-out test-set metrics",
        "",
        f"Positive class: **{result.labels.positive}**.",
        "",
        "| Metric | Value |",
        "|---|---|",
        f"| Accuracy | {_fmt(m.accuracy)} |",
        f"| Precision | {_fmt(m.precision)} |",
        f"| Recall | {_fmt(m.recall)} |",
        f"| F1 | {_fmt(m.f1)} |",
        f"| ROC-AUC | {_fmt(m.roc_auc)} |",
        "",
        *[f"- {note}" for note in m.notes],
        "",
        "Confusion matrix (rows = true group, columns = predicted group):",
        "",
        "| | " + " | ".join(m.confusion.columns) + " |",
        "|---|" + "---|" * len(m.confusion.columns),
        *[
            f"| {idx} | " + " | ".join(str(v) for v in row) + " |"
            for idx, row in zip(m.confusion.index, m.confusion.to_numpy(), strict=True)
        ],
        "",
        "## Cross-validation (training split only)",
        "",
        cv.reason,
    ]
    if cv.ran:
        lines += [
            "",
            f"Fold accuracies: {', '.join(f'{a:.3f}' for a in cv.accuracies)} "
            f"(mean {_fmt(cv.accuracy_mean)}, SD {_fmt(cv.accuracy_sd)}). "
            f"Mean fold ROC-AUC: {_fmt(cv.roc_auc_mean)}.",
        ]
    lines += [
        "",
        "## Correct labels vs permuted-label baseline",
        "",
        result.comparison_text(),
        "",
        "## Limitations",
        "",
        *[f"- {text}" for text in LIMITATIONS],
        "",
        "---",
        "",
        f"*{FULL_DISCLAIMER}*",
        "",
    ]
    return "\n".join(lines)


def write_ml_demo_outputs(result: MLDemoResult, output_dir: str | Path) -> dict[str, Path]:
    """Write labeled JSON, CSV, and Markdown outputs. No model weights are written.

    Args:
        result: Demo result.
        output_dir: Destination directory (for example ``outputs/ml_demo``).

    Returns:
        Mapping of output role to written path.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    header = [
        "Educational ML demonstration on synthetic data - NOT a clinical prediction model.",
        *output_statements(result.labels.positive),
        f"Classes: {result.labels.describe()}.",
    ]
    files = {
        "results_json": out / "ml_demo_results.json",
        "confusion_matrix": out / "confusion_matrix.csv",
        "test_predictions": out / "test_predictions.csv",
        "summary": out / "summary.md",
    }
    with files["results_json"].open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(result.to_dict(), handle, indent=2, sort_keys=True)
        handle.write("\n")
    confusion = result.metrics.confusion.rename_axis("true \\ predicted").reset_index()
    for key, frame in (("confusion_matrix", confusion), ("test_predictions", result.predictions)):
        text = labeled_csv_text(frame, index=False, extra_comments=header)
        files[key].write_text(text, encoding="utf-8", newline="\n")
    files["summary"].write_text(render_ml_summary(result), encoding="utf-8", newline="\n")
    return files
