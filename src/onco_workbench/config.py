"""Load and validate the YAML configuration.

Relative paths in the configuration are resolved against a project root. By default
the root is the parent of the directory containing the configuration file (config
files live in ``<root>/config/``), so the same file works on every operating system.

Implemented sections (``seed``, ``paths``, ``synthetic``, ``validation``,
``normalization``, ``comparison``, ``ranking``, ``qc``, ``figures``) are parsed
into typed, validated dataclasses; unknown keys are rejected. Sections for
features not yet built (currently ``ml_demo``) are kept as raw mappings in
:attr:`WorkbenchConfig.sections`.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

PROJECT_ROOT: Path = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH: Path = PROJECT_ROOT / "config" / "default.yaml"

_FORBIDDEN_ID_CHARACTERS = (",", "#", '"', "\n", "\r")


class ConfigError(ValueError):
    """Raised when a configuration file is missing, malformed, or has invalid values."""


def _fail(section: str, message: str) -> None:
    raise ConfigError(f"Invalid configuration [{section}]: {message}")


@dataclass(frozen=True)
class PathsConfig:
    """Absolute, resolved file-system locations used by the workbench."""

    raw_dir: Path
    synthetic_dir: Path
    processed_dir: Path
    outputs_dir: Path
    expression_file: Path
    metadata_file: Path
    ground_truth_file: Path
    manifest_file: Path


@dataclass(frozen=True)
class SyntheticConfig:
    """Parameters of the synthetic demo-data generator.

    All values describe SIMULATED data. Instances validate themselves on creation
    and raise :class:`ConfigError` with an explanatory message if a value is invalid.
    """

    n_samples: int
    n_genes: int
    group_labels: tuple[str, str]
    group_proportions: tuple[float, float]
    n_batches: int
    batch_shift_sd: float
    baseline_mean: float
    baseline_sd: float
    noise_sd: float
    n_signal_genes_up: int
    n_signal_genes_down: int
    signal_effect_range: tuple[float, float]
    missing_fraction: float
    gene_id_prefix: str
    sample_id_prefix: str
    float_decimals: int

    def __post_init__(self) -> None:
        section = "synthetic"
        if self.n_samples < 4:
            _fail(section, f"n_samples must be >= 4, got {self.n_samples}.")
        if self.n_genes < 1:
            _fail(section, f"n_genes must be >= 1, got {self.n_genes}.")
        if len(self.group_labels) != 2 or len(set(self.group_labels)) != 2:
            _fail(section, f"group_labels must be two distinct labels, got {self.group_labels}.")
        for label in self.group_labels:
            _check_identifier(section, "group label", label)
        if len(self.group_proportions) != 2:
            _fail(section, "group_proportions must contain exactly two values.")
        if not all(0.0 < p < 1.0 for p in self.group_proportions):
            _fail(section, f"group_proportions must be in (0, 1), got {self.group_proportions}.")
        if not math.isclose(sum(self.group_proportions), 1.0, abs_tol=1e-9):
            _fail(section, f"group_proportions must sum to 1, got {sum(self.group_proportions)}.")
        n_a, n_b = self.group_sizes
        if min(n_a, n_b) < 2:
            _fail(section, f"each group needs >= 2 samples; sizes would be {n_a} and {n_b}.")
        if not 1 <= self.n_batches <= min(n_a, n_b):
            _fail(
                section,
                f"n_batches must be between 1 and the smaller group size ({min(n_a, n_b)}), "
                f"got {self.n_batches}.",
            )
        for name in ("batch_shift_sd", "baseline_sd"):
            if getattr(self, name) < 0:
                _fail(section, f"{name} must be >= 0, got {getattr(self, name)}.")
        if self.noise_sd <= 0:
            _fail(section, f"noise_sd must be > 0, got {self.noise_sd}.")
        if self.n_signal_genes_up < 0 or self.n_signal_genes_down < 0:
            _fail(section, "n_signal_genes_up and n_signal_genes_down must be >= 0.")
        if self.n_signal_genes_up + self.n_signal_genes_down > self.n_genes:
            _fail(
                section,
                f"n_signal_genes_up + n_signal_genes_down "
                f"({self.n_signal_genes_up + self.n_signal_genes_down}) exceeds n_genes "
                f"({self.n_genes}).",
            )
        low, high = self.signal_effect_range
        if not 0 <= low <= high:
            _fail(section, f"signal_effect_range must satisfy 0 <= low <= high, got {low, high}.")
        if not 0.0 <= self.missing_fraction < 0.5:
            _fail(section, f"missing_fraction must be in [0, 0.5), got {self.missing_fraction}.")
        _check_identifier(section, "gene_id_prefix", self.gene_id_prefix)
        _check_identifier(section, "sample_id_prefix", self.sample_id_prefix)
        if not 1 <= self.float_decimals <= 8:
            _fail(section, f"float_decimals must be between 1 and 8, got {self.float_decimals}.")

    @property
    def group_sizes(self) -> tuple[int, int]:
        """Number of samples assigned to each of the two groups."""
        n_a = int(round(self.n_samples * self.group_proportions[0]))
        return n_a, self.n_samples - n_a


@dataclass(frozen=True)
class ValidationConfig:
    """Thresholds for input-data validation."""

    allowed_groups: tuple[str, ...] | None
    min_samples_per_group: int
    max_missing_fraction_per_gene: float
    max_missing_fraction_per_sample: float

    def __post_init__(self) -> None:
        section = "validation"
        if self.allowed_groups is not None and len(set(self.allowed_groups)) < 2:
            _fail(section, "allowed_groups must list at least two distinct groups, or be null.")
        if self.min_samples_per_group < 2:
            _fail(section, f"min_samples_per_group must be >= 2, got {self.min_samples_per_group}.")
        for name in ("max_missing_fraction_per_gene", "max_missing_fraction_per_sample"):
            value = getattr(self, name)
            if not 0.0 <= value <= 1.0:
                _fail(section, f"{name} must be in [0, 1], got {value}.")


@dataclass(frozen=True)
class NormalizationConfig:
    """Normalization switches (appropriate for synthetic continuous values only)."""

    median_center_samples: bool
    zscore_genes_for_plots: bool


SUPPORTED_FDR_METHODS: tuple[str, ...] = ("fdr_bh",)


@dataclass(frozen=True)
class ComparisonConfig:
    """Two-group comparison settings. Thresholds are display/flagging choices only."""

    group_a: str
    group_b: str
    min_non_missing_per_group: int
    fdr_method: str
    fdr_threshold: float
    effect_size_threshold: float

    def __post_init__(self) -> None:
        section = "comparison"
        _check_identifier(section, "group_a", self.group_a)
        _check_identifier(section, "group_b", self.group_b)
        if self.group_a == self.group_b:
            _fail(section, f"group_a and group_b must differ, both are {self.group_a!r}.")
        if self.min_non_missing_per_group < 2:
            _fail(
                section,
                f"min_non_missing_per_group must be >= 2, got {self.min_non_missing_per_group}.",
            )
        if self.fdr_method not in SUPPORTED_FDR_METHODS:
            _fail(
                section,
                f"fdr_method must be one of {', '.join(SUPPORTED_FDR_METHODS)}, "
                f"got {self.fdr_method!r}.",
            )
        if not 0.0 < self.fdr_threshold <= 1.0:
            _fail(section, f"fdr_threshold must be in (0, 1], got {self.fdr_threshold}.")
        if self.effect_size_threshold < 0:
            _fail(section, f"effect_size_threshold must be >= 0, got {self.effect_size_threshold}.")


@dataclass(frozen=True)
class RankingConfig:
    """Ranking-score settings: ``score = |d| * -log10(max(p_adj, p_floor))``."""

    p_floor: float
    top_n: int

    def __post_init__(self) -> None:
        if not 0.0 < self.p_floor < 1.0:
            _fail("ranking", f"p_floor must be in (0, 1), got {self.p_floor}.")
        if self.top_n < 1:
            _fail("ranking", f"top_n must be >= 1, got {self.top_n}.")


@dataclass(frozen=True)
class QCConfig:
    """Quality-control settings."""

    heatmap_max_samples: int
    pca_components: int

    def __post_init__(self) -> None:
        if self.heatmap_max_samples < 2:
            _fail("qc", f"heatmap_max_samples must be >= 2, got {self.heatmap_max_samples}.")
        if self.pca_components < 2:
            _fail("qc", f"pca_components must be >= 2, got {self.pca_components}.")


SUPPORTED_FIGURE_FORMATS: tuple[str, ...] = ("png", "svg", "pdf")


@dataclass(frozen=True)
class FiguresConfig:
    """Static figure export settings."""

    dpi: int
    format: str

    def __post_init__(self) -> None:
        if not 50 <= self.dpi <= 600:
            _fail("figures", f"dpi must be between 50 and 600, got {self.dpi}.")
        if self.format not in SUPPORTED_FIGURE_FORMATS:
            _fail(
                "figures",
                f"format must be one of {', '.join(SUPPORTED_FIGURE_FORMATS)}, "
                f"got {self.format!r}.",
            )


@dataclass(frozen=True)
class WorkbenchConfig:
    """Complete, validated configuration."""

    seed: int
    root: Path
    source: Path
    paths: PathsConfig
    synthetic: SyntheticConfig
    validation: ValidationConfig
    normalization: NormalizationConfig
    comparison: ComparisonConfig
    ranking: RankingConfig
    qc: QCConfig
    figures: FiguresConfig
    sections: Mapping[str, Any]

    def with_seed(self, seed: int) -> WorkbenchConfig:
        """Return a copy of this configuration with a different random seed.

        Args:
            seed: Non-negative integer seed.

        Returns:
            A new :class:`WorkbenchConfig`.
        """
        return dataclasses.replace(self, seed=_as_int("root", "seed", seed, minimum=0))

    def with_comparison(self, **changes: Any) -> WorkbenchConfig:
        """Return a copy with selected comparison settings replaced (and re-validated).

        Args:
            **changes: Field names of :class:`ComparisonConfig` and their new values.

        Returns:
            A new :class:`WorkbenchConfig`.

        Raises:
            ConfigError: If a new value is invalid.
        """
        return dataclasses.replace(self, comparison=dataclasses.replace(self.comparison, **changes))


def _check_identifier(section: str, name: str, value: str) -> None:
    if not isinstance(value, str) or not value.strip():
        _fail(section, f"{name} must be a non-empty string, got {value!r}.")
    if any(ch in value for ch in _FORBIDDEN_ID_CHARACTERS):
        _fail(section, f"{name} {value!r} must not contain commas, quotes, '#', or newlines.")


def _as_int(section: str, key: str, value: Any, minimum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int | float):
        _fail(section, f"{key} must be an integer, got {value!r}.")
    if isinstance(value, float) and not value.is_integer():
        _fail(section, f"{key} must be an integer, got {value!r}.")
    result = int(value)
    if minimum is not None and result < minimum:
        _fail(section, f"{key} must be >= {minimum}, got {result}.")
    return result


def _as_float(section: str, key: str, value: Any) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        _fail(section, f"{key} must be a number, got {value!r}.")
    result = float(value)
    if not math.isfinite(result):
        _fail(section, f"{key} must be finite, got {value!r}.")
    return result


def _as_pair(section: str, key: str, value: Any) -> tuple[Any, Any]:
    if not isinstance(value, Sequence) or isinstance(value, str) or len(value) != 2:
        _fail(section, f"{key} must be a list of exactly two values, got {value!r}.")
    return value[0], value[1]


def _section(raw: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    value = raw.get(name)
    if not isinstance(value, Mapping):
        raise ConfigError(f"Configuration is missing the required mapping section '{name}'.")
    return value


def _check_keys(section: str, values: Mapping[str, Any], expected: Sequence[str]) -> None:
    missing = [key for key in expected if key not in values]
    unknown = [key for key in values if key not in expected]
    if missing:
        _fail(section, f"missing required key(s): {', '.join(missing)}.")
    if unknown:
        _fail(section, f"unknown key(s): {', '.join(map(str, unknown))} (check for typos).")


def _parse_paths(values: Mapping[str, Any], root: Path) -> PathsConfig:
    names = [f.name for f in dataclasses.fields(PathsConfig)]
    _check_keys("paths", values, names)
    resolved: dict[str, Path] = {}
    for name in names:
        value = values[name]
        if not isinstance(value, str) or not value.strip():
            _fail("paths", f"{name} must be a non-empty path string, got {value!r}.")
        resolved[name] = (root / Path(value)).resolve()
    return PathsConfig(**resolved)


def _parse_synthetic(values: Mapping[str, Any]) -> SyntheticConfig:
    s = "synthetic"
    _check_keys(s, values, [f.name for f in dataclasses.fields(SyntheticConfig)])
    labels = _as_pair(s, "group_labels", values["group_labels"])
    proportions = _as_pair(s, "group_proportions", values["group_proportions"])
    effect = _as_pair(s, "signal_effect_range", values["signal_effect_range"])
    return SyntheticConfig(
        n_samples=_as_int(s, "n_samples", values["n_samples"]),
        n_genes=_as_int(s, "n_genes", values["n_genes"]),
        group_labels=(str(labels[0]), str(labels[1])),
        group_proportions=(
            _as_float(s, "group_proportions", proportions[0]),
            _as_float(s, "group_proportions", proportions[1]),
        ),
        n_batches=_as_int(s, "n_batches", values["n_batches"]),
        batch_shift_sd=_as_float(s, "batch_shift_sd", values["batch_shift_sd"]),
        baseline_mean=_as_float(s, "baseline_mean", values["baseline_mean"]),
        baseline_sd=_as_float(s, "baseline_sd", values["baseline_sd"]),
        noise_sd=_as_float(s, "noise_sd", values["noise_sd"]),
        n_signal_genes_up=_as_int(s, "n_signal_genes_up", values["n_signal_genes_up"]),
        n_signal_genes_down=_as_int(s, "n_signal_genes_down", values["n_signal_genes_down"]),
        signal_effect_range=(
            _as_float(s, "signal_effect_range", effect[0]),
            _as_float(s, "signal_effect_range", effect[1]),
        ),
        missing_fraction=_as_float(s, "missing_fraction", values["missing_fraction"]),
        gene_id_prefix=values["gene_id_prefix"],
        sample_id_prefix=values["sample_id_prefix"],
        float_decimals=_as_int(s, "float_decimals", values["float_decimals"]),
    )


def _parse_validation(values: Mapping[str, Any]) -> ValidationConfig:
    s = "validation"
    _check_keys(s, values, [f.name for f in dataclasses.fields(ValidationConfig)])
    allowed = values["allowed_groups"]
    if allowed is not None:
        if not isinstance(allowed, Sequence) or isinstance(allowed, str):
            _fail(s, f"allowed_groups must be a list of labels or null, got {allowed!r}.")
        allowed = tuple(str(label) for label in allowed)
    return ValidationConfig(
        allowed_groups=allowed,
        min_samples_per_group=_as_int(s, "min_samples_per_group", values["min_samples_per_group"]),
        max_missing_fraction_per_gene=_as_float(
            s, "max_missing_fraction_per_gene", values["max_missing_fraction_per_gene"]
        ),
        max_missing_fraction_per_sample=_as_float(
            s, "max_missing_fraction_per_sample", values["max_missing_fraction_per_sample"]
        ),
    )


def _as_bool(section: str, key: str, value: Any) -> bool:
    if not isinstance(value, bool):
        _fail(section, f"{key} must be true or false, got {value!r}.")
    return value


def _as_str(section: str, key: str, value: Any) -> str:
    if not isinstance(value, str) or not value.strip():
        _fail(section, f"{key} must be a non-empty string, got {value!r}.")
    return value


def _parse_normalization(values: Mapping[str, Any]) -> NormalizationConfig:
    s = "normalization"
    _check_keys(s, values, [f.name for f in dataclasses.fields(NormalizationConfig)])
    return NormalizationConfig(
        median_center_samples=_as_bool(s, "median_center_samples", values["median_center_samples"]),
        zscore_genes_for_plots=_as_bool(
            s, "zscore_genes_for_plots", values["zscore_genes_for_plots"]
        ),
    )


def _parse_comparison(values: Mapping[str, Any]) -> ComparisonConfig:
    s = "comparison"
    _check_keys(s, values, [f.name for f in dataclasses.fields(ComparisonConfig)])
    return ComparisonConfig(
        group_a=_as_str(s, "group_a", values["group_a"]),
        group_b=_as_str(s, "group_b", values["group_b"]),
        min_non_missing_per_group=_as_int(
            s, "min_non_missing_per_group", values["min_non_missing_per_group"]
        ),
        fdr_method=_as_str(s, "fdr_method", values["fdr_method"]),
        fdr_threshold=_as_float(s, "fdr_threshold", values["fdr_threshold"]),
        effect_size_threshold=_as_float(
            s, "effect_size_threshold", values["effect_size_threshold"]
        ),
    )


def _parse_ranking(values: Mapping[str, Any]) -> RankingConfig:
    s = "ranking"
    _check_keys(s, values, [f.name for f in dataclasses.fields(RankingConfig)])
    return RankingConfig(
        p_floor=_as_float(s, "p_floor", values["p_floor"]),
        top_n=_as_int(s, "top_n", values["top_n"]),
    )


def _parse_qc(values: Mapping[str, Any]) -> QCConfig:
    s = "qc"
    _check_keys(s, values, [f.name for f in dataclasses.fields(QCConfig)])
    return QCConfig(
        heatmap_max_samples=_as_int(s, "heatmap_max_samples", values["heatmap_max_samples"]),
        pca_components=_as_int(s, "pca_components", values["pca_components"]),
    )


def _parse_figures(values: Mapping[str, Any]) -> FiguresConfig:
    s = "figures"
    _check_keys(s, values, [f.name for f in dataclasses.fields(FiguresConfig)])
    return FiguresConfig(
        dpi=_as_int(s, "dpi", values["dpi"]),
        format=_as_str(s, "format", values["format"]),
    )


def load_config(
    path: str | Path | None = None, *, root: str | Path | None = None
) -> WorkbenchConfig:
    """Load, validate, and resolve a workbench configuration file.

    Args:
        path: YAML file to load. Defaults to ``config/default.yaml`` in the project.
        root: Directory that relative paths are resolved against. Defaults to the
            parent of the configuration file's directory.

    Returns:
        The validated :class:`WorkbenchConfig`.

    Raises:
        ConfigError: If the file is missing, is not valid YAML, or contains invalid values.
    """
    config_path = Path(path) if path is not None else DEFAULT_CONFIG_PATH
    config_path = config_path.resolve()
    if not config_path.is_file():
        raise ConfigError(f"Configuration file not found: {config_path}")
    try:
        raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ConfigError(f"Configuration file is not valid YAML: {config_path}\n{exc}") from exc
    if not isinstance(raw, Mapping):
        raise ConfigError(f"Configuration file must contain a YAML mapping: {config_path}")

    project_root = Path(root).resolve() if root is not None else config_path.parent.parent
    if "seed" not in raw:
        raise ConfigError("Configuration is missing the required key 'seed'.")

    typed = {
        "seed",
        "paths",
        "synthetic",
        "validation",
        "normalization",
        "comparison",
        "ranking",
        "qc",
        "figures",
    }
    return WorkbenchConfig(
        seed=_as_int("root", "seed", raw["seed"], minimum=0),
        root=project_root,
        source=config_path,
        paths=_parse_paths(_section(raw, "paths"), project_root),
        synthetic=_parse_synthetic(_section(raw, "synthetic")),
        validation=_parse_validation(_section(raw, "validation")),
        normalization=_parse_normalization(_section(raw, "normalization")),
        comparison=_parse_comparison(_section(raw, "comparison")),
        ranking=_parse_ranking(_section(raw, "ranking")),
        qc=_parse_qc(_section(raw, "qc")),
        figures=_parse_figures(_section(raw, "figures")),
        sections={key: value for key, value in raw.items() if key not in typed},
    )
