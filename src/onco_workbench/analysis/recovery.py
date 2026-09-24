"""Compare flagged genes with the synthetic ground truth.

This is a **workflow check**: on simulated data the planted differences are known,
so we can count how many were recovered. It verifies that the code behaves as
designed. It is **not** evidence of performance on any real data.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import pandas as pd

from onco_workbench.analysis.ranking import DIRECTION_HIGHER, DIRECTION_LOWER

_PLANTED_TO_OBSERVED = {"up": DIRECTION_HIGHER, "down": DIRECTION_LOWER}


@dataclass(frozen=True)
class RecoverySummary:
    """Counts comparing ``meets_thresholds`` with the planted signal genes."""

    n_genes: int
    n_planted: int
    n_flagged: int
    true_positives: int
    false_positives: int
    false_negatives: int
    direction_matches: int

    @property
    def sensitivity(self) -> float:
        """Share of planted genes that were flagged (NaN if none were planted)."""
        return self.true_positives / self.n_planted if self.n_planted else math.nan

    @property
    def false_discovery_proportion(self) -> float:
        """Share of flagged genes that were not planted (NaN if none were flagged)."""
        return self.false_positives / self.n_flagged if self.n_flagged else math.nan

    def as_rows(self) -> list[tuple[str, str]]:
        """Return ``(label, value)`` pairs for report tables."""

        def pct(value: float) -> str:
            return "n/a" if math.isnan(value) else f"{value:.1%}"

        return [
            ("Simulated signal genes (planted)", str(self.n_planted)),
            ("Genes meeting thresholds", str(self.n_flagged)),
            ("Planted genes recovered", f"{self.true_positives} ({pct(self.sensitivity)})"),
            ("Planted genes missed", str(self.false_negatives)),
            (
                "Flagged genes that were not planted",
                f"{self.false_positives} ({pct(self.false_discovery_proportion)} of flagged)",
            ),
            ("Recovered genes with the planted direction", str(self.direction_matches)),
        ]


def compare_with_ground_truth(results: pd.DataFrame, ground_truth: pd.DataFrame) -> RecoverySummary:
    """Count recovered, missed, and spurious genes against the planted ground truth.

    Args:
        results: Flagged results with ``gene_id``, ``meets_thresholds``, ``direction``.
        ground_truth: Table with ``gene_id``, ``is_signal``, and ``direction``.

    Returns:
        A :class:`RecoverySummary`.

    Raises:
        ValueError: If gene IDs in the two tables do not match.
    """
    truth = ground_truth.set_index("gene_id")
    if set(truth.index) != set(results["gene_id"]):
        raise ValueError("Ground-truth gene IDs do not match the analysed genes.")
    merged = results.set_index("gene_id").join(truth, rsuffix="_planted")
    planted = merged["is_signal"].astype(bool)
    flagged = merged["meets_thresholds"].astype(bool)
    true_positive = planted & flagged
    expected = merged["direction_planted"].map(_PLANTED_TO_OBSERVED)
    direction_ok = true_positive & (merged["direction"] == expected)
    return RecoverySummary(
        n_genes=len(merged),
        n_planted=int(planted.sum()),
        n_flagged=int(flagged.sum()),
        true_positives=int(true_positive.sum()),
        false_positives=int((flagged & ~planted).sum()),
        false_negatives=int((planted & ~flagged).sum()),
        direction_matches=int(direction_ok.sum()),
    )
