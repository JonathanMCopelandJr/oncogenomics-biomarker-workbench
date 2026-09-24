"""Seeded generator for SYNTHETIC transcriptomics-style demonstration data.

The data are simulated, not measured, and do not describe real patients or real
genes. The model for each value is::

    value[sample, gene] = baseline[gene]
                        + batch_offset[batch(sample), gene]
                        + noise[sample, gene]
                        + true_shift[gene] * (sample is in the second group)

followed by random masking of a small fraction of values (set to missing) and
rounding to ``float_decimals``. A configured number of "signal" genes get a
non-zero ``true_shift``. These planted differences are recorded in the ground-truth
table so that tests can check whether a workflow recovers them.

All randomness comes from one ``numpy.random.default_rng(seed)`` generator, used
in a fixed order. The same seed and configuration always produce the same data.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from onco_workbench.config import SyntheticConfig

SAMPLE_ID_COLUMN = "sample_id"
GROUP_COLUMN = "group"
BATCH_COLUMN = "batch"


@dataclass(frozen=True)
class SyntheticDataset:
    """A generated synthetic dataset.

    Attributes:
        expression: Samples as rows (index ``sample_id``), synthetic gene IDs as columns.
        metadata: One row per sample with ``sample_id``, ``group``, and ``batch``.
        ground_truth: One row per gene with ``gene_id``, ``is_signal``, ``direction``,
            and ``true_shift`` (the planted mean shift in the second group).
        seed: Seed that produced the dataset.
    """

    expression: pd.DataFrame
    metadata: pd.DataFrame
    ground_truth: pd.DataFrame
    seed: int


def make_ids(prefix: str, n: int, min_width: int) -> list[str]:
    """Create zero-padded identifiers such as ``SYN_G0001``.

    Args:
        prefix: Identifier prefix.
        n: Number of identifiers.
        min_width: Minimum number of digits.

    Returns:
        A list of ``n`` unique identifiers.
    """
    width = max(min_width, len(str(n)))
    return [f"{prefix}{i:0{width}d}" for i in range(1, n + 1)]


def generate_synthetic_dataset(config: SyntheticConfig, seed: int) -> SyntheticDataset:
    """Generate a reproducible synthetic dataset.

    Args:
        config: Validated generator parameters.
        seed: Non-negative integer seed for ``numpy.random.default_rng``.

    Returns:
        The generated :class:`SyntheticDataset`.

    Raises:
        ValueError: If ``seed`` is negative.
    """
    if seed < 0:
        raise ValueError(f"seed must be a non-negative integer, got {seed}.")
    rng = np.random.default_rng(seed)
    n, g = config.n_samples, config.n_genes
    label_a, label_b = config.group_labels
    n_a, n_b = config.group_sizes

    # 1. Group assignment: exact group sizes, random order.
    groups = rng.permutation(np.array([label_a] * n_a + [label_b] * n_b, dtype=object))

    # 2. Batch assignment, balanced within each group so batch is not confounded with group.
    batch_index = np.empty(n, dtype=int)
    for label in (label_a, label_b):
        members = np.flatnonzero(groups == label)
        batch_index[members] = rng.permutation(np.arange(members.size) % config.n_batches)

    # 3-5. Per-gene baseline, per-batch-per-gene offset, per-value noise.
    baseline = rng.normal(config.baseline_mean, config.baseline_sd, size=g)
    batch_offsets = rng.normal(0.0, config.batch_shift_sd, size=(config.n_batches, g))
    noise = rng.normal(0.0, config.noise_sd, size=(n, g))

    # 6. Planted signal genes: first n_up are shifted up in group B, the rest down.
    n_up, n_down = config.n_signal_genes_up, config.n_signal_genes_down
    signal_index = rng.choice(g, size=n_up + n_down, replace=False)
    low, high = config.signal_effect_range
    magnitudes = rng.uniform(low, high, size=n_up + n_down)
    true_shift = np.zeros(g)
    true_shift[signal_index[:n_up]] = magnitudes[:n_up]
    true_shift[signal_index[n_up:]] = -magnitudes[n_up:]
    true_shift = np.round(true_shift, config.float_decimals)

    in_group_b = (groups == label_b).astype(float)
    values = baseline[None, :] + batch_offsets[batch_index] + noise
    values = values + in_group_b[:, None] * true_shift[None, :]

    # 7. Missing values, then rounding (so in-memory data match what is written to disk).
    values[rng.random((n, g)) < config.missing_fraction] = np.nan
    values = np.round(values, config.float_decimals)

    sample_ids = make_ids(config.sample_id_prefix, n, min_width=3)
    gene_ids = make_ids(config.gene_id_prefix, g, min_width=4)

    expression = pd.DataFrame(
        values, index=pd.Index(sample_ids, name=SAMPLE_ID_COLUMN), columns=gene_ids
    )
    metadata = pd.DataFrame(
        {
            SAMPLE_ID_COLUMN: sample_ids,
            GROUP_COLUMN: groups.astype(str),
            BATCH_COLUMN: [f"Batch_{i + 1}" for i in batch_index],
        }
    )
    direction = np.where(true_shift > 0, "up", np.where(true_shift < 0, "down", "none"))
    ground_truth = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "is_signal": true_shift != 0,
            "direction": direction,
            "true_shift": true_shift,
        }
    )
    return SyntheticDataset(expression, metadata, ground_truth, seed)
