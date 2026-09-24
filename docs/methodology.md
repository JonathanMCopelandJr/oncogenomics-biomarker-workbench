# Methodology

> RESEARCH AND EDUCATION ONLY - NOT FOR CLINICAL USE. Applies to SYNTHETIC data only.

**Status:** Sections 1–2 are implemented and tested (Phase 2). Sections 3–6 are the
planned design and are finalized when their code is implemented.

## 1. Synthetic data (implemented: `onco_workbench.data.synthetic`)

For sample *i* and gene *j*:

```
x[i, j] = baseline[j] + batch_offset[batch(i), j] + noise[i, j] + true_shift[j] * B[i]
```

| Term | Distribution (default parameters) |
|---|---|
| `baseline[j]` | Normal(mean 8.0, sd 1.5), one per gene |
| `batch_offset[b, j]` | Normal(0, 0.3), one per batch per gene |
| `noise[i, j]` | Normal(0, 0.5), independent |
| `B[i]` | 1 if the sample is in `Group_B`, else 0 |
| `true_shift[j]` | 0 for most genes; for 20 genes Uniform(0.8, 2.0), for 20 genes the negative of Uniform(0.8, 2.0) |

After this, each value is independently set to missing with probability 0.005, and
values are rounded to 3 decimals.

- **Groups:** exact sizes (40/40 by default), assigned in random order.
- **Batches:** assigned in a balanced way *within* each group (20/20 per group), so
  batch is not confounded with group.
- **Randomness:** one `numpy.random.default_rng(seed)` generator, used in a fixed order
  (groups, batches, baselines, batch offsets, noise, signal genes, missingness).
- **Deliberate simplifications:** genes are independent (no correlation structure),
  the noise is homoscedastic, and there are no counts, library sizes, or outliers.
  These keep the ground truth unambiguous for testing and make the data unlike real
  RNA-seq data. See `docs/limitations_and_ethics.md`.

## 2. Validation (implemented: `onco_workbench.data.validation`)

`validate_dataset` collects **all** findings into a report instead of stopping at the
first one. **Errors** make analysis impossible or misleading: blank or duplicate IDs,
non-numeric or infinite values, genes or samples with no observed values, missing or
unexpected group labels, expression samples without metadata, and groups smaller than
the configured minimum. **Warnings** need review but do not block analysis: high
missingness, metadata rows without expression data, missing batch labels, and batch
fully confounded with group. The full rule list is in
[`data_dictionary.md`](data_dictionary.md#validation-rules).

## 3. Normalization (Phase 3)

- Optional per-sample median centering before group comparison.
- Per-gene z-scoring **only** for PCA and heatmaps.
- These steps suit simulated continuous values only. Real RNA-seq counts need
  count-aware methods, which are out of scope.

## 4. Group comparison (Phase 3)

For each gene: group means, mean difference (B minus A), Cohen's d (pooled SD), a Welch
two-sample t-test using non-missing values only, and a raw p-value. Benjamini-Hochberg
adjustment is applied across all tested genes.

## 5. Ranking (Phase 3)

`score = |d| x -log10(max(p_adj, p_floor))`. This is a documented heuristic for sorting
results, not a statistical test.

## 6. ML demonstration (Phase 3)

A stratified train/test split, and a `StandardScaler` plus `LogisticRegression` pipeline
fitted on the training data only. Cross-validation runs only when the sample size allows,
and a label-permutation sanity check is included. The demonstration is educational only.
