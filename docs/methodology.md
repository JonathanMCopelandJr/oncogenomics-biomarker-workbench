# Methodology

> RESEARCH AND EDUCATION ONLY - NOT FOR CLINICAL USE. Applies to SYNTHETIC data only.

**Status:** Sections 1–6 are implemented and tested. Section 6 (the ML demonstration)
is an optional, separate workflow (`obw ml-demo`).

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

## 3. QC and normalization (implemented: `analysis.qc`, `analysis.normalization`)

**QC summaries.** These are sample and gene counts; missing values overall, per gene,
and per sample; per-sample sum, mean, median, and SD of observed values; and samples
per group and batch.

**Normalization for the comparison** (`normalization.median_center_samples: true`).
Each value becomes `x - median(sample) + median(all sample medians)`. This removes
per-sample offsets while keeping the original scale. Missing values are ignored.

**For plots only** (PCA, correlation heatmap, top-gene heatmap):
1. Missing values are replaced with the gene mean.
2. Each gene is z-scored. If `zscore_genes_for_plots: false`, genes are only centered.

These transformed values are **never** used for statistical testing.

**PCA** uses scikit-learn's `PCA` (full SVD) on the prepared values after dropping
zero-variance genes. It shows scores colored by group and, when present, by batch.

**The correlation heatmap** shows Pearson correlation between samples on the prepared
values, for a deterministic, group-balanced subset of at most `qc.heatmap_max_samples`.
Each group contributes its first samples by ID.

These normalization steps suit simulated continuous values only. Real RNA-seq counts
need count-aware methods, which are out of scope.

## 4. Group comparison (implemented: `analysis.differential`, `analysis.multiple_testing`)

Every difference is **group B minus group A** (`comparison.group_b` minus
`comparison.group_a`). For each gene, using only observed values:

| Statistic | Formula |
|---|---|
| `mean_a`, `mean_b` | Group means |
| `mean_diff` | `mean_b - mean_a` |
| `cohens_d` | `mean_diff / pooled_sd`, where `pooled_sd = sqrt(((n_a-1)·var_a + (n_b-1)·var_b) / (n_a+n_b-2))` |
| `t_statistic` | Welch: `mean_diff / sqrt(var_a/n_a + var_b/n_b)` |
| `df` | Welch-Satterthwaite: `(var_a/n_a + var_b/n_b)² / ((var_a/n_a)²/(n_a-1) + (var_b/n_b)²/(n_b-1))` |
| `p_value` | Two-sided, from Student's t distribution with `df` degrees of freedom |
| `p_adj` | Benjamini-Hochberg (statsmodels `fdr_bh`) across all tested genes |

- A gene is tested only if both groups have at least
  `comparison.min_non_missing_per_group` observed values. Untested genes get missing
  statistics and are excluded from the BH count `m`.
- Undefined cases (zero variance in both groups) give missing values, never infinities.
- The test suite checks `t_statistic` and `p_value` against
  `scipy.stats.ttest_ind(equal_var=False)`, with and without missing values. It also
  checks the BH values against hand-computed examples.

## 5. Flagging and ranking (implemented: `analysis.ranking`, `analysis.recovery`)

- `meets_thresholds`: `p_adj <= comparison.fdr_threshold` **and**
  `|cohens_d| >= comparison.effect_size_threshold`. The boundaries are inclusive.
  Both thresholds are user-chosen display settings, not conclusions.
- `direction`: `higher_in_b`, `lower_in_b`, `no_difference`, or `not_tested`.
- `ranking_score = |d| × -log10(max(p_adj, ranking.p_floor))`. This is a documented
  sorting heuristic, not a statistical test. `rank` 1 is the highest score. Ties are
  broken by `gene_id`, and untested genes are placed last.
- **Ground-truth workflow check** (synthetic demo data only). Flagged genes are compared
  with the planted signal genes, counting recovered, missed, and spurious genes, plus
  direction agreement. This verifies that the code behaves as designed. High recovery is
  expected when planted shifts are large relative to the simulated noise. It is not
  evidence of performance on real data.

## 6. Optional ML demonstration (implemented: `ml.classifier_demo`; `obw ml-demo`)

> **Synthetic-data educational example only. NOT a clinical prediction model.** It
> distinguishes two arbitrary synthetic groups and predicts nothing about any person.

This is an educational demonstration of leakage-safe evaluation on synthetic data. All
preprocessing (mean imputation, standardization) is inside a scikit-learn `Pipeline`
fitted on the training split only. Cross-validation, when eligible, refits the whole
pipeline within each training fold and never uses the test set. No feature selection is
performed.

| Step | Detail (defaults from `config/default.yaml`, section `ml_demo`) |
|---|---|
| Classes | Positive class (1) = `comparison.group_b` (`Group_B`). Negative class (0) = `comparison.group_a` (`Group_A`). |
| Features | All 500 **raw** synthetic gene values. No feature selection, and nothing from the group comparison or ranking is used (no double-dipping). |
| Size check | Every class needs at least `min_samples_per_class` (10) samples, and at least 1 test and 2 training samples after the split. Otherwise `InsufficientDataError` is raised with the class name and count. |
| Split | `train_test_split(stratify=y, test_size=0.25, random_state=seed)`, giving 60 training samples (30 per class) and 20 test samples (10 per class) |
| Pipeline | `SimpleImputer(mean)`, then `StandardScaler`, then `LogisticRegression(C=1.0, max_iter=1000)`, fitted on training rows only. Test rows only pass through the already fitted pipeline. |
| Test metrics | Accuracy, and precision, recall, and F1 for the positive class, computed from the confusion counts. A metric whose denominator is zero is reported as **undefined** with a note, never as 0. ROC-AUC is reported only when the test set contains both classes. The confusion matrix has rows for the true group and columns for the predicted group. |
| Cross-validation | Runs only if every **training** class has at least `min_train_per_class_for_cv` (25) samples, which must be at least `cv_folds` (5). It uses `StratifiedKFold(5, shuffle=True, random_state=seed)` on the training split only and reports fold accuracy and ROC-AUC. Otherwise it is skipped and the reason is shown. |
| Permuted-label baseline | `n_label_permutations` (20) refits on training labels shuffled with `numpy.random.default_rng(seed)`, each scored on the unchanged test set. It is a **sanity-check baseline, not a biological benchmark**: it shows the scores a model gets when labels carry no information. |

**Why the ML preprocessing is separate from median centering.** The group comparison
(section 4) median-centers each sample. The ML pipeline instead uses raw values and
learns its own imputation and scaling from the training split. Keeping them separate
means every ML transformation is fitted inside the leakage-safe pipeline, and changing
the exploratory analysis never silently changes the model (or the reverse). The
baseline analysis and its preprocessing are unchanged by the ML demonstration.

**Tests** (`tests/test_ml_demo.py`) check:

- deterministic, stratified, disjoint splitting seeded from the configuration
- that the imputer and scaler statistics equal the training-set statistics
- **invariance**: adding +1000 to every test value (and blanking some) leaves the fitted
  imputer, scaler, model coefficients, and CV scores unchanged
- that the CV folds are built from training rows only
- metric values against hand counts and scikit-learn, plus the undefined-metric cases
- error handling for missing or too-small classes
- a deterministic permutation baseline
