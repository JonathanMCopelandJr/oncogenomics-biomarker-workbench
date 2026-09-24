# Methodology

> RESEARCH AND EDUCATION ONLY - NOT FOR CLINICAL USE. Applies to SYNTHETIC data only.

**Status:** Planned design (Phase 1). Each section is finalized, with exact formulas,
when its code is implemented and tested.

## 1. Synthetic data (Phase 2)

Continuous, log2-like values: a per-gene baseline mean, plus a small additive batch
offset, plus Gaussian noise. A configurable set of "signal" genes receives a mean shift
in `Group_B` (half up, half down). The planted signal genes are recorded in a separate
ground-truth file that is used only for testing.

## 2. Validation (Phase 2)

Checks the schema, unique sample and gene identifiers, allowed group values, numeric
expression values, missingness thresholds, minimum group sizes, and agreement between
matrix and metadata sample IDs. It also warns when batch is perfectly confounded with group.

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
