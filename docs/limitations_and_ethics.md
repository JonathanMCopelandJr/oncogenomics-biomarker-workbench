# Limitations and ethics

> **RESEARCH AND EDUCATION ONLY - NOT FOR CLINICAL USE.**

## Intended use

This repository is a portfolio and teaching project. It demonstrates software
engineering, data validation, basic statistics, and visualization practices on
**synthetic** data.

## Out of scope (explicitly not supported)

- Diagnosis, prognosis, screening, risk or survival prediction, or treatment selection.
- Claims about drug efficacy.
- Identifying or validating biomarkers.
- Any patient-level output or recommendation.
- Processing protected health information, controlled-access data, or private data.

## Data boundaries

- The project ships only data simulated by its own generator. Synthetic gene IDs
  (`SYN_G...`) and group labels (`Group_A`/`Group_B`) are arbitrary and carry no
  biological or clinical meaning.
- Users must not place patient data in the repository. `data/raw/` is ignored by Git.
- Any future public dataset would need a documented source, license, and explicit
  maintainer approval before it is added.

## Methodological limitations (known in advance)

- Simulated data have a planted, known structure. Recovering it shows the code works as
  designed. It does **not** show that the method works on real data.
- The Welch t-test on continuous values is a teaching simplification. Real RNA-seq
  analysis needs count-aware models, careful normalization, and design-aware handling
  of covariates.
- Benjamini-Hochberg control of the false discovery rate relies on assumptions about the
  dependence between tests. Correlated genes can affect its behavior.
- Small sample sizes make effect-size estimates unstable. Thresholds on volcano plots
  are display choices, not scientific conclusions.
- Batch effects confounded with group cannot be separated statistically. The generator
  balances batches by default, and validation warns about confounding.
- ML metrics on synthetic data are high **by construction** and must not be read as
  predictive performance on any real data.

## Communication standards

All outputs (dashboard, reports, exports) carry the DEMONSTRATION / SYNTHETIC label and
the nonclinical disclaimer. Results are described as "simulated differences detected by
the workflow", never as biomarkers or clinical findings.
