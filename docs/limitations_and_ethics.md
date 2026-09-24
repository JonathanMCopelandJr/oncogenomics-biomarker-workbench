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

## Feature status

- **Not implemented:** the educational machine-learning demonstration and the
  walkthrough notebook. The dashboard's Model Demonstration page shows the required
  warnings and says so explicitly. No model metrics are shown, simulated, or invented.
- **The dashboard is local-only.** It binds to `localhost`, has no authentication or
  upload feature, reads only the committed synthetic demo files, and turns off
  usage-statistics collection. It is not designed for deployment or multi-user use.

## Known technical limitations

- The demo data and analysis are verified byte-for-byte on Python 3.12 with the pinned
  NumPy 2.5.3. Python 3.11 uses NumPy 2.4.6, because newer NumPy releases require Python
  3.12 or later. It is checked only by CI. If a different NumPy version changed the
  random stream, the regeneration test would fail loudly.
- Figure PNG files may differ in bytes across platforms and matplotlib versions (font
  rendering). The CSV tables are the reproducible record.
- The package is meant to be used from a repository checkout (editable install). A
  built wheel does not include `config/`, `app/`, or the demo data.
- The synthetic generator is deliberately simple: independent genes, equal noise,
  continuous values, and large planted effects. Near-perfect recovery is therefore
  expected and is not a performance claim.
