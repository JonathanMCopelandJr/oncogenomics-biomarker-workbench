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

## Optional ML demonstration: limitations

`obw ml-demo` and the dashboard's Model Demonstration page train a logistic-regression
classifier on the synthetic data.

> **Synthetic-data educational example only. This is NOT a clinical prediction model.**
> No model here predicts anything about any person. The classifier distinguishes two
> arbitrary synthetic groups (`Group_A`, `Group_B`) that were simulated with deliberately
> planted differences. It does not classify cancer, predict any disease or outcome, or
> validate any biomarker.

- **High scores are expected by construction.** The synthetic groups differ in 40
  planted genes with large simulated effects, so near-perfect separation says only that
  the code works. It says nothing about real data.
- **The test set is tiny (20 samples).** Each misclassified sample changes accuracy by 5
  percentage points. Metrics from sets this small are unstable and must not be quoted as
  performance.
- **The data are simulated:** independent genes, uniform noise, continuous values, balanced
  batches. Real transcriptomics data are far more complex, and a model like this would
  need independent external validation that this project does not and cannot provide.
- **This is not a diagnostic, prognostic, or treatment tool.** The metrics are not
  evidence of clinical performance and must not inform any decision about any person's
  health.
- **The permuted-label baseline** shows the scores a model gets when labels carry no
  information. It is a sanity check, not a biological benchmark. Compare the correctly
  labeled result with it, not with an absolute standard.

The same wording appears on the dashboard page and in every ML output (`summary.md`,
`ml_demo_results.json`). It is defined once, in `onco_workbench.ml.classifier_demo`.

## Feature status

- **Walkthrough notebook:** it uses only the synthetic demo data and opens with the
  required disclaimer. It frames every apparent difference and ranked gene as an
  intentionally implanted synthetic signal, never as a biomarker, diagnostic result,
  clinical prediction, biological finding, or treatment insight.
- **The dashboard is local-only.** It binds to `localhost`, has no authentication or
  upload feature, reads only the committed synthetic demo files, and turns off
  usage-statistics collection. It is not designed for deployment or multi-user use.
- **Dashboard verification is local only.** The dashboard was visually opened and
  reviewed locally in a browser, and the screenshots in `docs/assets/screenshots/` were
  captured from that session on the synthetic demo data. This is a local verification.
  It is not a production deployment, a multi-user or security test, or an external
  usability study, and the screenshots show only simulated data.
- **Related work and the deep-learning roadmap are documentation only.**
  [`related_work.md`](related_work.md) links to the public list
  `hussius/deeplearning-biology` purely as a discovery index. This project does not
  copy, reproduce, vendor, depend on, benchmark against, or claim results from anything
  listed there. The deep-learning roadmap is **not implemented**. Any future module
  would use only this repository's synthetic data, remain educational and nonclinical,
  and would not use patient data, real clinical data, pretrained biological models,
  external model weights, or unverified external datasets.

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
