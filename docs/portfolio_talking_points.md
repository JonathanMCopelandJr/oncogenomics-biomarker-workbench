# Portfolio talking points

> Research and education project on **synthetic** data. It is not for clinical use.
> Every number below was produced by code in this repository during the Phase 5 audit.
> None of it is a claim about real biology or clinical performance.

## One-sentence summary

A reproducible Python workbench that validates, quality-checks, and compares two groups
in synthetic transcriptomics-style data. It provides a transparent statistics pipeline,
labeled reports, and a local Streamlit dashboard, with the nonclinical boundary built
into the code.

## Talking points

### 1. Testing an analysis with a planted answer

- The generator simulates 80 samples x 500 genes and plants a known shift in 40 genes
  (20 up, 20 down). The ground truth is stored separately.
- The pipeline recovers all 40 planted genes with 0 spurious flags at the default
  thresholds. **This is expected by construction**, because the planted shifts are large
  relative to the noise. I use it to show the code works, not to make a performance claim.
- A test regenerates the committed data from the seed and checks it byte for byte
  (SHA-256).

### 2. Transparent, verified statistics

- Welch's t-test is written out in vectorized NumPy (Welch-Satterthwaite df), along with
  Cohen's d (pooled SD) and NaN-aware Benjamini-Hochberg correction.
- Tests cross-check the p-values against `scipy.stats.ttest_ind(equal_var=False)`, with
  and without missing values, and check BH against hand-worked examples.
- The ranking score is documented as a sorting heuristic, not a test.

### 3. Designing the ethical boundary into the software

- One module holds all disclaimer text. It is reused by the CLI, the CSV headers, the
  figures, the Markdown report, the run manifest, and every dashboard page.
- Gene IDs are neutral (`SYN_G0001`) so no output can suggest a real gene. The group
  labels are arbitrary.
- Tests check that every dashboard page shows the banner and footer. The Model
  Demonstration page declares that it is not yet implemented, instead of showing
  invented metrics.

### 4. Data validation that helps the user

- 23 issue codes (18 errors, 5 warnings), all collected into one report, with offending
  IDs listed.
- Duplicate gene columns are preserved when loading (pandas would silently rename
  them) so that validation can report them.
- A test keeps the data dictionary in sync with the validation codes.

### 5. Architecture for testability

- The package is in a src layout. The CLI, the Streamlit pages, and the dashboard helpers
  are thin layers over the same pipeline functions.
- The dashboard's data layer has no Streamlit dependency and is unit tested. The pages
  are smoke tested with Streamlit's `AppTest`, including their interactions.
- Figures use matplotlib's object-oriented API (no global `pyplot` state) and a
  colorblind-safety-validated palette, with marker shapes as a second cue.

### 6. Reproducibility and engineering hygiene

- Seeded generation, pinned direct dependencies, `pathlib` everywhere, LF line-ending
  normalization, and run manifests with versions, git state, parameters, and checksums.
- The Phase 5 audit reinstalled the project from scratch in a new environment. The full
  suite passed there, and the regenerated data matched the committed checksums.
- CI runs lint, tests, validation, and an end-to-end run on Linux, Windows, and macOS,
  with read-only permissions and no secrets.

## Honest limitations to mention

- The statistics are a teaching simplification. They are not appropriate for real
  RNA-seq counts.
- The ML demonstration is planned but not built.
- Python 3.11 is verified by CI only. The local audit used Python 3.12.
