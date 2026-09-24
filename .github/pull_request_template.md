## Summary

<!-- What does this change do, and why? Keep changes small and reviewable. -->

## Checklist

- [ ] `python -m ruff format --check src tests app` and `python -m ruff check src tests app` pass
- [ ] `python -m pytest` passes locally
- [ ] New or changed behavior is covered by tests
- [ ] Relevant docs are updated (`README.md`, `docs/`, `CHANGELOG.md`)
- [ ] **No clinical, diagnostic, prognostic, treatment, or biomarker-validation claims** are added
- [ ] **No real patient data, restricted data, secrets, or credentials** are added
- [ ] Outputs and UI still carry the DEMONSTRATION / SYNTHETIC label and the disclaimer
- [ ] If generator output changed, the demo data were regenerated with `obw generate-data`
      and the change is explained here
