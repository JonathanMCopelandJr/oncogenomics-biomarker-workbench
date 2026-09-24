# data/synthetic/

> **DEMONSTRATION / SYNTHETIC DATA.** Every value in this directory is simulated by
> this repository's generator. It is **not** patient data and does not describe real
> genes, real samples, or real biology.

**Status:** Placeholder. The generator and the demo files are added in Phase 2.

Planned files (all regenerable from the seed in `config/default.yaml`):

| File | Contents |
|---|---|
| `demo_expression.csv` | Samples as rows, synthetic gene IDs (`SYN_G0001`, ...) as columns, continuous log2-like values |
| `demo_metadata.csv` | `sample_id`, synthetic `group` (`Group_A` / `Group_B`), synthetic `batch` |
| `demo_ground_truth.csv` | Which synthetic genes were simulated with a group difference, used only to test that the workflow recovers them |

Field definitions will be documented in [`docs/data_dictionary.md`](../../docs/data_dictionary.md).

License of the generated demo data: CC0-1.0 (public-domain dedication), pending
confirmation in Phase 2.
