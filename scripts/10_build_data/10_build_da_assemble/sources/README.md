# Step 10 (build_da_assemble) raw sources

Raw inputs for `../differential_analysis_results.R` that are **not** regenerated inside this
pipeline — external artifacts vendored here so the step runs offline. Same convention as
`../../09_build_da/sources/` and `../../06_generate_qc_norm/sources/`.

Unlike those two, everything here is small and text, so it is committed rather than fetched
by `make sources`.

## `contrast_converter.txt` — for `../differential_analysis_results.R`

The 33 acute-bout contrasts, as a two-column tab-delimited table:

| column | meaning |
|---|---|
| `contrast` | the model-matrix expression, e.g. `group_timepointADUEndur.during_20_min - group_timepointADUEndur.pre_exercise - ...`. This is the join key — it is what step 09 writes into each DA table's `contrast` column. |
| `contrast_short` | the human-readable form, e.g. `Endur.during_20_min - Control.during_20_min (delta-delta)`. |

**Provenance:** copied verbatim from
`MotrpacHumanPreSuspensionAnalysis/data-raw/contrast_converter.txt`
(commit `99a97e2`, 2026-01-30). Not consortium-gated: it names contrasts, not samples.

**This file is authored, not derived** — which is why it is vendored rather than rebuilt
from the contrasts present in the freeze. Two things depend on its exact contents:

- `contrast_short` is hand-written. Nothing in the freeze carries it, and every other
  contrast column the package exposes (`contrast_type`, `contrast_category`,
  `randomGroupCode`, `Timepoint`) is parsed out of it by regex in the builder. A typo here
  silently re-categorises contrasts.
- **Row order is significant.** It becomes `contrast_order`, and it sets the factor level
  order of `contrast` / `contrast_short` / `contrast_category` that every `*_DA` object
  inherits. Re-sorting this file changes the objects.

The builder fails loudly if a DA table carries a contrast this file does not list, so a
freeze refit that adds or renames a contrast will not pass unnoticed.
