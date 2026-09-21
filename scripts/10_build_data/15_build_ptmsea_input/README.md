# Step 15 — build_ptmsea_input

Formats the phosphoproteomics differential-analysis results as a **PTM-SEA input `.gct`**:
one row per confidently localized phosphosite, one column per selected contrast, values =
the DA **z-statistics**. The other three statistics (`logFC`, `p_value`, `adj_p_value`) and
the feature annotation ride along as row descriptors.

Requires steps 06, 08 and 10 — **06** because `confident_site` is written into the prot-ph
`metadata_features` freeze there, and a freeze built before that column was kept will fail
this step's gate (see *Where `confident_site` comes from* below).

Output, in `scripts/10_build_data/data/`:

| File | What it is |
|---|---|
| `PTMSEA_INPUT_muscle.gct` | **the deliverable** — 16,469 sites × 6 contrasts |
| `PTMSEA_INPUT_adipose.gct` | **the deliverable** — 7,287 sites × 2 contrasts |
| `PTMSEA_INPUT.rda` | the same two GCT objects as a named list, for the tests |

Adipose has two contrasts, not six, because it was only sampled at `post_3.5_4_hr`; muscle
carries `post_15_30_45_min`, `post_3.5_4_hr` and `post_24_hr` × endurance/resistance.

## ⚠ The enrichment itself is NOT in this repo

**The `.gct` is an input file, not a result.** PTM-SEA is run from the **PTM-SEA / ssGSEA2.0
Docker image maintained by the Broad Institute**, scored against **PTMsigDB**. **That run is
out of scope for `motrpac-human-presuspension-repro`**, and this step's contract ends the moment a valid `.gct`
is on disk.

Concretely, this repo:

- does **not** pull, build, or run the Broad container, and has no `docker/` service for it
  (`make docker-*` is this pipeline's own R environment, nothing to do with PTM-SEA);
- does **not** read PTM-SEA enrichment scores anywhere in Stage 1, 2 or 3 — no other step
  consumes `PTMSEA_INPUT*`, and nothing is uploaded to the bucket (Stage 2 syncs
  `staging/freeze` only, and these files are in `scripts/10_build_data/data/`);
- makes **no reproducibility claim** over the container version, the PTMsigDB release, or the
  parameters it is run with. Those are pinned outside this pipeline, so a PTM-SEA result is
  an **external artifact** — reproducing it is not something `make data` can promise.

The handoff is manual and deliberate: build the `.gct` here, hand it to the container, and
record the container version / PTMsigDB release / parameters wherever that result is
published. If PTM-SEA is ever brought in-repo it needs its own stage with its own pinning,
not an extension of this step.

## What is in this folder

| File | Kind | Notes |
|------|------|-------|
| `build.sh` | driver | runs the builder, then the tests unless `RUN_TESTS=0` |
| `PTMSEA_INPUT.R` | adapted | loads the DA + QC globals, attaches `confident_site`, loops tissues, writes the `.gct` |
| `preprocess_PTMSEA.R` | **vendored** | `preprocess_PTMSEA()` |
| `ptmsea_tests.R` | tests | scope, confident-site filter, GCT structure, file round-trip |

Vendored **code** sits alongside the stem, the same convention as steps 12 and 13. Upstream
has no `data-raw` object generator for this at all — `preprocess_PTMSEA()` is a package
*function* whose caller writes the file — so `PTMSEA_INPUT.R` is the caller this pipeline
did not previously have, and the function is vendored rather than reimplemented.

## Vendored code — provenance

Copied from `MotrpacHumanPreSuspensionAnalysis/data-raw/preprocess_PTMSEA.R`. (The copy at
`MotrpacHumanPreSuspension/R/preprocess_PTMSEA.R` is the same file with the
`MotrpacHumanPreSuspensionAnalysis::` prefixes already dropped and a `TODO` about
`id`/`feature_id`; both are accounted for below.) Every change is marked in place with the
upstream original kept commented directly above it:

1. `check_package_installation(pkg = "cmapR")` dropped — `PTMSEA_INPUT.R` fails up front if
   cmapR is missing, exactly as `CAMERA_RESULTS.R` does for TMSig.
2. `MotrpacHumanPreSuspensionAnalysis::{MUSCLE,ADIPOSE}_PROT_PH_DA` become the pipeline's
   local globals of the same name (step 10), loaded by the stem. Same for the `*_PROT_PH_QC`
   objects (step 08). This is the dependency the stage removes.
3. The `id` → `feature_id` rename is now conditional. Upstream carries
   `TODO: @Chris -> Remove this chunk once id/feature_id is changed`; the locally built QC
   objects already key `feature_metadata` on `feature_id`, so applying it unconditionally
   errors. Guarded rather than deleted, so the file still works against an unchanged package
   object.
4. `GCT()` → `cmapR::GCT()`, and `pivot_wider`/`right_join`/`starts_with` namespaced —
   upstream relies on the Analysis package's `@importFrom`, and nothing attaches those here.

Modelling and reshaping are otherwise unchanged: the same `pivot_wider`, the same
`right_join`, the same confident-site filter, the same z-score matrix with tissue-prefixed
column names.

## Where `confident_site` comes from

`preprocess_PTMSEA()` keeps only fully localized phosphosites —
`prot_da_wide[prot_da_wide$confident_site, ]` — and reads that column off
`*_PROT_PH_QC$feature_metadata`.

This step **reads** that column; it does not reconstruct it. That was not true until step 06
was fixed. `.annotate_prot_ph()` used to `select()` down to
`assay / feature_id / entrez_gene / gene_symbol / ensembl_gene / uniprot / flanking_sequence`
before writing `metadata_features`, so `confident_site` was dropped on the way into the
freeze — even though it was sitting on the raw `ratio-results` `rdesc` the whole time and
needed no recomputation.

**This was not a local-regeneration regression.** The shipped
`MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_QC$feature_metadata` has the same 7 columns and
no `confident_site` either, which means upstream's own `preprocess_PTMSEA()` does not run
against the shipped object: the filter subscripts with `NULL` and silently returns **zero
rows**. Keeping the column fixes that too.

The column is now carried at both per-tissue tiers:

| Tier | Object | Scope |
|---|---|---|
| freeze | `*_prot-ph_metadata_features_v2.0.txt` (step 06) | per tissue, exact |
| `.rda` | `*_PROT_PH_QC$feature_metadata` (step 08) | per tissue, exact — **what this step filters on** |

Counts: muscle 16,786 of 18,548 confident; adipose 15,147 of 21,022.

Regenerating the prot-ph freeze for this changed **only** the two `metadata_features` files —
the `qc-norm` matrices came back byte-identical, so the DA in steps 09/10 is untouched.

### `HUMAN_FEATURE_TO_GENE` does not carry it

`confident_site` is a **per-tissue measurement**, but `HUMAN_FEATURE_TO_GENE` is keyed on
`(assay, feature_id)` and has **no tissue column**. 7,865 phosphosites appear in both prot-ph
freeze files and **859 of them disagree** (782 confident in muscle only, 77 in adipose only).
Carried raw, those 859 survive de-duplication as two rows sharing one key, which quietly turns
the map into a one-to-many join for every consumer — step 12 joins DA results against it, so
those features would be double-counted. Collapsing them instead (a site confident in the map
only where it is confident in every tissue that measured it) removes the duplicate keys but
publishes a value that is not the measurement for either tissue.

Step 07 therefore drops the column (map `.txt` `v2.2`). The exact per-tissue value stays on
`*_PROT_PH_QC$feature_metadata`, which is what this step filters on, so the `.gct` is
unaffected — verified byte-identical output. Anything needing site localization must read the
QC object, not the map.

`flanking_sequence` needs none of this and is carried: it disagrees on **0** of the 7,865
shared sites, since the sequence around a residue is a property of the protein rather than the
tissue. That is why it could be added to the map without touching the key.

### The gate

`PTMSEA_INPUT.R` refuses to build unless `confident_site` is present, logical, and free of
`NA`. Both failure modes of `prot_da_wide[prot_da_wide$confident_site, ]` are silent:

- **column absent** → `NULL` subscript → every row dropped → an empty GCT;
- **value `NA`** → an **all-`NA` row** that survives into the matrix instead of being removed,
  which PTM-SEA reads as real data.

`ptmsea_tests.R` guards the same thing from the other side, re-deriving the expected row set
from the raw `ratio-results` file independently of the QC object the builder reads.

## Build parameters

Parameters, not prompts:

```bash
STEPS="15" make data                                      # defaults: both tissues, EE-CON + RE-CON
PTMSEA_TISSUES="muscle" STEPS="15" make data              # one tissue
PTMSEA_CONTRAST_CATEGORY="EE-CON" STEPS="15" make data    # one contrast category
RUN_TESTS=0 STEPS="15" make data                          # skip the tests
```

`PTMSEA_CONTRAST_TYPE` (default `exercise_with_controls`) is the third knob. Upstream's own
default is muscle only; this stem builds both tissues, since the pipeline regenerates both.

## Tests

`ptmsea_tests.R` re-derives the expected row and column sets from the inputs rather than from
the builder, so a filtering or reshaping mistake fails instead of being mirrored:

- every selected contrast is a column, tissue-prefixed, and nothing else is;
- every confident DA site is a row, and **no non-confident site is kept** (the regression
  guard for the all-`NA`-row failure above) — `confident_site` re-read from the raw file by
  glob, not through the catalog lookup the builder uses;
- the GCT is structurally sound: numeric matrix, unique rids, `rid`/`cid` aligned to the
  matrix, `rdesc` carrying `feature_id` and **no** `z.std` columns (values must not also be
  descriptors);
- the written file re-parses, and round-trips to within `5e-5` — `cmapR::write_gct` rounds to
  4 decimal places, so that bound is asserted against the actual write precision rather than
  a loose tolerance;
- **NA stays NA** through the round trip: a contrast a site was not tested in must not become
  a `0`, which PTM-SEA would score as a real null result;
- the filename is not dimension-stamped (`appenddim = FALSE`), so a handoff script can name
  the file.

There is no diff-vs-package tier here, by design rather than by gap: upstream ships a
function, not a data object, so there is nothing to diff against.
