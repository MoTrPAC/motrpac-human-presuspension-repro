# Step 12 — build_camera_results

Pre-ranked CAMERA (Correlation Adjusted MEan RAnk) molecular-signature analysis of the
differential-analysis results. For each tissue × ome and each contrast, it tests whether the
z-statistics of the features in a gene / kinase / metabolite set are shifted relative to
everything else, and reports a direction, a two-sample t-statistic, its standard-Normal
equivalent, and a BH-adjusted p-value.

Output: `scripts/10_build_data/data/CAMERA_RESULTS.rda` — one object, ~1.0M rows × 19 cols.
Requires steps 03, 04, 07 and 10.

## What is in this folder

| File | Kind | Notes |
|------|------|-------|
| `build.sh` | driver | runs the builder, then the tests unless `RUN_TESTS=0` |
| `CAMERA_RESULTS.R` | adapted | builds the nested `DA_list`, calls `run_cameraPR()`, saves |
| `run_cameraPR.R` | **vendored** | `run_cameraPR()` + `.prepare_DA_results()`, `.create_index()`, `.prepare_sets()`, `.prepare_DA_results_and_sets()` |
| `camera_tests.R` | tests | scope, joins, statistics, BH grouping, diff-vs-package |

Vendored **code** sits alongside the stem, the same convention as
`09_build_da/da_common.R`. Upstream's
`data-raw/CAMERA_RESULTS.R` is a single line — `CAMERA_RESULTS <-
MotrpacHumanPreSuspensionAnalysis::run_cameraPR()` — so all the substance is in the function,
which is why it is vendored rather than reimplemented.

## Vendored code — provenance

Copied from `MotrpacHumanPreSuspensionAnalysis/R/run_cameraPR.R`. Every change is marked in
place with the upstream original kept commented directly above it, and they are nearly all
the same kind: `MotrpacHumanPreSuspensionAnalysis::<OBJECT>` becomes the pipeline's local
global of the same name — `MOLECULAR_SIGNATURES` (step 03), `SET_TO_ID` (step 04),
`HUMAN_FEATURE_TO_GENE` (step 07), `CONTRAST_CONVERTER` (step 10). `CAMERA_RESULTS.R` loads
those four and the vendored functions close over them.

Two structural changes:

- `check_package_installation(pkg = "TMSig")` is dropped; `CAMERA_RESULTS.R` fails up front if
  TMSig is missing.
- `.prepare_DA_results()` no longer falls back to `load_differential_analysis()` when
  `DA_list` is `NULL`. That fallback reads the objects out of the installed analysis package,
  which is the dependency this stage removes, so a missing `DA_list` is now an error.

Modelling and reshaping are otherwise unchanged.

## Scope

Five omes, which is `run_cameraPR()`'s own default and matches the released object:
`transcript-rna-seq`, `prot-pr`, `prot-ph`, `prot-ol`, `metab` — 11 tissue × ome datasets.

The two clinical DA objects drop out of `selected_omes` on their own (clinical chemistry is
9 analytes, not an enrichment target). Epigen is assembled by step 10, so it does reach the
output dir; it is filtered out at the file glob instead, because letting `selected_omes` drop
it would mean loading 20.4M rows into `DA_list` first only to hand them over and discard them.

Metabolomics is tested as **one dataset per tissue**, not per platform: the step-10 `METAB`
objects stack every platform under `assay = "metab"` with the platform in its own column, and
upstream tests them together.

Databases are every collection in `MOLECULAR_SIGNATURES` except `PTMSIGDB`, again the
upstream default — 12 of them. PTMSigDB sets carry a `;u`/`;d` direction suffix that only the
PTMsigDB branch of `.prepare_sets()` handles, and the released object has no PTMSIGDB rows
either.

## Phosphoproteomics depends on a step-07 fix

`run_cameraPR()` keys prot-ph enrichment on the **flanking sequence**, not the `feature_id`,
because PhosphoSitePlus kinase sets (the `PSP` database) are defined over flanking sequences.

Step 07's `select()` used to drop `flanking_sequence` from `HUMAN_FEATURE_TO_GENE` — a gap
`feature_to_gene_tests.R` recorded as a SKIP with "revisit". Left as it was, this step would
have failed for prot-ph. The column is present in the prot-ph `metadata_features` freeze
files all along, so step 07 now keeps it (`any_of()`, since only prot-ph has it), and both
that test and `camera_tests.R` assert its presence rather than skipping.

## Comparison against the released package

`camera_tests.R` reports a diff against the shipped `CAMERA_RESULTS`. It is **INFO, never
FAIL** — this object is computed from locally refit DA, so it moves with the freeze.

As built against freeze `v2.0`: 1,016,739 rows against the released 1,018,860 (0.2% apart),
same 19 columns in the same order, the same 5 assays and the same 12 databases.
