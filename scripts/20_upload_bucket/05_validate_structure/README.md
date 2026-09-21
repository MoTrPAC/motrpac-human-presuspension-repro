# Step 05 — validate_structure

Checks the staging bucket against the manifest in `lib/required_structure.R` and the
column schemas in `lib/expected_columns.R`.

Six checks, all reported in `logs/upload_validate_report.tsv`:

| Check | Fails the run? |
|-------|----------------|
| `validate:required` — every required manifest row resolves to a file whose header carries the required columns and none of the forbidden ones | yes |
| `validate:optional` — same, for optional rows (`removed-samples`, the `resources/` caches) | no |
| `validate:coverage` — every file in the bucket is claimed by a manifest row or by `CARRIED_FORWARD_RE` | no |
| `validate:qc-report` — each subdir in `QC_REPORT_SUBDIRS` still has an HTML QC report | no |
| `validate:versions` — every versioned file in the bucket carries the version `config/file_versions.json` names | yes |
| `validate:forbidden` — no file present, at any version and on any row, carries a column its schema forbids (case-insensitive). This is what stops a removed-samples file with participant identifiers, which is an optional row with no required columns | yes |

Per-row detail: `logs/structure_validation_<ts>.tsv`.

Writes nothing to GCS. Lists the bucket itself rather than reusing step 01's snapshot,
which step 04 will have invalidated.

## When it skips

**SKIP (77) when step 03 still lists files waiting to be uploaded.** On a dry run the
bucket is untouched, so validating it reports on the *previous* release — every DA table
in the bucket today predates the v2.0 schema, so a dry run would end in 45 failures that
say nothing about the release being built. `FORCE_VALIDATE=1` validates the bucket as it
stands, which is how to check a bucket someone else uploaded to.

## Matching

The manifest describes a file; this step composes its name and matches it as a path
prefix:

```
<gcs_subdir>/<category_dir>/human-precovid-sed-adu_<tissue_code>_<ome>_<data_category>_<data_details>_v
```

`category_dir` is the data category, except that `imputed` files live in `qc-norm/`
beside the matrix they were imputed from. Where a stem appears at more than one version,
the highest wins — step 03 is what reports the duplicate.

Upstream instead greps the whole listing once per manifest field and keeps whatever
survives all five, unanchored. That is loose enough to cross rows — `da` appears inside
`metadata`, `samples` inside `removed-samples`, which is why upstream has to wrap
`data_details` in underscores — and it means a file sitting in the wrong subdir still
matches. Here it does not, and the coverage check reports it.

## Headers, not files

Column checks read the first line by range request (`gsutil cat -r 0-1048575`) rather
than downloading the file. A qc-norm matrix is hundreds of MB and only its header is
under test.

Column names are unquoted except in the methylcap DA tables, which are external MALAX
GLMM output copied into the freeze verbatim and were written with `quote = TRUE`.
`read.csv()` strips quotes for free; splitting a raw line does not, so the step strips
them explicitly. Without that, every methylcap column reads as `"feature_id"` with the
quotes attached and the whole schema looks absent.

`qc-norm` and `imputed` matrices additionally fail on a header of `feature_id` alone —
that is a matrix that lost its samples.

A `.txt.gz` file would come back as gzip bytes, so the schema check is skipped for
compressed files and the row says so. Nothing ships compressed today; teaching the reader
to decompress can wait until something does.

## Which subdirs should have a QC report

`QC_REPORT_SUBDIRS` in `lib/required_structure.R`, not "every subdir in the manifest".
Neither metabolomics subdir has an HTML report — not in the staging bucket, not in
production v1.3, not in the freeze — so checking all five would warn on those two on every
run forever, which trains you to ignore the check. Add a subdir there when a report for it
starts shipping.

## Adapted from upstream

`google_cloud_bucket_checks/05_validate_structure.R`.

- **No value comparison.** Upstream downloads every matched qc-norm and metadata file,
  caches it, and runs `compare_qc_norm()` against the installed
  `MotrpacHumanPreSuspensionData`. Stage 2 does not depend on that package — Stage 3 is
  what updates it, and Stage 1's own tests already compare the freeze against it — so this
  step checks structure and schema and leaves values to the stage that owns them.
  Consequently `qc_norm_visualization_helpers.R` is not ported and
  `SKIP_TISSUE_NOT_IN_PKG` has nothing to configure.
- **Anchored matching** (above).
- **Coverage is checked.** Upstream validates the manifest against the bucket but never
  the bucket against the manifest, so a file nothing expects goes unnoticed.
- **The `clinical-chemistry` special case is gone.** Upstream hardcodes
  `if (ome == "clinical-chemistry") tissue_code = "t02-plasma"` because its manifest
  appends clinical rows without one. The two clinical assays are ordinary omes in
  `OME_TISSUE_CODE` now and carry their own tissue code.
- **Headers by range request**, not whole-file downloads.

## Tests

`validate_tests.R`, run from `build.sh` unless `RUN_TESTS=0`. Asserts that no composed
prefix contains an `NA` (an ome missing from `.qc_details()` or `.gcs_subdir()` would
compose one silently), that no prefix is a prefix of another — which under `startsWith`
matching would let one row claim another's file — that every manifest row resolves to a
column schema, and, when a freeze is present, that every required row matches a freeze
file and every freeze file is claimed by a row.

That last tier runs against the freeze, not the bucket: the only snapshot available to it
is step 01's, taken *before* step 04 uploads. Against the freeze it is a different check
anyway — the manifest against what Stage 1 actually wrote.
