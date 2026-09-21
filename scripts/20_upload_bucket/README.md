# Stage 2 — upload to bucket (step by step)

Takes the freeze Stage 1 wrote to `staging/freeze/` and brings the staging bucket into
line with it. **The folder sequence is the run order** — `20_upload_bucket.sh` runs each
`NN_*/build.sh` in turn, the same shape `10_build_data.sh` uses for Stage 1.

Each subfolder exposes one `build.sh` that sources `scripts/lib/common.sh` and does one
step. Its exit code is the step status: **0 = PASS**, **77 = SKIP**, anything else = FAIL.

Promotion to production is deliberately not here — see `make promote`.

## Steps

| Step | Does | Writes to GCS |
|------|------|---------------|
| `01_snapshot_buckets` | MD5 + size inventory of the production and staging buckets | no |
| `02_seed_staging` | mirrors production into staging at the top of a release cycle | **yes** (SKIP unless `SEED_STAGING=1`) |
| `03_diff_freeze_vs_staging` | classifies every freeze artifact against staging | no |
| `04_upload_to_staging` | uploads the ADDED/MODIFIED/REPLACE files, verifies, retires superseded versions | **yes** (SKIP unless `APPLY=1`) |
| `05_validate_structure` | checks the bucket against the manifest and column schemas | no |

Only steps 02 and 04 write, and both are opt-in. A bare `make upload` inventories,
diffs, prints an upload plan and stops — the same posture `make promote` takes toward
production.

## Run

```bash
make upload                          # dry run: snapshot, diff, plan (writes nothing to GCS)
APPLY=1 make upload                  # actually upload, then validate
STEPS="01 03" make upload            # only steps whose folder name matches a token
STOP_ON_FAIL=false make upload       # keep going past failures
FREEZE_DIR=/path/to/freeze make upload   # read the freeze from another checkout
```

Per-step logs: `logs/upload_<step>.log`. Summary: `logs/upload_bucket_report.tsv`
(`status` ∈ PASS/SKIP/FAIL). Each step also writes its own report TSV under `logs/`.

When the diff carries `REPLACE` rows, step 04 also writes `logs/upload_replace_README.md`
— the only record that those objects were overwritten in place, since the bucket itself
keeps no trace of it. See that step's README.

### Env knobs

| Variable | Effect |
|----------|--------|
| `APPLY=1` | step 04 uploads. Without it, it writes `logs/upload_plan_<ts>.tsv` and SKIPs. |
| `SEED_STAGING=1` | step 02 rsyncs production into staging. Once per cycle — see that step's README. |
| `ALLOW_CONFLICT=1` | step 03 downgrades CONFLICT to REPLACE (overwrite in place). Read step 03's README first. |
| `FORCE_VALIDATE=1` | step 05 validates the bucket even with uploads still pending. |
| `FREEZE_DIR` | where the freeze lives. Defaults to `${PRECOVID_ROOT}/staging/freeze`. |
| `STEPS`, `STOP_ON_FAIL` | as in Stage 1. |

## Change types

Step 03 assigns exactly one to every file on either side. Step 04 uploads the first
four and touches nothing else.

| | Meaning | Step 04 |
|---|---|---|
| `ADDED` | freeze stem not in staging | upload |
| `MODIFIED` | **different** version, **different** md5 | upload, then remove the superseded object |
| `REVERSIONED` | **different** version, **same** md5 | upload, then remove the superseded object; WARN |
| `REPLACE` | a CONFLICT, under `ALLOW_CONFLICT=1` | overwrite in place |
| `UNCHANGED` | same version, same md5 | — |
| `CONFLICT` | same version, **different** md5 | fails the run |
| `DUPLICATE` | staging holds two versions of one stem | reported only |
| `CARRIED` | in staging, not regenerated, expected | reported only |
| `ORPHANED` | in staging, not regenerated, not expected | reported only |

Nothing is ever deleted except the specific object a verified upload supersedes.

## Shared resources (stage level)

- `lib/bucket_helpers.R` — bucket and freeze inventories, BIC name parsing, the version
  ranking, and the paths. Reads `CURRENT_VERSION`/`NEW_VERSION`/the two bucket URLs out of
  `config/pipeline.env` by parsing the file, not the environment, for the reason
  `10_build_data/lib/file_versions.R` gives about `NEW_VERSION`: those values are pinned
  there so a step run outside the driver cannot be pointed somewhere else.
- `lib/required_structure.R` — the manifest, generated from the preflight
  `OME_TISSUE_CODE` object. Also holds `CARRIED_FORWARD_RE`, the list of bucket files
  Stage 1 does not produce.
- `lib/expected_columns.R` — per-ome column schemas, required and forbidden.
- `check_required_inputs.sh` — the pre-run gate: tooling, bucket read access, the freeze,
  and the manifest's input object.
- `data/snapshots/`, `data/diffs/` — step outputs, the same way Stage 1 steps write into
  `scripts/10_build_data/data/`. Gitignored.

## Relationship to `google_cloud_bucket_checks`

This stage is the port of
`MotrpacHumanPreSuspensionAnalysis/data-raw/google_cloud_bucket_checks/`, whose steps
01–05 map one-to-one onto the folders here. The stub this replaced shelled out to that
directory's `run_validation_pipeline.R`; it no longer does, and nothing here reads that
directory.

Porting it in resolves the TODO the stub carried — upstream's `config.R` pins
`PRODUCTION_BUCKET`, `STAGING_BUCKET`, `CURRENT_VERSION` and `NEW_VERSION` a second time,
so the release a run reads and the release it writes were decided in two files that had
to be kept in agreement by hand. They are now read from `config/pipeline.env` only.

Each step's README lists what it changed from its upstream counterpart. The changes that
alter behaviour rather than plumbing:

- **The diff is against staging, not a production snapshot** (step 03). Upstream seeds
  staging from production and then diffs against the production snapshot, which is
  equivalent only on the first run of a cycle. The staging bucket has since moved ahead of
  production v1.3, so diffing against production would call every already-uploaded file
  MODIFIED and try to retire a production-versioned object that staging does not have.
- **A content change with no version bump fails** (step 03). Upstream skips any file whose
  basename already exists in staging, taking the name as proof of the contents. See that
  step's README for what this catches in the current freeze.
- **The manifest covers metab and clinical** (`lib/required_structure.R`). Upstream's
  `da_details` lookup lists only the five non-metab omes, so 34 of the 48 DA files in the
  bucket had no manifest entry and were never validated.
- **Validation reads headers, not files** (step 05). Upstream downloads every matched file
  to compare values against the installed data package. Stage 2 does not depend on that
  package — Stage 3 is what updates it, and Stage 1's own tests already compare the freeze
  against it — so this checks structure and schema and leaves values to the stage that
  owns them.

## Tests

Steps 03 and 05 ship `diff_tests.R` and `validate_tests.R`, run from their `build.sh`
unless `RUN_TESTS=0`, against `10_build_data/lib/test_helpers.R` — the same harness and
the same PASS/FAIL/SKIP/INFO vocabulary Stage 1's suites use.

They cover the logic that decides what gets uploaded and what gets **deleted**, all of
which is pure and needs no bucket: version ranking (`v1.10` must outrank `v1.9`), join-key
derivation, BIC name parsing, base64→hex md5 conversion, the full
`classify_freeze_file()` truth table, and the composed file names.

Both suites are scoped to what Stage 1 does not already assert. Stage 1 runs first over
the same `OME_TISSUE_CODE` catalog and checks every freeze file's columns, so what is left
to these is the part it has no view of — chiefly the comparison between a freeze file and
the bucket object it would replace, which is step 03's `CONFLICT`. Each suite says at its
own `== manifest ==` and `== column schemas ==` headers what it leaves to Stage 1 and why.

The tiers that want real file names read the freeze or step 01's snapshot, and SKIP
without it, so the suites run on a fresh clone.

## Still to do

- Step 05 has no equivalent of upstream's value comparison. If Stage 2 should assert on
  values as well as structure, the comparison belongs against the freeze it just
  uploaded, not against an installed package.
- `02_seed_staging` opens a new staging bucket by mirroring production, but the bucket
  name itself is still edited by hand in `config/pipeline.env`.
- `docs/dependency_graph/build_depgraph.R` expands steps, stems, tests and `lib/` for
  `scripts/10_build_data` only — every path is hardcoded to it, from when Stage 2 was a
  single stub node. Stage 2 is now step-structured the same way, so `make depgraph` draws
  it as one driver node with no steps. Generalizing that expansion over the Makefile's
  stage list would pick both stages up.
