# Step 04 — upload_to_staging

Applies the step-03 diff to the staging bucket: uploads the `ADDED`, `MODIFIED`,
`REVERSIONED` and `REPLACE` rows, verifies each upload by reading the object's md5 back,
and only then removes the version it supersedes.

**Dry run by default.** Without `APPLY=1` it writes the plan to
`logs/upload_plan_<ts>.tsv` and exits 77, so the driver records that nothing was written.

```bash
APPLY=1 STEPS=04 make upload
```

Requires step 03, and write access to `STAGING_BUCKET` (probed before the R script runs).

Output: `logs/upload_log_<ts>.tsv` — one row per file with both md5s, the verification
result, and whether the superseded object was removed. Summary in
`logs/upload_apply_report.tsv`.

`logs/upload_replace_README.md` — written whenever the diff carries `REPLACE` rows, on the
dry run as well as the real one. Every other thing this step does stays legible from the
bucket afterwards: an `ADDED` object is a path that did not exist, and `MODIFIED` /
`REVERSIONED` move to a new version and take the old object with them. A `REPLACE` leaves
no trace — same path, same version, different bytes — so nothing in the bucket
distinguishes an object that was overwritten from one that was never touched, and any
consumer pinned to that version silently reads different content than before. This file is
the only record of it, and carries both md5s per object so you can tell which of the two
you are holding. The dry run writes it too, so the list can be read *before* the run that
discards the old bytes.

## Why it is opt-in

This is the only step in the pipeline that writes to a bucket other people read. Stage 3
is plan-only for the same reason, and `make promote` refuses to touch production at all.
A dry run gives you the complete list of what would change before anything does.

## Order of operations

Per file: **copy → verify → retire**. The superseded object is removed only after its
replacement is on the bucket and its md5 matches the local file, so a failed or truncated
upload leaves the previous version in place. A verification failure stops the run at the
end, after the log is written, with every superseded version still present.

`REPLACE` files overwrite their own path, so there is nothing to retire — the guard is
both the `change_type` check and an explicit refusal to delete the object just written.

`MODIFIED` and `REVERSIONED` are the two types that retire anything, and the script keeps
them in one `SUPERSEDING` vector rather than testing for `MODIFIED` in four places. Both
moved the version, so both land at a new path and leave the old object behind; a
`REVERSIONED` file skipped here would come back as a `DUPLICATE` next cycle. That the
bytes happen to be identical changes nothing about which object has to go.

## Adapted from upstream

`google_cloud_bucket_checks/04_upload_to_staging.R`.

- **Dry run by default.** Upstream uploads as soon as it is sourced;
  `run_validation_pipeline.R` sources it unconditionally as step 4 of 5.
- **The superseded path comes from the diff.** Upstream derives it by substituting
  `PRODUCTION_BUCKET` for `STAGING_BUCKET` in the production path, which only resolves to
  a real object while staging is still a copy of production. Step 03 already knows which
  staging object it matched; step 04 removes that one.
- **`REPLACE` and `REVERSIONED` are handled.** New change types, so upstream has no notion
  of either.
- Both scripts convert gsutil's base64 md5 to hex before comparing; that conversion now
  lives once in `lib/bucket_helpers.R` rather than inline.
