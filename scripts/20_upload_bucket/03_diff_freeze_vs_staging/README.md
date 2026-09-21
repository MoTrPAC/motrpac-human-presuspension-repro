# Step 03 — diff_freeze_vs_staging

Classifies every artifact in the Stage 1 freeze against what the staging bucket already
holds, and every staging file the freeze does not produce. Writes
`scripts/20_upload_bucket/data/diffs/diff_<ts>.tsv` and points `latest.tsv` at it. That
file is step 04's entire input.

Reads no bucket — the staging side comes from step 01's snapshot. Writes nothing to GCS.
Requires step 01.

Columns: `change_type`, `rel_path`, `local_path`, `gcs_path`, `staging_gcs_path`,
`join_key`, `ome`, `tissue_code`, `data_category`, `data_details`, `old_version`,
`new_version`, `old_md5`, `new_md5`.

## The join key

The path relative to the bucket root with its `_v<major>.<minor>` suffix removed, so
`v1.2` and `v2.0` of one file match each other. Every file, `resources/` included, is
versioned from v2.0; a pre-v2.0 incumbent with no version token keys on its full path.

Upstream keys on the **basename** alone. Keying on the path means a file that moves
between top-level subdirs shows up as an add plus an orphan instead of silently matching
its old location — which is what the clinical split looks like this cycle:
`clinical_chemistry/` is retired and its two assays are filed under `proteomics/`
(`prot-clinical`) and `metabolomics-targeted/` (`metab-t-clinical`). A move that large
should be visible in the diff.

## CONFLICT

Same version string, different bytes. **Fails the run by default.**

Upstream skips any file whose exact basename already exists in staging — "definitively
unchanged", on the reasoning that the version is in the name. But
`config/file_versions.json` is hand-maintained, so a file regenerated with new content and
left at its old version is exactly the mistake it can make, and skipping it means the
bucket keeps the stale bytes under a number that now means something else.

As of the v2.0 freeze this fires on **44 files**. One example:
`transcriptomics/metadata/human-precovid-sed-adu_t06-muscle_transcript-rna-seq_metadata_samples_v1.2.txt`
has 485 samples in the freeze and 484 in the bucket — one sample was added, at the same
version. Under upstream's rule that file is never uploaded and the bucket never gets the
sample.

The fix is one of two things, and the step will not choose for you:

- **Bump them** in `config/file_versions.json` and re-run Stage 1. They then come back as
  `MODIFIED` and the old version is retired properly.
- **Overwrite in place** with `ALLOW_CONFLICT=1`, which reclassifies them as `REPLACE`.
  Only for a staging bucket, and only after looking at what actually differs — the diff
  TSV carries both md5s, and `logs/upload_03_*.log` lists every path.

An object with no md5 (a composite upload) at a matching version is also `CONFLICT`: not
being able to tell is not the same as unchanged.

Nothing in Stage 1 catches this. Its one bucket-facing check, `diff_freeze_vs_bucket()`,
correlates qc-norm matrices per feature against a 0.5 floor and tolerates version drift by
design, so a file whose content changed under an unchanged version passes it.

## REVERSIONED

Different version string, **same** bytes. Uploaded like any other version bump, and
reported `WARN`.

It is the mirror image of `CONFLICT`. Both mean the freeze and `config/file_versions.json`
disagree about whether a file changed this cycle, and both are mistakes a hand-maintained
map can make; they differ in what the mistake costs. `CONFLICT` would publish new bytes
under a version number that already means something else, which is not recoverable by
looking at the bucket afterwards, so it stops the run. `REVERSIONED` only moves identical
content to a new number, so nothing is lost either way and the run continues.

It is still worth seeing, because shipping one is not free: every consumer pinning the old
version has to move for a file whose contents did not change, and the old object is
retired to make it happen. The usual cause is a version bumped in the map for a file that
turned out not to need regenerating, and the usual fix is to revert that bump and re-run
Stage 1. Keeping the bump is a legitimate choice too, which is why this warns rather than
failing.

Before this type existed these files were `MODIFIED` and invisible — a bump with no
content change looked exactly like a bump with one. Step 04 treats the two identically
(upload, verify, retire the superseded object); the split is in the reporting only.

An object with no md5 at a moved version stays `MODIFIED`: not being able to tell the
bytes are identical is not the same as knowing they are, the same guard `UNCHANGED` uses.

## CARRIED vs ORPHANED

Both are files in staging the freeze does not produce. `CARRIED` matches
`CARRIED_FORWARD_RE` in `lib/required_structure.R` — the 7 `*_qc-report_*.html` reports
and `resources/motrpac_human-precovid_1kg_pca` (`.csv` through v1.4, `_v2.0.csv` after), all
produced outside this pipeline.
`ORPHANED` is everything else, and means either a stale artifact from a retired layout or
something new that belongs in the manifest.

Neither is ever deleted. This step's only effect on the bucket is the plan it hands step 04.

## Version drift

Each freeze file's name is stamped by `freeze_version()` at build time, so it agrees with
`config/file_versions.json` by construction. It stops agreeing when the map is edited
without rebuilding — and then the version this step plans to upload is not the version the
map decided on. Reported as a warning, since the fix is a rebuild and that is Stage 1's
business. A file with no version token (only a pre-v2.0 incumbent) is exempt.

## Adapted from upstream

`google_cloud_bucket_checks/03_diff_local_vs_staging.R`.

- Diffs against the **staging** snapshot, not the production one. Upstream's flow makes
  those identical (it seeds staging from production first); ours does not, because staging
  has moved ahead of production v1.3 across several cycles. Diffing against production
  would call every already-uploaded file `MODIFIED` and then try to retire an object at a
  production-versioned path that staging does not have.
- Same-version content differences are `CONFLICT`, not a silent skip.
- The join key is the relative path, not the basename.
- `CARRIED` / `ORPHANED` / `DUPLICATE` are new. Upstream drops staging-only files from the
  output entirely, so a stale or misfiled object in the bucket is invisible.
- `REVERSIONED` is new. Upstream compares basenames only, so a version bump is a
  non-match either way and it cannot tell a bump that carried a content change from one
  that did not.
- The local side is `staging/freeze/`, which mirrors the bucket layout already. Upstream
  points `LOCAL_UPDATED_DIR` at `data/tmp/freeze` in the precovid repo.

## Tests

`diff_tests.R`, run from `build.sh` unless `RUN_TESTS=0`. The classification lives in
`classify_freeze_file()` in `lib/bucket_helpers.R` rather than inline here precisely so it
can be tested without a bucket — it is the decision that tells step 04 which object to
delete. The suite covers its full truth table (including the NA-md5 case and every
`ALLOW_CONFLICT` interaction), the version ranking that decides which staging object is
the incumbent, join-key derivation for versioned and unversioned names, BIC name parsing,
base64→hex md5 conversion, and the manifest's `required` flags.
