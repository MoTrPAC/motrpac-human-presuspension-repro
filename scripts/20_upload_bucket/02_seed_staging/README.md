# Step 02 — seed_staging

Mirrors `PRODUCTION_BUCKET` into `STAGING_BUCKET` with `gsutil rsync -r`, so a release
cycle starts from a copy of what is live.

**SKIP (77) unless `SEED_STAGING=1`.** This is a once-per-cycle step, and the staging
bucket named in `config/pipeline.env` is already seeded.

```bash
SEED_STAGING=1 STEPS=02 make upload
```

Requires write access to `STAGING_BUCKET`.

## Why it is opt-in

`rsync` without `-d` only adds and overwrites. Run it mid-cycle and every production file
that step 04 has already superseded and removed comes back, leaving staging with both the
old and the new version of those files — which step 03 then reports as `DUPLICATE` and
step 05 resolves by taking the higher version, so the damage is quiet rather than loud.

Seeding also invalidates step 01's staging snapshot, so this step deletes
`latest_staging.tsv` on success and step 03 will refuse to run until step 01 is re-run.

## Adapted from upstream

`google_cloud_bucket_checks/02_copy_to_staging.R`.

Same `rsync`, same single-threaded insistence — **do not add `-m`**; parallel rsync
between GCS buckets fails intermittently on large transfers, and a partial seed is worse
than a slow one. Upstream carries the same note.

Upstream's post-copy file count check is gone: it compared against the production
snapshot it had just copied from, which re-derives what `rsync`'s own exit code already
says. Step 01 re-run after seeding gives the same assurance with a real inventory.

Ported from R to shell — it was `system(paste(GSUTIL, "rsync", ...))` in a script whose
only other work was that count check.
