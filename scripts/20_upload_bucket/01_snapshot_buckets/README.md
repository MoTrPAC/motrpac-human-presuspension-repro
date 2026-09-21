# Step 01 — snapshot_buckets

Timestamped MD5 + size inventory of both buckets, written to
`scripts/20_upload_bucket/data/snapshots/`:

```
snapshot_production_<ts>.tsv   latest_production.tsv
snapshot_staging_<ts>.tsv      latest_staging.tsv
```

Columns: `gcs_path`, `size_bytes`, `md5` (hex), `rel_path`, `snapshot_timestamp`,
`bucket`. Step 03 reads `latest_staging.tsv`; the production snapshot is the baseline
record of what was live when the cycle started and nothing consumes it.

Read-only against GCS. Requires read access to both buckets.

## Adapted from upstream

`google_cloud_bucket_checks/01_snapshot_production.R`.

**Staging is inventoried too.** Upstream snapshots production only, because in its flow
staging is a fresh copy of production and the two are the same thing. They are not the
same thing here: `gs://pre-cawg/staging_20260720` has moved ahead of production v1.3
across several cycles, and it is the bucket step 03 diffs against.

**One listing call, not one per file.** Upstream runs `gsutil stat` once per object to
get its md5 — 227 round trips against staging, ~220 against production. `gsutil ls -L -r`
returns `Content-Length` and `Hash (md5)` for every object in one pass. Same columns out.

**md5 is stored as hex.** gsutil reports base64; `tools::md5sum()` on the local freeze
returns hex. Converting once here means every later comparison is a string equality.

## Objects with no md5

Composite objects (uploaded in parallel chunks) carry only a crc32c. Their `md5` is `NA`,
the step warns with a count, and step 03 cannot content-compare them — a same-version
match with an `NA` md5 is classified `CONFLICT` rather than `UNCHANGED`, because "we
cannot tell" is not "unchanged". Neither bucket currently holds any.
