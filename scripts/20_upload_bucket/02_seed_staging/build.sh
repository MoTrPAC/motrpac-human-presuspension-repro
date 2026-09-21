#!/usr/bin/env bash
# Step 02 — seed the staging bucket from production.
#
# Mirrors PRODUCTION_BUCKET into STAGING_BUCKET with `gsutil rsync -r`, so the release
# cycle starts from a copy of what is live. This runs ONCE, at the top of a cycle.
#
# SKIP (77) unless SEED_STAGING=1, because re-running it mid-cycle undoes work: rsync
# without -d only adds and overwrites, so every production file that step 04 has already
# superseded and removed from staging comes back, and the bucket ends up carrying both
# the old and the new version of those files. The staging bucket named in
# config/pipeline.env is already seeded; leave this step skipped unless you are opening a
# new one.
#
# Requires write access to STAGING_BUCKET.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"

if [[ "${SEED_STAGING:-0}" != "1" ]]; then
  warn "not seeding: set SEED_STAGING=1 to mirror ${PRODUCTION_BUCKET} -> ${STAGING_BUCKET}"
  exit 77
fi

gcs_can_write "${STAGING_BUCKET}" \
  || die "no write access to ${STAGING_BUCKET} — run preflight / re-auth before seeding"

log "seeding ${STAGING_BUCKET} from ${PRODUCTION_BUCKET} (single-threaded; several minutes)"

# Do not add -m here. Parallel rsync between GCS buckets fails intermittently on large
# transfers — rate limiting or connection contention — and a partial seed is worse than a
# slow one. This is the same note upstream 02_copy_to_staging.R carries.
"${GSUTIL}" rsync -r "${PRODUCTION_BUCKET}/" "${STAGING_BUCKET}/" \
  || die "gsutil rsync failed"

# The seed invalidates step 01's staging snapshot, so require a re-snapshot rather than
# letting step 03 diff against a listing taken before the copy.
rm -f "${PIPELINE_ROOT}/scripts/20_upload_bucket/data/snapshots/latest_staging.tsv"
ok "staging seeded — re-run step 01 to snapshot the seeded bucket"
