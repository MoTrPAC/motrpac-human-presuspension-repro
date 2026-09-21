#!/usr/bin/env bash
# Step 04 — upload the ADDED, MODIFIED and REVERSIONED files to the staging bucket.
#
# DRY RUN BY DEFAULT. Without APPLY=1 it writes the plan to logs/upload_plan_<ts>.tsv and
# exits 77 (SKIP), so the driver records that nothing was written. With APPLY=1 it copies
# each file, verifies its md5 by reading the object back, and only then removes the
# version it supersedes.
#
# This is the one step that writes to a bucket other people read, which is why it is
# opt-in — the same posture `make promote` takes toward production.
#
#   APPLY=1 STEPS=04 make upload
#
# Requires step 03 and write access to STAGING_BUCKET.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export STAGE_LIB="${HERE}/../lib"

if [[ "${APPLY:-0}" == "1" ]]; then
  gcs_can_write "${STAGING_BUCKET}" \
    || die "no write access to ${STAGING_BUCKET} — run preflight / re-auth before uploading"
fi

"${RSCRIPT}" "${HERE}/upload_to_staging.R"
