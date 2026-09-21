#!/usr/bin/env bash
# Step 01 — snapshot the production and staging buckets.
#
# Writes a timestamped MD5 + size inventory of each to
# scripts/20_upload_bucket/data/snapshots/, and points latest_production.tsv /
# latest_staging.tsv at the newest pair. Every later step reads the staging snapshot
# rather than re-listing the bucket.
#
# Read-only against GCS. Requires read access to both buckets.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export STAGE_LIB="${HERE}/../lib"
"${RSCRIPT}" "${HERE}/snapshot_buckets.R"
