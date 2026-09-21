#!/usr/bin/env bash
# Step 03 — diff the Stage 1 freeze against the staging bucket.
#
# Classifies every freeze artifact as ADDED / MODIFIED / REVERSIONED / UNCHANGED / CONFLICT,
# and every staging file the freeze does not produce as CARRIED / ORPHANED / DUPLICATE. Writes
# scripts/20_upload_bucket/data/diffs/diff_<ts>.tsv and points latest.tsv at it — that
# file is step 04's entire input.
#
# Reads no bucket: the staging side comes from step 01's snapshot. Writes nothing to GCS.
# Fails on CONFLICT (content changed, version did not).
#
# Requires step 01. Reads ${PRECOVID_ROOT}/staging/freeze — set FREEZE_DIR to point at a
# freeze built in another checkout.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export STAGE_LIB="${HERE}/../lib"
# The tests run whether or not the stem does. common.sh sets -e, so a bare invocation
# would skip them on any non-zero exit — and the non-zero exit here is CONFLICT, the
# verdict classify_freeze_file() produces and the tests exist to cover. Gating them on a
# clean run would leave them unrun in exactly the state that needs them. They read the
# freeze and step 01's snapshot, never this run's result, so running them after a failed
# stem is meaningful. The stem's code is propagated afterwards.
rc=0
"${RSCRIPT}" "${HERE}/diff_freeze_vs_staging.R" || rc=$?
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/diff_tests.R"
fi
exit ${rc}
