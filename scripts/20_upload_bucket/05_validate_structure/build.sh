#!/usr/bin/env bash
# Step 05 — validate the staging bucket against the manifest.
#
# Checks that every required ome x tissue x data_category file is present, that its header
# carries the required columns and none of the forbidden ones, that no file in the bucket
# is unaccounted for, and that each subdir still has its HTML QC report.
#
# Reads headers by range request, so it does not download the matrices. Writes nothing to
# GCS. Fails if any REQUIRED manifest row fails.
#
# SKIPs (77) when step 03 still lists files waiting to be uploaded — on a dry run the
# bucket holds the previous release, and validating it would report that release's schema
# as this release's failures. FORCE_VALIDATE=1 validates the bucket as it stands.
#
# Lists the bucket itself rather than reusing step 01's snapshot, which step 04 will have
# invalidated. Output: logs/structure_validation_<ts>.tsv.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export STAGE_LIB="${HERE}/../lib"
# The tests run whether or not the validator does. common.sh sets -e, so a bare invocation
# would skip them on the 77 above — which is the DEFAULT posture, since step 04 uploads
# nothing without APPLY=1. That would leave the manifest and the composed names unchecked
# on every dry run. They read the manifest and the freeze, never the bucket or this run's
# result, so they are meaningful after a skipped validation. The code is propagated
# afterwards, so the driver still records 77 as SKIP.
rc=0
"${RSCRIPT}" "${HERE}/validate_structure.R" || rc=$?
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/validate_tests.R"
fi
exit ${rc}
