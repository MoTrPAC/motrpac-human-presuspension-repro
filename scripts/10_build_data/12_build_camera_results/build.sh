#!/usr/bin/env bash
# Step 12 — CAMERA_RESULTS: pre-ranked CAMERA molecular-signature analysis of the DA results.
#
# Tests, for each tissue x ome and contrast, whether the z-statistics of the features in a
# gene / kinase / metabolite set are shifted relative to the rest. Reads the step-10 *_DA
# objects plus MOLECULAR_SIGNATURES (03), SET_TO_ID (04) and HUMAN_FEATURE_TO_GENE (07).
# run_cameraPR() is vendored alongside; see CAMERA_RESULTS.R for what is adapted.
#
# Requires steps 03, 04, 07 and 10. Output: scripts/10_build_data/data/CAMERA_RESULTS.rda.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${RSCRIPT}" "${HERE}/CAMERA_RESULTS.R"
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/camera_tests.R"
fi
