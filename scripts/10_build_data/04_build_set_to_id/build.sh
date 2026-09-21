#!/usr/bin/env bash
# Step 04 — build SET_TO_ID from local MOLECULAR_SIGNATURES.rda.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${RSCRIPT}" "${HERE}/SET_TO_ID.R"
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/set_to_id_tests.R"
fi
