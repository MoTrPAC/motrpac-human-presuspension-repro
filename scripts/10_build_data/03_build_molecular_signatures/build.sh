#!/usr/bin/env bash
# Step 03 — build MOLECULAR_SIGNATURES from the 7 preflight GMTs.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${RSCRIPT}" "${HERE}/MOLECULAR_SIGNATURES.R"
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/molecular_signatures_tests.R"
fi
