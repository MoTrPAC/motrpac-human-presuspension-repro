#!/usr/bin/env bash
# Step 08 — package the *_QC data objects (+ blood_transcript_1/2) from the local
# staging/freeze qc-norm + metadata + pheno. Requires steps 02 and 06 to have run: the
# qc-norm and metadata come from 06, and lib/qc_helpers.R loads step 02's pheno.rda on source.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${RSCRIPT}" "${HERE}/qc_norm_results.R"
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/qc_object_tests.R"
fi
