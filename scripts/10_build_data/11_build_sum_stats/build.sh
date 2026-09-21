#!/usr/bin/env bash
# Step 11 — the *_SUM_STATS data objects: the public aggregate layer (per group x timepoint
# Count/Mean/SD), built from the step-08 *_QC objects and filtered to the features present in
# the step-10 assembled *_DA objects. Requires steps 08 and 10 to have run — not the step-09
# freeze, which for metabolomics is in a different feature namespace and matches nothing; see
# all_group_stats.R's header for what that cost.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${RSCRIPT}" "${HERE}/all_group_stats.R"
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/sum_stats_tests.R"
fi
