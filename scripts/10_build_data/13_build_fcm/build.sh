#!/usr/bin/env bash
# Step 13 — FCM: fuzzy c-means clustering of the acute DA results, plus enrichment of the
# clusters.
#
# Clusters features by the shape of their z-score trajectory across the exercise_with_controls
# contrasts (FCM_CLUSTERS), then asks which molecular signatures follow each cluster centroid —
# by CAMERA-PR on the membership probabilities (FCM_CAMERA) and by ORA on the hard-assigned
# clusters (FCM_ORA). Reads the step-10 *_DA objects plus MOLECULAR_SIGNATURES (03), SET_TO_ID
# (04), HUMAN_FEATURE_TO_GENE (07) and CONTRAST_CONVERTER (10). run_cmeans(),
# run_cluster_cameraPR() and run_cluster_ORA() are vendored alongside; see
# FCM_clustering_results.R for what is adapted.
#
# The number of clusters is a build parameter, not a prompt: set FCM_K_ADIPOSE / FCM_K_BLOOD /
# FCM_K_MUSCLE (defaults 13 / 12 / 12) after reading the sweep fcm_diagnostics.R writes to
# scripts/10_build_data/13_build_fcm/fcm_diagnostics/. Skip the sweep with FCM_DIAG=0; widen it with
# FCM_DIAG_KMIN / FCM_DIAG_KMAX / FCM_DIAG_REPEATS.
#
# Requires steps 03, 04, 07 and 10. Output:
# scripts/10_build_data/data/{FCM_CLUSTERS,FCM_CAMERA,FCM_ORA}.rda
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${RSCRIPT}" "${HERE}/FCM_clustering_results.R"
if [[ "${FCM_DIAG:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/fcm_diagnostics.R"
fi
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/fcm_tests.R"
fi
