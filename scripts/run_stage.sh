#!/usr/bin/env bash
# One stage of the DAG, by name. This is the single definition of "what does stage X
# run and which stamp does it set" — the Makefile calls it, and so does the SLURM
# chain submitter, so a batch run and a `make` run execute the same thing.
#
#   bash scripts/run_stage.sh preflight
#   bash scripts/run_stage.sh data
#
# It does not resolve dependencies: that is the Makefile's job locally, and
# --dependency=afterok's job on the cluster. It runs the one stage it was given and
# touches its stamp on success.
#
# Stage-internal execution (local vs. per-step SLURM submission) is decided by the
# stage driver via run_step, not here.

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/common.sh"

usage() {
  cat <<'EOF'
usage: run_stage.sh <stage>

stages:
  preflight                  Stage 0 — permissions + external-asset checks
  data-objects               Stage 0 — build the leaf (deps=-) data objects
  data                       Stage 1 — regenerate data objects
  upload                     Stage 2 — snapshot, diff, plan/apply the bucket upload
  update-relevant-packages   Stage 3 — plan the package updates
EOF
}

stage="${1:-}"
[[ -n "${stage}" ]] || { usage >&2; die "no stage given"; }

STAMP_DIR="${PIPELINE_ROOT}/.stamps"
mkdir -p "${STAMP_DIR}"

case "${stage}" in
  preflight)
    bash "${PIPELINE_ROOT}/scripts/check_existing_outputs.sh"
    bash "${PIPELINE_ROOT}/scripts/00_preflight/00_preflight.sh"
    ;;
  data-objects)
    bash "${PIPELINE_ROOT}/scripts/check_existing_outputs.sh"
    bash "${PIPELINE_ROOT}/scripts/00_preflight/data-raw/build_data_objects.sh"
    ;;
  data)
    bash "${PIPELINE_ROOT}/scripts/10_build_data/10_build_data.sh"
    ;;
  upload)
    bash "${PIPELINE_ROOT}/scripts/20_upload_bucket/20_upload_bucket.sh"
    ;;
  update-relevant-packages)
    bash "${PIPELINE_ROOT}/scripts/30_update_relevant_packages/30_update_relevant_packages.sh"
    ;;
  *)
    usage >&2; die "unknown stage: ${stage}"
    ;;
esac

# data-objects is out of band (nothing depends on it), so it gets no stamp.
[[ "${stage}" == "data-objects" ]] || touch "${STAMP_DIR}/${stage}"
