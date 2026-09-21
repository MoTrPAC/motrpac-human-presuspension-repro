#!/usr/bin/env bash
# Convenience wrapper: run the full DAG with a timestamped combined log.
# Usage:
#   ./run_all.sh                        # normal run (outputs overwritten in place)
#   CLEAN=1 ./run_all.sh                # delete ALL existing outputs first (via
#                                       # scripts/check_existing_outputs.sh), then
#                                       # rebuild from scratch. Clears stale leftovers
#                                       # such as a prior-version freeze file lingering
#                                       # beside the current one. Because it also wipes
#                                       # the 00_preflight leaf objects (which `make all`
#                                       # does not rebuild), it re-runs `make data-objects`
#                                       # after the wipe.
#   CLEAN=1 CLEAN_LOGS=1 ./run_all.sh   # also wipe logs/ (old run reports + per-object logs)
set -euo pipefail
cd "$(dirname "$0")"

# The clean-up runs BEFORE the tee log below is opened, so this run's own log survives.
if [[ "${CLEAN:-0}" == "1" ]]; then
  if [[ "${CLEAN_LOGS:-0}" == "1" ]]; then
    echo "CLEAN=1 CLEAN_LOGS=1 — deleting existing outputs + logs + stamps for a from-scratch rebuild"
    bash scripts/check_existing_outputs.sh delete logs
  else
    echo "CLEAN=1 — deleting existing outputs + stamps for a from-scratch rebuild"
    bash scripts/check_existing_outputs.sh delete
  fi
  rm -rf .stamps                       # force every stage to re-run
fi

mkdir -p logs
ts="$(date +'%Y%m%d_%H%M%S')"

run() {
  [[ "${CLEAN:-0}" == "1" ]] && make data-objects   # rebuild the leaf objects wiped above
  make all
}

run 2>&1 | tee "logs/run_all_${ts}.log"
