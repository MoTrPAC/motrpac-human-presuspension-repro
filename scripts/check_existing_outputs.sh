#!/usr/bin/env bash
# Existing-output helper — one script, two modes, over a single list of the
# pipeline's expected-output locations (the `OUTPUTS` array below is the source of
# truth; edit it there and both modes follow).
#
# Why it exists: the build stages write their outputs IN PLACE and overwrite. So a
# file that a *given* run happens not to produce (a dropped ome, a renamed version)
# is not overwritten — it just lingers and can be mistaken for current output. This
# helper either warns about that up front, or clears the slate for a clean rebuild.
#
# The managed output locations:
#   staging/freeze                     qc-norm / DA freeze files (Stage 1)
#   staging/raw-files                  cached bucket downloads (raw inputs + test refs)
#   staging/scion                      step-16 SCION network tables
#   scripts/10_build_data/data         Stage-1 .rda data objects
#   scripts/00_preflight/data          preflight leaf .rds objects
#   scripts/00_preflight/output        preflight gmt.gz + json
#   .../13_build_fcm/fcm_diagnostics   step-13 cluster-count sweep (tsv/png/pdf)
# `logs/` is NOT in that list; it is only touched by the opt-in `logs` argument.
#
# ── Modes ─────────────────────────────────────────────────────────────────────
#   check   (default)  Warn (never fail) if any managed location already holds
#                      files, then print how to clean. The Makefile runs this at the
#                      start of every build via the `precheck` target, so you see the
#                      warning automatically before anything is overwritten.
#   delete             Remove every file under each managed location (directories are
#                      kept, contents wiped). Prints a per-location count. Add a second
#                      argument `logs` to ALSO wipe logs/ (run reports + per-object
#                      logs). This is what `CLEAN=1 ./run_all.sh` calls.
#
# ── Examples ──────────────────────────────────────────────────────────────────
#   check_existing_outputs.sh                 # warn about existing outputs (what `make` does)
#   check_existing_outputs.sh delete          # wipe every managed location, keep logs/
#   check_existing_outputs.sh delete logs     # wipe outputs AND logs/
#
#   CLEAN=1 ./run_all.sh                       # from-scratch rebuild (calls `delete`)
#   CLEAN=1 CLEAN_LOGS=1 ./run_all.sh          # from-scratch rebuild, also clears logs/
#
#   NODES-style targeted cleanup is not offered here — to wipe just one location,
#   delete it by hand, e.g.  rm -rf staging/freeze
#
# Exit status is always 0 (check never blocks; delete reports and returns success).
#
# Usage: check_existing_outputs.sh [delete [logs]]
set -uo pipefail
MODE="${1:-check}"
DELETE_LOGS=0; case "${2:-}" in logs|--logs) DELETE_LOGS=1;; esac
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
c_yel=$'\033[33m'; c_rst=$'\033[0m'

# label -> path (relative to repo root)  ::  the expected-output locations
declare -a OUTPUTS=(
  "freeze (qc-norm/DA)::staging/freeze"
  "bucket-download cache::staging/raw-files"
  "SCION networks::staging/scion"
  "stage-1 data objects::scripts/10_build_data/data"
  "preflight data objects::scripts/00_preflight/data"
  "preflight output (gmt/json)::scripts/00_preflight/output"
  "FCM diagnostics::scripts/10_build_data/13_build_fcm/fcm_diagnostics"
)

if [[ "${MODE}" == "delete" ]]; then
  total=0
  for entry in "${OUTPUTS[@]}"; do
    label="${entry%%::*}"; rel="${entry##*::}"; d="${ROOT}/${rel}"
    [[ -d "$d" ]] || continue
    n=$(find "$d" -type f 2>/dev/null | wc -l | tr -d ' ')
    (( n > 0 )) || continue
    find "$d" -mindepth 1 -delete 2>/dev/null || rm -rf "${d:?}"/* 2>/dev/null || true
    printf '  deleted %5s file(s) from %s\n' "$n" "$rel"
    total=$((total + n))
  done
  printf '  Removed %s stale output file(s) across %s location(s).\n' "$total" "${#OUTPUTS[@]}"
  if (( DELETE_LOGS )); then
    d="${ROOT}/logs"
    if [[ -d "$d" ]]; then
      n=$(find "$d" -type f 2>/dev/null | wc -l | tr -d ' ')
      find "$d" -mindepth 1 -delete 2>/dev/null || rm -rf "${d:?}"/* 2>/dev/null || true
      printf '  deleted %5s log file(s) from logs/\n' "$n"
    fi
  fi
  exit 0
fi

found=0
for entry in "${OUTPUTS[@]}"; do
  label="${entry%%::*}"; rel="${entry##*::}"; d="${ROOT}/${rel}"
  [[ -d "$d" ]] || continue
  n=$(find "$d" -type f 2>/dev/null | wc -l | tr -d ' ')
  (( n > 0 )) || continue
  printf '  %s⚠%s %-24s %4s existing file(s) in %s\n' "$c_yel" "$c_rst" "$label" "$n" "$rel"
  found=1
done

if (( found )); then
  printf '\n  %sExisting outputs will be OVERWRITTEN in place by this build.%s\n' "$c_yel" "$c_rst"
  printf '  To regenerate cleanly from scratch, run:  %sCLEAN=1 ./run_all.sh%s\n' "$c_yel" "$c_rst"
  printf '  (or delete a folder by hand, e.g. rm -rf staging/freeze)\n\n'
else
  printf '  No existing outputs found — clean build.\n'
fi
exit 0
