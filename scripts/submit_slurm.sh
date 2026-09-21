#!/usr/bin/env bash
# Submit the whole DAG to SLURM as one chain of dependent jobs, and return.
#
#   bash scripts/submit_slurm.sh              # submit preflight -> data -> {upload, packages}
#   bash scripts/submit_slurm.sh --stages "data upload"
#   bash scripts/submit_slurm.sh --dry-run    # print the batch scripts, submit nothing
#   bash scripts/submit_slurm.sh status       # squeue for the last chain submitted
#   bash scripts/submit_slurm.sh cancel       # scancel that chain
#
# Two granularities exist and this is the coarse one. `EXECUTOR=slurm make data`
# submits each STEP as its own job and blocks on the login node until the stage
# finishes — right for a single expensive step you want to size precisely, wrong for
# an overnight run because it needs your shell to stay alive. This submits each
# STAGE as one job and exits: the whole pipeline then runs unattended, ordered by
# --dependency=afterok, and each stage's steps run inside its own allocation. Because
# the generated job scripts set EXECUTOR=local, a stage job never spawns a second
# generation of jobs from a compute node.
#
# Resources per stage come from the "jobs" entries in config/slurm.json keyed by
# stage name (preflight, data, upload, update-relevant-packages) — a stage job must
# be sized for its heaviest step, since every step shares its one allocation.
#
# Env passes through to the jobs (--export=ALL), so the usual knobs work:
#   APPLY=1 bash scripts/submit_slurm.sh          # stage 2 actually uploads
#   STEPS="06 09" bash scripts/submit_slurm.sh --stages data

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/common.sh"

# An explicit submit command submits, whatever "enabled" says in the config.
export EXECUTOR=slurm

# The stage DAG: submission order, and each stage's parent. upload and
# update-relevant-packages both hang off data and therefore run concurrently.
STAGE_ORDER=(preflight data upload update-relevant-packages)
stage_parent() {
  case "$1" in
    preflight)                ;;
    data)                     printf 'preflight' ;;
    upload)                   printf 'data' ;;
    update-relevant-packages) printf 'data' ;;
  esac
}

STAGES="${STAGE_ORDER[*]}"
CHAIN_DIR="$(slurm_log_dir)"
CHAIN_LATEST="${CHAIN_DIR}/chain_latest.tsv"

# ---- status / cancel subcommands -------------------------------------------
chain_ids() {
  [[ -f "${CHAIN_LATEST}" ]] || die "no chain recorded yet (${CHAIN_LATEST#"${PIPELINE_ROOT}"/})"
  awk 'NR > 1 && $2 != "DRYRUN" {printf "%s%s", sep, $2; sep=","}' "${CHAIN_LATEST}"
}

case "${1:-submit}" in
  status)
    ids="$(chain_ids)"
    [[ -n "${ids}" ]] || die "chain record holds no job ids"
    log "chain: ${ids}"
    squeue -j "${ids}" 2>/dev/null || true
    have sacct && sacct -X -j "${ids}" -o JobID,JobName%28,State,Elapsed,MaxRSS,ExitCode --units=G 2>/dev/null || true
    exit 0
    ;;
  cancel)
    ids="$(chain_ids)"
    warn "cancelling ${ids}"
    scancel "${ids//,/ }"
    exit 0
    ;;
  submit) shift ;;
esac

while (( $# )); do
  case "$1" in
    --stages) STAGES="$2"; shift 2 ;;
    --dry-run) export SLURM_DRY_RUN=1; shift ;;
    -h|--help) sed -n '2,26p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

selected() { local s; for s in ${STAGES}; do [[ "$1" == "${s}" ]] && return 0; done; return 1; }

slurm_require_tools data
mkdir -p "${CHAIN_DIR}/rc"
ts="$(date +'%Y%m%d_%H%M%S')"
chain_tsv="${CHAIN_DIR}/chain_${ts}.tsv"
printf 'stage\tjob_id\tdepends_on\tlog\n' > "${chain_tsv}"

log "submitting chain (stages: ${STAGES})"
[[ "${SLURM_DRY_RUN:-0}" == "1" ]] && warn "DRY RUN — nothing will be submitted"

declare -a JOB_STAGES=() JOB_IDS=()
job_id_of() {
  local i
  for i in "${!JOB_STAGES[@]}"; do
    [[ "${JOB_STAGES[$i]}" == "$1" ]] && { printf '%s' "${JOB_IDS[$i]}"; return 0; }
  done
  return 1
}

# The parent a job actually waits on is its nearest SELECTED ancestor. Submitting
# `--stages "data upload"` from a checkout whose preflight already passed must not
# invent a dependency on a job that was never submitted.
resolve_dep() {
  local p="$1"
  while [[ -n "${p}" ]]; do
    if selected "${p}"; then job_id_of "${p}" && return 0; fi
    p="$(stage_parent "${p}")"
  done
  return 1
}

for stage in "${STAGE_ORDER[@]}"; do
  selected "${stage}" || continue

  dep="$(resolve_dep "$(stage_parent "${stage}")" || true)"
  logfile="${CHAIN_DIR}/${stage}_${ts}.out"
  rcfile="${CHAIN_DIR}/rc/${stage}_${ts}.rc"

  # slurm_submit's `die` cannot escape a command substitution, so a failed sbatch
  # arrives as an empty id. Stop here rather than submitting the rest of the chain
  # with a broken dependency.
  jid="$(slurm_submit "${stage}" "${logfile}" "${rcfile}" "${dep}" -- \
          bash "${PIPELINE_ROOT}/scripts/run_stage.sh" "${stage}")"
  [[ -n "${jid}" ]] || die "submission failed at stage ${stage}; chain so far: ${chain_tsv#"${PIPELINE_ROOT}"/}"

  JOB_STAGES+=("${stage}"); JOB_IDS+=("${jid}")
  printf '%s\t%s\t%s\t%s\n' "${stage}" "${jid}" "${dep:--}" "${logfile#"${PIPELINE_ROOT}"/}" >> "${chain_tsv}"
  ok "${stage} -> job ${jid}${dep:+ (after ${dep})} — cpus=$(slurm_cfg "${stage}" cpus_per_task) mem=$(slurm_cfg "${stage}" mem) time=$(slurm_cfg "${stage}" time)"
done

cp "${chain_tsv}" "${CHAIN_LATEST}"
log "chain recorded: ${chain_tsv#"${PIPELINE_ROOT}"/}"
if [[ "${SLURM_DRY_RUN:-0}" != "1" ]]; then
  log "watch:  make slurm-status     cancel:  make slurm-cancel"
fi
