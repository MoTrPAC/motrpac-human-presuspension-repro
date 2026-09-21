#!/usr/bin/env bash
# SLURM submission layer. Sourced by scripts/lib/common.sh, so every stage script
# already has it. Nothing here runs at source time except reading two env vars.
#
# The unit of submission is the thing the drivers already treat as a unit: a step's
# build.sh, or a whole stage. Drivers call run_step instead of `bash build.sh`, and
# that single call site decides local-vs-cluster. With EXECUTOR=local (the default)
# run_step is a `bash build.sh >log 2>&1` and the pipeline behaves exactly as it did
# before this file existed; with EXECUTOR=slurm the same command is wrapped in a
# generated sbatch script, submitted, and waited on. Per-step PASS/SKIP/FAIL
# accounting, exit codes and log paths are identical either way, which is the whole
# point — a cluster run and a laptop run produce the same report.
#
# Resources come from config/slurm.json, never from this file. See docs/SLURM.md.
#
# Env knobs (all optional):
#   EXECUTOR=slurm|local          override config/slurm.json's "enabled"
#   SLURM_CONFIG=/path/slurm.json use a different config file
#   PRECOVID_SLURM_<FIELD>=...    override one resolved field for every job, e.g.
#                                 PRECOVID_SLURM_ACCOUNT, PRECOVID_SLURM_PARTITION,
#                                 PRECOVID_SLURM_TIME, PRECOVID_SLURM_MEM,
#                                 PRECOVID_SLURM_CPUS_PER_TASK
#   SLURM_DRY_RUN=1               generate + print the sbatch script, submit nothing
#                                 (works with no Slurm installed — this is what makes
#                                 the wiring testable off-cluster)
#   SLURM_JSON_BACKEND=R          force the Rscript config reader instead of jq
#                                 (the two must agree; this is how that is checked)

SLURM_CONFIG="${SLURM_CONFIG:-${PIPELINE_ROOT}/config/slurm.json}"
SLURM_DRY_RUN="${SLURM_DRY_RUN:-0}"

# ---- Reading the config ----------------------------------------------------
# One accessor, two backends. jq if it is on PATH (it usually is on a cluster and
# costs milliseconds); otherwise Rscript + jsonlite, which this pipeline already
# requires. Both implement the same contract: walk a key path, print nothing for a
# missing/null value, print one line per element for an array, one line otherwise.
# Keeping the accessor this dumb is what lets the merge logic below live in bash
# instead of being written twice.
_slurm_json_get() {
  [[ -f "${SLURM_CONFIG}" ]] || return 0
  if [[ "${SLURM_JSON_BACKEND:-auto}" != "R" ]] && have jq; then
    jq -r --args '
      getpath($ARGS.positional)
      | if . == null then empty
        elif type == "array" then .[]
        elif type == "object" then empty
        else . end' "$@" < "${SLURM_CONFIG}" 2>/dev/null
  else
    "${RSCRIPT}" -e '
      a <- commandArgs(trailingOnly = TRUE)
      x <- tryCatch(jsonlite::fromJSON(a[1], simplifyVector = TRUE), error = function(e) NULL)
      for (k in a[-1]) {
        if (is.null(x) || !is.list(x) || !k %in% names(x)) { x <- NULL; break }
        x <- x[[k]]
      }
      if (!is.null(x) && !is.list(x) && length(x)) cat(as.character(x), sep = "\n")
    ' "${SLURM_CONFIG}" "$@" 2>/dev/null
  fi
}

# jq prints JSON booleans as true/false, R prints them as TRUE/FALSE. Normalize
# once here so no caller has to know which backend answered.
_slurm_json_bool() { _slurm_json_get "$@" | tr '[:upper:]' '[:lower:]'; }

# slurm_cfg JOB_KEY FIELD [DEFAULT]
# Resolution order, first non-empty wins:
#   1. PRECOVID_SLURM_<FIELD>   — one env var overrides every job in the run
#   2. jobs.<JOB_KEY>.<FIELD>   — this step's or stage's entry
#   3. defaults.<FIELD>         — the config-wide default
#   4. DEFAULT argument         — last resort, so a config missing a key still submits
slurm_cfg() {
  local job_key="$1" field="$2" fallback="${3:-}" v env_name
  env_name="PRECOVID_SLURM_$(printf '%s' "${field}" | tr '[:lower:]' '[:upper:]')"
  v="${!env_name:-}"
  [[ -n "${v}" ]] && { printf '%s' "${v}"; return 0; }
  v="$(_slurm_json_get jobs "${job_key}" "${field}")"
  [[ -n "${v}" ]] && { printf '%s' "${v}"; return 0; }
  v="$(_slurm_json_get defaults "${field}")"
  [[ -n "${v}" ]] && { printf '%s' "${v}"; return 0; }
  printf '%s' "${fallback}"
}

# Multi-valued fields (modules, pre_commands, sbatch_args) merge rather than
# override: the job's list is appended to the config-wide one, because a step that
# needs one extra module should not have to restate the common ones.
slurm_cfg_list() {
  local job_key="$1" field="$2"
  case "${field}" in
    modules|pre_commands) _slurm_json_get environment "${field}" ;;
    *)                    _slurm_json_get defaults "${field}" ;;
  esac
  _slurm_json_get jobs "${job_key}" "${field}"
}

# ---- Are we submitting? ----------------------------------------------------
# EXECUTOR wins over the config so a one-off `EXECUTOR=local make data` works on a
# cluster checkout, and `EXECUTOR=slurm` works without editing the config at all.
slurm_enabled() {
  case "${EXECUTOR:-}" in
    slurm) return 0 ;;
    local) return 1 ;;
  esac
  [[ "$(_slurm_json_bool enabled)" == "true" ]]
}

# Fail loudly at the first submission rather than silently running a 40-hour job on
# a login node because sbatch was missing.
slurm_require_tools() {
  [[ "${SLURM_DRY_RUN}" == "1" ]] && return 0
  local t
  for t in sbatch squeue scancel; do
    have "${t}" || die "EXECUTOR=slurm but ${t} is not on PATH — is this a submit host?"
  done
  local acct; acct="$(slurm_cfg "$1" account)"
  [[ -n "${acct}" && "${acct}" != "CHANGE_ME" ]] \
    || die "config/slurm.json: defaults.account is unset — put your SCG account there (see docs/SLURM.md)."
}

slurm_log_dir() {
  local d; d="$(_slurm_json_get logs dir)"
  d="${d:-logs/slurm}"
  [[ "${d}" == /* ]] || d="${PIPELINE_ROOT}/${d}"
  printf '%s' "${d}"
}

# Slurm job names cannot usefully carry the "/" in a step key.
_slurm_job_name() { printf 'precovid-%s' "${1//\//-}"; }

# ---- Writing the batch script ----------------------------------------------
# The generated script is a complete, self-contained record of how the job ran:
# every #SBATCH directive, every module, the exact command. It is kept under
# logs/slurm/scripts/ so a failed job can be re-run by hand with `sbatch <script>`
# and so the resources a given run actually asked for are recoverable afterwards.
#
# slurm_write_script JOB_KEY LOG_FILE RC_FILE CMD...   -> prints the script path
slurm_write_script() {
  local job_key="$1" logfile="$2" rcfile="$3"; shift 3
  local sdir; sdir="$(slurm_log_dir)/scripts"
  mkdir -p "${sdir}"
  local script="${sdir}/$(_slurm_job_name "${job_key}")_$(date +'%Y%m%d_%H%M%S')_$$.sbatch"

  local account partition qos time nodes ntasks cpus mem constraint mail_type mail_user
  account="$(slurm_cfg   "${job_key}" account)"
  partition="$(slurm_cfg "${job_key}" partition)"
  qos="$(slurm_cfg       "${job_key}" qos)"
  time="$(slurm_cfg      "${job_key}" time      "04:00:00")"
  nodes="$(slurm_cfg     "${job_key}" nodes     "1")"
  ntasks="$(slurm_cfg    "${job_key}" ntasks    "1")"
  cpus="$(slurm_cfg      "${job_key}" cpus_per_task "1")"
  mem="$(slurm_cfg       "${job_key}" mem       "16G")"
  constraint="$(slurm_cfg "${job_key}" constraint)"
  mail_type="$(slurm_cfg "${job_key}" mail_type)"
  mail_user="$(slurm_cfg "${job_key}" mail_user)"

  {
    printf '#!/usr/bin/env bash\n'
    printf '#SBATCH --job-name=%s\n'      "$(_slurm_job_name "${job_key}")"
    printf '#SBATCH --output=%s\n'        "${logfile}"
    printf '#SBATCH --error=%s\n'         "${logfile}"
    printf '#SBATCH --chdir=%s\n'         "${PIPELINE_ROOT}"
    [[ -n "${account}"    ]] && printf '#SBATCH --account=%s\n'       "${account}"
    [[ -n "${partition}"  ]] && printf '#SBATCH --partition=%s\n'     "${partition}"
    [[ -n "${qos}"        ]] && printf '#SBATCH --qos=%s\n'           "${qos}"
    [[ -n "${time}"       ]] && printf '#SBATCH --time=%s\n'          "${time}"
    [[ -n "${nodes}"      ]] && printf '#SBATCH --nodes=%s\n'         "${nodes}"
    [[ -n "${ntasks}"     ]] && printf '#SBATCH --ntasks=%s\n'        "${ntasks}"
    [[ -n "${cpus}"       ]] && printf '#SBATCH --cpus-per-task=%s\n' "${cpus}"
    [[ -n "${mem}"        ]] && printf '#SBATCH --mem=%s\n'           "${mem}"
    [[ -n "${constraint}" ]] && printf '#SBATCH --constraint=%s\n'    "${constraint}"
    [[ -n "${mail_type}"  ]] && printf '#SBATCH --mail-type=%s\n'     "${mail_type}"
    [[ -n "${mail_user}"  ]] && printf '#SBATCH --mail-user=%s\n'     "${mail_user}"
    local extra; while IFS= read -r extra; do
      [[ -n "${extra}" ]] && printf '#SBATCH %s\n' "${extra}"
    done < <(slurm_cfg_list "${job_key}" sbatch_args)

    printf '\nset -euo pipefail\n'
    printf 'cd %q\n' "${PIPELINE_ROOT}"

    local m; while IFS= read -r m; do
      [[ -n "${m}" ]] && printf 'module load %q\n' "${m}"
    done < <(slurm_cfg_list "${job_key}" modules)

    local pc; while IFS= read -r pc; do
      [[ -n "${pc}" ]] && printf '%s\n' "${pc}"
    done < <(slurm_cfg_list "${job_key}" pre_commands)

    local tmpdir; tmpdir="$(_slurm_json_get environment tmpdir)"
    [[ -n "${tmpdir}" ]] && printf 'export TMPDIR=%q\n' "${tmpdir}"

    # A job never submits jobs. Step-granularity submission happens on the submit
    # host; stage-granularity submission (scripts/submit_slurm.sh) puts a whole
    # stage inside one allocation and its steps must then run there, not spawn a
    # second generation of jobs from a compute node.
    printf 'export EXECUTOR=local\n'

    # The two core-count knobs in pipeline.env are written as ${VAR:-default}, so an
    # exported value wins. Binding them to the allocation is the only way a job that
    # asked for 16 CPUs actually uses 16 — the defaults are tuned for a 12-core Mac,
    # and a fork-based worker pool that ignores its cgroup will oversubscribe the
    # node. Turn this off with environment.bind_r_cores=false if a step should keep
    # the pipeline.env values.
    if [[ "$(_slurm_json_bool environment bind_r_cores)" != "false" ]]; then
      printf 'export PARALLEL_CORES="${PARALLEL_CORES:-${SLURM_CPUS_PER_TASK:-1}}"\n'
      printf 'export VARIANCEPARTITION_PARALLEL_CORES="${VARIANCEPARTITION_PARALLEL_CORES:-${SLURM_CPUS_PER_TASK:-1}}"\n'
    fi

    printf '\necho "[slurm] job ${SLURM_JOB_ID:-?} on $(hostname) — %s"\n' "${job_key}"
    printf 'echo "[slurm] cpus=${SLURM_CPUS_PER_TASK:-?} mem=%s time=%s"\n' "${mem}" "${time}"

    # The exit code is written where the waiter can read it without depending on the
    # cluster's accounting database. sacct is the fallback, not the source of truth:
    # it lags, it can be disabled, and it cannot report a job that Slurm killed
    # before the payload ran. A file the job itself wrote is unambiguous.
    printf 'rc=0\n'
    printf '%q ' "$@"; printf '|| rc=$?\n'
    printf 'printf "%%s\\n" "${rc}" > %q\n' "${rcfile}"
    printf 'echo "[slurm] exit ${rc}"\n'
    printf 'exit "${rc}"\n'
  } > "${script}"

  chmod +x "${script}"
  printf '%s' "${script}"
}

# ---- Submitting and waiting ------------------------------------------------
# slurm_submit JOB_KEY LOG_FILE RC_FILE [DEPENDENCY] -- CMD...  -> prints job id
slurm_submit() {
  local job_key="$1" logfile="$2" rcfile="$3" dependency="${4:-}"; shift 4
  [[ "${1:-}" == "--" ]] && shift

  # Everything this function says goes to stderr: its stdout is the job id, and a
  # caller captures it.
  local script; script="$(slurm_write_script "${job_key}" "${logfile}" "${rcfile}" "$@")"
  if [[ "${SLURM_DRY_RUN}" == "1" ]]; then
    log "DRY RUN — would submit ${job_key}${dependency:+ (after ${dependency})}:" >&2
    sed 's/^/    /' "${script}" >&2
    printf 'DRYRUN'
    return 0
  fi

  local args=(--parsable)
  [[ -n "${dependency}" ]] && args+=(--dependency="afterok:${dependency}" --kill-on-invalid-dep=yes)
  local export_mode; export_mode="$(_slurm_json_get environment export)"
  [[ -n "${export_mode}" ]] && args+=(--export="${export_mode}")

  local jid
  jid="$(sbatch "${args[@]}" "${script}" | tr -d '[:space:]')" \
    || die "sbatch failed for ${job_key} (script: ${script#"${PIPELINE_ROOT}"/})"
  [[ -n "${jid}" ]] || die "sbatch returned no job id for ${job_key}"
  printf '%s' "${jid}"
}

# slurm_wait JOB_ID RC_FILE -> returns the job's exit code
#
# Polls squeue rather than using `sbatch --wait` so the state transitions
# (PENDING -> RUNNING -> gone) are visible in the stage log while a multi-hour step
# runs, and so the same submit path serves both the blocking driver and the
# fire-and-forget chain submitter.
slurm_wait() {
  local jid="$1" rcfile="$2"
  local poll; poll="$(_slurm_json_get poll_seconds)"; poll="${poll:-30}"
  local state last_state=""

  # An interrupted driver must not leave the job running: the next run would then
  # have two jobs writing the same freeze files.
  trap 'warn "interrupted — cancelling job '"${jid}"'"; scancel '"${jid}"' 2>/dev/null || true' INT TERM

  # squeue exiting non-zero is not the same as squeue reporting nothing. It fails
  # transiently when the scheduler is loaded, and permanently once a finished job
  # ages past MinJobAge ("Invalid job id"). Treating either as "job gone" the first
  # time would abandon a running job; treating both as fatal would hang on a job that
  # merely aged out. So tolerate a few failures, then fall through to the rc file,
  # which answers both cases correctly.
  local out fails=0
  while true; do
    if out="$(squeue -h -j "${jid}" -o '%T' 2>/dev/null)"; then
      fails=0
      state="$(printf '%s' "${out}" | head -n1)"
      [[ -z "${state}" ]] && break
      [[ "${state}" != "${last_state}" ]] && { log "job ${jid}: ${state}"; last_state="${state}"; }
    else
      fails=$((fails + 1))
      (( fails >= 3 )) && break
      warn "squeue failed for job ${jid} (${fails}/3) — retrying"
    fi
    sleep "${poll}"
  done
  trap - INT TERM

  # squeue drops a job a few seconds before its stdout is flushed and before sacct
  # sees it, so give the rc file a moment to land.
  local i
  for i in 1 2 3 4 5 6; do
    [[ -f "${rcfile}" ]] && break
    sleep 2
  done

  if [[ -f "${rcfile}" ]]; then
    local rc; rc="$(tr -dc '0-9' < "${rcfile}")"
    return "${rc:-1}"
  fi

  # No rc file: the payload never finished. Slurm killed it (TIMEOUT, OUT_OF_MEMORY,
  # NODE_FAIL, cancelled) or the node died. Report what accounting knows, if it knows.
  local sacct_state=""
  if have sacct; then
    sacct_state="$(sacct -n -X -j "${jid}" -o State%30 2>/dev/null | head -n1 | tr -d ' ')"
  fi
  err "job ${jid} produced no exit code${sacct_state:+ — sacct state: ${sacct_state}}"
  return 1
}

# Elapsed time and peak RSS for the run report, so the next edit to slurm.json is
# informed by what the job used rather than by a guess. Best effort: silent when
# accounting is unavailable.
slurm_usage() {
  local jid="$1"
  have sacct || return 0
  sacct -n -j "${jid}" -o Elapsed,MaxRSS,State --units=G 2>/dev/null \
    | awk 'NF && $2 != "" {printf "elapsed=%s maxrss=%s state=%s", $1, $2, $3; exit}'
}

# ---- The one call site the drivers use -------------------------------------
# run_step JOB_KEY LOG_FILE CMD...   -> the command's exit code
#
# Local: run CMD with output to LOG_FILE. Slurm: submit CMD as a job whose stdout is
# LOG_FILE and block until it finishes. Both leave the same artefact behind, which is
# why the drivers' PASS/SKIP/FAIL logic needs no branch of its own.
run_step() {
  local job_key="$1" logfile="$2"; shift 2
  local rc

  # errexit is a global shell option, not a function-local one: a bare `set -e` in
  # here would switch it back on in the CALLER, whose next non-zero return would then
  # kill the driver. A step exiting 77 (SKIP) is normal in this pipeline, so that
  # would turn every skipped step into an aborted stage. Save what the caller had,
  # run with errexit off so a failing step returns its code, restore on the way out.
  local had_errexit=0
  [[ $- == *e* ]] && had_errexit=1
  set +e

  if slurm_enabled; then
    slurm_require_tools "${job_key}"
    local rcdir; rcdir="$(slurm_log_dir)/rc"
    mkdir -p "${rcdir}" "$(dirname "${logfile}")"
    local rcfile="${rcdir}/$(_slurm_job_name "${job_key}").$$.rc"
    rm -f "${rcfile}"

    # slurm_submit runs in a command substitution, so its `die` can only kill that
    # subshell — an empty job id is how a failed submission arrives here.
    local jid; jid="$(slurm_submit "${job_key}" "${logfile}" "${rcfile}" "" -- "$@")"
    if [[ -z "${jid}" ]]; then
      if (( had_errexit )); then set -e; fi
      die "submission failed for ${job_key} (see the sbatch error above)"
    fi

    if [[ "${jid}" == "DRYRUN" ]]; then
      log "DRY RUN — not waiting on ${job_key}"
      rc=0
    else
      log "submitted ${job_key} as job ${jid} (cpus=$(slurm_cfg "${job_key}" cpus_per_task), mem=$(slurm_cfg "${job_key}" mem), time=$(slurm_cfg "${job_key}" time))"
      slurm_wait "${jid}" "${rcfile}"
      rc=$?
      local usage; usage="$(slurm_usage "${jid}")"
      [[ -n "${usage}" ]] && log "job ${jid} ${usage}"
      rm -f "${rcfile}"
    fi
  else
    "$@" > "${logfile}" 2>&1
    rc=$?
  fi

  if (( had_errexit )); then set -e; fi
  return "${rc}"
}
