# Running the pipeline on SLURM (SCG)

The pipeline runs unchanged on a laptop. This document is about the other mode: every
script executed as a batch job on Stanford's SCG cluster, with resources and account
read from `config/slurm.json`.

Nothing about the DAG changes. The stages, the steps, the per-step `PASS/SKIP/FAIL`
report, the log paths and the exit codes are identical whether a step ran here or on a
compute node — that is the design constraint the whole layer is built around, because a
cluster run whose report cannot be compared to a local run is not a reproduction.

## The one switch

```bash
make data EXECUTOR=slurm     # every STEP of stage 1 becomes its own job
make data                    # unchanged: everything runs right here
```

`EXECUTOR` overrides `"enabled"` in the config, so a cluster checkout can set
`"enabled": true` once and forget it, and still run something locally with
`EXECUTOR=local make data`.

## Two granularities, and when each is right

|  | `EXECUTOR=slurm make data` | `make slurm` |
|---|---|---|
| unit of submission | one job per **step** | one job per **stage** |
| resources | per step, from `jobs."10_build_data/09_build_da"` | per stage, from `jobs.data` |
| ordering | the driver waits for each job in turn | `--dependency=afterok` between stages |
| your shell | must stay alive for the whole stage | free the moment the jobs are queued |
| good for | sizing or re-running one expensive step | the overnight full run |

**Per step.** The stage driver submits a step, waits for it, records the result, submits
the next. Each step gets exactly the CPUs, memory and walltime it needs, which matters
here because the range is wide — `05_stage_metabolomics_cvs` is a file copy and
`09_build_da` fits mixed models for hours on 16 cores. The cost is a login-node shell
babysitting the queue for as long as the stage takes.

**Per stage.** `make slurm` submits `preflight → data → {upload, update-relevant-packages}`
as four dependent jobs and returns. `upload` and `update-relevant-packages` both hang off
`data`, so they run concurrently. Each stage job runs its own steps *inside its one
allocation*, which means the allocation must be sized for the heaviest step in the stage —
stage 1 asks for 16 CPUs and 128 GB for 48 hours because step 09 needs it, and the file
copies inherit that. Wasteful, and the right trade when the alternative is a shell that
has to survive two days.

```bash
make slurm                       # submit the chain, print job ids, return
make slurm STAGES="data upload"  # part of it (dependencies fall back to the nearest selected ancestor)
make slurm DRY_RUN=1             # print the generated batch scripts, submit nothing
make slurm-status                # squeue + sacct for the chain last submitted
make slurm-cancel                # scancel that chain
APPLY=1 make slurm               # stage 2 actually writes to the staging bucket
```

A job never submits jobs: every generated script sets `EXECUTOR=local`, so a stage job's
steps run inside its allocation instead of spawning a second generation of jobs from a
compute node.

## config/slurm.json

Every field resolves in this order, first non-empty winning:

1. `PRECOVID_SLURM_<FIELD>` in the environment — one override for the whole run, e.g.
   `PRECOVID_SLURM_PARTITION=interactive make data EXECUTOR=slurm`
2. `jobs.<job key>.<field>` — that step's or stage's entry
3. `defaults.<field>` — the config-wide default

A **job key** is either a stage name (`preflight`, `data`, `upload`,
`update-relevant-packages`, `data-objects`) or a step path
(`10_build_data/09_build_da`, `20_upload_bucket/04_upload_to_staging`). A step with no
entry simply gets `defaults`, so adding a step to a stage never breaks submission.

The env-var layer is namespaced `PRECOVID_SLURM_*` rather than `SLURM_*` on purpose:
Slurm sets its own `SLURM_*` variables inside every job, and a config layer that read
those would change its answers depending on whether it was consulted from a login node
or a compute node.

```jsonc
{
  "enabled": false,          // true on a cluster checkout; EXECUTOR overrides it
  "defaults": {
    "account": "CHANGE_ME",  // your SCG account — preflight FAILS while this is unset
    "partition": "batch",
    "time": "04:00:00",
    "cpus_per_task": 4,
    "mem": "32G",
    "mail_type": "",         // "FAIL" to be told when a job dies
    "mail_user": "",
    "sbatch_args": []        // raw directives, e.g. ["--exclusive"]
  },
  "environment": {
    "modules": [],           // e.g. ["R/4.4.1", "gcloud"] — module load, in order
    "pre_commands": [],      // shell lines run after the modules, before the step
    "export": "ALL",         // sbatch --export; ALL passes STEPS=, APPLY=, RERUN_ATAC= through
    "bind_r_cores": true
  },
  "jobs": { "10_build_data/09_build_da": { "cpus_per_task": 16, "mem": "192G" } }
}
```

`modules` and `pre_commands` **merge** rather than override: a step's list is appended to
the config-wide one, so a step needing one extra module does not have to restate the
common ones. Scalars override.

### Before the first run

1. Set `defaults.account`. Preflight fails while it reads `CHANGE_ME`, and it also warns
   if `sacctmgr` says you have no association with that account — an account you are not
   associated with is accepted by `sbatch` and then rejected by the scheduler, so the job
   simply never starts and nothing says why.
2. Set `environment.modules` to whatever SCG needs for R 4.4.1 and the Google Cloud SDK.
   The batch script has no interactive shell profile behind it, so an `Rscript` that
   works when you log in may not exist inside a job.
3. Check that `gcloud` credentials are visible from a compute node. Stage 2 and parts of
   stage 1 read consortium buckets; on SCG `~/.config/gcloud` is on shared storage, so
   they are, but this is the first thing to check if a job fails on bucket access alone.
4. `make slurm DRY_RUN=1` and read the generated scripts.

### bind_r_cores

The two core-count knobs in `config/pipeline.env` are tuned for a 12-core Mac
(`PARALLEL_CORES=10`, `VARIANCEPARTITION_PARALLEL_CORES=8`). Left alone, a job that asked
for 16 CPUs would use 8 of them, and a job that asked for 2 would try to fork 8 workers
into a 2-CPU cgroup. With `bind_r_cores: true` (the default) the generated script exports
both from `SLURM_CPUS_PER_TASK`, so the allocation and the worker pool agree. Both are
written as `${VAR:-default}` in `pipeline.env`, so an explicit
`PARALLEL_CORES=4 make data EXECUTOR=slurm` still wins.

Because `variancePartition` runs on forked workers over a parent that reaches ~4 GB,
**`mem` and `cpus_per_task` have to move together** in this file. Raising cores without
raising memory is how a stage-1 job gets OOM-killed at hour six.

The resource numbers shipped in `config/slurm.json` are starting points, not
measurements. After a real run, `make slurm-status` prints `Elapsed` and `MaxRSS` per
job, and per-step jobs log the same line into the stage log — tune from those.

## Where things land

```
logs/slurm/scripts/<job>_<ts>_<pid>.sbatch  the generated batch script, kept
logs/slurm/rc/<job>.<pid>.rc           the exit code a single job wrote for itself
logs/slurm/rc/<stage>_<ts>.rc          the same, in chain mode, which keys on the stage
logs/slurm/<stage>_<ts>.out            stage-job stdout (chain mode)
logs/slurm/chain_<ts>.tsv              stage -> job id -> dependency, per submission
logs/slurm/chain_latest.tsv            what slurm-status and slurm-cancel read
logs/data_<step>.log                   unchanged: the step's own log, wherever it ran
```

The generated script is a complete record of how a job ran — every directive, every
module, the exact command — so a failed job can be re-run by hand with
`sbatch logs/slurm/scripts/<that file>`, and the resources a past run asked for stay
recoverable. `logs/` is gitignored and `make clean` clears it.

**The exit code comes from a file the job wrote, not from `sacct`.** Accounting lags,
can be disabled, and cannot report a job Slurm killed before the payload started. The
`.rc` file is unambiguous when it exists; when it does not, the job died without running
to completion (`TIMEOUT`, `OUT_OF_MEMORY`, `NODE_FAIL`, cancelled) and the waiter reports
whatever `sacct` knows as the reason.

Interrupting a driver that is waiting on a job cancels that job. Leaving it running would
put two jobs on the same freeze files the next time the stage was started.

## Implementation

| File | Role |
|---|---|
| `config/slurm.json` | account, partitions, per-job resources, modules |
| `scripts/lib/slurm.sh` | config resolution, script generation, submit/wait; defines `run_step` |
| `scripts/run_stage.sh` | the single definition of what each stage runs — used by the Makefile *and* by batch jobs |
| `scripts/submit_slurm.sh` | the chained stage submission, plus `status` and `cancel` |

`run_step JOB_KEY LOG_FILE CMD...` is the only place the local/cluster decision is made.
The stage drivers call it where they used to call `bash build.sh >log 2>&1`; local mode
*is* that command. Everything downstream — the report rows, the `STOP_ON_FAIL` handling,
the `77 = SKIP` convention — needs no branch of its own.

The config is read by `jq` when it is on `PATH` and by `Rscript` + `jsonlite` otherwise
(this pipeline already requires both R and that package), which keeps submission working
on a cluster with no `jq` module. The two backends have to agree; `SLURM_JSON_BACKEND=R`
forces the fallback so that can be checked.

### Not done here

Submission granularity stops at the step. Step 06 runs 13 qc-norm stems and step 09 runs
one DA stem per ome-tissue inside a single job, and several of those are independent —
an obvious array-job win. It is not built because the stems are ordered (the imputed prot
stems consume the qc-norm the same step just wrote) and that ordering lives inside the
step, not in a manifest a submitter could read. Giving those two steps a large allocation
recovers most of the parallelism through R's own worker pools.
