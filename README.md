# motrpac-human-presuspension-repro

End-to-end reproduction pipeline for the MoTrPAC Human Pre-Suspension results, run as a single
dependency graph. Scoped to two repos: `MotrpacHumanPreSuspensionData` (data objects) and
`MotrpacHumanPreSuspensionAnalysis` (summary/DA/enrichment + the GCS upload pipeline it wraps).

## How this repo fits with the others

The Pre-Suspension human work is split across four repositories. Individual-level molecular
and phenotypic data cannot be distributed publicly, so they live in a separate, access-gated
package, and everything that *can* be released publicly (aggregate results, all analysis
code) is kept clear of them.

```
                MoTrPAC BIC — consortium GCS buckets
                     (raw assay data, gated)
                                │
                         motrpac-human-presuspension-repro
        normalizes omics data, applies statistical models,
        builds every data object, versions it, uploads it,
                and carries it into both packages
                                │
              ┌─────────────────┴─────────────────┐
              ▼                                   ▼
 MotrpacHumanPreSuspensionData      MotrpacHumanPreSuspensionAnalysis
 subject-level data — access-gated  aggregate results — public
              └─────────────────┬─────────────────┘
                                ▼
               motrpac-human-presuspension-acute
               manuscript figure code + QC vignettes
                                ▼
                          manuscripts
```

| Repository | What it holds | Access |
|---|---|---|
| [`motrpac-human-presuspension-repro`](https://github.com/MoTrPAC/motrpac-human-presuspension-repro) | the end-to-end rebuild pipeline and its pinned software environment | code; a full run needs consortium bucket access |
| [`MotrpacHumanPreSuspensionData`](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionData) | subject-level molecular and phenotypic data objects | formal data-access request to the consortium |
| [`MotrpacHumanPreSuspensionAnalysis`](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis) | differential analysis, group summary statistics, enrichment, clustering, feature-to-gene map, plotting functions | public |
| [`motrpac-human-presuspension-acute`](https://github.com/MoTrPAC/motrpac-human-presuspension-acute) | per-manuscript figure code and QC vignettes | code public; some panels need Data access |

## The DAG
```
00_preflight ─▶ 10_build_data ─▶ 20_upload_bucket ─▶ 30_update_relevant_packages
```
The `Makefile` is the graph; each stage is a shell script in `scripts/` that delegates to the R code
already living in the two repos. `config/pipeline.env` is the single source of truth (versions,
buckets, repo paths). 

There is one run mode: every node runs. Consortium access to the gated bucket is a hard requirement
and preflight fails without it.

That graph is also drawn, at a finer grain than the four boxes above. `make depgraph` parses the
Makefile, the step drivers, the R stems and `config/file_versions.json` and writes an interactive map
to [`docs/dependency_graph/`](docs/dependency_graph/README.md): every make target, build step,
generator stem, freeze file and external input, wired by what produces and consumes what and sized by
blast radius. It reads source text and config only — no pipeline run, no `staging/` read, no
consortium access — so it works on a fresh clone, and it derives the stage list from the Makefile
rather than a hardcoded copy, so it follows a renamed stage on its own. `docs/dependency_graph/app.R`
is a Shiny explorer over the same tables: pick a node and trace what it depends on, or what breaks if
it changes.

## Stages
| Stage | Script | Does |
|-------|--------|------|
| 0 | `scripts/00_preflight/00_preflight.sh` | tools, gcloud auth, sibling repos, bucket access, external assets; fails if the consortium-gated bucket is unreadable. **Fully implemented.** Ships a `data-raw/` that builds the leaf (`deps=-`) data objects from vendored sources (`make data-objects`). |
| 1 | `scripts/10_build_data/10_build_data.sh` | regenerate `.rda` objects step by step — runs each numbered subfolder's `build.sh` in order (`01_stage_clinical` … `16_build_scion`); per-step PASS/SKIP/FAIL report; `STEPS=`, `STOP_ON_FAIL` knobs. `docs/data_objects.tsv` is the object inventory (no longer drives the build). Adapted steps: molecular signatures, set-to-id, qc-norm, feature-to-gene, `*_QC` objects, DA assembly, summary stats, `CAMERA_RESULTS`, `FCM_*`, `UTORONTO_TFs`, `PTMSEA_INPUT`, SCION networks (step 16, not a data object — on by default, `RUN_SCION=FALSE` opts out). |
| 2 | `scripts/20_upload_bucket/20_upload_bucket.sh` | bring the staging bucket in line with the Stage 1 freeze, step by step — runs each numbered subfolder's `build.sh` in order (`01_snapshot_buckets` … `05_validate_structure`); same PASS/SKIP/FAIL report and `STEPS=`/`STOP_ON_FAIL` knobs as Stage 1. Ports `google_cloud_bucket_checks/` in, so the buckets and versions are read from `config/pipeline.env` alone. **Writes to GCS only under `APPLY=1` (upload) or `SEED_STAGING=1` (seed); a bare `make upload` snapshots, diffs and plans.** |
| 3 | `scripts/30_update_relevant_packages/30_update_relevant_packages.sh` | carry the regenerated data into `MotrpacHumanPreSuspension{Data,Analysis}`. **Audit implemented** — validates both package trees, checks the dependency graph for cycles, inventories code no longer needed, and fails the stage on a hard cycle or a structural ERROR (`STRICT=0` to downgrade, `AUDIT_ONLY=1` to stop there). Reports land in `logs/package_audit/`. **Carry implemented** — routes every built object, assembles a **test package** per repo under `staging/test-packages/`, re-documents it, bumps version and NEWS, and runs its test suite. Neither package checkout is written to; promoting a test package into its repo stays manual. |

## Access, and what is public

What specifically needs MoTrPAC consortium credentials: regenerating qc-norm and individual-level DA,
epigenomics, clinical chemistry — anything that reads the private GCS buckets. That is why preflight
refuses an unreadable gated bucket outright rather than producing a partial result. Request access
from the consortium.

Individual-level molecular and phenotypic data are not distributed. What *is* publicly released is
the aggregate layer — group and timepoint summary statistics (`*_SUM_STATS`), enrichment and
clustering built from summary DA (`CAMERA_RESULTS`, `FCM_*`), and the acute (`ADU_BAS`) analyses —
shipped through the `MotrpacHumanPreSuspensionAnalysis` package rather than rebuilt by this pipeline.

## Usage
```bash
make all           # full chain from a fresh clone (needs consortium access + gcloud auth):
                   #   preflight -> sources -> data-objects -> data -> packages + upload plan
make preflight     # tools, auth and bucket access only
make sources       # fetch the staged inputs git can't hold (skips what is present)
make env           # record the run environment (see below)
make env-diff      # how does this machine differ from that record?
make depgraph      # redraw the pipeline dependency graph -> docs/dependency_graph/
make slurm         # run the whole chain on SCG as batch jobs (see below)
make docker-build  # containerized environment instead (see below)
make help

# `make all` runs every step: ATAC is refit, SCION is inferred, and networks already on
# disk are re-inferred. Opt out of the expensive ones rather than into them.
# SCION is the hours-to-days step — submit it (make slurm, or EXECUTOR=slurm) rather than
# run it inline. SCION_PERMUTATIONS stays 0: edges are cut by a stricter weight cutoff,
# because the networks are too dense for a permutation null to be the boundary.
RERUN_ATAC=FALSE RUN_SCION=FALSE make all   # the opt-outs: stage ATAC, skip SCION
```

## What is not in git

The repo carries its own vendored inputs, so a clone is self-contained except for seven files that
exceed GitHub's 100 MB per-file limit: the three MethylCap beta-value matrices, the three MethylCap
DA tables and the Ensembl v105 TxDb. `RERUN_ATAC=FALSE` stages ATAC from the release rather than
refitting it and adds four more — the ATAC qc-norm and DA tables for `t05-pbmc` and
`t06-muscle`. Stage 1 will not start without them —
`10_build_data/check_required_inputs.sh` fails first with the list.

One command supplies all seven:

```bash
make sources                 # download everything that is missing
make sources DRY_RUN=1       # preview: what would transfer
make sources FORCE=1         # re-fetch even when present
make sources GROUPS=methylcap_da
```

That is `config/copy_from_source.sh`, out of band like `make env` — no stage depends on it. The six
MethylCap files come from the bucket's `epigenomics/` tiers; the TxDb comes from `resources/`, which
step 06's `stage_ensembl_txdb` stem publishes into the freeze so Stage 2 uploads it with everything
else. Downloading it beats rebuilding it, so the Ensembl builder
(`build_ensembl_v105_cache.R`) is only a fallback for when `resources/` does not carry it yet — and it
is the one path that needs a live Ensembl connection.

The RefMet/KEGG snapshot under `00_preflight/data-raw/sources/refmet/` plays the same offline role
for the metabolomics stem — it replaces a live Metabolomics Workbench POST and a live KEGG fetch that
otherwise moved `refmet_name`/`refmet_id`/`kegg_id` with the databases — but it is small enough to
live in git, so it has no `make sources` group and a clone already has it. Step 06's
`stage_refmet_map` stem publishes it to `resources/` alongside the TxDb.

Everything is skipped when already present at the right size, so a second `make sources` costs
seconds and an interrupted transfer resumes rather than restarting. `GROUPS=` are
`methylcap_qc_norm`, `methylcap_da`, `atac_qc_norm`, `atac_da` and `ensembl_v105`. The two ATAC
groups are fetched only while `RERUN_ATAC=FALSE` — under the default `TRUE` the pipeline fits those
tables itself and needs neither.

Also absent from git, by design: `staging/` and `logs/`, and the built data objects under
`scripts/**/data/` and `scripts/00_preflight/output/` — those are what the pipeline produces, so
committing them would let a stale copy pass for current output.

## Running on SLURM (SCG)

The same DAG, submitted to a cluster instead of run here. `config/slurm.json` holds the
account, the partitions and the per-job resources; nothing else changes — same stages,
same steps, same `PASS/SKIP/FAIL` report, same log paths, same exit codes, so a cluster
run and a laptop run produce comparable output.

```bash
make slurm                    # one job per STAGE, chained by --dependency=afterok, returns
make slurm DRY_RUN=1          # print the generated batch scripts, submit nothing
make slurm-status             # squeue + sacct for the chain last submitted
make slurm-cancel             # scancel it
make data EXECUTOR=slurm      # the other granularity: one job per STEP, sized per step,
                              # blocking until the stage finishes
```

`make slurm` is for the overnight run — it queues four dependent jobs and hands your
shell back, and `upload` and `update-relevant-packages` run concurrently once `data`
succeeds. `EXECUTOR=slurm` on a stage target is for sizing or re-running one expensive
step, since each step then gets its own allocation.

Set `defaults.account` in `config/slurm.json` before the first run — preflight fails
while it reads `CHANGE_ME` — and set `environment.modules` to whatever SCG needs for R
and the Cloud SDK, since a batch job has no login profile behind it.
[`docs/SLURM.md`](docs/SLURM.md) covers the config, the two granularities and how job
exit codes are recovered.

## Running in Docker

The container is a snapshot of this repo's R environment: R 4.4.1, the Google Cloud SDK, and the
~436-package dependency closure **pinned to the exact versions installed on the machine that
generated `docs/ENVIRONMENT.md`**. It is built *from* `docs/environment/package_versions.tsv`, the
same file `make env` writes, so the image and the provenance record cannot drift apart — re-run
`make env`, rebuild, and the image follows.

```bash
make docker-build                                # once; see "Build time" below
docker compose run --rm precovid make preflight  # start here
make docker-shell                                # or poke around interactively
make docker-verify                               # what actually landed vs. the reference
```

### Scope

**Only this repo is mounted.** The parent directory is not, so no other checkout on your machine is
visible to the container:

| Mounted | At | Why |
|---|---|---|
| this repo | `/github/motrpac-human-presuspension-repro` | the working directory |
| `~/.config/gcloud` | `/root/.config/gcloud` | consortium bucket access; run `gcloud auth login` **on the host** and the container reuses those credentials |
| `precovid-local-lib` volume | `/local-lib` | the only writable R library, so anything installed at run time survives between runs |

**The image is the software environment, not the data.** It neither contains nor substitutes for
consortium access, and preflight still fails without it. Stages 1–3 delegate to R code in the
sibling repos (`MotrpacHumanPreSuspension{Data,Analysis}`). Each of those is a commented-out mount in
`docker-compose.yml` — uncomment individually, and only for a run that needs it. When
`MotrpacHumanPreSuspension{Data,Analysis}` are mounted, the entrypoint installs them into
`/local-lib` automatically and reinstalls only when their version or commit changes.

So with the default mounts, `make docker-preflight` passes the tooling and R-package checks
and **fails on `repo:*`** — that is the scope boundary doing its job, not a broken image. What it
proves is the part the container is responsible for: every tool on `PATH`, every R package
importable. Read `logs/preflight_report.tsv` row by row rather than the
exit code.

### How the versions are pinned

Each install stage runs twice. A bulk pass installs the whole set from one repository so the
resolver can sort out dependency order; a pin pass then moves any package whose version does not
match the record onto the recorded one, with `dependencies = FALSE`. Pinning only what actually
drifted is what keeps this tractable — pinning all 436 in dependency order is not.

- **CRAN** — `remotes::install_version()` against the recorded version. (`make env` normalizes every
  version separator to `.`, so CRAN's `1.4-8` is recorded as `1.4.8`; both `remotes` and the
  verification script compare as `package_version` objects, which treats the two as equal.)
- **Bioconductor** — your library spans several releases at once, so no single release reproduces
  it. The installer searches Bioconductor 3.19 through 3.23 for whichever one publishes the exact
  recorded version and installs from there.
- **GitHub** — the eight packages installed from a commit are pinned to that commit. For these the
  SHA is the only version information that exists, so it is the one pin that cannot be
  reconstructed from a version number.

Anything that cannot be resolved is reported as drift rather than silently substituted — see
`docs/DOCKER_DIFF.md` and `/opt/reference/build_report.tsv` inside the image.

**One thing worth knowing about the reference library:** it works, but it cannot be rebuilt from
scratch in the obvious order. `ggtree 3.12.0` resolves `ggplot2:::check_linewidth` while
byte-compiling, and `ggplot2 4.0` removed that function — so installing the recorded `ggtree` on top
of the recorded `ggplot2 4.0.2` fails, taking `enrichplot`, `ChIPseeker` and `clusterProfiler` down
with it. On the reference machine they work only because they were compiled back when `ggplot2` was
still 3.x and the upgrade happened around them afterwards. The installer reproduces that history
deliberately: it holds `ggplot2` at 3.5.2 for the Bioconductor stage (the one release with both
`check_linewidth` and the `is_ggplot` export the recorded `patchwork` needs), moves the stragglers'
CRAN dependencies to their recorded versions, retries them, then restores `ggplot2` to the recorded
4.0.2. The end state matches the
reference package for package. `HOLDBACKS` in `docker/install_packages.R` is where such cases are
recorded, each with the reason.

### Build time and architecture

The Dockerfile is architecture-agnostic. Pinning exact versions means most packages are fetched from
the CRAN archive and compiled from source, which costs the binary advantage x86 would otherwise have.
A measured `linux/arm64` build on a 12-core Apple-silicon machine takes about **70 minutes** (32 min
Bioconductor, 22 min CRAN, 14 min GitHub); a slower or more constrained host takes correspondingly
longer.

- **Native (default).** On Apple silicon this builds `linux/arm64`, which runs at full speed
  afterwards. Posit's package manager publishes binaries for x86 only, so every package compiles.
- **`make docker-build PLATFORM=linux/amd64`.** Binaries apply only to the minority of packages
  already at the pinned version in the snapshot, so the build is modestly faster but the result runs
  under emulation. Turn on Docker Desktop's Rosetta setting first; without it, emulation is 5–10×
  slower rather than ~1.5×.

Build once, either way — the layers are cached, and a rebuild after `make env` re-runs only the R
install stages.

**Check the image, not the build log.** `make docker-verify` writes `docs/DOCKER_DIFF.md` from the
image's own library; a good build reports **0 differ**, with the only `missing` rows being the two
packages installed from mounted repos at run time. The build log is not a substitute: BuildKit
truncates a long step's output, and a transient network failure during a pin pass leaves packages at
the wrong version without failing the build.

If the GitHub stage hits API rate limits, export a token and build with it as a secret:
`docker build --secret id=github_pat,env=GITHUB_PAT -f docker/Dockerfile .`

### If you would rather not use Docker

[`docs/ENVIRONMENT.md`](docs/ENVIRONMENT.md) records the machine, date and exact software the released
results were generated with: OS, R and Bioconductor releases, every tool version, the commit of each
repository, and the version and install source of every R package in the dependency closure (also as
a table in [`docs/environment/package_versions.tsv`](docs/environment/package_versions.tsv)). Install
R 4.4.1 and match it, then **capture your own environment and report how it differs**:

```bash
make env        # record this machine -> docs/ENVIRONMENT.md
make env-diff   # compare it against the reference -> docs/ENV_DIFF.md
```

`make env-diff` writes `docs/ENV_DIFF.md` and `docs/environment/env_diff.tsv`, classifying every
package as `match`, `differs`, `missing` or `extra`. Include that report when you report results
that diverge — a `differs` row on a directly-declared package is the first thing worth checking.
It is the same script the container runs (`make docker-verify`), so the two reports are comparable.


## Layout
```
config/     pipeline.env (source of truth) — values every stage reads
            file_versions.json — the version each file will carry in the next staging
            bucket upload; hand-maintained, NOT a mirror of the current bucket. A
            [tissue_code][ome][data_category][data_details] tree — the load_qc() nesting
            one level deeper — whose path IS the filename, with a why_bumped note per
            file. write_with_path_name() resolves each file's _v<ver> from it, so
            bumping a file for a release is an edit here
            capture_environment.{sh,R}  provenance capture, run on demand via `make env`
            slurm.json — SCG account, partitions and per-job resources; read only when
            a run submits to SLURM (docs/SLURM.md)
docker/     STATUS.md  what is built, what is verified, how to resume (start here)
            Dockerfile + install_packages.R (builds the image FROM the capture above)
            verify_environment.R  reference-vs-current diff, used by both `make docker-verify`
                                  and `make env-diff`
            entrypoint.sh, install_local_pkgs.sh, config.json  run-time wiring
docker-compose.yml  mounts (this repo, gcloud creds, writable library) + env
scripts/    one subfolder per stage (00_preflight/, 10_build_data/, …) + lib/common.sh
            run_stage.sh    what each stage runs, defined once — used by the Makefile
                            and by SLURM batch jobs alike
            submit_slurm.sh the chained stage submission (`make slurm`) + status/cancel
            lib/slurm.sh    run_step: the one place local-vs-cluster is decided
            00_preflight/config/    external_assets.tsv + gmt_files.tsv (Stage 0 manifests)
            00_preflight/data-raw/  builders for leaf data objects (+ vendored raw sources)
            00_preflight/data/      built .rds data objects only (R-package standard)
            00_preflight/output/    gmt.gz + quantid_bucket_files.json (quant-id catalog)
docs/       data_objects.tsv (object inventory)
            dependency_graph/  interactive pipeline map + node/edge tables (`make depgraph`)
            ENVIRONMENT.md + environment/  generated run-provenance record (`make env`)
logs/       per-stage reports, per-object logs
```

## Status
All four stages are implemented; the table at the top of this file says what each one does. What
none of them do is write outside this repo without being told to: Stage 2 uploads to GCS only under
`APPLY=1`, Stage 3 builds test packages under `staging/` and never touches either package checkout,
and `make promote` is manual and never runs as part of a stage.

Three Stage 1 steps stage their object rather than rebuilding it — `01_stage_clinical`,
`02_stage_pheno` and `05_stage_metabolomics_cvs` — so their output is carried from upstream, not
reproduced here. Each says so in its own `build.sh`.

## License
The pipeline code and documentation are MIT licensed (see `LICENSE`). Vendored third-party and
consortium files under `sources/` and `prebuilt/` keep their original terms; each source directory's
README records them.
