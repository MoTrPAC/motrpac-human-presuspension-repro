# motrpac-human-presuspension-repro — the pipeline DAG.
# Each stage is a shell script under scripts/; this Makefile encodes the order
# preflight -> data -> upload -> update-relevant-packages and provides per-stage
# entry points.
#
#   make all                      # run the whole chain from a fresh clone:
#                                 #   preflight -> sources -> data-objects -> data -> packages + upload
#   make preflight                # Stage 0 only
#   make data                     # Stage 1 (implies preflight)
#   make upload                   # Stage 2 (implies data) — dry run; APPLY=1 to write to GCS
#   make update-relevant-packages # Stage 3 (implies data) — audit the two packages, plan the carry
#   make promote                  # gated: rsync staging -> production (manual, never automatic)
#   make slurm                    # submit the whole chain to SLURM and return (see docs/SLURM.md)
#   make slurm-status             # squeue/sacct for the chain last submitted
#   make slurm-cancel             # scancel that chain
#   make sources                  # on demand: fetch/build the staged inputs too big to keep in git
#   make env                      # on demand: record OS/date/package versions -> docs/ENVIRONMENT.md
#   make env-diff                 # on demand: how does this machine differ from that record?
#   make depgraph                 # on demand: redraw the pipeline dependency graph -> docs/dependency_graph/
#   make docker-*                 # build/run the containerized environment (see README)
#   make clean                    # remove stamps + logs
#
# Stamps under .stamps/ record that a stage completed, and gate `promote`.
#
# They do NOT make a re-run resume. The stage targets are .PHONY, and a recipe-less
# rule on a phony prerequisite is always out of date, so `$(STAMPS)/preflight:
# preflight` and `$(STAMPS)/data: data` resolve by running the stage every time —
# `make upload` runs preflight and data whether or not their stamps exist. Those two
# stamps carry the DAG and nothing else; to skip work inside a stage, use STEPS=.
#
# .stamps/upload is different, because nothing here knows how to build it. Only
# run_stage.sh writes it, on a successful upload, so `promote: $(STAMPS)/upload`
# cannot be satisfied any other way — without a completed upload in this checkout,
# make refuses with "No rule to make target `.stamps/upload'". That is the interlock
# on the one target that writes to production, and `make clean` re-arms it.
#
# Every stage target runs `scripts/run_stage.sh <stage>`, which is the single
# definition of what a stage does and which stamp it sets. The SLURM chain submitter
# runs the same script inside a batch job, so `make data` and a submitted data job
# are the same work.
#
# Two ways to use the cluster (docs/SLURM.md):
#   make data EXECUTOR=slurm   submit each STEP as its own job, block until done
#   make slurm                 submit each STAGE as one dependent job, return now

SHELL := /bin/bash
SCRIPTS := scripts
STAMPS  := .stamps

.PHONY: all precheck preflight sources env env-diff depgraph data-objects data upload \
        update-relevant-packages promote clean help \
        slurm slurm-status slurm-cancel \
        docker-build docker-shell docker-preflight docker-verify docker-reinstall

# Containerized environment. PLATFORM=linux/amd64 builds the x86 image, which uses
# precompiled P3M binaries instead of compiling from source.
COMPOSE   := docker compose
DC_RUN    := $(COMPOSE) run --rm precovid
REF_TSV   := docs/environment/package_versions.tsv

# Serial on purpose: make runs prerequisites left to right, and that order is the fresh-clone
# order. preflight fails fast on missing tools or bucket access; sources fetches the
# gitignored inputs; data-objects builds the leaf objects Stage 1's input gate requires.
# sources and data-objects stay out of band for the per-stage targets.
.NOTPARALLEL:
all: preflight sources data-objects update-relevant-packages upload ## fresh clone to upload: preflight, sources, data objects, data, packages, upload plan

$(STAMPS):
	@mkdir -p $(STAMPS)

precheck: ## warn if existing outputs would be overwritten (delete freeze folder to regen from scratch)
	@bash $(SCRIPTS)/check_existing_outputs.sh

# The stage targets carry the DAG and the stamps; run_stage.sh carries what each
# stage runs (precheck included) and touches its own stamp.
preflight: | $(STAMPS) ## Stage 0 — permissions + external-asset checks
	@bash $(SCRIPTS)/run_stage.sh preflight

$(STAMPS)/preflight: preflight

data-objects: ## Stage 0 — build the leaf (deps=-) data objects from vendored sources
	@bash $(SCRIPTS)/run_stage.sh data-objects

data: $(STAMPS)/preflight ## Stage 1 — regenerate data objects (EXECUTOR=slurm to submit each step)
	@bash $(SCRIPTS)/run_stage.sh data

$(STAMPS)/data: data

upload: $(STAMPS)/data ## Stage 2 — snapshot, diff and plan the bucket upload (APPLY=1 to write; no promote)
	@bash $(SCRIPTS)/run_stage.sh upload

update-relevant-packages: $(STAMPS)/data ## Stage 3 — audit the data + analysis packages, then plan the carry (writes nothing outside logs/)
	@bash $(SCRIPTS)/run_stage.sh update-relevant-packages

promote: $(STAMPS)/upload ## gated: promote staging -> production (manual)
	@echo "Promotion is manual. Review logs/structure_validation_*.tsv, then run the"
	@echo "rsync documented in google_cloud_bucket_checks/README.md after both package PRs merge."

# ---- SLURM ------------------------------------------------------------------
# Submit-and-return: one job per stage, ordered by --dependency=afterok, each stage
# running its steps inside its own allocation. Nothing here blocks, so the pipeline
# survives losing your shell. Resources and account live in config/slurm.json.
# The per-step alternative — `make data EXECUTOR=slurm` — needs no target of its own;
# it is the same `make data`, deciding differently at the one call site in run_step.

slurm: ## submit the whole chain to SLURM and return (DRY_RUN=1, STAGES="data upload")
	@bash $(SCRIPTS)/submit_slurm.sh $(if $(DRY_RUN),--dry-run) $(if $(STAGES),--stages "$(STAGES)")

slurm-status: ## squeue + sacct for the chain last submitted
	@bash $(SCRIPTS)/submit_slurm.sh status

slurm-cancel: ## scancel every job in the chain last submitted
	@bash $(SCRIPTS)/submit_slurm.sh cancel

# Out of band, like `env`: no stage depends on this. Run it after a fresh clone —
# the seven files it provides are gitignored for exceeding GitHub's 100 MB file
# limit. Skips whatever is already present, so re-running is cheap.
sources: ## fetch/build the staged inputs too big for git (DRY_RUN=1, FORCE=1, GROUPS=)
	@bash config/copy_from_source.sh $(if $(DRY_RUN),--dry-run) $(if $(FORCE),--force) $(GROUPS)

# Out of band: no stage depends on this and it depends on no stage. Run it when
# the environment changes or before committing regenerated outputs.
env: ## record OS, run date, tool + R package versions to docs/ENVIRONMENT.md
	@bash config/capture_environment.sh

env-diff: ## compare THIS machine's R library against docs/environment/package_versions.tsv
	@Rscript docker/verify_environment.R \
	  $(REF_TSV) docs/environment/env_diff.tsv docs/ENV_DIFF.md

# Out of band, like `env`. Parses source text and config only — it does not run
# the pipeline and does not read staging/, so it is safe on a fresh clone.
# Re-run it after changing the Makefile, a build.sh, a stem, or file_versions.json.
depgraph: ## redraw the pipeline dependency graph -> docs/dependency_graph/
	@Rscript docs/dependency_graph/build_depgraph.R

# ---- Containerized environment ---------------------------------------------
# Out of band, like `env`: no pipeline stage depends on these.

docker-build: ## build the image (PLATFORM=linux/amd64 for the fast binary build)
	@DOCKER_DEFAULT_PLATFORM=$(PLATFORM) $(COMPOSE) build $(if $(NO_CACHE),--no-cache)

docker-shell: ## interactive shell inside the container
	@$(DC_RUN) bash

docker-preflight: ## run Stage 0 inside the container
	@$(DC_RUN) make preflight

# SKIP_LOCAL_PKG_INSTALL: verification only inspects the image's own library, so
# there is no reason to spend minutes installing the mounted 600 MB data package.
docker-verify: ## write docs/environment/docker_diff.tsv + docs/DOCKER_DIFF.md from the image
	@$(COMPOSE) run --rm -e SKIP_LOCAL_PKG_INSTALL=1 precovid \
	  Rscript /usr/local/bin/verify_environment.R \
	  $(REF_TSV) docs/environment/docker_diff.tsv docs/DOCKER_DIFF.md

docker-reinstall: ## force-reinstall the mounted source packages (only if you enabled those mounts)
	@$(COMPOSE) run --rm -e FORCE=1 precovid true

clean: ## remove stamps and logs
	@rm -rf $(STAMPS) logs/*

help: ## list targets
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
	 awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-24s\033[0m %s\n",$$1,$$2}'
