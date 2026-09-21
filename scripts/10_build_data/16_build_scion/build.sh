#!/usr/bin/env bash
# Step 16 — SCION regulatory networks: prot-ph TF regulators over DE muscle transcript
# targets, and blood metabolite regulators over the same targets, per exercise group.
#
# ON by default, like every step: `make all` runs the whole graph and a step is opted out
# of rather than into. The output is not a data object and nothing downstream in this
# pipeline reads it: Figure 6B is a Cytoscape layout of a hand-merged export of these
# tables, and ED8A's TFEB target list comes out of that same session. This step regenerates
# what those manual steps started from. It is also the one step here measured in hours to
# days — a random-forest fit per cluster per group.
#
# Knobs (config/pipeline.env):
#   RUN_SCION=TRUE            run it at all; FALSE opts out and the step exits 77
#   SCION_PERMUTATIONS=0      null replicates after the observed network (0 = observed only;
#                             edges are cut by a weight cutoff, not a permutation null)
#   SCION_CORES=11            cores handed to run_SCION()
#   SCION_FORCE=TRUE          re-infer networks already on disk; FALSE resumes
#
# Submit it rather than running it inline: EXECUTOR=slurm keys this step in config/slurm.json.
#
# Requires steps 06 (prot-ph imputed freeze matrix), 07, 08, 10 and 14.
# Output: staging/scion/<network>/<group>/ and staging/scion/scion_manifest.tsv
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../lib/common.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ "$(printf '%s' "${RUN_SCION:-TRUE}" | tr '[:lower:]' '[:upper:]')" != "TRUE" ]]; then
  log "RUN_SCION=${RUN_SCION:-TRUE} — opted out; unset RUN_SCION to infer the networks"
  exit 77
fi

"${RSCRIPT}" "${HERE}/SCION_NETWORKS.R"
if [[ "${RUN_TESTS:-1}" == "1" ]]; then
  "${RSCRIPT}" "${HERE}/scion_tests.R"
fi
