#!/usr/bin/env Rscript
# Stage 1 QC-norm stem: removed-samples metadata (staging only — nothing is computed here).
#
# The removed-samples files are release artifacts in their own right: the bucket-checks
# validator expects a metadata/removed-samples entry per ome x tissue alongside qc-norm and
# the other metadata (see data-raw/google_cloud_bucket_checks/required_structure.R, where the
# row exists but is marked required = FALSE). No stem generates them — they are produced
# upstream, downloaded from the production bucket, and vendored under
# scripts/00_preflight/data-raw/sources/removed_samples/, where OUTLIERS.R reads them to build
# the OUTLIERS object.
#
# This stem is what puts them into the freeze, so the rest of the pipeline and the freeze
# tests see the same file set the release does. Their filenames already carry the full BIC
# convention, so they are copied verbatim — the same treatment the vendored MethylCap
# beta-value matrices get in generate_methylcap_qc_norm.R.
#
# It lives here rather than in OUTLIERS.R so that .freeze_subdir_for_ome() can be reused from
# lib/qc_helpers.R: preflight cannot source that file, because qc_helpers.R loads OUTLIERS.rds
# at source time and OUTLIERS.R is what builds it.
#
# Output: staging/freeze/<ome-group>/metadata/
#   human-precovid-sed-adu_<tissue_code>_<ome>_metadata_removed-samples_v<x.y>.txt
suppressWarnings(suppressMessages({ library(dplyr) }))
ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
source(file.path(ROOT, "scripts", "10_build_data", "lib", "qc_helpers.R"))

stage_removed_samples = function(repo_local_dir){
  src_dir = file.path(.QC_REPO, "scripts", "00_preflight", "data-raw", "sources", "removed_samples")
  files = list.files(src_dir, pattern = "removed-samples.*\\.txt$", full.names = TRUE)
  if (length(files) == 0)
    stop("no removed-samples files in ", src_dir)

  copied = 0L; unmapped = character(0)
  for (path in files) {
    # human-precovid-sed-adu_<tissue_code>_<ome>_metadata_removed-samples_v<x.y>.txt
    token = sub("_metadata_removed-samples.*$", "", sub("^human-precovid-sed-adu_", "", basename(path)))
    ome = sub("^[^_]+_", "", token)
    subdir = .freeze_subdir_for_ome(ome)
    if (is.na(subdir)) { unmapped = c(unmapped, basename(path)); next }

    dest_dir = file.path(repo_local_dir, "freeze", subdir, "metadata")
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    if (file.copy(path, file.path(dest_dir, basename(path)), overwrite = TRUE))
      copied = copied + 1L
  }

  if (length(unmapped) > 0)
    warning("removed-samples: no freeze subdir for ", length(unmapped), " file(s): ",
            paste(utils::head(unmapped, 3), collapse = ", "))
  message("removed-samples: staged ", copied, " of ", length(files), " file(s)")
  return(invisible(copied))
}

out_base <- .STAGING
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
stage_removed_samples(out_base)
message("removed-samples staged under ", file.path(out_base, "freeze"))
