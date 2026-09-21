#!/usr/bin/env Rscript
# Stage 1 QC-norm stem: the Ensembl v105 TxDb as a freeze resource (staging only —
# nothing is computed here).
#
# The TxDb is 175 MB, past GitHub's per-file limit, so unlike the rest of sources/ it
# cannot be vendored in git. Publishing it under resources/ makes the bucket its
# distribution channel: config/copy_from_source.sh downloads it from there on a fresh
# clone instead of rebuilding it from Ensembl, which is slow and flaky.
#
# Copied verbatim, and versioned from v2.0 on. The published name carries a _v<version>
# suffix like every other file in the collection; the local cache under sources/ keeps
# its bare name, since that copy is an input to this repo rather than a release
# artifact. The version is a release stamp only — the Ensembl release this TxDb is built
# from is in the filename already, and a rebuild of the same Ensembl release is not a
# content change, so the entry in config/file_versions.json stays put until it is.
#
# Output: staging/freeze/resources/txdb_hsapiens_ensembl_v105_v<version>.sqlite
ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}

# freeze_version(): per-file version from config/file_versions.json. Sourced directly
# rather than via lib/qc_helpers.R, which this staging-only stem has no other use for.
source(file.path(ROOT, "scripts", "10_build_data", "lib", "file_versions.R"))

txdb <- file.path(ROOT, "scripts", "00_preflight", "data-raw", "sources", "ensembl_v105",
                  "txdb_hsapiens_ensembl_v105.sqlite")
if (!file.exists(txdb))
  stop("missing ", txdb, " — run `make sources` to fetch or build the Ensembl v105 cache")

resources_dir <- file.path(ROOT, "staging", "freeze", "resources")
dir.create(resources_dir, showWarnings = FALSE, recursive = TRUE)
stem <- sub("[.]sqlite$", "", basename(txdb))
out  <- file.path(resources_dir, paste0(stem, "_v", freeze_version(stem, "sqlite"), ".sqlite"))

# Skip the 175 MB copy when the freeze already holds an identical file.
if (file.exists(out) && file.size(out) == file.size(txdb)) {
  message("ensembl v105 TxDb: already staged -> ", out)
} else {
  if (!file.copy(txdb, out, overwrite = TRUE))
    stop("failed to copy ", txdb, " -> ", out)
  message("ensembl v105 TxDb: staged ", format(file.size(out), big.mark = ","), " bytes -> ", out)
}
