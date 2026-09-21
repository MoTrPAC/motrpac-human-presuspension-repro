#!/usr/bin/env Rscript
# Stage 1 QC-norm stem: the frozen RefMet mapping as a freeze resource (staging
# only — nothing is computed here).
#
# generate_metab_qc_norm folds RefMet into each platform's metadata_features file,
# but only the columns the freeze needs (refmet_name, refmet_id, kegg_id) and only
# for features that survived QC. Publishing the snapshot itself under resources/
# gives the collection the whole mapping: every name the study asked RefMet about,
# what RefMet standardized it to, and the class hierarchy and cross-references
# (HMDB, ChEBI, LIPID MAPS, PubChem, KEGG) that come with it. That is also what
# makes the metabolomics annotation auditable — it is the exact table the freeze
# was built against, not whatever RefMet returns today.
#
# Copied verbatim, and versioned from v2.0 on — same treatment as
# resources/txdb_hsapiens_ensembl_v105.sqlite. The _v<ver> suffix is a release stamp,
# not a statement about the snapshot: the freeze version tracks the study release, so
# the snapshot's own build date still travels beside it in the provenance JSON, and
# re-fetching RefMet without a content change leaves the version where it is.
#
# Output: staging/freeze/resources/motrpac_human-precovid_refmet-map_v<version>.txt
#         staging/freeze/resources/motrpac_human-precovid_refmet-map_provenance_v<version>.json
ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}

# freeze_version(): per-file version from config/file_versions.json. Sourced directly
# rather than via lib/qc_helpers.R, which this staging-only stem has no other use for.
source(file.path(ROOT, "scripts", "10_build_data", "lib", "file_versions.R"))

refmet_dir <- file.path(ROOT, "scripts", "00_preflight", "data-raw", "sources", "refmet")

# source -> published stem and extension. The version is not written here: it comes from
# the resources entry in config/file_versions.json, keyed by the stem.
staged <- list(
  "refmet_name_map.txt"  = c(stem = "motrpac_human-precovid_refmet-map",            ext = "txt"),
  "refmet_snapshot.json" = c(stem = "motrpac_human-precovid_refmet-map_provenance", ext = "json")
)

missing <- names(staged)[!file.exists(file.path(refmet_dir, names(staged)))]
if (length(missing))
  stop("missing RefMet snapshot file(s): ", paste(missing, collapse = ", "),
       "\n  in ", refmet_dir,
       "\n  build them: Rscript ", file.path(refmet_dir, "build_refmet_cache.R"))

resources_dir <- file.path(ROOT, "staging", "freeze", "resources")
dir.create(resources_dir, showWarnings = FALSE, recursive = TRUE)

for (src_name in names(staged)) {
  stem <- staged[[src_name]][["stem"]]
  ext  <- staged[[src_name]][["ext"]]
  src <- file.path(refmet_dir, src_name)
  out <- file.path(resources_dir,
                   paste0(stem, "_v", freeze_version(stem, ext), ".", ext))
  if (!file.copy(src, out, overwrite = TRUE))
    stop("failed to copy ", src, " -> ", out)
  message("refmet snapshot: staged ", format(file.size(out), big.mark = ","), " bytes -> ", out)
}
