#!/usr/bin/env Rscript
# Step 05 side output: the METABOLOMICS_CVS table as a freeze resources .txt.
#
# The bucket carries this table under resources/ alongside the feature-to-gene map, so
# the freeze needs it there too — otherwise resources/ is short one file relative to the
# reference and Stage 2 has nothing local to diff against. Reads the .rda that build.sh
# has just staged into data/; it does not rebuild the object (step 05 is still a stage,
# not a local build — see build.sh's TODO).
#
# Format follows gs://pre-cawg/staging_20260720/resources/, i.e. the 9 columns the object
# already carries, in the object's own order:
#   tissue_assay_sites, tissue_code, assay, site, feature_id, feature_cv, refmet_name,
#   tissue, lowest_CV
# Row order is the object's, which is the bucket's.
#
# Written UNQUOTED, matching write_with_path_name()'s convention for every other freeze
# file.

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
# freeze_version(): per-file version from config/file_versions.json
source(file.path(ROOT, "scripts", "10_build_data", "lib", "file_versions.R"))

rda <- file.path(ROOT, "scripts", "10_build_data", "data", "METABOLOMICS_CVS.rda")
if (!file.exists(rda))
  stop("missing ", rda, " — build.sh stages it from the Analysis package before this runs")
e <- new.env()
load(rda, envir = e)
METABOLOMICS_CVS <- e$METABOLOMICS_CVS

resources_dir <- file.path(ROOT, "staging", "freeze", "resources")
dir.create(resources_dir, showWarnings = FALSE, recursive = TRUE)

stem <- "motrpac_human-precovid_metabolite-cv"
out  <- file.path(resources_dir, paste0(stem, "_v", freeze_version(stem), ".txt"))
utils::write.table(as.data.frame(METABOLOMICS_CVS), file = out,
                   row.names = FALSE, sep = "\t", quote = FALSE, na = "NA")
message("METABOLOMICS_CVS: ", nrow(METABOLOMICS_CVS), " rows -> ", out)
