#!/usr/bin/env Rscript
# Leaf data object: COVARIATES_FILE  (deps=-, non-gated)
# Vendored + adapted from .../data-raw/COVARIATES_FILE.R
# Change: read the vendored sources/covariates_pre_cawg.csv instead of resolving
# precovid_repo_path from ~/config.json; saveRDS instead of usethis::use_data.
source(file.path(Sys.getenv("PREFLIGHT_DR"), "lib", "helpers.R"))

src <- file.path(.dataraw_dir(), "sources", "covariates_pre_cawg.csv")
if (!file.exists(src))
  stop("missing source: ", src, " (copy library/covariates_pre_cawg.csv here)")

COVARIATES_FILE <- read.csv(src)
saveRDS(COVARIATES_FILE, file.path(.out_dir(), "COVARIATES_FILE.rds"))
message("COVARIATES_FILE: ", nrow(COVARIATES_FILE), " rows -> ", .out_dir())
