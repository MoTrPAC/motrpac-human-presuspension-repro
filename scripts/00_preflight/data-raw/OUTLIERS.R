#!/usr/bin/env Rscript
# Leaf data object: OUTLIERS  (deps=-)
# Vendored + adapted from .../data-raw/OUTLIERS.R
# Changes: reads the removed-samples files vendored under
# sources/removed_samples/ (downloaded from the production bucket) instead
# of a live gsutil fetch, so it builds offline; tissue/ome are derived from the
# filename via the vendored .find_tissue/.find_ome; saveRDS not usethis::use_data.
# NOTE: the raw files hold vial labels (consortium-gated); see the folder README.
suppressWarnings(suppressMessages(library(dplyr)))
source(file.path(Sys.getenv("PREFLIGHT_DR"), "lib", "helpers.R"))

src_dir <- file.path(.dataraw_dir(), "sources", "removed_samples")
files <- list.files(src_dir, pattern = "removed-samples.*\\.txt$", full.names = TRUE)
if (length(files) == 0) stop("no removed-samples files in ", src_dir)

outlier_list <- lapply(files, function(path) {
  x <- read.csv(path, sep = "\t", check.names = FALSE)
  tissue <- .find_tissue(basename(path)); ome <- .find_ome(basename(path))
  if ("sample" %in% names(x)) x <- dplyr::rename(x, vialLabel = sample)
  if ("Sample" %in% names(x)) x <- dplyr::rename(x, vialLabel = Sample)
  if ("PC" %in% names(x))     x <- dplyr::rename(x, reason = PC)
  x %>% dplyr::select(vialLabel, reason) %>% dplyr::mutate(tissue = tissue, ome = ome)
})
OUTLIERS <- do.call(rbind, outlier_list)
# Manually identified outliers (blood rna-seq 11263010401, the six ATAC mix-ups) are rows in
# the removed-samples files themselves, appended by append_manual_outliers.R, so the object
# and the file the freeze ships list the same samples. Nothing is added here.

saveRDS(OUTLIERS, file.path(.out_dir(), "OUTLIERS.rds"))
message("OUTLIERS: ", nrow(OUTLIERS), " rows (from ", length(files), " vendored files) -> ", .out_dir())
