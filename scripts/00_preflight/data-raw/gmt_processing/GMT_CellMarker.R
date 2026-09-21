#!/usr/bin/env Rscript
# Leaf data object: cellmarker.v2024.symbols.gmt.gz  (deps=-, non-gated)
# Vendored + adapted from .../data-raw/gmt_processing/GMT_CellMarker.R
# Source downloaded 2025-01-23 from https://maayanlab.cloud/Enrichr/#libraries
suppressWarnings(suppressMessages(library(dplyr)))
source(file.path(Sys.getenv("PREFLIGHT_DR"), "lib", "helpers.R"))

src <- file.path(.dataraw_dir(), "gmt_processing", "sources", "CellMarker_2024.txt.gz")
cellmarker <- TMSig::readGMT(src, check = FALSE)
cellmarker <- cellmarker[grepl("human", names(cellmarker), ignore.case = TRUE)]
names(cellmarker) <- paste0("CELLMARKER_", names(cellmarker))

out <- .gmt_out(cellmarker, "cellmarker.v2024.symbols.gmt")
message("CellMarker: ", length(cellmarker), " sets -> ", out)
