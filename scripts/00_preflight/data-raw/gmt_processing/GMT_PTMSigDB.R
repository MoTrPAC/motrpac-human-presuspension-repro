#!/usr/bin/env Rscript
# Leaf data object: ptmsigdb.v2.0.flanking.gmt.gz  (deps=-, non-gated)
# Vendored + adapted from .../data-raw/gmt_processing/GMT_PTMSigDB.R
# Changes: local source path; vendored .writeGMT; fixed the .writeGMT(path=file)
# and gzip(file) bugs (undefined `file`). Raw source is the PTMsigDB v2.0 flanking
# GMT (ptm.sig.db.all.flanking.human.v2.0.0.gmt.gz).
suppressWarnings(suppressMessages({library(dplyr); library(TMSig)}))
source(file.path(Sys.getenv("PREFLIGHT_DR"), "lib", "helpers.R"))

src <- file.path(.dataraw_dir(), "gmt_processing", "sources",
                 "ptm.sig.db.all.flanking.human.v2.0.0.gmt.gz")
ptmsigdb <- readGMT(src) %>%
  lapply(function(xi) {
    xi <- xi[grepl("-p;[ud]$", xi)]                       # phosphosites only
    xi <- sub("-p", "", xi)                               # drop phospho indicator
    xi <- gsub("_", "-", xi)
    xi <- sub("(^.{7})(.{1})(.*$)", "\\1\\L\\2\\E\\3", xi, perl = TRUE)
    xi
  })
ptmsigdb <- ptmsigdb[lengths(ptmsigdb) > 0]
names(ptmsigdb) <- paste0("PTMSIGDB_", names(ptmsigdb))

out <- .gmt_out(ptmsigdb, "ptmsigdb.v2.0.flanking.gmt")
message("PTMSigDB: ", length(ptmsigdb), " sets -> ", out)
