#!/usr/bin/env Rscript
# Leaf data object: SPLICING_DA  (deps=-)
# Vendored + adapted from .../data-raw/SPLICING_DA.R
# Reads differential_splicing.rds (differential alternative-splicing results,
# FDR < 0.05) supplied by Zidong Zhang (zidong.zhang@mssm.edu). The analytical
# methods (alignment, isoform quantification, splicing detection, testing) are
# described in "Exercise modulation of the alternative splicing landscape in
# human tissues" and are not rederived here.
#
# Unlike the other leaf objects this raw file is NOT downloadable from a bucket;
# it is supplied by Zidong Zhang (MSSM) and vendored at
# sources/differential_splicing.rds.
source(file.path(Sys.getenv("PREFLIGHT_DR"), "lib", "helpers.R"))

src <- file.path(.dataraw_dir(), "sources", "differential_splicing.rds")
if (!file.exists(src))
  stop("missing source: ", src,
       "\n  Obtain differential_splicing.rds from Zidong Zhang (MSSM) and place it there.")

SPLICING_DA <- readRDS(src)
saveRDS(SPLICING_DA, file.path(.out_dir(), "SPLICING_DA.rds"))
message("SPLICING_DA: ", length(SPLICING_DA), " elements -> ", .out_dir())
