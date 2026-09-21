#!/usr/bin/env Rscript
# Leaf data object: phosphositeplus.v6.7.1.1.flanking.gmt.gz  (deps=-, non-gated)
# Vendored + adapted from .../data-raw/gmt_processing/GMT_PSP_kinase.R
# Changes: local source path; vendored .writeGMT; fixed gzip(file) bug.
# NOTICE: PhosphoSitePlus data is for non-commercial use only (v6.7.1.1).
suppressWarnings(suppressMessages(library(dplyr)))
source(file.path(Sys.getenv("PREFLIGHT_DR"), "lib", "helpers.R"))

src <- file.path(.dataraw_dir(), "gmt_processing", "sources", "Kinase_Substrate_Dataset.gz")
kinase_sets <- src %>%
  gzfile() %>%
  read.csv(skip = 2, sep = "\t") %>%
  filter(KIN_ORGANISM == "human", KIN_ORGANISM == SUB_ORGANISM) %>%
  mutate(human_flanking = toupper(`SITE_...7_AA`),
         human_flanking = gsub("_", "-", human_flanking),
         human_flanking = sub("(^.{7})(.{1})(.*$)", "\\1\\L\\2\\E\\3",
                              human_flanking, perl = TRUE)) %>%
  select(human_flanking, kinase = GENE) %>%
  distinct() %>%
  unstack()

names(kinase_sets) <- paste0("PSP_", names(kinase_sets))

out <- .gmt_out(kinase_sets, "phosphositeplus.v6.7.1.1.flanking.gmt")
message("PSP kinase: ", length(kinase_sets), " sets -> ", out)
