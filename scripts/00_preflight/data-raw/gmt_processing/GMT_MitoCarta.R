#!/usr/bin/env Rscript
# Leaf data object: mitocarta3.0.symbols.gmt.gz  (deps=-, non-gated)
# Vendored + adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/gmt_processing/GMT_MitoCarta.R
# Changes: local source path; vendored .writeGMT; output via .gmt_out(); fixed
# the upstream gzip(file) bug (undefined `file`).
suppressWarnings(suppressMessages(library(dplyr)))
source(file.path(Sys.getenv("PREFLIGHT_DR"), "lib", "helpers.R"))

src <- file.path(.dataraw_dir(), "gmt_processing", "sources", "Human.MitoCarta3.0.xls")
mitocarta <- readxl::read_xls(src, sheet = "C MitoPathways") %>%
  select(MitoPathway, Genes) %>%
  filter(!is.na(MitoPathway)) %>%
  mutate(Genes = strsplit(Genes, split = ", ")) %>%
  {structure(.$Genes, names = .$MitoPathway)}

names(mitocarta) <- paste0("MITOCARTA_", names(mitocarta))

out <- .gmt_out(mitocarta, "mitocarta3.0.symbols.gmt")
message("MitoCarta: ", length(mitocarta), " sets -> ", out)
