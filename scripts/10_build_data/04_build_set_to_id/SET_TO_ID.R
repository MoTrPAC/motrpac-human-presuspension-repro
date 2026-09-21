#!/usr/bin/env Rscript
# Stage 1 data object: SET_TO_ID  (deps: MOLECULAR_SIGNATURES)
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/SET_TO_ID.R.
#
# Local-.rda build: upstream did library(...) and used MOLECULAR_SIGNATURES by
# name from the INSTALLED package; here it is load()ed from the local stage-1
# MOLECULAR_SIGNATURES.rda instead. Output is save()d as a local .rda. The stray
# library(MotrpacHumanPreSuspensionData) is dropped — SET_TO_ID uses no Data object.
#
# NOTE: this data/SET_TO_ID.rda is staging — it will eventually be organized as
# the SET_TO_ID data object in the Analysis package.
suppressWarnings(suppressMessages(library(dplyr)))

# Root from the PRECOVID_ROOT config env (config/pipeline.env); fall back to an
# upward search for config/pipeline.env so standalone runs work.
ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
}
data_dir <- file.path(ROOT, "scripts", "10_build_data", "data")
molsig <- file.path(data_dir, "MOLECULAR_SIGNATURES.rda")
if (!file.exists(molsig)) stop("missing MOLECULAR_SIGNATURES.rda — build it first: ", molsig)
load(molsig)  # -> MOLECULAR_SIGNATURES

## Helper functions ----
# Split strings on `split`, then keep the leading components whose total length
# does not exceed `n`.
.cutstr <- function(x, split = "", n = Inf) {
  x <- strsplit(x, split = split)
  x <- vapply(x, function(xi) {
    keep <- cumsum(nchar(xi)) + nchar(split) * (seq_along(xi) - 1L) <= n
    xi <- paste(xi[keep], collapse = split)
    return(xi)
  }, character(1L))
  return(x)
}
.cutstr <- Vectorize(.cutstr)

# Number of digits used for each ID is a function of the number of sets
n_digits <- floor(log10(sum(lengths(MOLECULAR_SIGNATURES)))) + 1L

SET_TO_ID <- data.frame(
  database = rep(names(MOLECULAR_SIGNATURES), lengths(MOLECULAR_SIGNATURES)),
  set = unlist(lapply(MOLECULAR_SIGNATURES, names))
) %>%
  mutate(
    database = ifelse(grepl("^PTMSIGDB", database),
                      sub("(^PTMSIGDB_[^_]+).*", "\\1", set),
                      database),
    database = factor(database, levels = unique(database))
  ) %>%
  # PTMSigDB last, then CellMarker, to keep set_id stable
  arrange(grepl("^CELLMARKER", database), grepl("^PTMSIGDB", database), set) %>%
  mutate(set_id = sprintf(paste0("%0", n_digits, "d"), seq_len(n())),
         set_short = ifelse(database %in% c("MITOCARTA", "PSP", "REFMET", "CELLMARKER"),
                            sub("^[^_]+_", "", set), set),
         set_short = ifelse(database == "CELLMARKER", sub(" Human$", "", set_short), set_short),
         temp = .cutstr(set_short,
                        split = ifelse(database %in% c("MITOCARTA", "PSP", "REFMET", "CELLMARKER"),
                                       " ", "_"), n = 50L),
         set_short = ifelse(nchar(set_short) > nchar(temp) + 5L + n_digits + 4L,
                            sprintf("%s...(%s)", temp, set_id), set_short)) %>%
  select(-temp) %>%
  relocate(set_id, .before = set) %>%
  `rownames<-`(NULL) %>%
  mutate(collection = case_when(
    database %in% c("MITOCARTA", "REFMET", "PSP", "CELLMARKER") ~ database,
    grepl("^PTMSIGDB", database) ~ "PTMSIGDB",
    grepl("^GO", database) ~ "C5",
    TRUE ~ "C2"
  ),
  collection = factor(collection, levels = unique(collection))) %>%
  relocate(collection, .before = everything()) %>%
  mutate(across(.cols = everything(), .fns = ~ structure(.x, names = NULL)))

save(SET_TO_ID, file = file.path(data_dir, "SET_TO_ID.rda"), version = 3, compress = TRUE)
message("SET_TO_ID: ", nrow(SET_TO_ID), " rows -> ", file.path(data_dir, "SET_TO_ID.rda"))
