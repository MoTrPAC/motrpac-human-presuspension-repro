# Shared RefMet definitions: how the pipeline names metabolites, and how it reads
# the frozen RefMet/KEGG snapshot that sits next to this file.
#
# Both sides of the snapshot need this code and must agree exactly:
#   - build_refmet_cache.R  derives the query-name universe it sends to RefMet
#   - lib/qc_helpers.R      sources this so the metabolomics stem resolves those
#                           same names against the snapshot at run time
# If the two derived names differently the snapshot would look complete while
# silently missing keys, so the derivation lives here once and only here.
#
# Deliberately dependency-free (base R + dplyr/stringr) and free of pipeline data
# objects, so the Stage 0 builder can source it without pulling in Stage 1's
# preflight/pheno objects.

suppressWarnings(suppressMessages(library(dplyr)))

REFMET_MAP_FILE  <- "refmet_name_map.txt"
REFMET_KEGG_FILE <- "kegg_compound_list.txt"
REFMET_META_FILE <- "refmet_snapshot.json"

# The RefMet batch endpoint. The only URL the snapshot is built from; nothing in a
# pipeline run may contact it.
REFMET_ENDPOINT <- "https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_new_minID.php"

# RefMet writes "-" for a field it has no value for, and the batch endpoint mixes
# "-" with empty strings for the same meaning. Both collapse to this placeholder
# before any join, and `na_if(x, REFMET_MISSING)` turns it into NA afterwards.
REFMET_MISSING <- "-"

# Where the snapshot files sit. A caller that already knows the repo root (both do)
# sets REFMET_DIR before sourcing; otherwise fall back to walking up for
# config/pipeline.env, so `Rscript refmet_lib.R`-style standalone use still works.
if (!exists("REFMET_DIR")) {
  REFMET_DIR <- local({
    d <- Sys.getenv("PRECOVID_ROOT")
    if (!nzchar(d)) {
      d <- normalizePath(getwd())
      while (!file.exists(file.path(d, "config", "pipeline.env")) && dirname(d) != d) d <- dirname(d)
    }
    file.path(normalizePath(d), "scripts", "00_preflight", "data-raw", "sources", "refmet")
  })
}

.refmet_dir <- function() REFMET_DIR

# --- name derivation ---------------------------------------------------------

# Manual corrections that align the submitted refmet_name with the current RefMet
# standard. Two correction sources are merged into one pass: feature_id-based
# overrides for cases where the stored refmet_name is wrong/missing, and a
# name-to-name map covering capitalization errors, typos, slash-delimited
# ambiguities, and lab-internal naming conventions. Some labs submitted
# annotations using outdated LIPID MAPS names or lab-internal aliases that
# predate the current RefMet standard.
# metab: data frame with columns feature_id and refmet_name
# Returns: metab with refmet_name corrected in place
refmet_fix_names <- function(metab) {
  # feature_id-based overrides
  metab <- metab %>%
    dplyr::mutate(
      refmet_name = dplyr::case_when(
        feature_id == "13-HODE" ~ "13-HODE",
        feature_id == "13 HODE" ~ "13-HODE", #new for 1.4
        feature_id == "FA(20:1)" ~ "Eicosenoic acid", #new for 1.4
        feature_id == "FA(18:2)" ~ "Linoleic acid", #new for 1.4
        feature_id == "FA(22:6)" ~ "Docosahexaenoic acid", #new for 1.4
        TRUE ~ refmet_name
      )
    )

  refmet_name_map <- c(
    "cholesterol sulfate"              = "Cholesterol sulfate",
    "oleamide"                         = "Oleamide",
    "tridecylamine"                    = "Tridecylamine",
    "PGE2 thanolamide"                 = "PGE2 ethanolamide",
    "12(13)-EpMOE"                     = "12(13)-EpOME",
    "Docosatetraenoic aicd"            = "Docosatetraenoic acid",
    "Tetracosenoic aicd"               = "Tetracosenoic acid",
    "Prostagladin"                     = "Prostaglandin", #new for 1.4

    "Chenodeoxycholic acid\\Deoxycholic acid"   = "Deoxycholic acid",
    "Glycocholic acid\\Glycohyocholic acid"     = "Glycocholic acid",

    "N-Lauroylglycine"                 = "NAGly 12:0",
    "N-linoleoylglycine"               = "NAGly 18:2(9Z,12Z)",
    "N-Oleoyl glycine"                 = "NAGly 18:1(9Z)",
    "N-Undecanoylglycine"              = "NAGly 11:0",
    "N-Myristoylglycine"               = "NAGly 14:0",

    "CAR 12:0-OH"                      = "CAR 12:0;OH",
    "CAR 14:0-OH"                      = "CAR 14:0;OH",
    "CAR 14:1-OH"                      = "CAR 14:1;OH",
    "CAR 4:0-OH"                       = "CAR 4:0;OH",

    "DG 16:0_16:0_0:0"                 = "DG 16:0_16:0",
    "DG 16:0_16:1_0:0"                 = "DG 16:0_16:1",
    "DG 16:0_18:1_0:0"                 = "DG 16:0_18:1",
    "DG 18:2_18:2_0:0"                 = "DG 18:2_18:2",
    "DG 18:2_18:3_0:0"                 = "DG 18:2_18:3",

    "13-OxoODE(13-KODE)"               = "13-Oxo-ODE",
    "2-Arachidonoyl Glycerol (2AG)"    = "MG 0:0/20:4/0:0",
    "5,6-DiHET"                        = "5,6-DiHETE",
    "8(9)-DiHET"                       = "8,9-DiHETE",
    "CoA(15:0)_and_CoA(C14:1-OH)"     = "Pentadecanoyl-CoA/Hydroxytetradecenoyl-CoA",
    "CoA(2:0-COOH)_and_CoA(4:0-OH)"   = "Malonyl-CoA/Hydroxybutyryl-CoA",
    "Linoleoyl Ethanolamide (LEA)"     = "Linoleoyl-EA",
    "Oleoyl Ethanolamide (OEA)"        = "Oleoyl-EA",
    "PC(O-33:2)>PC(O-15:0/18:2)"      = "PC O-16:1/20:4",
    "PC(O-36:5)<PC(O-16:1/20:4)"      = "PC O-16:1/20:4",
    "PE(36:4)>(16:0_20:4)"            = "PE 16:0_20:4",
    "PE(38:4)>(PE(18:0_20:4)"         = "PE 18:0_20:4",
    "Stearoyl Ethanolamide (ceramid)"  = "Stearoyl-EA"
  )

  metab <- metab %>%
    dplyr::mutate(refmet_name = dplyr::recode(refmet_name, !!!refmet_name_map, .default = refmet_name))

  return(metab)
}

# The exact string RefMet is asked about, given a (fixed) refmet_name. Untargeted
# platforms append an LC suffix (_hp_a, _rp_b, ...) that RefMet does not know, so
# it is stripped before lookup. This is the snapshot's key column.
# refmet_name: character vector
# Returns: character vector of the same length
refmet_lookup_key <- function(refmet_name) {
  tosearch <- "_hp_|_rp_|_rn_|_in_|_lp_|_ln_"
  key <- dplyr::if_else(
    grepl(tosearch, refmet_name),
    gsub("(.*)(_\\w{2}_\\w{1})", "\\1", refmet_name),
    refmet_name
  )
  trimws(key)
}

# --- snapshot readers --------------------------------------------------------

# Restore the runtime's missing-value placeholder. The snapshot is written with
# empty cells so it reads as an ordinary TSV; the joins downstream expect "-".
.refmet_restore_missing <- function(df) {
  df[is.na(df)] <- REFMET_MISSING
  df[df == ""]  <- REFMET_MISSING
  df
}

.refmet_read_tsv <- function(path) {
  # colClasses/na.strings keep every field a literal string: RefMet ships names
  # containing quotes and "NA"-like tokens, and the runtime join is string-based.
  read.csv(path, sep = "\t", check.names = FALSE, colClasses = "character",
           na.strings = character(0), quote = "", comment.char = "")
}

.refmet_require <- function(path) {
  if (!file.exists(path))
    stop("RefMet snapshot not found: ", path,
         "\n  fetch it: make sources",
         "\n  or rebuild it: Rscript ", file.path(.refmet_dir(), "build_refmet_cache.R"))
  path
}

# The frozen RefMet batch-endpoint response, keyed by `Input name` (the value
# refmet_lookup_key() produces). Every column RefMet returned is kept.
refmet_map <- function() {
  .refmet_restore_missing(.refmet_read_tsv(.refmet_require(file.path(.refmet_dir(), REFMET_MAP_FILE))))
}

# The frozen KEGGREST::keggList("compound") response: kegg_entry ("cpd:C00031")
# and the ";"-delimited synonym string KEGG returns as the value.
kegg_compound_list <- function() {
  .refmet_read_tsv(.refmet_require(file.path(.refmet_dir(), REFMET_KEGG_FILE)))
}

# Build date / endpoint / row counts for the snapshot, as written by the builder.
refmet_snapshot_meta <- function() {
  jsonlite::fromJSON(.refmet_require(file.path(.refmet_dir(), REFMET_META_FILE)))
}
