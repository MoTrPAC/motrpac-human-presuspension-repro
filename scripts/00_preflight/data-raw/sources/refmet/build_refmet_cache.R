#!/usr/bin/env Rscript
# One-time builder for the local RefMet / KEGG annotation snapshot.
#
# This is the ONLY script that contacts Metabolomics Workbench or KEGG. It builds
# two artifacts into this folder so the metabolomics stem can annotate features
# OFFLINE and deterministically (both were live calls inside
# .build_metab_refmet_map(), which made refmet_name / refmet_id / kegg_id depend on
# the database state on the day the pipeline happened to run):
#
#   1. refmet_name_map.txt      — the RefMet batch-endpoint response for every
#        metabolite name the pipeline asks about, keyed by `Input name`. All 14
#        returned columns are kept (class hierarchy, formula, HMDB/ChEBI/LIPID
#        MAPS/PubChem/KEGG cross-references, RefMet_ID).
#   2. kegg_compound_list.txt   — KEGGREST::keggList("compound") verbatim, used to
#        fill residual kegg_id values RefMet has no cross-reference for.
#
# Plus refmet_snapshot.json, recording the build date, endpoint and row counts.
#
# The query universe comes from the raw metabolomics metabolite-metadata files
# (staging/raw-files/metabolomics_qc_norm), read straight off disk rather than
# through the BIC parser: the parser only ever subsets rows, so the raw files are a
# superset of what reaches .build_metab_refmet_map(). Names are put through
# refmet_fix_names() + refmet_lookup_key() from refmet_lib.R — the same code the
# stem uses at run time, so a name the stem will look up cannot be missing here.
#
# A refresh keeps every key already in the snapshot and adds whatever the raw files
# now contain, so re-running with only some platforms downloaded never drops
# coverage. Both artifacts are skipped if they already exist; set FORCE=1 to rebuild.
# See README.md for provenance.

OUT_DIR <- dirname(normalizePath(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
if (is.na(OUT_DIR) || !nzchar(OUT_DIR)) OUT_DIR <- getwd()

REFMET_DIR <- OUT_DIR
suppressWarnings(suppressMessages(library(dplyr)))
source(file.path(OUT_DIR, "refmet_lib.R"))

MAP_PATH  <- file.path(OUT_DIR, REFMET_MAP_FILE)
KEGG_PATH <- file.path(OUT_DIR, REFMET_KEGG_FILE)
META_PATH <- file.path(OUT_DIR, REFMET_META_FILE)

FORCE <- !Sys.getenv("FORCE") %in% c("", "0", "false", "FALSE")

# Raw metabolomics inputs. Defaults to this repo's staging scratch; override with
# METAB_RAW_DIR to build from a copy elsewhere.
RAW_DIR <- Sys.getenv("METAB_RAW_DIR")
if (!nzchar(RAW_DIR)) {
  root <- Sys.getenv("PRECOVID_ROOT")
  if (!nzchar(root)) {
    root <- normalizePath(getwd())
    while (!file.exists(file.path(root, "config", "pipeline.env")) && dirname(root) != root)
      root <- dirname(root)
  }
  RAW_DIR <- file.path(root, "staging", "raw-files", "metabolomics_qc_norm")
}

MAX_TRIES <- 5
with_retry <- function(what, expr) {
  for (i in seq_len(MAX_TRIES)) {
    res <- tryCatch(force(expr), error = function(e) e)
    if (!inherits(res, "error")) return(res)
    message(sprintf("[%s] attempt %d/%d failed: %s", what, i, MAX_TRIES, conditionMessage(res)))
    Sys.sleep(5)
  }
  stop(sprintf("%s: exhausted %d attempts (endpoint unreachable?)", what, MAX_TRIES))
}

# Write a snapshot table: tab-delimited, unquoted, missing cells left empty. The
# readers in refmet_lib.R restore "-" before joining, so the file on disk stays a
# plain TSV for anyone reading it out of resources/.
write_snapshot <- function(df, path) {
  df[] <- lapply(df, function(x) { x <- as.character(x); x[is.na(x) | x == REFMET_MISSING] <- ""; x })
  utils::write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
}

# --- 1. query-name universe --------------------------------------------------
# Only the `_named-` metabolite-metadata files: .build_metab_refmet_map() keeps
# is_named == TRUE rows, and the `_unnamed-` files are exactly the rows it drops.
collect_query_names <- function() {
  if (!dir.exists(RAW_DIR))
    stop("raw metabolomics inputs not found: ", RAW_DIR,
         "\n  download them first (step 06 metab stem with redownload=TRUE), or set METAB_RAW_DIR")
  files <- list.files(RAW_DIR, pattern = "_named-metadata-metabolites.*\\.txt$",
                      full.names = TRUE, recursive = TRUE)
  if (!length(files))
    stop("no `_named-metadata-metabolites` files under ", RAW_DIR)
  message("Reading ", length(files), " metabolite-metadata files ...")

  annot <- dplyr::bind_rows(lapply(files, function(f) {
    d <- read.csv(f, sep = "\t", check.names = FALSE, colClasses = "character")
    if (!all(c("metabolite_name", "refmet_name") %in% names(d))) {
      message("  skipping (no metabolite_name/refmet_name): ", basename(f))
      return(NULL)
    }
    d[, c("metabolite_name", "refmet_name")]
  }))

  keys <- annot %>%
    dplyr::rename(feature_id = metabolite_name) %>%
    refmet_fix_names() %>%
    dplyr::mutate(lookup_refmet = refmet_lookup_key(refmet_name)) %>%
    dplyr::pull(lookup_refmet)

  keys <- unique(keys[!is.na(keys) & nzchar(keys) & keys != REFMET_MISSING])
  message("  ", nrow(annot), " annotated features -> ", length(keys), " distinct RefMet queries")
  sort(keys)
}

# --- 2. RefMet batch response ------------------------------------------------
# One POST with every name, newline-delimited — the same request the stem used to
# make at run time, so the frozen response is what a live run would have received.
query_refmet <- function(keys) {
  mets <- paste(keys, collapse = "\n")
  h <- curl::new_handle()
  curl::handle_setform(h, metabolite_name = mets)
  req <- with_retry("refmet batch POST", curl::curl_fetch_memory(REFMET_ENDPOINT, handle = h))
  if (req$status_code != 200)
    stop("RefMet returned HTTP ", req$status_code, " for the batch query")

  # check.names = FALSE keeps RefMet's own header ("Input name", "Standardized
  # name", ...) so the snapshot reads as the endpoint's table, not a mangled copy.
  res <- utils::read.csv(text = rawToChar(req$content), header = TRUE,
                         colClasses = "character", na.strings = character(0),
                         check.names = FALSE, quote = "", comment.char = "", sep = "\t")
  if (!"Input name" %in% names(res))
    stop("unexpected RefMet response columns: ", paste(names(res), collapse = ", "))
  unique(res)
}

if (file.exists(MAP_PATH) && !FORCE) {
  message("RefMet map already present, skipping build: ", MAP_PATH)
  refmet_tab <- .refmet_read_tsv(MAP_PATH)
} else {
  keys <- collect_query_names()

  # Union with what the snapshot already covers, so refreshing from a partial set
  # of downloaded platforms cannot silently shrink coverage.
  if (file.exists(MAP_PATH)) {
    prev <- .refmet_read_tsv(MAP_PATH)$`Input name`
    added <- setdiff(keys, prev)
    keys <- sort(union(keys, prev))
    message("  refresh: ", length(prev), " already covered, ", length(added), " new, ",
            length(keys), " total")
  }

  message("Querying RefMet for ", length(keys), " names ...")
  refmet_tab <- query_refmet(keys)

  missing <- setdiff(keys, refmet_tab$`Input name`)
  if (length(missing))
    stop("RefMet did not echo ", length(missing), " of the ", length(keys),
         " submitted names (e.g. ", paste(utils::head(missing, 3), collapse = " | "),
         ") — the snapshot would be incomplete; re-run when the endpoint is healthy")

  write_snapshot(refmet_tab, MAP_PATH)
  resolved <- sum(!refmet_tab$`Standardized name` %in% c("", REFMET_MISSING))
  message(sprintf("  saved %s (%d rows, %d resolved to a RefMet name)",
                  MAP_PATH, nrow(refmet_tab), resolved))
}

# --- 3. KEGG compound list ---------------------------------------------------
# Stored verbatim (entry + the ";"-delimited synonym string); the name-to-id
# parsing stays in the stem so the snapshot is a plain copy of what KEGG serves.
if (file.exists(KEGG_PATH) && !FORCE) {
  message("KEGG compound list already present, skipping fetch: ", KEGG_PATH)
  kegg_tab <- .refmet_read_tsv(KEGG_PATH)
} else {
  message("Fetching the KEGG compound list ...")
  all_compounds <- with_retry("keggList compound", KEGGREST::keggList("compound"))
  kegg_tab <- data.frame(kegg_entry = names(all_compounds),
                         kegg_names = as.character(all_compounds),
                         stringsAsFactors = FALSE)
  write_snapshot(kegg_tab, KEGG_PATH)
  message(sprintf("  saved %s (%d compounds)", KEGG_PATH, nrow(kegg_tab)))
}

# --- 4. provenance -----------------------------------------------------------
meta <- list(
  built_on          = format(Sys.Date()),
  refmet_endpoint   = REFMET_ENDPOINT,
  kegg_source       = "KEGGREST::keggList(\"compound\")",
  refmet_queries    = nrow(refmet_tab),
  refmet_resolved   = sum(!refmet_tab$`Standardized name` %in% c("", REFMET_MISSING)),
  kegg_compounds    = nrow(kegg_tab)
)
writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE), META_PATH)
message("RefMet snapshot build complete — ", META_PATH)
