#!/usr/bin/env Rscript
# Per-file freeze versions, resolved from config/file_versions.json.
#
# What the map holds: the intended version for every file in the NEXT staging-bucket
# upload — the version this pipeline stamps on a file as it writes it. The map is
# hand-maintained, not a snapshot of the bucket. Its entries were seeded once from a
# listing of gs://pre-cawg/staging_20260720 and are edited by hand from there on, so
# they may legitimately differ from what that bucket carries today.
#
# Why this exists: the collection is deliberately NOT uniformly versioned. A file's
# version is the last collection release in which its CONTENT changed, so a file
# regenerated without a content change keeps its existing version while a changed one
# moves to NEW_VERSION. As seeded on 2026-07-29 that meant 54 files at v1.2, 51 at
# v1.3 and 121 at v1.4, with prot-*, ATAC and transcriptomics qc-norm at 1.2 and most
# DA tables at 1.3. The incoming release is v2.0, and the map has since been bumped
# for it: every metadata_features file, both clinical assays (prot-clinical and
# metab-t-clinical, all four files each), the feature-to-gene mapping and the prot-ol
# and transcript-rna-seq qc-norm matrices are now 2.0, while everything else holds the
# version of the release that last changed it. No
# single NEW_VERSION can name all of them, and spreading the versions across ~20
# hardcoded `version = "1.4"` call sites (which is what this replaces) put them out
# of reach of any single edit.
#
# Bumping a file for a release is therefore an edit to config/file_versions.json, not
# to a call site: that map is where the version is decided. NEW_VERSION in
# config/pipeline.env only supplies the fallback for a file the map does not list.
#
# HOW THE MAP IS SHAPED. `files` is a tree, keyed the way a freeze filename is built:
#
#   files[[tissue_code]][[ome]][[data_category]][[data_details]]
#     -> list(version =, ext =, why_bumped =)
#
# the same nesting load_qc()/load_differential_analysis() return, one level deeper. A
# filename is its path: the stem_prefix, then the four keys joined by "_", then
# _v<version>.<ext>. Every file whose name stops at data_category — the prot qc-report
# HTMLs — is keyed under the reserved data_details "_none" (.FV_NO_DETAILS below), so
# the tree stays a uniform four levels deep; no real data_details starts with "_".
# `resources` is a flat sibling map for the bucket's resources/ tier, whose filenames
# follow no tissue/ome convention.
#
# Two views come out of that one tree:
#   freeze_version(stem, ext)  the lookup every writer uses — unchanged by the reshape,
#                              because the tree is flattened back to stem.ext keys here
#   freeze_files()             the inventory, one row per file, with tissue_code / ome /
#                              data_category / data_details already split out
#
# freeze_files() is what replaced the string-splitting that consumers used to do on the
# flat keys: the tree already holds those fields, so nothing needs to re-derive them
# from a filename and nothing has to special-case the entries that have no
# data_details token.
#
# why_bumped is carried through to freeze_files() but nothing here acts on it — it is a
# hand-written note recording what content change earned the current version, and lives
# in the map so a bump and its reason travel together.
#
# Standalone-safe: sourced by lib/qc_helpers.R, by
# 07_build_human_feature_to_gene/HUMAN_FEATURE_TO_GENE.R (which deliberately does
# not source qc_helpers.R) and by docs/dependency_graph/build_depgraph.R, so it
# resolves the repo root itself and pulls in nothing but jsonlite.

.FV_ROOT <- local({
  r <- Sys.getenv("PRECOVID_ROOT")
  if (nzchar(r)) return(normalizePath(r))
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "config", "pipeline.env")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config", "pipeline.env")))
    stop("cannot find motrpac-human-presuspension-repro root; set PRECOVID_ROOT")
  d
})

# The data_details key standing for "this filename has no data_details token".
.FV_NO_DETAILS <- "_none"

# Every hand-editable field is validated on the way in. The map is edited by hand, and
# a mistyped entry that fell through to the NEW_VERSION fallback would ship a file at
# the wrong version silently; `where` names the tree path so the fix is obvious.
.fv_entry <- function(e, where){
  if (!is.list(e) || !all(c("version", "ext", "why_bumped") %in% names(e)))
    stop("config/file_versions.json: entry at ", where,
         " needs version, ext and why_bumped")
  bad <- vapply(c("version", "ext", "why_bumped"),
                function(f) !is.character(e[[f]]) || length(e[[f]]) != 1L, logical(1))
  if (any(bad))
    stop("config/file_versions.json: ", where, " has non-string ",
         paste(names(bad)[bad], collapse = ", "))
  if (!nzchar(e$version) || !nzchar(e$ext))
    stop("config/file_versions.json: ", where, " has an empty version or ext")
  e
}

# One row per file, flattened out of the tree: the four tree levels become columns and
# the path becomes the composed filename. A resources/ entry has no tissue_code, ome or
# category; a tree entry keyed "_none" has no data_details.
.fv_rows <- function(fv){
  row <- function(stem, e, tissue_code = NA_character_, ome = NA_character_,
                  data_category = NA_character_, data_details = NA_character_)
    data.frame(key = paste0(stem, ".", e$ext), stem = stem, ext = e$ext,
               version = e$version, tissue_code = tissue_code, ome = ome,
               data_category = data_category, data_details = data_details,
               why_bumped = e$why_bumped, stringsAsFactors = FALSE)

  out <- list()
  for (tissue in names(fv$files))
    for (ome in names(fv$files[[tissue]]))
      for (category in names(fv$files[[tissue]][[ome]]))
        for (details in names(fv$files[[tissue]][[ome]][[category]])) {
          where <- paste("files", tissue, ome, category, details, sep = " / ")
          e <- .fv_entry(fv$files[[tissue]][[ome]][[category]][[details]], where)
          named <- !identical(details, .FV_NO_DETAILS)
          out[[length(out) + 1L]] <- row(
            paste(c(fv$stem_prefix, tissue, ome, category, if (named) details),
                  collapse = "_"),
            e, tissue, ome, category, if (named) details else NA_character_)
        }
  for (stem in names(fv$resources))
    out[[length(out) + 1L]] <- row(
      stem, .fv_entry(fv$resources[[stem]], paste("resources", stem, sep = " / ")))

  d <- do.call(rbind, out)
  dup <- d$key[duplicated(d$key)]
  if (length(dup))
    stop("config/file_versions.json: two entries compose the same filename: ",
         paste(unique(dup), collapse = ", "))
  d[order(d$key), , drop = FALSE]
}

.FV_JSON <- local({
  f <- file.path(.FV_ROOT, "config", "file_versions.json")
  if (!file.exists(f))
    stop("missing version map: ", f, "\n  restore it from git — it is hand-maintained, ",
         "not regenerated from the staging bucket")
  fv <- jsonlite::fromJSON(f, simplifyVector = FALSE)
  if (!is.character(fv$stem_prefix) || length(fv$stem_prefix) != 1L)
    stop("config/file_versions.json: stem_prefix must be a single string")
  fv
})

# The constant leading field of every BIC-named freeze file. Read from the map rather
# than written here so the tree really is the whole filename rule; qc_helpers.R's
# write_with_path_name() composes stems from this same value.
FREEZE_STEM_PREFIX <- .FV_JSON$stem_prefix

.FV_MAP <- .fv_rows(.FV_JSON)

# version_counts is the map's own self-check on hand edits, and until now nothing
# compared it to the entries it claims to tally. A warning rather than an error: a
# stale count means the tally was not updated alongside a bump, which is worth saying
# out loud but is not a reason to stop a build from writing correctly-versioned files.
local({
  declared <- .FV_JSON$version_counts
  if (!length(declared)) return(invisible(NULL))
  actual <- table(.FV_MAP$version)
  count <- function(x, v) if (v %in% names(x)) as.integer(x[[v]]) else 0L
  vers <- sort(union(names(declared), names(actual)))
  off <- vers[vapply(vers, function(v) count(declared, v) != count(actual, v), logical(1))]
  if (length(off))
    warning("config/file_versions.json: version_counts disagrees with the entries — ",
            paste(sprintf("v%s says %d, found %d", off,
                          vapply(off, function(v) count(declared, v), integer(1)),
                          vapply(off, function(v) count(actual, v), integer(1))),
                  collapse = "; "), call. = FALSE)
})

# Lookup by composed filename. A list, not a named vector: `[[` on a missing name
# returns NULL here, which is what the fallback below tests for.
.FILE_VERSIONS <- stats::setNames(as.list(.FV_MAP$version), .FV_MAP$key)

# The freeze inventory: one row per file, with the tree path already split into
# tissue_code / ome / data_category / data_details (all NA for a resources/ entry, as
# is data_details for a file whose name stops at the category). Callers that need to
# reason about the collection as a whole — the DA completeness tests, the dependency
# graph — read this instead of re-parsing filenames.
freeze_files <- function() .FV_MAP

# NEW_VERSION from config/pipeline.env, for stems the map does not cover.
# Read from the file only, never from the environment. NEW_VERSION is pinned in
# pipeline.env (see its header), so honouring an inherited value would let a stem
# run outside a driver — where nothing sources pipeline.env to overwrite it — name
# freeze files after a release the config never chose. Parsing the file is what
# makes the pinned value authoritative on that path too.
.fv_new_version <- function(){
  f <- file.path(.FV_ROOT, "config", "pipeline.env")
  ln <- grep("(^|export )NEW_VERSION=", readLines(f), value = TRUE)
  if (!length(ln)) return("")
  gsub('^"|"$', "", sub(".*NEW_VERSION=", "", ln[length(ln)]))
}

# Version for a composed file stem, e.g.
#   freeze_version("human-precovid-sed-adu_t06-muscle_transcript-rna-seq_qc-norm_log-cpm")
#
# A stem absent from the map falls back to NEW_VERSION: the map enumerates the files
# the next upload is expected to contain, so a miss means this pipeline emits something
# the map does not plan for, and a genuinely new artifact belongs at the version this
# run writes. The fallback is announced rather than silent — a miss is usually a renamed
# data_details, and silently versioning it NEW_VERSION would hide the rename. If the
# file is meant to ship, add it to the map rather than leaning on the fallback.
freeze_version <- function(stem, ext = "txt"){
  v <- .FILE_VERSIONS[[paste0(stem, ".", ext)]]
  if (!is.null(v) && nzchar(v)) return(v)
  fallback <- .fv_new_version()
  if (!nzchar(fallback))
    stop("no version for '", stem, ".", ext, "' in config/file_versions.json ",
         "and NEW_VERSION is not set in config/pipeline.env")
  message("file_versions: '", stem, ".", ext, "' is not in config/file_versions.json — ",
          "using NEW_VERSION ", fallback, "; add an entry if this file is here to stay")
  fallback
}
