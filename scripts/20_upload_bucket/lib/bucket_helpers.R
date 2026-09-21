#!/usr/bin/env Rscript
# Shared helpers for the Stage 2 steps: bucket and freeze inventories, BIC file-name
# parsing, and the paths every step reads and writes.
#
# Standalone-safe in the same sense as 10_build_data/lib/file_versions.R: it resolves the
# repo root itself and pulls in nothing but base R plus base64enc, so a step script can be
# run on its own with Rscript.

.BH_ROOT <- local({
  r <- Sys.getenv("PRECOVID_ROOT")
  if (nzchar(r)) return(normalizePath(r))
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "config", "pipeline.env")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config", "pipeline.env")))
    stop("cannot find motrpac-human-presuspension-repro root; set PRECOVID_ROOT")
  d
})

# Values pinned in config/pipeline.env. Parsed from the file rather than read from the
# environment for the reason file_versions.R gives about NEW_VERSION: the buckets a run
# reads and writes are decided in that file and nowhere else, so a step invoked outside
# the driver — where nothing has sourced pipeline.env — cannot be pointed at a different
# bucket by an inherited variable.
.pipeline_env <- function(key) {
  f <- file.path(.BH_ROOT, "config", "pipeline.env")
  ln <- grep(paste0("(^|export )", key, "="), readLines(f), value = TRUE)
  if (!length(ln)) return("")
  v <- gsub('^"|"$', "", sub(paste0(".*", key, "="), "", ln[length(ln)]))
  # pipeline.env interpolates ${CURRENT_VERSION} / ${NEW_VERSION} into the bucket paths.
  for (k in c("CURRENT_VERSION", "NEW_VERSION")) {
    if (grepl(paste0("\\$\\{", k, "\\}"), v)) v <- gsub(paste0("\\$\\{", k, "\\}"), .pipeline_env(k), v)
  }
  v
}

CURRENT_VERSION   <- .pipeline_env("CURRENT_VERSION")
NEW_VERSION       <- .pipeline_env("NEW_VERSION")
PRODUCTION_BUCKET <- .pipeline_env("PRODUCTION_BUCKET")
STAGING_BUCKET    <- .pipeline_env("STAGING_BUCKET")
GSUTIL            <- Sys.getenv("GSUTIL", unset = "gsutil")

# BIC file-name header shared by every uploaded artifact.
FILE_HEADER <- "human-precovid-sed-adu"

# Stage 1 writes the freeze here. FREEZE_DIR overrides it, which is what lets a step run
# against a freeze built in another checkout — staging/ is gitignored, so a worktree has
# no freeze of its own.
FREEZE_DIR <- local({
  d <- Sys.getenv("FREEZE_DIR")
  if (nzchar(d)) normalizePath(d, mustWork = FALSE) else file.path(.BH_ROOT, "staging", "freeze")
})

# Step outputs (snapshots, diffs) live in the stage's own data/ dir, the same way
# 10_build_data steps write their .rda into scripts/10_build_data/data/. Run reports and
# per-step logs go to the repo-level logs/ dir like every other stage.
STAGE_DIR    <- file.path(.BH_ROOT, "scripts", "20_upload_bucket")
STAGE_DATA   <- file.path(STAGE_DIR, "data")
SNAPSHOT_DIR <- file.path(STAGE_DATA, "snapshots")
DIFF_DIR     <- file.path(STAGE_DATA, "diffs")
LOG_DIR      <- file.path(.BH_ROOT, "logs")

for (d in c(SNAPSHOT_DIR, DIFF_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Extensions that are upload artifacts.
#
# .sqlite and .json are here for resources/: step 06's stage_ensembl_txdb publishes the
# offline Ensembl v105 TxDb and stage_refmet_map the RefMet snapshot plus its provenance
# JSON, both so the bucket becomes their distribution point (`make sources` fetches the
# TxDb back out of resources/ rather than rebuilding it against a live Ensembl). From v2.0
# they carry a _v<version> token like every other file; the join key still accepts a
# pre-v2.0 incumbent with none.
UPLOAD_EXT_RE <- "\\.(txt|txt\\.gz|csv|html|sqlite|json)$"

.timestamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

# ---- Bucket inventory --------------------------------------------------------------
# One `gsutil ls -L -r` pass over the whole bucket instead of a `gsutil stat` per file.
# Upstream 01_snapshot_production.R shells out once per object, which is ~227 round trips
# against this bucket; -L returns the same Content-Length and Hash (md5) fields for every
# object in a single call. Same columns out, minutes faster.
gcs_inventory <- function(bucket) {
  raw <- system(paste(GSUTIL, "ls -L -r", shQuote(paste0(bucket, "/**"))), intern = TRUE)

  # Each object is a `gs://...:` header line followed by indented key/value fields.
  starts <- grep("^gs://.*:$", raw)
  if (!length(starts)) {
    return(data.frame(gcs_path = character(0), size_bytes = character(0),
                      md5 = character(0), stringsAsFactors = FALSE))
  }
  ends <- c(starts[-1] - 1L, length(raw))

  rows <- lapply(seq_along(starts), function(i) {
    block <- raw[starts[i]:ends[i]]
    path  <- sub(":$", "", block[1])
    field <- function(label) {
      ln <- grep(label, block, fixed = TRUE, value = TRUE)
      if (!length(ln)) return(NA_character_)
      trimws(sub("^[^:]*:", "", ln[1]))
    }
    data.frame(
      gcs_path   = path,
      size_bytes = field("Content-Length:"),
      md5        = md5_base64_to_hex(field("Hash (md5):")),
      stringsAsFactors = FALSE
    )
  })

  inv <- do.call(rbind, rows)
  inv[grepl(UPLOAD_EXT_RE, inv$gcs_path), , drop = FALSE]
}

# gsutil reports md5 base64-encoded; tools::md5sum returns hex. Everything downstream
# compares hex, so the conversion happens once, here.
md5_base64_to_hex <- function(b64) {
  if (is.na(b64) || !nzchar(b64)) return(NA_character_)
  tryCatch({
    raw <- base64enc::base64decode(b64)
    # base64decode does not reject malformed input — it skips what it cannot read and
    # returns whatever it got, so a truncated or garbled Hash (md5) line would otherwise
    # become a plausible-looking hex string that silently never matches. An md5 is 16
    # bytes; anything else is not one.
    if (length(raw) != 16L) return(NA_character_)
    paste(sprintf("%02x", as.integer(raw)), collapse = "")
  }, error = function(e) NA_character_)
}

# ---- Freeze inventory --------------------------------------------------------------

freeze_inventory <- function(dir = FREEZE_DIR) {
  if (!dir.exists(dir))
    stop("freeze dir not found: ", dir, "\n  run `make data` first, or set FREEZE_DIR")

  files <- list.files(dir, recursive = TRUE, full.names = TRUE)
  files <- files[grepl(UPLOAD_EXT_RE, files)]
  if (!length(files))
    stop("no upload artifacts under ", dir, " — is this a freeze dir?")

  data.frame(
    local_path = files,
    rel_path   = sub(paste0("^", normalizePath(dir), "/?"), "", normalizePath(files)),
    size_bytes = as.character(file.size(files)),
    md5        = as.character(tools::md5sum(files)),
    stringsAsFactors = FALSE
  )
}

# ---- File naming -------------------------------------------------------------------
# BIC pattern: human-precovid-sed-adu_{tissue_code}_{ome}_{data_category}_{data_details}_v{ver}.txt
# resources/ files do not follow it (motrpac-mappings-human-feature-to-gene_v2.0.txt,
# motrpac_human-precovid_1kg_pca_v2.0.csv), so every field comes back NA for them and callers
# fall back to the path.

parse_bic_name <- function(path) {
  base  <- basename(path)
  parts <- strsplit(sub(paste0(UPLOAD_EXT_RE, "$"), "", base), "_")[[1]]
  na <- list(tissue_code = NA_character_, ome = NA_character_,
             data_category = NA_character_, data_details = NA_character_,
             version = strip_to_version(base))
  if (length(parts) < 6 || parts[1] != FILE_HEADER) return(na)

  list(
    tissue_code   = parts[2],
    ome           = parts[3],
    data_category = parts[4],
    data_details  = paste(parts[5:(length(parts) - 1)], collapse = "_"),
    version       = sub("^v", "", parts[length(parts)])
  )
}

# The trailing _v<major>.<minor> of a file name, or NA when it carries no version.
strip_to_version <- function(path) {
  m <- regmatches(basename(path), regexpr("_v[0-9]+\\.[0-9]+(?=(\\.[a-z]+)+$)",
                                          basename(path), perl = TRUE))
  if (!length(m)) return(NA_character_)
  sub("^_v", "", m)
}

# Sortable rank for a "<major>.<minor>" version string; unversioned files rank lowest.
# Used wherever the highest version of a stem has to win — a plain string sort would put
# v1.10 below v1.9.
.ver_rank <- function(v) {
  p <- do.call(rbind, lapply(strsplit(ifelse(is.na(v), "0.0", v), ".", fixed = TRUE),
                             function(x) as.numeric(c(x, "0")[1:2])))
  p[, 1] * 1000 + p[, 2]
}

# The join key: the path with its version suffix removed. Upstream 03_diff_local_vs_staging.R
# keys on the BASENAME alone; this keys on the path relative to the bucket root, so a file
# that moves between top-level subdirs is reported as an add plus an orphan rather than
# silently matching its old location. That is what the clinical split looks like this cycle
# — clinical_chemistry/ retired, prot-clinical and metab-t-clinical filed under
# proteomics/ and metabolomics-targeted/ — and it should be visible in the diff, not hidden.
#
# The version token and the extension come off together, so an unversioned path keeps its
# extension and a versioned one does not: "foo.txt" keys as "foo.txt", "foo_v2.0.txt" as
# "foo". A file that gains a version token therefore does not pair with its own unversioned
# predecessor — the new name is an ADDED and the old one an ORPHANED, and Stage 2 deletes
# only what it just superseded, so the predecessor stays in the bucket.
#
# v2.0 gives every file a version token, and one incumbent was unversioned in v1.3 and so
# hit this: resources/motrpac-mappings-human-feature-to-gene.txt, whose replacement is
# motrpac-mappings-human-feature-to-gene_v2.0.txt. It was deleted from the bucket by hand.
# A future unversioned-to-versioned rename needs the same manual deletion, or this function
# has to strip an optional version token while keeping the extension so the two pair up.
join_key <- function(rel_path) {
  sub("_v[0-9]+\\.[0-9]+(\\.[a-z]+)+$", "", rel_path)
}

# Extension of a relative path, kept so ADDED rows can be reassembled into a file name.
path_ext <- function(rel_path) {
  m <- regmatches(rel_path, regexpr(UPLOAD_EXT_RE, rel_path))
  if (!length(m)) NA_character_ else m
}

# ---- Change classification ---------------------------------------------------------
# The decision step 04 acts on, including which object it DELETES, so it lives here as a
# pure function rather than inline in the step — that is what makes it testable without a
# bucket. See 03_diff_freeze_vs_staging/README.md for what each type means.
#
# old_version / old_md5 are NA when the stem is not in the bucket at all.
classify_freeze_file <- function(new_version, old_version, new_md5, old_md5,
                                 allow_conflict = FALSE) {
  if (is.na(old_version) && is.na(old_md5)) return("ADDED")
  if (!identical(new_version, old_version)) {
    # REVERSIONED is the mirror image of CONFLICT: same bytes under a new number, where
    # CONFLICT is new bytes under the same number. Both mean the freeze and
    # config/file_versions.json disagree about whether this file changed this cycle, and
    # both are mistakes the hand-maintained map can make — the difference is that this one
    # is not dangerous. The upload is a pure rename: no consumer reads different bytes, it
    # just reads the same ones at a new version, and everything pinning the old number has
    # to move for no reason. So it uploads and retires the old object exactly as MODIFIED
    # does, but is counted separately and reported WARN, because the usual fix is to revert
    # the bump in the map rather than to ship it.
    if (!is.na(old_md5) && identical(new_md5, old_md5)) return("REVERSIONED")
    return("MODIFIED")
  }
  # Same version. Equal bytes is the only way to conclude "unchanged" — an NA md5 (a
  # composite object carries only a crc32c) means we cannot tell, and cannot tell is not
  # the same as unchanged. That same guard is why an md5-less object at a moved version
  # comes back MODIFIED rather than REVERSIONED above.
  if (!is.na(old_md5) && identical(new_md5, old_md5)) return("UNCHANGED")
  if (allow_conflict) "REPLACE" else "CONFLICT"
}

# EQUIVALENT is to CONFLICT what REVERSIONED is to MODIFIED — the md5s disagree and the
# content does not. Recomputing a matrix under a different BLAS or summation order moves
# the last bits of a double, and write.table prints enough significant digits to show it,
# so a file nothing meaningful changed in comes back as new bytes under its old version.
# Comparing the numbers rather than the bytes separates that from a real edit.
#
# Called only on a file that already classified CONFLICT or REPLACE, because it downloads
# the incumbent: the md5 comparison above needs no bucket read, and this one does.
#
# Everything it cannot positively establish is FALSE, so the caller keeps the CONFLICT:
# a non-.txt file, a download or parse failure, a changed shape or column set, a changed
# column type, any difference in a non-numeric column, or a table with no numeric column
# at all — where differing bytes can only mean a real difference.
freeze_numerically_equivalent <- function(local_path, staging_gcs, tolerance = 1e-9) {
  if (is.na(local_path) || is.na(staging_gcs)) return(FALSE)
  if (!identical(tolower(sub("^\\.", "", path_ext(local_path))), "txt")) return(FALSE)

  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  rc <- suppressWarnings(system2(GSUTIL, c("cp", shQuote(staging_gcs), shQuote(tmp)),
                                 stdout = FALSE, stderr = FALSE))
  if (!identical(as.integer(rc), 0L) || !file.exists(tmp)) return(FALSE)

  read_tsv <- function(p) tryCatch(read.csv(p, sep = "\t", check.names = FALSE),
                                   error = function(e) NULL)
  a <- read_tsv(local_path); b <- read_tsv(tmp)
  if (is.null(a) || is.null(b)) return(FALSE)
  if (!identical(dim(a), dim(b)) || !identical(colnames(a), colnames(b))) return(FALSE)

  num <- vapply(a, is.numeric, logical(1))
  if (!identical(num, vapply(b, is.numeric, logical(1)))) return(FALSE)
  if (!any(num)) return(FALSE)
  if (any(!num) && !identical(a[!num], b[!num])) return(FALSE)

  isTRUE(all.equal(as.matrix(a[num]), as.matrix(b[num]), tolerance = tolerance))
}

# ---- Reporting ---------------------------------------------------------------------
# Same three-column report shape scripts/lib/common.sh writes, so a Stage 2 R step's
# report reads like a shell stage's.

report_write <- function(df, path) {
  write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
  message("wrote ", path)
}

# The written record of an in-place overwrite.
#
# Every other thing step 04 does is legible from the bucket afterwards: an ADDED object is
# a path that did not exist, a MODIFIED or REVERSIONED one moved to a new version and took
# its predecessor with it. A REPLACE leaves no trace at all. The path is the same, the
# version is the same, only the bytes changed, so six months from now nothing in the bucket
# distinguishes an object that was overwritten from one that was never touched -- and
# anything pinning that version silently reads different content than it did before.
#
# Hence this file. It is the only place that overwrite is recorded, so it carries both md5s
# per object: the pair is what lets someone confirm which of the two versions they have.
write_replace_readme <- function(rows, path, applied, bucket, verified = NULL) {
  if (nrow(rows) == 0) return(invisible(NULL))

  verb <- if (applied) "were overwritten" else "would be overwritten"
  out <- c(
    "# Objects overwritten in place -- staging bucket record",
    "",
    paste0("Bucket:    ", bucket),
    paste0("Written:   ", .timestamp(), if (applied) "" else "  (DRY RUN -- nothing was written)"),
    paste0("Objects:   ", nrow(rows)),
    "",
    "## What this file records",
    "",
    paste0("These objects ", verb, " at their EXISTING version. Their content changed; ",
           "their version did not."),
    "",
    "Step 03 classifies that as CONFLICT and stops the run, because publishing new bytes",
    "under a version number that already means something else cannot be untangled from the",
    "bucket afterwards. This run set ALLOW_CONFLICT=1, which reclassifies them as REPLACE",
    "and overwrites in place -- a deliberate choice that the freeze is right and the version",
    "is not moving.",
    "",
    "The consequence to be aware of: nothing in the bucket marks these objects as changed.",
    "Any consumer pinned to this version gets different content than it did before, with no",
    "version bump to signal it. That is what this file exists to record.",
    "",
    "## Objects",
    ""
  )

  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    out <- c(out,
      paste0("### ", r$rel_path),
      "",
      paste0("- version:  v", r$new_version, " (unchanged)"),
      paste0("- object:   ", r$gcs_path),
      paste0("- md5 before: ", r$old_md5),
      paste0("- md5 after:  ", r$new_md5),
      if (!is.null(verified)) paste0("- md5-verified after upload: ", verified[i]) else NULL,
      ""
    )
  }

  out <- c(out,
    "## Checking one of these",
    "",
    "    gsutil stat <object>            # Hash (md5) is base64; the values above are hex",
    "",
    "An object matching `md5 after` is the replacement. One matching `md5 before` predates",
    "this run.",
    ""
  )

  writeLines(out, path)
  message("wrote ", path)
  invisible(path)
}
