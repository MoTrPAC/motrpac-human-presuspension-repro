# Shared helpers for the Stage 3 carry steps.
#
# Every step is a standalone Rscript that reads a TSV written by the previous one
# and writes its own, so a step can be re-run in isolation and the intermediate
# state is inspectable rather than held in memory.

# Rscript does not set a script-directory variable, so recover it from the
# --file= argument the front end passes through.
script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}

# ---- Argument handling -----------------------------------------------------

#' Pull `--name value` pairs out of the command line
#'
#' @param args character; `commandArgs(trailingOnly = TRUE)`.
#' @param names character; the option names to look for, without the `--`.
#' @return a named list; missing options are NULL.
parse_opts <- function(args, names) {
  out <- setNames(vector("list", length(names)), names)
  for (n in names) {
    flag <- paste0("--", n)
    i <- which(args == flag)
    if (length(i)) out[[n]] <- args[i[1] + 1L]
  }
  out
}

# ---- TSV I/O ---------------------------------------------------------------
# read.csv/write.table per the repo convention, with check.names off because
# object names and column names here are data, not syntactic R names.

read_tsv <- function(path) {
  read.csv(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
           quote = "", comment.char = "")
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(path)
}

# ---- Reporting -------------------------------------------------------------

new_report <- function() {
  e <- new.env(parent = emptyenv())
  e$rows <- list()
  e$fails <- 0L
  e$warns <- 0L
  e
}

#' Record one check
#'
#' @param rpt the accumulator from `new_report()`.
#' @param status character; one of "PASS", "WARN", "FAIL".
#' @param check,detail character; what was checked and what was found.
record <- function(rpt, status, check, detail = "") {
  status <- match.arg(status, c("PASS", "WARN", "FAIL"))
  rpt$rows[[length(rpt$rows) + 1L]] <-
    data.frame(status = status, check = check, detail = detail,
               stringsAsFactors = FALSE)
  if (status == "FAIL") rpt$fails <- rpt$fails + 1L
  if (status == "WARN") rpt$warns <- rpt$warns + 1L
  message(sprintf("[%-4s] %s%s", status, check,
                  if (nzchar(detail)) paste0(" — ", detail) else ""))
  invisible(rpt)
}

report_frame <- function(rpt) {
  if (!length(rpt$rows)) {
    return(data.frame(status = character(), check = character(),
                      detail = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, rpt$rows)
}

# ---- Git -------------------------------------------------------------------
# Read-only interrogation of a package checkout. Every function returns NA
# rather than erroring when the path is not a git repo, so the caller decides
# whether that is fatal.

git_out <- function(repo, ...) {
  res <- suppressWarnings(system2("git", c("-C", shQuote(repo), ...),
                                  stdout = TRUE, stderr = FALSE))
  if (!is.null(attr(res, "status")) && attr(res, "status") != 0) return(NA_character_)
  if (!length(res)) return("")
  paste(res, collapse = "\n")
}

is_git_repo <- function(repo) identical(git_out(repo, "rev-parse", "--is-inside-work-tree"), "true")
git_branch  <- function(repo) git_out(repo, "rev-parse", "--abbrev-ref", "HEAD")
git_sha     <- function(repo) git_out(repo, "rev-parse", "--short", "HEAD")
git_dirty   <- function(repo) {
  s <- git_out(repo, "status", "--porcelain")
  if (is.na(s) || !nzchar(s)) character() else strsplit(s, "\n", fixed = TRUE)[[1]]
}

# ---- DESCRIPTION -----------------------------------------------------------

read_description <- function(pkg_root) {
  dcf <- file.path(pkg_root, "DESCRIPTION")
  if (!file.exists(dcf)) return(NULL)
  as.list(read.dcf(dcf)[1, ])
}

#' Rewrite a single DESCRIPTION field in place, preserving every other line
#'
#' `write.dcf()` would reflow continuation lines across the whole file, which
#' turns a one-field edit into an unreviewable diff. This edits the field's own
#' line and leaves the rest of the file byte-for-byte.
set_description_field <- function(pkg_root, field, value) {
  path <- file.path(pkg_root, "DESCRIPTION")
  lines <- readLines(path, warn = FALSE)
  i <- grep(paste0("^", field, ":"), lines)
  if (!length(i)) stop("DESCRIPTION has no ", field, " field: ", path)
  lines[i[1]] <- paste0(field, ": ", value)
  writeLines(lines, path)
  invisible(value)
}

# ---- Data objects ----------------------------------------------------------

#' Load the single object out of an .rda, by name
load_rda_object <- function(path, object) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!object %in% ls(e, all.names = TRUE)) {
    stop("'", object, "' not found in ", path, "; contains: ",
         paste(ls(e, all.names = TRUE), collapse = ", "))
  }
  get(object, envir = e)
}

#' Compare a rebuilt object against the one a package currently ships
#'
#' Byte comparison is not usable here: the pipeline writes gzip and the packages
#' store bzip2, so every object differs on disk whether or not its content moved.
#'
#' @return a list with `verdict` and `detail`. `verdict` is one of IDENTICAL,
#'   SCHEMA-CHANGE, VALUES-DIFFER, ADDED.
compare_object <- function(new, old) {
  if (is.null(old)) return(list(verdict = "ADDED", detail = describe_object(new)))
  if (identical(new, old)) return(list(verdict = "IDENTICAL", detail = ""))

  if (is.data.frame(new) && is.data.frame(old)) {
    gained <- setdiff(names(new), names(old))
    lost   <- setdiff(names(old), names(new))
    detail <- sprintf("%s -> %s", dim_str(old), dim_str(new))
    if (length(gained)) detail <- paste0(detail, "; +cols: ", paste(gained, collapse = ","))
    if (length(lost))   detail <- paste0(detail, "; -cols: ", paste(lost, collapse = ","))
    verdict <- if (length(gained) || length(lost)) "SCHEMA-CHANGE" else "VALUES-DIFFER"
    return(list(verdict = verdict, detail = detail))
  }

  if (is.list(new) && is.list(old)) {
    gained <- setdiff(names(new), names(old))
    lost   <- setdiff(names(old), names(new))
    detail <- sprintf("list %d -> %d elements", length(old), length(new))
    if (length(gained) || length(lost)) {
      detail <- paste0(detail, "; +", length(gained), " -", length(lost), " names")
      return(list(verdict = "SCHEMA-CHANGE", detail = detail))
    }
    # A QC object is list(qc_norm, feature_metadata, sample_metadata). Its element names do
    # not move when a component gains or loses a COLUMN, so comparing names alone reported
    # VALUES-DIFFER and NEWS -- which reads +cols/-cols off this detail -- said nothing about
    # the change. Recurse one level into the data-frame components, naming each one so two
    # components changing differently stay distinguishable.
    parts <- character(0)
    for (nm in names(new)) {
      a <- new[[nm]]; b <- old[[nm]]
      if (!is.data.frame(a) || !is.data.frame(b)) next
      g <- setdiff(names(a), names(b)); l <- setdiff(names(b), names(a))
      if (!length(g) && !length(l)) next
      part <- nm
      if (length(g)) part <- paste0(part, ": +cols: ", paste(g, collapse = ","))
      if (length(l)) part <- paste0(part, ": -cols: ", paste(l, collapse = ","))
      parts <- c(parts, part)
    }
    if (length(parts))
      return(list(verdict = "SCHEMA-CHANGE",
                  detail = paste0(detail, "; ", paste(parts, collapse = "; "))))
    return(list(verdict = "VALUES-DIFFER", detail = detail))
  }

  list(verdict = "VALUES-DIFFER",
       detail = sprintf("%s -> %s", describe_object(old), describe_object(new)))
}

dim_str <- function(x) if (is.null(dim(x))) paste0("length ", length(x)) else paste(dim(x), collapse = "x")

describe_object <- function(x) {
  paste0(paste(class(x), collapse = "/"), " ", dim_str(x))
}

#' Write an object into a package `data/` directory the way the package stores it
#'
#' Both packages hold their `data/` bzip2-compressed; the pipeline writes gzip.
#' Copying without resaving grows the package and earns an R CMD check NOTE
#' recommending exactly this. Serialization stays at version 3, which both
#' packages already use uniformly.
save_package_object <- function(value, object, data_dir) {
  path <- file.path(data_dir, paste0(object, ".rda"))
  e <- new.env(parent = emptyenv())
  assign(object, value, envir = e)
  save(list = object, envir = e, file = path, compress = "bzip2", version = 3)
  path
}
