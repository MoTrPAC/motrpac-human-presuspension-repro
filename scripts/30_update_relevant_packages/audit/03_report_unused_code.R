# Stage 3, step 3 — inventory of code that is no longer needed.
#
# For every function defined in R/ and every object in data/, count call sites in
# three widening scopes and classify from the counts:
#
#   own package    R/, tests/, vignettes/ of the defining package
#   consumers      the repos passed via --consumer (analysis repos, this pipeline)
#   exported?      whether NAMESPACE exports the name
#
# Classification, in priority order:
#   DEPRECATED    roxygen or a comment marks it deprecated/superseded/legacy
#   LIVE          referenced outside its own definition file, anywhere
#   API-ONLY      exported, but no reference anywhere we can see. Removing it is
#                 an API break for users outside these repos, so it is never
#                 reported as removable on its own.
#   DEAD          not exported and referenced nowhere outside its definition file
#
# Only DEAD is a removal candidate without further judgement. That conservatism
# is the point: a false DEAD verdict deletes working code, and the search scope
# here is the filesystem, not the world.
#
# Usage:
#   Rscript 03_report_unused_code.R <pkg_root> [<pkg_root> ...] \
#       [--consumer <repo>]... --out <dir>

script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
source(file.path(script_dir(), "lib", "pkg_introspect.R"))

# ---- Arguments -------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
out_dir <- "."
consumers <- character(0)

# Flags and positional roots may appear in any order.
keep <- character(0)
while (length(args)) {
  if (args[1] == "--out") {
    out_dir <- args[2]; args <- args[-(1:2)]
  } else if (args[1] == "--consumer") {
    consumers <- c(consumers, args[2]); args <- args[-(1:2)]
  } else {
    keep <- c(keep, args[1]); args <- args[-1]
  }
}
pkg_roots <- normalizePath(keep, mustWork = FALSE)
consumers <- normalizePath(consumers, mustWork = FALSE)
consumers <- consumers[dir.exists(consumers)]
if (!length(pkg_roots)) stop("usage: 03_report_unused_code.R <pkg_root> ... [--consumer <repo>]... --out <dir>")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Source corpora --------------------------------------------------------

# Read every text file under `root` once. Grepping in memory beats spawning a
# process per name — there are hundreds of names and a handful of repos.
# Paths that hold a COPY of a package rather than code that uses one. A vendored
# snapshot contains every definition verbatim, so counting it as consumer evidence
# marks the entire package live and the audit reports nothing. Worktrees are
# duplicate checkouts of the repo being audited, with the same effect.
VENDORED <- "(^|/)(staging/pkgs-[^/]*|\\.claude/worktrees|renv|packrat|node_modules)(/|$)"

corpus <- function(root, subdirs = NULL) {
  base <- if (is.null(subdirs)) root else file.path(root, subdirs)
  base <- base[dir.exists(base)]
  if (!length(base)) return(list())
  files <- unlist(lapply(base, list.files, recursive = TRUE, full.names = TRUE))
  files <- files[grepl("\\.(R|r|Rmd|rmd|qmd|Rnw|sh|json|yml|yaml|md)$", files)]
  files <- files[!grepl(VENDORED, files)]
  info <- file.info(files)
  files <- files[!is.na(info$size) & !info$isdir & info$size < 5e6]
  out <- lapply(files, function(f) tryCatch(readLines(f, warn = FALSE),
                                            error = function(e) character(0)))
  names(out) <- files
  out
}

# Hit count for `name` bounded by R-identifier characters.
#
# \b cannot be used here. Most internals are dot-prefixed, and in `pkg:::.find_ome`
# the character before the dot is ':' — both non-word, so \b matches nothing and
# every dot-prefixed helper silently reports zero hits. Explicit lookaround on the
# R identifier set (letters, digits, dot, underscore) is what actually works.
count_hits <- function(name, corp, exclude = character(0), qualify = NULL) {
  if (!length(corp)) return(0L)
  esc <- gsub("([.\\\\])", "\\\\\\1", name)
  # A package-private name cannot be reached from another repo except through
  # `Pkg:::name`. Without this, an identically-named private helper in a sibling
  # package (both zzz.R files define .block_set_env) counts as external usage and
  # keeps genuinely unused code alive in the report.
  pat <- if (is.null(qualify)) {
    paste0("(?<![A-Za-z0-9._])", esc, "(?![A-Za-z0-9._])")
  } else {
    paste0(qualify, ":::?", esc, "(?![A-Za-z0-9._])")
  }
  keep <- setdiff(names(corp), exclude)
  sum(vapply(keep, function(f) {
    length(grep(pat, corp[[f]], perl = TRUE))
  }, integer(1)))
}

message("consumer repos: ", if (length(consumers)) paste(basename(consumers), collapse = ", ") else "(none)")
consumer_corpora <- lapply(consumers, corpus)
names(consumer_corpora) <- basename(consumers)

# ---- Per package -----------------------------------------------------------

deprecation_marks <- "deprecat|supersed|legacy|no longer used|obsolete|do not use|TODO:? *remove"

all_rows <- list()

for (pkg_root in pkg_roots) {
  if (!dir.exists(pkg_root)) {
    message("[skip] not a directory: ", pkg_root)
    next
  }
  desc <- read_description(pkg_root)
  pkg <- desc[["Package"]]
  lazy_data <- "LazyData" %in% names(desc) &&
    grepl("true|yes", desc[["LazyData"]], ignore.case = TRUE)
  ns <- parse_namespace(pkg_root)
  # Does this package reach its own data objects by constructed name? load_qc(),
  # load_clinical_data(), load_summary_stats() and load_differential_analysis()
  # all enumerate data() and eval(parse(text = nm)). Where that is true, NO name
  # search can decide whether a data object is used — the call site never spells
  # the name out. Saying so beats emitting 70 confident false verdicts.
  r_text_all <- unlist(lapply(pkg_r_files(pkg_root), readLines, warn = FALSE))
  dynamic_data <- any(grepl("eval\\(parse\\(text|utils::data\\(|\\bdata\\(package\\s*=",
                            r_text_all))
  defs <- pkg_function_defs(pkg_root)
  data_objs <- pkg_data_objects(pkg_root)

  own_r <- corpus(pkg_root, "R")
  # tests/vignettes/inst CONSUME the package. data-raw/ PRODUCES it, so a hit
  # there is the generator writing the object, not evidence anyone reads it —
  # counting it as usage would mark every shipped object live by construction.
  own_other <- corpus(pkg_root, c("tests", "vignettes", "inst"))
  own_dataraw <- corpus(pkg_root, "data-raw")
  # sandbox/ is scratch space; a hit there does not keep code alive.
  own_sandbox <- corpus(pkg_root, "sandbox")

  # Intra-package function calls come from the parse tree, not from grep. A helper
  # is very often called by the exported function directly above it in the SAME
  # file, so excluding the defining file to avoid counting the definition would
  # bury every one of those call sites and report live helpers as dead.
  use_counts <- pkg_symbol_use_counts(pkg_root)
  # S3 methods are reached by dispatch, never by name, so NAMESPACE registration
  # is their only call site. `S3method(print,motrdat)` registers `print.motrdat`.
  s3_methods <- gsub(",", ".", ns$s3methods, fixed = TRUE)
  entry_points <- c(ns$exports, s3_methods,
                    ".onLoad", ".onUnload", ".onAttach", ".onDetach")

  # Sibling packages count as consumers of each other.
  sibling_corpora <- consumer_corpora
  for (other in setdiff(pkg_roots, pkg_root)) {
    if (dir.exists(other)) {
      sibling_corpora[[basename(other)]] <- corpus(other, c("R", "tests", "vignettes", "data-raw"))
    }
  }

  message(sprintf("%s: %d functions, %d data objects", pkg, nrow(defs), length(data_objs)))

  classify <- function(name, kind, def_file) {
    in_own_r <- if (name %in% names(use_counts)) use_counts[[name]] else 0L
    in_own_other <- count_hits(name, own_other)
    in_dataraw <- count_hits(name, own_dataraw)
    in_sandbox <- count_hits(name, own_sandbox)
    exported <- name %in% entry_points
    # Exported names and data objects are reachable by bare name after library();
    # private helpers only through Pkg:::name.
    qualify <- if (exported || identical(kind, "data")) NULL else pkg
    in_consumers <- sum(vapply(sibling_corpora,
                               function(cp) count_hits(name, cp, qualify = qualify),
                               integer(1)))
    # Bare mentions of a private helper in another repo. Not a legal call, so it
    # does not make the helper live — but it is not nothing either: either the
    # caller is broken, or the helper was meant to be exported.
    bare_consumers <- if (is.null(qualify)) in_consumers else
      sum(vapply(sibling_corpora, function(cp) count_hits(name, cp), integer(1)))

    # Deprecation marks: look at the definition file only, near the definition.
    deprecated <- FALSE
    if (!is.na(def_file) && nzchar(def_file) && file.exists(def_file)) {
      txt <- tryCatch(readLines(def_file, warn = FALSE), error = function(e) character(0))
      deprecated <- any(grepl(deprecation_marks, txt, ignore.case = TRUE, perl = TRUE))
    }

    status <- if (deprecated) {
      "DEPRECATED"
    } else if (in_own_r + in_own_other + in_consumers > 0L) {
      "LIVE"
    } else if (identical(kind, "data") && dynamic_data) {
      "DYNAMIC-ACCESS"
    } else if (in_dataraw > 0L) {
      # Reached only by the generators. Not dead — but it ships to users who can
      # never reach it, so it belongs in data-raw/ rather than in R/.
      "BUILD-ONLY"
    } else if (!exported && bare_consumers > 0L) {
      "CALLED-UNEXPORTED"
    } else if (identical(kind, "data") && lazy_data) {
      # Under LazyData every documented object is reachable as `pkg::obj` by any
      # user. Shipping data IS the product of a data package, so "nothing in these
      # repos greps the name" says nothing about whether it is needed.
      "SHIPPED-DATA"
    } else if (exported) {
      "API-ONLY"
    } else {
      "DEAD"
    }
    data.frame(
      package = pkg, name = name, kind = kind, exported = exported,
      status = status,
      hits_own_R = in_own_r, hits_own_tests_vignettes = in_own_other,
      hits_data_raw = in_dataraw, hits_sandbox = in_sandbox,
      hits_consumers = in_consumers, hits_consumers_bare = bare_consumers,
      defined_in = if (is.na(def_file)) NA_character_ else
        sub(paste0("^", pkg_root, "/?"), "", def_file),
      stringsAsFactors = FALSE)
  }

  fn_rows <- lapply(seq_len(nrow(defs)), function(i) {
    classify(defs$name[i], "function", defs$file[i])
  })
  data_rows <- lapply(data_objs, function(o) {
    classify(o, "data", NA_character_)
  })
  all_rows <- c(all_rows, fn_rows, data_rows)

  # Duplicate definitions of one name are a maintenance hazard regardless of use.
  dup <- defs$name[duplicated(defs$name)]
  if (length(dup)) {
    message(sprintf("  [WARN] %s: %d name(s) defined more than once: %s",
                    pkg, length(unique(dup)), paste(unique(dup), collapse = ", ")))
  }
}

inv <- do.call(rbind, all_rows)
if (is.null(inv)) {
  inv <- data.frame(package = character(0), name = character(0), kind = character(0),
                    exported = logical(0), status = character(0))
}
inv <- inv[order(inv$package,
                 factor(inv$status, levels = c("DEAD", "DEPRECATED", "API-ONLY", "LIVE")),
                 inv$name), ]
write_tsv(inv, file.path(out_dir, "unused_code_inventory.tsv"))

# ---- Commented-out code blocks ---------------------------------------------
# Runs of >= 5 consecutive comment lines whose text parses as code rather than
# prose. Roxygen (#') is documentation and is excluded.

block_rows <- list()
for (pkg_root in pkg_roots) {
  if (!dir.exists(pkg_root)) next
  pkg <- read_description(pkg_root)[["Package"]]
  for (f in pkg_r_files(pkg_root)) {
    txt <- readLines(f, warn = FALSE)
    is_code_comment <- grepl("^\\s*#(?!')", txt, perl = TRUE) &
      grepl("(<-|\\(|\\)|=|\\{|\\})", sub("^\\s*#+", "", txt))
    r <- rle(is_code_comment)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    sel <- which(r$values & r$lengths >= 5L)
    for (k in sel) {
      block_rows[[length(block_rows) + 1L]] <- data.frame(
        package = pkg,
        file = sub(paste0("^", pkg_root, "/?"), "", f),
        start = starts[k], end = ends[k], n_lines = r$lengths[k],
        first_line = trimws(txt[starts[k]]),
        stringsAsFactors = FALSE)
    }
  }
}
blocks <- do.call(rbind, block_rows)
if (is.null(blocks)) {
  blocks <- data.frame(package = character(0), file = character(0), start = integer(0),
                       end = integer(0), n_lines = integer(0), first_line = character(0))
}
blocks <- blocks[order(-blocks$n_lines), ]
write_tsv(blocks, file.path(out_dir, "unused_code_commented_blocks.tsv"))

# ---- Summary ---------------------------------------------------------------

summary_tab <- as.data.frame(table(inv$package, inv$status), stringsAsFactors = FALSE)
names(summary_tab) <- c("package", "status", "n")
summary_tab <- summary_tab[summary_tab$n > 0L, ]
write_tsv(summary_tab, file.path(out_dir, "unused_code_summary.tsv"))

for (i in seq_len(nrow(summary_tab))) {
  message(sprintf("  %-40s %-12s %d", summary_tab$package[i],
                  summary_tab$status[i], summary_tab$n[i]))
}
message(sprintf("commented-out code blocks (>=5 lines): %d", nrow(blocks)))
message("wrote inventory to ", file.path(out_dir, "unused_code_inventory.tsv"))

# Reporting dead code is never itself a failure — the stage decides what to do.
quit(status = 0L)
