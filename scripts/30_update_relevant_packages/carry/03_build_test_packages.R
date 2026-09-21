# Stage 3 carry, step 3 — build a test package per repo and carry the objects in.
#
# The real package checkouts are never written to. Each test package is a copy of
# its source tree with a fresh data/ assembled from the routing table, so the
# carry can be documented, checked and tested end to end and then thrown away.
#
# Two things happen here that a plain copy would get wrong.
#
# Compression. Both packages store data/ bzip2-compressed; the pipeline writes
# gzip (`save(..., compress = TRUE)`). Copying the files across would grow both
# packages and earn the R CMD check NOTE that recommends exactly this resave. So
# every object is written out through save(compress = "bzip2"), not copied.
#
# Comparison. Because of that same compression difference, every rebuilt object
# differs from its packaged counterpart at the byte level whether or not its
# content moved — the raw diff reports 106 changes where two of them are pure
# recompression. Verdicts here are computed on the loaded objects.
#
# Usage:
#   Rscript 03_build_test_packages.R --routing <tsv> --data-pkg <dir> \
#     --analysis-pkg <dir> --out-root <dir> --out <dir>
#
# Writes:
#   <out-root>/<Package>/          the test package
#   <out>/carry_manifest.tsv       one row per object, with a content verdict
#   <out>/carry_report.tsv

.here <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
source(file.path(.here, "lib", "carry_helpers.R"))

opt <- parse_opts(commandArgs(trailingOnly = TRUE),
                  c("routing", "data-pkg", "analysis-pkg", "out-root", "out"))
for (need in c("routing", "data-pkg", "analysis-pkg", "out-root", "out")) {
  if (is.null(opt[[need]])) stop("missing required option --", need)
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
rpt <- new_report()

routing <- read_tsv(opt$routing)
src_pkg <- c(data = opt$`data-pkg`, analysis = opt$`analysis-pkg`)

# ---- Build the test package trees -------------------------------------------
# Everything except data/ is copied verbatim; data/ is assembled object by
# object below so that what the test package ships is exactly what the routing
# table says, with nothing inherited by accident.
#
# .git is excluded deliberately: a test package is a build product, not a
# checkout, and copying the history would invite committing into it.

SKIP_TOP <- c(".git", ".Rproj.user", "data", "sandbox", ".Rhistory")

test_pkg <- character()
for (dest in names(src_pkg)) {
  from <- src_pkg[[dest]]
  pkg  <- basename(from)
  to   <- file.path(opt$`out-root`, pkg)

  # out-root comes from TEST_PKG_ROOT, which a caller can set to anything, and the unlink
  # below is recursive. Everything else in this stage is safe because it only ever writes
  # under staging/; this is the one line that would take a caller at their word and delete
  # a real checkout, so the target has to prove it is a build directory first. A git
  # working tree never is one.
  if (identical(normalizePath(to, mustWork = FALSE), normalizePath(from, mustWork = FALSE)))
    stop("--out-root would build ", pkg, " on top of its own source: ", to)
  if (dir.exists(file.path(to, ".git")))
    stop("refusing to delete ", to, " — it is a git checkout, not a test-package build ",
         "directory. Check TEST_PKG_ROOT.")

  unlink(to, recursive = TRUE)
  dir.create(to, recursive = TRUE, showWarnings = FALSE)
  for (entry in setdiff(list.files(from, all.files = TRUE, no.. = TRUE), SKIP_TOP)) {
    file.copy(file.path(from, entry), to, recursive = TRUE, copy.date = TRUE)
  }
  dir.create(file.path(to, "data"), showWarnings = FALSE)
  test_pkg[dest] <- to
  record(rpt, "PASS", paste0("tree:", pkg), paste("built at", to))
}

# ---- Carry the objects -------------------------------------------------------

rows <- list()
add_row <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)

carry <- routing[routing$action %in% c("replace", "add"), , drop = FALSE]
carry <- carry[order(carry$destination, carry$object), ]

t0 <- Sys.time()
for (i in seq_len(nrow(carry))) {
  object <- carry$object[i]
  dest   <- carry$destination[i]
  path   <- carry$source_path[i]

  value <- if (identical(carry$source_format[i], "rds")) readRDS(path) else load_rda_object(path, object)

  old_path <- file.path(src_pkg[[dest]], "data", paste0(object, ".rda"))
  old <- if (file.exists(old_path)) load_rda_object(old_path, object) else NULL
  cmp <- compare_object(value, old)

  out_path <- save_package_object(value, object, file.path(test_pkg[[dest]], "data"))

  add_row(object = object, destination = dest, action = carry$action[i],
          source_stage = carry$source_stage[i], verdict = cmp$verdict, detail = cmp$detail,
          bytes_before = if (file.exists(old_path)) file.size(old_path) else NA_integer_,
          bytes_after = file.size(out_path))

  rm(value, old); gc(FALSE)
  if (i %% 25L == 0L) {
    message(sprintf("  ... %d/%d objects (%.1f min elapsed)", i, nrow(carry),
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

# ---- Carry the orphans verbatim ---------------------------------------------
# An object nothing rebuilt keeps exactly the bytes the package already ships.
# Resaving it would produce a diff with no change behind it.

orph <- routing[routing$action == "orphan-intended", , drop = FALSE]
for (i in seq_len(nrow(orph))) {
  dest <- orph$destination[i]
  from <- file.path(src_pkg[[dest]], "data", paste0(orph$object[i], ".rda"))
  to   <- file.path(test_pkg[[dest]], "data", paste0(orph$object[i], ".rda"))
  file.copy(from, to, overwrite = TRUE, copy.date = TRUE)
  add_row(object = orph$object[i], destination = dest, action = "orphan-intended",
          source_stage = NA_character_, verdict = "CARRIED-FORWARD", detail = orph$reason[i],
          bytes_before = file.size(from), bytes_after = file.size(to))
}

# ---- Withdraw the dropped orphans -------------------------------------------
# Nothing is copied: data/ was assembled from scratch, so an object is dropped by
# not writing it. The row is still recorded, because a release that silently
# ships one fewer object is the thing worth being able to see.

drop <- routing[routing$action == "orphan-dropped", , drop = FALSE]
for (i in seq_len(nrow(drop))) {
  dest <- drop$destination[i]
  from <- file.path(src_pkg[[dest]], "data", paste0(drop$object[i], ".rda"))
  add_row(object = drop$object[i], destination = dest, action = "orphan-dropped",
          source_stage = NA_character_, verdict = "REMOVED", detail = drop$reason[i],
          bytes_before = if (file.exists(from)) file.size(from) else NA_integer_,
          bytes_after = NA_integer_)
}

manifest <- do.call(rbind, rows)
write_tsv(manifest, file.path(opt$out, "carry_manifest.tsv"))

# ---- Verify the result -------------------------------------------------------

for (dest in names(test_pkg)) {
  pkg <- basename(test_pkg[[dest]])
  want <- sort(manifest$object[manifest$destination == dest & manifest$verdict != "REMOVED"])
  have <- sort(sub("\\.rda$", "", list.files(file.path(test_pkg[[dest]], "data"), pattern = "\\.rda$")))
  if (identical(want, have)) {
    record(rpt, "PASS", paste0("data:", pkg),
           sprintf("%d object(s) in data/, matching the routing table", length(have)))
  } else {
    record(rpt, "FAIL", paste0("data:", pkg),
           sprintf("routed %d but data/ holds %d; missing: %s; unexpected: %s",
                   length(want), length(have),
                   paste(setdiff(want, have), collapse = ","),
                   paste(setdiff(have, want), collapse = ",")))
  }

  # Every object must load back under its own name, or LazyData will not build.
  bad <- character()
  for (o in have) {
    ok <- tryCatch({ load_rda_object(file.path(test_pkg[[dest]], "data", paste0(o, ".rda")), o); TRUE },
                   error = function(e) FALSE)
    if (!ok) bad <- c(bad, o)
  }
  if (length(bad)) {
    record(rpt, "FAIL", paste0("loadable:", pkg),
           sprintf("%d object(s) do not load under their own name: %s",
                   length(bad), paste(bad, collapse = ", ")))
  } else {
    record(rpt, "PASS", paste0("loadable:", pkg), "every object loads under its own name")
  }

  comp <- tools::checkRdaFiles(file.path(test_pkg[[dest]], "data"))
  wrong <- rownames(comp)[comp$compress != "bzip2"]
  if (length(wrong)) {
    record(rpt, "WARN", paste0("compress:", pkg),
           sprintf("%d file(s) not bzip2 (carried-forward orphans keep their original compression)",
                   length(wrong)))
  } else {
    record(rpt, "PASS", paste0("compress:", pkg), "every data file is bzip2")
  }
}

for (v in names(table(manifest$verdict))) {
  record(rpt, "PASS", paste0("verdict:", v), sprintf("%d object(s)", table(manifest$verdict)[[v]]))
}

write_tsv(report_frame(rpt), file.path(opt$out, "carry_report.tsv"))
message(sprintf("\n%d check(s): %d FAIL, %d WARN — %.1f min",
                length(rpt$rows), rpt$fails, rpt$warns,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (rpt$fails > 0L) quit(status = 1L)
