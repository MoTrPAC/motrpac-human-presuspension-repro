# Stage 3 carry, step 6 — run each package's own test suite against the carried
# payload, after updating the expectations the new payload legitimately breaks.
#
# Two test edits are made, and the difference between them matters. One is an
# expectation that pinned a row count which the rebuild changed — the test was
# right about the old object and has to be re-pinned. The other is a new
# assertion: the DA tables lost two columns, and nothing pinned their absence, so
# the removal could silently reverse. An expectation that is only updated, never
# added, lets a payload change erase its own evidence.
#
# Tests run through pkgload::load_all() on the test package. Cross-package calls
# still resolve against whatever is in .libPaths(), so a test in the analysis
# package that reaches into the data package exercises the *installed* data
# package, not the test one. That is recorded per suite rather than worked
# around, because installing an 800 MB pair to close the gap costs more than the
# gap does.
#
# Usage:
#   Rscript 06_check_and_test.R --out-root <dir> --out <dir> [--check 1]

.here <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
source(file.path(.here, "lib", "carry_helpers.R"))

opt <- parse_opts(commandArgs(trailingOnly = TRUE), c("out-root", "out", "check"))
for (need in c("out-root", "out")) if (is.null(opt[[need]])) stop("missing required option --", need)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
rpt <- new_report()

PKG <- c(data = "MotrpacHumanPreSuspensionData",
         analysis = "MotrpacHumanPreSuspensionAnalysis")
root <- function(dest) file.path(opt$`out-root`, PKG[[dest]])

# ---- Re-pin the expectations the new payload changes -------------------------

tf <- file.path(root("analysis"), "tests", "testthat", "test-data_objects.R")
if (!file.exists(tf)) {
  record(rpt, "FAIL", "test:data_objects", "test-data_objects.R not found")
} else {
  lines <- readLines(tf, warn = FALSE)
  outliers <- load_rda_object(file.path(root("analysis"), "data", "OUTLIERS.rda"), "OUTLIERS")

  for (spec in list(list(pat = "^\\s*expect_equal\\(nrow\\(OUTLIERS\\), [0-9]+L\\)\\s*$",
                         new = sprintf("  expect_equal(nrow(OUTLIERS), %dL)", nrow(outliers)),
                         what = "nrow"),
                    list(pat = "^\\s*expect_equal\\(ncol\\(OUTLIERS\\), [0-9]+L\\)\\s*$",
                         new = sprintf("  expect_equal(ncol(OUTLIERS), %dL)", ncol(outliers)),
                         what = "ncol"))) {
    i <- grep(spec$pat, lines)
    if (!length(i)) {
      record(rpt, "WARN", paste0("test:OUTLIERS:", spec$what), "expectation not found; left alone")
    } else if (identical(trimws(lines[i[1]]), trimws(spec$new))) {
      record(rpt, "PASS", paste0("test:OUTLIERS:", spec$what), "already matches the carried object")
    } else {
      record(rpt, "PASS", paste0("test:OUTLIERS:", spec$what),
             sprintf("%s -> %s", trimws(lines[i[1]]), trimws(spec$new)))
      lines[i[1]] <- spec$new
    }
  }
  writeLines(lines, tf)
}

# A regression test for the removal itself. Without it, the two columns could
# come back and every existing assertion would still pass.
da_test <- file.path(root("analysis"), "tests", "testthat", "test-load_differential_analysis.R")
if (file.exists(da_test)) {
  lines <- readLines(da_test, warn = FALSE)
  # Matched on the bare column name only: the DA tables now carry CI.L_calculated, and a
  # fixed "CI.L" substring would see that as the regression test already being present.
  if (any(grepl("CI\\.L(?!_calculated)", lines, perl = TRUE))) {
    record(rpt, "PASS", "test:CI-columns", "regression test already present")
  } else {
    writeLines(c(lines, "",
      "test_that(\"differential analysis results no longer carry the CI.L/CI.R columns\", {",
      "  # Dropped for v2.0. Pinned so the removal cannot silently reverse.",
      "  expect_false(any(c(\"CI.L\", \"CI.R\") %in% colnames(ADIPOSE_PROT_PR_DA)))",
      "  expect_false(any(c(\"CI.L\", \"CI.R\") %in% colnames(MUSCLE_TRNSCRPT_DA)))",
      "})"), da_test)
    record(rpt, "PASS", "test:CI-columns", "added a regression test for the removed columns")
  }
}

# ---- Run the suites ----------------------------------------------------------

if (!requireNamespace("testthat", quietly = TRUE) || !requireNamespace("pkgload", quietly = TRUE)) {
  record(rpt, "FAIL", "test:deps", "testthat and pkgload are required")
} else {
  all_results <- list()
  for (dest in names(PKG)) {
    p <- root(dest)
    res <- tryCatch(
      testthat::test_local(p, reporter = testthat::SilentReporter$new(), stop_on_failure = FALSE),
      error = function(e) e)

    if (inherits(res, "error")) {
      record(rpt, "FAIL", paste0("test:", PKG[[dest]]), conditionMessage(res))
      next
    }

    df <- as.data.frame(res)
    df$package <- PKG[[dest]]
    all_results[[dest]] <- df[, c("package", "file", "test", "nb", "failed", "skipped",
                                  "error", "warning")]

    n_fail <- sum(df$failed) + sum(df$error)
    if (n_fail > 0L) {
      bad <- df[df$failed > 0 | df$error, ]
      record(rpt, "FAIL", paste0("test:", PKG[[dest]]),
             sprintf("%d failing expectation(s) across %d test(s): %s",
                     sum(df$failed), nrow(bad), paste(bad$test, collapse = "; ")))
    } else {
      record(rpt, "PASS", paste0("test:", PKG[[dest]]),
             sprintf("%d expectation(s) in %d test(s) passed, %d skipped",
                     sum(df$nb), nrow(df), sum(df$skipped)))
    }
  }
  if (length(all_results)) {
    write_tsv(do.call(rbind, all_results), file.path(opt$out, "test_results.tsv"))
  }
}

# ---- Optional R CMD check ----------------------------------------------------
# Off by default: building the tarball copies ~800 MB of data twice, and the
# structural findings it would surface are already covered by the Stage 3 audit.

if (!is.null(opt$check) && opt$check %in% c("1", "TRUE", "true")) {
  for (dest in names(PKG)) {
    p <- root(dest)
    lib <- tempfile("checklib"); dir.create(lib)
    args <- c("CMD", "check", "--no-manual", "--no-build-vignettes", "--no-vignettes",
              "--no-examples", "--no-tests", "--output", dirname(p), shQuote(p))
    out <- suppressWarnings(system2(file.path(R.home("bin"), "R"), args,
                                    stdout = TRUE, stderr = TRUE))
    status <- attr(out, "status")
    errs <- grep("^(ERROR|WARNING)", out, value = TRUE)
    if (!is.null(status) && status != 0) {
      record(rpt, "FAIL", paste0("check:", PKG[[dest]]),
             paste(utils::head(errs, 3), collapse = " | "))
    } else {
      record(rpt, "PASS", paste0("check:", PKG[[dest]]),
             sprintf("%d NOTE(s), no ERROR or WARNING", sum(grepl("^NOTE", out))))
    }
    writeLines(out, file.path(opt$out, paste0("check_", PKG[[dest]], ".log")))
  }
}

write_tsv(report_frame(rpt), file.path(opt$out, "check_test_report.tsv"))
message(sprintf("\n%d check(s): %d FAIL, %d WARN", length(rpt$rows), rpt$fails, rpt$warns))
if (rpt$fails > 0L) quit(status = 1L)
