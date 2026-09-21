# Stage 3 carry, step 1 — validate the two package checkouts and the build we are
# about to carry into them.
#
# Nothing is written to either package repo, here or anywhere else in the carry:
# step 3 builds *test packages* under staging/ from these trees. That changes what
# each gate is for but not whether it is worth having — the test package inherits
# the branch, the version and the uncommitted edits of whatever it was copied
# from, so an unrecorded source state produces an unattributable result.
#
# Usage:
#   Rscript 01_validate_targets.R --data-pkg <dir> --analysis-pkg <dir> \
#     --build-data <dir> --preflight-data <dir> --out <dir>
#
# Writes:
#   <out>/validate_targets.tsv   one row per check
#   <out>/source_state.tsv       branch, sha, version and dirty-file count per package
# Exits non-zero on any FAIL.

.here <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
source(file.path(.here, "lib", "carry_helpers.R"))

opt <- parse_opts(commandArgs(trailingOnly = TRUE),
                  c("data-pkg", "analysis-pkg", "build-data", "preflight-data", "out"))
for (need in c("data-pkg", "analysis-pkg", "build-data", "out")) {
  if (is.null(opt[[need]])) stop("missing required option --", need)
}
out_dir <- opt$out
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rpt <- new_report()

# The versions the v1.3 comparison baseline was pinned to. A checkout older than its own baseline would mean
# the packages moved backwards since the comparison ran, and every verdict in
# reports/ would describe a tree that no longer exists.
BASELINE <- list(MotrpacHumanPreSuspensionData     = "0.0.1.100",
                 MotrpacHumanPreSuspensionAnalysis = "0.2.1")

state <- list()

for (repo in c(opt$`data-pkg`, opt$`analysis-pkg`)) {
  pkg <- basename(repo)

  if (!dir.exists(repo)) {
    record(rpt, "FAIL", paste0("repo:", pkg), paste("not found at", repo))
    next
  }

  desc <- read_description(repo)
  if (is.null(desc)) {
    record(rpt, "FAIL", paste0("repo:", pkg), "no DESCRIPTION — not an R package")
    next
  }
  if (!identical(desc$Package, pkg)) {
    record(rpt, "FAIL", paste0("repo:", pkg),
           paste0("DESCRIPTION declares Package: ", desc$Package))
    next
  }

  # --- git state -----------------------------------------------------------
  if (!is_git_repo(repo)) {
    record(rpt, "FAIL", paste0("git:", pkg), "not a git checkout; the carry would be unattributable")
    next
  }
  branch <- git_branch(repo)
  sha    <- git_sha(repo)
  dirty  <- git_dirty(repo)

  record(rpt, "PASS", paste0("git:", pkg), sprintf("branch %s at %s", branch, sha))

  # Writing into a package repo from a default branch is the one thing this
  # stage must never do. It does not write at all, so this is a WARN — but it
  # stays a check, because it is the gate that matters if the carry is ever
  # pointed at the repos directly.
  if (branch %in% c("main", "master")) {
    record(rpt, "WARN", paste0("git:branch:", pkg),
           sprintf("on %s — a direct carry would have to branch first", branch))
  } else if (identical(branch, "HEAD")) {
    record(rpt, "FAIL", paste0("git:branch:", pkg), "detached HEAD")
  } else {
    record(rpt, "PASS", paste0("git:branch:", pkg), sprintf("on working branch %s", branch))
  }

  # Uncommitted edits are carried into the test package along with everything
  # else. That is usually what you want — the Analysis package's in-flight
  # roxygen edit belongs to this change — but it has to be on the record.
  if (length(dirty)) {
    record(rpt, "WARN", paste0("git:clean:", pkg),
           sprintf("%d uncommitted path(s) will be carried into the test package: %s",
                   length(dirty), paste(sub("^...", "", dirty), collapse = ", ")))
  } else {
    record(rpt, "PASS", paste0("git:clean:", pkg), "working tree clean")
  }

  # --- version -------------------------------------------------------------
  ver <- desc$Version
  base <- BASELINE[[pkg]]
  if (is.null(base)) {
    record(rpt, "WARN", paste0("version:", pkg), sprintf("%s — no baseline pin recorded", ver))
  } else if (package_version(ver) < package_version(base)) {
    record(rpt, "FAIL", paste0("version:", pkg),
           sprintf("%s is older than the v1.3 comparison baseline %s", ver, base))
  } else {
    record(rpt, "PASS", paste0("version:", pkg),
           sprintf("%s, at or ahead of the v1.3 baseline %s", ver, base))
  }

  state[[pkg]] <- data.frame(package = pkg, path = repo, branch = branch, sha = sha,
                             version = ver, dirty_paths = length(dirty),
                             stringsAsFactors = FALSE)
}

# --- the build being carried -------------------------------------------------

n_rda <- length(list.files(opt$`build-data`, pattern = "\\.rda$"))
if (!dir.exists(opt$`build-data`) || n_rda == 0L) {
  record(rpt, "FAIL", "build:stage1",
         paste0("no .rda under ", opt$`build-data`, " — run `make data` first, or point ",
                "--build-data at the checkout that holds the build"))
} else {
  record(rpt, "PASS", "build:stage1", sprintf("%d .rda in %s", n_rda, opt$`build-data`))
}

if (!is.null(opt$`preflight-data`)) {
  n_rds <- length(list.files(opt$`preflight-data`, pattern = "\\.rds$"))
  if (n_rds == 0L) {
    record(rpt, "WARN", "build:stage0",
           paste0("no .rds under ", opt$`preflight-data`,
                  " — the leaf objects will be left at their packaged values"))
  } else {
    record(rpt, "PASS", "build:stage0", sprintf("%d .rds in %s", n_rds, opt$`preflight-data`))
  }
}

# --- write -------------------------------------------------------------------

write_tsv(report_frame(rpt), file.path(out_dir, "validate_targets.tsv"))
if (length(state)) write_tsv(do.call(rbind, state), file.path(out_dir, "source_state.tsv"))

message(sprintf("\n%d check(s): %d FAIL, %d WARN",
                length(rpt$rows), rpt$fails, rpt$warns))
if (rpt$fails > 0L) quit(status = 1L)
