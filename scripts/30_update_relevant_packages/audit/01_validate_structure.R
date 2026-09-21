# Stage 3, step 1 — validate the structure of the two packages this pipeline feeds.
#
# Static checks only: nothing is installed, loaded or attached, so this runs
# against a package that does not currently install. That is deliberate — a
# structural audit is most useful exactly when R CMD check will not complete.
#
# Usage:
#   Rscript 01_validate_structure.R <pkg_root> [<pkg_root> ...] --out <dir>
#
# Writes, per package:
#   <out>/structure_<pkg>.tsv   one row per finding
# and returns a non-zero exit status if any ERROR-severity finding is present.

# Rscript does not set a script-directory variable, so recover it from the
# --file= argument the front end passes through.
script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
source(file.path(script_dir(), "lib", "pkg_introspect.R"))

# ---- Argument handling -----------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
out_dir <- "."
if ("--out" %in% args) {
  i <- which(args == "--out")
  out_dir <- args[i + 1L]
  args <- args[-c(i, i + 1L)]
}
pkg_roots <- args
if (!length(pkg_roots)) stop("usage: 01_validate_structure.R <pkg_root> ... --out <dir>")

# ---- Finding accumulator ---------------------------------------------------

findings <- new.env(parent = emptyenv())
findings$rows <- list()

add <- function(severity, check, detail, where = NA_character_) {
  findings$rows[[length(findings$rows) + 1L]] <-
    data.frame(severity = severity, check = check, detail = detail,
               where = where, stringsAsFactors = FALSE)
  invisible(NULL)
}

# ---- The checks ------------------------------------------------------------

# DESCRIPTION fields that must be real, not the roxygen skeleton's placeholder text.
check_description <- function(pkg_root, desc) {
  placeholder <- "What the package does"
  for (field in c("Title", "Description")) {
    if (!field %in% names(desc)) {
      add("ERROR", "DESCRIPTION:field", paste0("missing required field ", field), "DESCRIPTION")
      next
    }
    if (grepl(placeholder, desc[[field]], fixed = TRUE)) {
      add("WARNING", "DESCRIPTION:placeholder",
          paste0(field, " still holds the package.skeleton placeholder text"), "DESCRIPTION")
    }
  }
  # roxygen2 writes RoxygenNote; Config/roxygen2/version is roxygen 8's field. Both
  # present with different values means two roxygen versions have written this tree.
  rn <- if ("RoxygenNote" %in% names(desc)) desc[["RoxygenNote"]] else NA_character_
  cv <- if ("Config/roxygen2/version" %in% names(desc)) desc[["Config/roxygen2/version"]] else NA_character_
  if (!is.na(rn) && !is.na(cv) && !identical(rn, cv)) {
    add("NOTE", "DESCRIPTION:roxygen-version",
        sprintf("RoxygenNote (%s) and Config/roxygen2/version (%s) disagree", rn, cv),
        "DESCRIPTION")
  }
}

# Every data/*.rda needs a \docType{data} man page, and every data man page needs
# a matching object. Both directions are R CMD check findings.
check_data_docs <- function(pkg_root) {
  objs <- pkg_data_objects(pkg_root)
  doc <- pkg_data_doc_aliases(pkg_root)
  documented <- doc$aliases

  undocumented <- setdiff(objs, documented)
  for (o in undocumented) {
    add("WARNING", "data:undocumented",
        paste0("data/", o, ".rda has no \\docType{data} man page"),
        paste0("data/", o, ".rda"))
  }
  # An @rdname topic name is not itself an object, so it is not an orphan.
  orphaned <- setdiff(documented, c(objs, doc$topics))
  for (o in orphaned) {
    add("WARNING", "data:orphan-doc",
        paste0("documented as data but no object in data/: ", o),
        "man/")
  }
  add("NOTE", "data:inventory",
      sprintf("%d objects in data/, %d documented, %d undocumented, %d orphan docs (%d @rdname topics excluded)",
              length(objs), length(intersect(objs, documented)),
              length(undocumented), length(orphaned),
              length(setdiff(doc$topics, objs))),
      "data/")
}

# An Imports entry nothing references is dead weight that still gates installation.
check_declared_deps <- function(pkg_root, desc, ns) {
  imports <- description_deps(desc, "Imports")
  depends <- description_deps(desc, "Depends")
  suggests <- description_deps(desc, "Suggests")

  ns_pkgs <- unique(c(ns$imports, ns$import_from$pkg))
  qualified <- pkg_qualified_calls(pkg_root)
  attached <- pkg_attached_calls(pkg_root)
  called <- sort(unique(qualified$pkg))

  used_in_r <- function(p) p %in% called || p %in% attached || p %in% ns_pkgs

  for (p in imports) {
    if (!used_in_r(p)) {
      add("NOTE", "deps:declared-unused",
          paste0("Imports: ", p, " — no `", p, "::` call, no library()/require() and no NAMESPACE import in R/"),
          "DESCRIPTION")
    }
  }
  for (p in depends) {
    if (!used_in_r(p)) {
      add("NOTE", "deps:declared-unused",
          paste0("Depends: ", p, " — attached but never referenced in R/"), "DESCRIPTION")
    }
  }

  # The reverse: a `pkg::` call whose package is declared nowhere is an R CMD
  # check ERROR, because the namespace will not be available at run time.
  # A package qualifying its own namespace (`Pkg::obj` inside Pkg) is legal and is
  # the idiomatic way to reach a LazyData object, so it is not a missing dependency.
  #
  # Only `base` is exempt from declaration. stats, utils, grDevices, methods and
  # the rest ship with R but are still ordinary namespaces: importing from one
  # without an Imports entry is the R CMD check ERROR "Namespace dependency not
  # required". Whitelisting them here is what hid the missing grDevices.
  declared <- c(imports, depends, suggests, description_deps(desc, "LinkingTo"),
                desc[["Package"]], "base")
  for (p in setdiff(called, declared)) {
    where <- qualified[qualified$pkg == p, ][1, ]
    add("ERROR", "deps:undeclared",
        sprintf("R/ calls %s::, but %s is in no DESCRIPTION dependency field (first at %s:%d)",
                p, p, basename(where$file), where$line),
        "DESCRIPTION")
  }

  # NAMESPACE imports a package DESCRIPTION never declares. R CMD check reports
  # this as "Namespace dependency not required", and it is easy to introduce: a
  # roxygen @importFrom regenerates NAMESPACE without touching DESCRIPTION, so the
  # two drift apart with nothing in R/ to show for it.
  for (p in setdiff(ns_pkgs, declared)) {
    add("ERROR", "deps:namespace-not-declared",
        paste0("NAMESPACE imports from ", p, " but DESCRIPTION declares it nowhere"),
        "NAMESPACE")
  }

  # A Suggests used in R/ without a requireNamespace guard fails for any user who
  # did not install the optional package.
  # Only requireNamespace()/loadNamespace() actually guard — library() inside a
  # package body does not make an optional dependency optional.
  guards <- pkg_attached_calls(pkg_root, c("requireNamespace", "loadNamespace"))
  for (p in intersect(suggests, called)) {
    if (p %in% guards) next
    where <- qualified[qualified$pkg == p, ][1, ]
    add("WARNING", "deps:unguarded-suggests",
        sprintf("Suggests: %s is called at %s:%d with no requireNamespace() guard",
                p, basename(where$file), where$line),
        "R/")
  }
}

# Exports must resolve to something the package actually defines.
check_namespace <- function(pkg_root, ns, defs) {
  if (!ns$generated) {
    add("NOTE", "namespace:hand-written",
        "NAMESPACE has no roxygen2 header — it is hand-maintained and can drift from R/",
        "NAMESPACE")
  }
  data_objs <- pkg_data_objects(pkg_root)
  # A re-export (`@importFrom magrittr %>%` + `@export`) is defined in the
  # upstream package, not in R/, so importFrom names count as known.
  reexports <- ns$import_from$fn
  known <- unique(c(defs$name, data_objs, reexports))
  for (e in setdiff(ns$exports, known)) {
    add("ERROR", "namespace:export-undefined",
        paste0("export(", e, ") but no definition in R/, no object in data/, and no importFrom to re-export"),
        "NAMESPACE")
  }
}

# man/*.Rd that roxygen would not regenerate: no roxygen block anywhere claims
# the topic, so `devtools::document()` leaves it behind forever.
check_stale_rd <- function(pkg_root) {
  man_dir <- file.path(pkg_root, "man")
  if (!dir.exists(man_dir)) return(invisible(NULL))
  rds <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)
  for (rd in rds) {
    txt <- readLines(rd, warn = FALSE)
    if (!any(grepl("Generated by roxygen2", txt, fixed = TRUE))) {
      add("NOTE", "man:not-roxygen-generated",
          paste0(basename(rd), " was not written by roxygen2 — document() will not update it"),
          file.path("man", basename(rd)))
    }
  }
}

# Directories that ship in the tarball unless .Rbuildignore excludes them.
check_build_hygiene <- function(pkg_root) {
  entries <- setdiff(list.files(pkg_root, all.files = TRUE, no.. = TRUE), c(".git"))
  standard <- c("R", "man", "data", "DESCRIPTION", "NAMESPACE", "tests",
                "vignettes", "inst", "src", "LICENSE", "LICENSE.md", ".Rbuildignore",
                ".gitignore", ".gitattributes", "NEWS.md", "README.md")
  for (e in setdiff(entries, standard)) {
    if (!rbuildignore_covers(pkg_root, e)) {
      sev <- if (e %in% c("sandbox", ".Rhistory", ".Rproj.user", ".DS_Store",
                          "data-raw", ".claude", ".quarto")) "WARNING" else "NOTE"
      add(sev, "build:not-ignored",
          paste0(e, " is not matched by .Rbuildignore — it ships in the tarball"),
          e)
    }
  }
}

# ---- Drive ------------------------------------------------------------------

exit_status <- 0L
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (pkg_root in pkg_roots) {
  pkg_root <- normalizePath(pkg_root, mustWork = FALSE)
  if (!dir.exists(pkg_root)) {
    message("[FAIL] not a directory: ", pkg_root)
    exit_status <- 1L
    next
  }
  findings$rows <- list()
  desc <- read_description(pkg_root)
  pkg_name <- desc[["Package"]]
  ns <- parse_namespace(pkg_root)
  defs <- pkg_function_defs(pkg_root)

  check_description(pkg_root, desc)
  check_data_docs(pkg_root)
  check_declared_deps(pkg_root, desc, ns)
  check_namespace(pkg_root, ns, defs)
  check_stale_rd(pkg_root)
  check_build_hygiene(pkg_root)

  res <- do.call(rbind, findings$rows)
  if (is.null(res)) {
    res <- data.frame(severity = character(0), check = character(0),
                      detail = character(0), where = character(0))
  }
  res <- res[order(factor(res$severity, levels = c("ERROR", "WARNING", "NOTE", "STYLE")),
                   res$check), ]
  path <- file.path(out_dir, paste0("structure_", pkg_name, ".tsv"))
  write_tsv(res, path)

  n_err <- sum(res$severity == "ERROR")
  n_warn <- sum(res$severity == "WARNING")
  message(sprintf("[%s] %s: %d ERROR, %d WARNING, %d NOTE  -> %s",
                  if (n_err) "FAIL" else " OK ", pkg_name, n_err, n_warn,
                  sum(res$severity == "NOTE"), path))
  if (n_err > 0L) exit_status <- 1L
}

quit(status = exit_status)
