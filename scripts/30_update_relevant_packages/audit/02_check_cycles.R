# Stage 3, step 2 — circular dependency and circular reference check.
#
# Three independent graphs are checked, because a cycle in each one breaks
# something different:
#
#   package graph   DESCRIPTION Depends/Imports/LinkingTo between the audited
#                   packages. A cycle here means neither package installs from a
#                   clean library: each needs the other present first. This is
#                   the fatal one.
#   build graph     data-raw/ generators. A cycle means neither package's data
#                   objects can be regenerated from source without the other
#                   already built — a bootstrap deadlock rather than an install
#                   failure.
#   call graph      functions within one package calling each other in a loop.
#                   Mutual recursion is legal R and sometimes intended, so these
#                   are reported, not failed on.
#
# Suggests is deliberately NOT an edge in the package graph: R permits a cycle
# through Suggests precisely because the dependency is optional and resolved
# lazily. A cycle that runs through Suggests in one direction is reported as SOFT.
#
# Usage:
#   Rscript 02_check_cycles.R <pkg_root> [<pkg_root> ...] --out <dir>
#
# Exit: 1 if a hard package-level cycle exists, else 0.

script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
source(file.path(script_dir(), "lib", "pkg_introspect.R"))

args <- commandArgs(trailingOnly = TRUE)
out_dir <- "."
if ("--out" %in% args) {
  i <- which(args == "--out")
  out_dir <- args[i + 1L]
  args <- args[-c(i, i + 1L)]
}
pkg_roots <- normalizePath(args, mustWork = FALSE)
if (!length(pkg_roots)) stop("usage: 02_check_cycles.R <pkg_root> ... --out <dir>")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

roots <- pkg_roots[dir.exists(pkg_roots)]
if (!length(roots)) stop("none of the given package roots exist")

descs <- lapply(roots, read_description)
names(descs) <- vapply(descs, function(d) d[["Package"]], character(1))
names(roots) <- names(descs)
audited <- names(descs)

message("auditing: ", paste(audited, collapse = ", "))

# ---- 1. Package-level dependency graph -------------------------------------

edges <- list()
for (p in audited) {
  d <- descs[[p]]
  for (field in c("Depends", "Imports", "LinkingTo")) {
    for (q in intersect(description_deps(d, field), audited)) {
      if (q == p) next
      edges[[length(edges) + 1L]] <-
        data.frame(from = p, to = q, field = field, kind = "hard",
                   stringsAsFactors = FALSE)
    }
  }
  for (q in intersect(description_deps(d, "Suggests"), audited)) {
    if (q == p) next
    edges[[length(edges) + 1L]] <-
      data.frame(from = p, to = q, field = "Suggests", kind = "soft",
                 stringsAsFactors = FALSE)
  }
}
pkg_edges <- do.call(rbind, edges)
if (is.null(pkg_edges)) {
  pkg_edges <- data.frame(from = character(0), to = character(0),
                          field = character(0), kind = character(0))
}
write_tsv(pkg_edges, file.path(out_dir, "cycles_package_edges.tsv"))

hard <- pkg_edges[pkg_edges$kind == "hard", , drop = FALSE]
hard_cycles <- find_cycles(data.frame(caller = hard$from, callee = hard$to))
all_cycles <- find_cycles(data.frame(caller = pkg_edges$from, callee = pkg_edges$to))

# ---- 2. Cross-package reference capacity -----------------------------------
# For every ordered pair, where and how does one package name the other? This is
# what turns "there is an edge" into "here is the line to change".

mention_rows <- list()
for (p in audited) {
  for (q in setdiff(audited, p)) {
    m <- pkg_mentions(roots[[p]], q)
    if (!nrow(m)) next
    m$from <- p
    m$to <- q
    mention_rows[[length(mention_rows) + 1L]] <- m
  }
}
mentions <- do.call(rbind, mention_rows)
if (is.null(mentions)) {
  mentions <- data.frame(from = character(0), to = character(0), file = character(0),
                         line = integer(0), capacity = character(0), text = character(0))
} else {
  mentions <- mentions[, c("from", "to", "file", "line", "capacity", "text")]
}
write_tsv(mentions, file.path(out_dir, "cycles_cross_references.tsv"))

# A runtime edge is one that exists in shipped code: an R/ call or a NAMESPACE
# import. Those are what actually force a load order.
runtime_edge <- function(from, to) {
  sel <- mentions$from == from & mentions$to == to &
    mentions$capacity %in% c("R/ runtime call", "R/ runtime via string", "NAMESPACE")
  sum(sel)
}

declared_hard <- function(from, to) {
  any(pkg_edges$from == from & pkg_edges$to == to & pkg_edges$kind == "hard")
}
declared_any <- function(from, to) {
  any(pkg_edges$from == from & pkg_edges$to == to)
}

# ---- 3. data-raw build-order graph -----------------------------------------
# An edge from P to Q means P's generators need Q available to run.

build_edges <- list()
for (p in audited) {
  dr <- file.path(roots[[p]], "data-raw")
  if (!dir.exists(dr)) next
  files <- list.files(dr, recursive = TRUE, full.names = TRUE,
                      pattern = "\\.[Rr]$")
  if (!length(files)) next
  txt <- unlist(lapply(files, function(f) {
    tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
  }))
  for (q in setdiff(audited, p)) {
    hits <- grep(q, txt, fixed = TRUE, value = TRUE)
    hits <- hits[!grepl("^\\s*#", hits)]
    if (length(hits)) {
      build_edges[[length(build_edges) + 1L]] <-
        data.frame(from = p, to = q, n_refs = length(hits),
                   example = hits[1], stringsAsFactors = FALSE)
    }
  }
}
build_graph <- do.call(rbind, build_edges)
if (is.null(build_graph)) {
  build_graph <- data.frame(from = character(0), to = character(0),
                            n_refs = integer(0), example = character(0))
}
write_tsv(build_graph, file.path(out_dir, "cycles_build_order.tsv"))
build_cycles <- find_cycles(data.frame(caller = build_graph$from,
                                       callee = build_graph$to))

# ---- 4. Function-level call cycles, per package ----------------------------

fn_rows <- list()
for (p in audited) {
  g <- pkg_call_graph(roots[[p]])
  comps <- find_cycles(g)
  defs <- pkg_function_defs(roots[[p]])
  for (comp in comps) {
    loc <- defs[defs$name %in% comp, , drop = FALSE]
    fn_rows[[length(fn_rows) + 1L]] <- data.frame(
      package = p,
      cycle = paste(comp, collapse = " -> "),
      size = length(comp),
      files = paste(unique(basename(loc$file)), collapse = ", "),
      stringsAsFactors = FALSE)
  }
  for (nm in pkg_self_recursive(roots[[p]])) {
    loc <- defs[defs$name == nm, , drop = FALSE]
    fn_rows[[length(fn_rows) + 1L]] <- data.frame(
      package = p, cycle = paste0(nm, " -> ", nm), size = 1L,
      files = paste(unique(basename(loc$file)), collapse = ", "),
      stringsAsFactors = FALSE)
  }
}
fn_cycles <- do.call(rbind, fn_rows)
if (is.null(fn_cycles)) {
  fn_cycles <- data.frame(package = character(0), cycle = character(0),
                          size = integer(0), files = character(0))
}
write_tsv(fn_cycles, file.path(out_dir, "cycles_function_level.tsv"))

# ---- 5. Verdict -------------------------------------------------------------

verdict <- list()
say <- function(status, scope, detail) {
  verdict[[length(verdict) + 1L]] <<-
    data.frame(status = status, scope = scope, detail = detail, stringsAsFactors = FALSE)
  message(sprintf("[%s] %s — %s", status, scope, detail))
}

if (length(hard_cycles)) {
  for (comp in hard_cycles) {
    say("FAIL", "package-graph",
        paste0("hard dependency cycle (Depends/Imports): ",
               paste(comp, collapse = " <-> "),
               " — neither installs into a clean library"))
  }
} else if (length(all_cycles)) {
  for (comp in all_cycles) {
    say("WARN", "package-graph",
        paste0("soft cycle through Suggests: ", paste(comp, collapse = " <-> "),
               " — R tolerates this; installation still resolves"))
  }
} else {
  say("PASS", "package-graph", "no dependency cycle among the audited packages")
}

# A declared hard Imports edge with no runtime reference is a cycle risk carried
# for nothing: dropping it removes the edge without changing behaviour.
for (i in seq_len(nrow(hard))) {
  n <- runtime_edge(hard$from[i], hard$to[i])
  if (n == 0L) {
    say("WARN", "package-graph",
        sprintf("%s declares %s: %s but never calls it from R/ or NAMESPACE — the edge is removable",
                hard$from[i], hard$field[i], hard$to[i]))
  }
}

# The dangerous asymmetry: shipped R/ code calls into a package that DESCRIPTION
# never declares. R CMD check cannot see it, no resolver installs it, and if the
# declared graph already runs the other way this is a runtime cycle that the
# metadata reports as acyclic.
for (p in audited) {
  for (q in setdiff(audited, p)) {
    n <- runtime_edge(p, q)
    if (n > 0L && !declared_any(p, q)) {
      if (declared_hard(q, p)) {
        say("FAIL", "package-graph",
            sprintf(paste0("de-facto cycle: %s declares Imports: %s, and %s calls back into %s ",
                           "at %d undeclared runtime site(s). DESCRIPTION looks acyclic; the code is not."),
                    q, p, p, q, n))
      } else {
        say("WARN", "package-graph",
            sprintf(paste0("%s calls %s at %d runtime site(s) but declares it in no ",
                           "DESCRIPTION field — invisible to R CMD check and to every resolver"),
                    p, q, n))
      }
    }
  }
}

# A hard Imports edge that no Remotes entry can resolve. Neither package is on
# CRAN, so an unlisted GitHub-only dependency makes a fresh install impossible.
for (i in seq_len(nrow(hard))) {
  from <- hard$from[i]; to <- hard$to[i]
  rem <- if ("Remotes" %in% names(descs[[from]])) descs[[from]][["Remotes"]] else ""
  if (is.na(rem)) rem <- ""
  if (!grepl(to, rem, fixed = TRUE)) {
    say("FAIL", "package-graph",
        sprintf(paste0("%s declares %s: %s but lists no Remotes entry for it — neither package ",
                       "is on CRAN, so a clean install of %s cannot resolve it"),
                from, hard$field[i], to, from))
  }
}

if (length(build_cycles)) {
  for (comp in build_cycles) {
    say("FAIL", "build-graph",
        paste0("data-raw bootstrap deadlock: ", paste(comp, collapse = " <-> "),
               " — neither package's data objects regenerate from a clean state"))
  }
} else {
  say("PASS", "build-graph", "data-raw generators alone have an acyclic build order")
}

# The bootstrap question is not answered by either graph alone. Read both as one
# "to produce X you must already have Y" relation:
#   install edge  P Imports Q          -> installing P needs Q
#   build edge    P's data-raw uses Q  -> regenerating P's data/ needs Q
# A cycle across the union is a bootstrap deadlock even when each graph is
# separately acyclic — which is the shape this pair actually has.
bootstrap <- rbind(
  data.frame(caller = hard$from, callee = hard$to,
             why = paste0("Imports: ", hard$to), stringsAsFactors = FALSE),
  data.frame(caller = build_graph$from, callee = build_graph$to,
             why = paste0("data-raw uses ", build_graph$to), stringsAsFactors = FALSE))
write_tsv(bootstrap, file.path(out_dir, "cycles_bootstrap_edges.tsv"))
boot_cycles <- find_cycles(bootstrap[, c("caller", "callee")])

if (length(boot_cycles)) {
  for (comp in boot_cycles) {
    detail <- vapply(comp, function(p) {
      e <- bootstrap[bootstrap$caller == p & bootstrap$callee %in% comp, ]
      paste0(p, " needs ", paste(unique(e$why), collapse = " + "))
    }, character(1))
    say("FAIL", "bootstrap",
        paste0("install+build deadlock across ", paste(comp, collapse = " <-> "),
               ": ", paste(detail, collapse = "; "),
               ". Bootstrapping requires a previously built copy of one of them."))
  }
} else {
  say("PASS", "bootstrap", "install and data-raw graphs are jointly acyclic")
}

if (nrow(fn_cycles)) {
  mutual <- fn_cycles[fn_cycles$size > 1L, , drop = FALSE]
  selfrec <- fn_cycles[fn_cycles$size == 1L, , drop = FALSE]
  if (nrow(mutual)) {
    say("WARN", "call-graph",
        sprintf("%d mutually recursive function group(s); see cycles_function_level.tsv",
                nrow(mutual)))
  }
  if (nrow(selfrec)) {
    say("NOTE", "call-graph",
        sprintf("%d self-recursive function(s) — verify each has a base case", nrow(selfrec)))
  }
} else {
  say("PASS", "call-graph", "no circular function references within either package")
}

verdict_df <- do.call(rbind, verdict)
write_tsv(verdict_df, file.path(out_dir, "cycles_verdict.tsv"))

quit(status = if (any(verdict_df$status == "FAIL")) 1L else 0L)
