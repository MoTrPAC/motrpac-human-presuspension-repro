#!/usr/bin/env Rscript
# =============================================================================
# build_depgraph.R — the motrpac-human-presuspension-repro dependency graph.
#
# A "targets-style" schematic of THIS repo: which Make target runs which stage
# driver, which driver runs which step, which step runs which R stem, what each
# stem reads and writes, and what the widest-blast-radius nodes are.
#
# It is the pipeline-shaped sibling of the package-shaped graph in
# notes/PreCovid Comp/dependency_graph/, which maps the three MoTrPAC R packages
# (data-raw/*.R -> .rda objects -> NAMESPACE exports -> figure scripts). That
# topology does not exist here: this repo is a Make/shell pipeline whose currency
# is BIC-named freeze files, so the node and edge taxonomies are rebuilt from the
# ground up while the rendering layer (visNetwork + the two CSVs + app.R) is the
# same.
#
# Strategy = hybrid, same as the package graph:
#
#   * AUTO-DERIVED (reproducible) --------------------------------------------
#       - make targets   : Makefile target/prereq/recipe parse, incl. the
#                          $(STAMPS)/x -> x stamp indirection
#       - stage drivers  : whatever the Makefile recipes actually invoke, so a
#                          renamed stage (30_figures -> 30_update_relevant_packages)
#                          is picked up without editing this script
#       - steps          : the STEP_DIRS array in 10_build_data.sh
#       - stems          : the ALL_STEMS arrays in each build.sh, plus the
#                          "${HERE}/<name>.R" single-stem form
#       - preflight      : the BUILDERS / PREBUILT_COPY arrays in
#                          00_preflight/data-raw/build_data_objects.sh
#       - freeze outputs : every key in config/file_versions.json, decomposed into
#                          tissue_code / ome / data_category / data_details and
#                          grouped into one node per (ome, category, details)
#       - gates          : the label::path arrays in check_required_inputs.sh
#       - external inputs: 00_preflight/config/external_assets.tsv + gmt_files.tsv,
#                          whose consumed_by / builder_script columns already name
#                          the consuming script
#       - defines/calls  : the intra-repo function call graph, recovered by parsing
#                          every .R under scripts/ with R's own parser
#                          (getParseData — AST tokens, NO code execution)
#       - sourced_by     : source(...) targets, matched on basename
#       - ledger         : docs/data_objects.tsv (object / node / step / generator /
#                          repo / deps) is hand-maintained upstream truth for the
#                          object-level DAG; its deps column is read as edges rather
#                          than re-derived
#
#   * CURATED (what is not cleanly derivable) --------------------------------
#       - OME_FAMILY: which stem family owns an ome whose stem never names it as a
#         string literal. Only the metabolomics stems need this — they loop over
#         omes pulled out of OME_TISSUE_CODE at run time, so no literal exists to
#         grep. Every other ome is auto-attributed from a literal in its stem.
#       - FREEZE_CONSUMERS: which later steps read which freeze category. The reads
#         go through load_qc_local() / list.files() over staging/freeze, so the
#         filename never appears in the consuming source either.
#
# It parses SOURCE TEXT AND CONFIG ONLY. It does not run the pipeline, does not
# read staging/freeze, and does not need the gitignored build outputs
# (scripts/00_preflight/data/*.rds, staging/**) to be present — so it works on a
# fresh clone.
#
# Output (written next to this script):
#   depgraph_nodes.csv          node table
#   depgraph_edges.csv          edge table
#   depgraph_freeze_files.csv   the per-file freeze inventory behind the grouped
#                               freeze nodes (one row per file_versions.json key)
#   depgraph.html               Layer 1 — the full graph
#   depgraph_summary.html       Layer 2 — collapsed to stage super-nodes
#
# Re-run any time the pipeline changes:
#   make depgraph                       # or:
#   Rscript docs/dependency_graph/build_depgraph.R [--root /path/to/motrpac-human-presuspension-repro]
# =============================================================================

suppressWarnings(suppressMessages({
  .have <- vapply(c("visNetwork", "htmlwidgets", "igraph", "dplyr", "jsonlite"),
                  requireNamespace, logical(1), quietly = TRUE)
}))
if (!all(.have)) {
  stop("Missing R packages: ", paste(names(.have)[!.have], collapse = ", "),
       "\nInstall with: install.packages(c(",
       paste(sprintf('"%s"', names(.have)[!.have]), collapse = ", "), "))")
}
suppressWarnings(suppressMessages({
  library(dplyr)
  library(visNetwork)
}))

# ---- locate the repo --------------------------------------------------------
# This script lives at <repo>/docs/dependency_graph/, so the root is two levels up.
# --root overrides that, which is how you point it at another checkout.
.args <- commandArgs(trailingOnly = FALSE)
.this <- sub("^--file=", "", grep("^--file=", .args, value = TRUE)[1])
SCRIPT_DIR <- if (!is.na(.this)) normalizePath(dirname(.this)) else getwd()

.cli <- commandArgs(trailingOnly = TRUE)
.root_arg <- if (length(.cli) && .cli[1] == "--root") .cli[2] else NA_character_
ROOT <- if (!is.na(.root_arg)) {
  normalizePath(.root_arg)
} else {
  normalizePath(file.path(SCRIPT_DIR, "..", ".."))
}
if (!file.exists(file.path(ROOT, "config", "pipeline.env")))
  stop("not a motrpac-human-presuspension-repro checkout (no config/pipeline.env): ", ROOT)
message("repo root  : ", ROOT)
message("output dir : ", SCRIPT_DIR)

# =============================================================================
# helpers
# =============================================================================
rd <- function(...) tryCatch(readLines(file.path(ROOT, ...), warn = FALSE),
                             error = function(e) character(0))
has <- function(...) file.exists(file.path(ROOT, ...))
rel <- function(p) sub(paste0("^", ROOT, "/?"), "", p)

# Capture group 1 of the first match on each line, dropping non-matches.
grab <- function(x, pat, group = 2L) {
  m <- regmatches(x, regexec(pat, x, perl = TRUE))
  out <- vapply(m, function(z) if (length(z) >= group) z[group] else NA_character_, "")
  out[!is.na(out)]
}

# Every element of a shell array literal:  NAME=( a b \n c )  ->  c("a","b","c").
# Handles the multi-line form the drivers use and strips comments and quotes.
sh_array <- function(lines, name) {
  start <- grep(paste0("^\\s*(declare -a )?", name, "=\\("), lines)
  if (!length(start)) return(character(0))
  start <- start[1]
  end <- start
  while (end <= length(lines) && !grepl("\\)", sub(paste0(".*", name, "=\\("), "", lines[end])) &&
         !grepl("^\\s*\\)", lines[end])) end <- end + 1
  # widen until the closing paren is seen (the first line may itself close it)
  while (end <= length(lines) && !grepl("\\)", lines[end])) end <- end + 1
  body <- lines[start:min(end, length(lines))]
  body <- sub(paste0(".*", name, "=\\("), "", body)
  body <- sub("\\).*$", "", body)
  body <- sub("#.*$", "", body)
  toks <- unlist(strsplit(paste(body, collapse = " "), "[[:space:]]+"))
  toks <- gsub('^"|"$|^\'|\'$', "", toks)
  toks[nzchar(toks)]
}

# Expand the handful of Make variables the recipes use.
mk_expand <- function(x) {
  x <- gsub("\\$\\(SCRIPTS\\)", "scripts", x)
  x <- gsub("\\$\\(STAMPS\\)", ".stamps", x)
  x <- gsub("\\$\\(MAKEFILE_LIST\\)", "Makefile", x)
  x
}

# =============================================================================
# node / edge accumulators
#
# Nodes carry: id, label, type, stage, ome, origin, detail, n_downstream,
# hotspot, level. Edges carry: from, to, type. Both are written to CSV and both
# are what app.R reads, so the column names are part of the contract.
# =============================================================================
NODES <- list(); EDGES <- list()

add_node <- function(id, label, type, stage = "stage1", ome = "shared",
                     origin = "auto", detail = "") {
  if (!is.null(NODES[[id]])) return(invisible(id))
  NODES[[id]] <<- data.frame(id = id, label = label, type = type, stage = stage,
                             ome = ome, origin = origin, detail = detail,
                             stringsAsFactors = FALSE)
  invisible(id)
}
add_edge <- function(from, to, type) {
  EDGES[[length(EDGES) + 1L]] <<- data.frame(from = from, to = to, type = type,
                                             stringsAsFactors = FALSE)
  invisible(NULL)
}

# =============================================================================
# 1. Makefile — targets, ordering, and the scripts each target runs
# =============================================================================
mk <- rd("Makefile")

# Target lines: `name: prereqs ## help`. Skip pattern rules and variable
# assignments; `.PHONY`/`.stamps` bookkeeping lines are filtered below.
tgt_idx <- grep("^[a-zA-Z_][a-zA-Z0-9_.-]*:", mk)
mk_targets <- list()
for (i in tgt_idx) {
  ln    <- mk[i]
  name  <- sub(":.*$", "", ln)
  rhs   <- sub("^[^:]*:", "", ln)
  help  <- if (grepl("##", rhs)) trimws(sub("^.*##\\s*", "", rhs)) else ""
  prereq <- trimws(sub("##.*$", "", rhs))
  prereq <- gsub("\\|", " ", prereq)                       # order-only marker
  prereq <- unlist(strsplit(mk_expand(prereq), "[[:space:]]+"))
  prereq <- prereq[nzchar(prereq) & prereq != ".stamps"]
  prereq <- sub("^\\.stamps/", "", prereq)                 # stamp -> the target
  # recipe: the indented lines until the next non-indented, non-blank line
  j <- i + 1L; recipe <- character(0)
  while (j <= length(mk) && (grepl("^\t", mk[j]) || !nzchar(trimws(mk[j])))) {
    if (grepl("^\t", mk[j])) recipe <- c(recipe, mk[j]); j <- j + 1L
  }
  mk_targets[[name]] <- list(name = name, help = help, prereq = unique(prereq),
                             recipe = mk_expand(recipe))
}
# `$(STAMPS)/preflight: preflight` lines only restate the stamp indirection.
mk_targets <- mk_targets[!grepl("^\\.stamps", names(mk_targets))]
mk_targets[["help"]] <- NULL

# Scripts a recipe invokes: `bash <path>` / `Rscript <path>` / `$(RSCRIPT) <path>`.
recipe_scripts <- function(recipe) {
  pat <- "(?:bash|sh|Rscript|\\$\\(RSCRIPT\\))\\s+([A-Za-z0-9_./-]+\\.(?:sh|R))"
  unique(grab(recipe, paste0(".*", pat), 2L))
}

for (t in mk_targets) {
  id <- paste0("make:", t$name)
  add_node(id, paste("make", t$name), "make_target", "orchestration",
           detail = t$help)
}
for (t in mk_targets) {
  for (p in t$prereq)
    if (!is.null(mk_targets[[p]])) add_edge(paste0("make:", p), paste0("make:", t$name), "orders")
}
message("make       : ", length(mk_targets), " targets")

# =============================================================================
# 2. Stage drivers, steps, stems — everything the Makefile actually runs
# =============================================================================
# A driver is any script a Make recipe invokes. Deriving the set this way rather
# than hardcoding it means a renamed or added stage shows up without touching
# this script.
STAGE_OF <- function(path) {
  if (grepl("00_preflight", path)) "stage0"
  else if (grepl("10_build_data", path)) "stage1"
  else if (grepl("20_upload", path)) "stage2"
  else if (grepl("scripts/30_", path)) "stage3"
  else "orchestration"
}

# One indirection to see through: the stage targets no longer name their stage
# script, they all run `scripts/run_stage.sh <stage>`. That script is the single
# definition of what a stage runs, because the SLURM chain submitter runs the same
# script inside a batch job — but taking the recipes literally would collapse four
# stages onto one dispatcher node. So parse its case arms into stage -> scripts and
# expand. A target with no arm keeps whatever its recipe named, so this stays quiet
# if run_stage.sh is ever dropped.
RUN_STAGE <- "scripts/run_stage.sh"
run_stage_map <- function() {
  lines <- rd(RUN_STAGE)
  out <- list(); arm <- NULL
  for (ln in lines) {
    hit <- grab(ln, "^\\s*([a-z0-9-]+)\\)\\s*$", 2L)
    if (length(hit)) { arm <- hit[1]; next }
    if (grepl("^\\s*;;", ln)) { arm <- NULL; next }
    if (is.null(arm)) next
    sc <- grab(ln, "\\$\\{PIPELINE_ROOT\\}/(scripts/[A-Za-z0-9_./-]+\\.sh)", 2L)
    if (length(sc)) out[[arm]] <- union(out[[arm]], sc)
  }
  out
}
RUN_STAGE_MAP <- run_stage_map()

drivers <- character(0)
for (t in mk_targets) {
  for (s in recipe_scripts(t$recipe)) {
    expanded <- if (identical(s, RUN_STAGE) && !is.null(RUN_STAGE_MAP[[t$name]]))
      RUN_STAGE_MAP[[t$name]] else s
    for (sc in expanded) {
      if (!has(sc)) next
      id <- paste0("sh:", sc)
      add_node(id, basename(sc), "stage_driver", STAGE_OF(sc),
               detail = paste("run by make", t$name))
      add_edge(paste0("make:", t$name), id, "runs")
      drivers <- union(drivers, sc)
    }
  }
}
message("drivers    : ", length(drivers), " (", paste(basename(drivers), collapse = ", "), ")")

# ---- Stage 1 steps ----------------------------------------------------------
BUILD_DATA <- "scripts/10_build_data/10_build_data.sh"
step_dirs <- if (has(BUILD_DATA)) sh_array(rd(BUILD_DATA), "STEP_DIRS") else character(0)
STEP_INDEX <- stats::setNames(seq_along(step_dirs), step_dirs)
message("steps      : ", length(step_dirs))

# The ome universe, taken from the freeze inventory (section 4 groups these; the
# stem scan below needs the vocabulary first).
#
# The one place this script runs repo code rather than parsing it: lib/file_versions.R
# is sourced for freeze_files(), which walks the
# [[tissue_code]][[ome]][[data_category]][[data_details]] tree in
# config/file_versions.json. It reads that one config file and defines functions —
# no pipeline, no build outputs — so the "config only" contract above still holds,
# and the alternative (re-splitting composed filenames on "_" here) would be a second
# definition of the naming rule to keep in step with the map.
Sys.setenv(PRECOVID_ROOT = ROOT)
source(file.path(ROOT, "scripts", "10_build_data", "lib", "file_versions.R"))
FREEZE <- freeze_files()
names(FREEZE)[match(c("data_category", "data_details"), names(FREEZE))] <- c("category", "details")
OMES <- sort(unique(stats::na.omit(FREEZE$ome)))

# Which omes a stem names as string literals. Cheap, exact, and it covers every
# stem except the metabolomics ones (they enumerate omes out of OME_TISSUE_CODE
# at run time, so no literal exists — see OME_FAMILY below).
stem_omes <- function(path) {
  txt <- tryCatch(readLines(file.path(ROOT, path), warn = FALSE), error = function(e) character(0))
  pd <- tryCatch(utils::getParseData(parse(text = txt, keep.source = TRUE)),
                 error = function(e) NULL)
  if (is.null(pd)) return(character(0))
  lit <- unique(gsub('^["\']|["\']$', "", pd$text[pd$token == "STR_CONST"]))
  intersect(OMES, lit)
}

STEPS <- list()
for (sd in step_dirs) {
  bs <- file.path("scripts/10_build_data", sd, "build.sh")
  step_id <- paste0("step:", sd)
  bl <- rd(bs)
  # exit 77 means SKIP, which covers two different things. A step that never runs
  # anything is a placeholder; one that exits 77 behind a conditional and calls
  # Rscript below it is adapted and runs unless opted out of (step 16). Both report
  # SKIP to 10_build_data.sh; only the first is unbuilt.
  skips <- any(grepl("^\\s*exit 77", bl))
  runs_r <- any(grepl("\\$\\{RSCRIPT\\}", bl))
  stub <- skips && !runs_r
  gated <- skips && runs_r
  add_node(step_id, sd, "step", "stage1",
           origin = if (stub) "stub" else "auto",
           detail = if (stub) "not yet adapted (exit 77 SKIP)"
                    else if (gated) "adapted; on by default (exit 77 SKIP when opted out)"
                    else "")
  add_edge(paste0("sh:", BUILD_DATA), step_id, "runs")

  # Stems: the ALL_STEMS array (steps 06/09) or the "${HERE}/<name>.R" form.
  stems <- sh_array(bl, "ALL_STEMS")
  if (!length(stems)) {
    inline <- grab(bl, '.*\\$\\{(?:HERE|GEN_DIR)\\}/([A-Za-z0-9_.-]+)\\.R', 2L)
    stems <- unique(inline)
  }
  tests <- stems[grepl("_tests$", stems)]
  stems <- stems[!grepl("_tests$", stems)]

  for (st in stems) {
    p <- file.path("scripts/10_build_data", sd, paste0(st, ".R"))
    if (!has(p)) next
    so <- stem_omes(p)
    sid <- paste0("stem:", sd, "/", st)
    add_node(sid, st, "stem", "stage1",
             ome = if (length(so) == 1) so else "shared",
             detail = if (length(so)) paste(so, collapse = ", ") else "")
    add_edge(step_id, sid, "runs")
    STEPS[[sid]] <- list(step = sd, stem = st, path = p, omes = so)
  }
  for (tt in tests) {
    p <- file.path("scripts/10_build_data", sd, paste0(tt, ".R"))
    if (!has(p)) next
    tid <- paste0("test:", sd, "/", tt)
    add_node(tid, tt, "test_script", "stage1")
    add_edge(step_id, tid, "tested_by")
  }
  # Steps 01/02/05 stage prebuilt objects in rather than building them; the copy
  # source is an external input, so record it as one.
  if (any(grepl("DATA_PKG_REPO", bl))) {
    add_node("src:DATA_PKG_REPO", "MotrpacHumanPreSuspensionData", "source", "external",
             detail = "sibling package checkout — objects staged verbatim")
    add_edge("src:DATA_PKG_REPO", step_id, "consumed_by")
  }
  if (any(grepl("sources/", bl))) {
    vid <- paste0("src:vendored/", sd)
    add_node(vid, paste0(sd, "/sources"), "source", "external",
             detail = "vendored raw source under the step's sources/")
    add_edge(vid, step_id, "consumed_by")
  }
}

# Test scripts that live beside a step but are not named in an ALL_STEMS array
# (single-stem steps invoke them directly by filename).
for (sd in step_dirs) {
  d <- file.path(ROOT, "scripts/10_build_data", sd)
  if (!dir.exists(d)) next
  for (f in list.files(d, pattern = "_tests\\.R$")) {
    tid <- paste0("test:", sd, "/", sub("\\.R$", "", f))
    if (is.null(NODES[[tid]])) {
      add_node(tid, sub("\\.R$", "", f), "test_script", "stage1")
      add_edge(paste0("step:", sd), tid, "tested_by")
    }
  }
}

# ---- Stage 0 preflight builders --------------------------------------------
BDO <- "scripts/00_preflight/data-raw/build_data_objects.sh"
pf_builders <- if (has(BDO)) sh_array(rd(BDO), "BUILDERS") else character(0)
pf_prebuilt <- if (has(BDO)) sh_array(rd(BDO), "PREBUILT_COPY") else character(0)
for (b in pf_builders) {
  p <- file.path("scripts/00_preflight/data-raw", b)
  bid <- paste0("pf:", b)
  add_node(bid, basename(b), "preflight_builder", "stage0",
           origin = if (has(p)) "auto" else "missing")
  add_edge(paste0("sh:", BDO), bid, "runs")
}
for (b in pf_prebuilt) {
  bid <- paste0("pfout:", b)
  add_node(bid, b, "preflight_object", "stage0", origin = "vendored",
           detail = "download-only GMT, vendored prebuilt and copied verbatim")
  add_edge(paste0("sh:", BDO), bid, "produces")
}
message("preflight  : ", length(pf_builders), " builders, ", length(pf_prebuilt), " prebuilt copies")

# =============================================================================
# 3. docs/data_objects.tsv — the curated object ledger
#
# The ledger is the hand-maintained record of every data object the collection
# ships: which node group it belongs to, which step and generator build it, and
# which node groups it depends on. Its deps column is the object-level DAG, so it
# is read as edges rather than re-derived from source.
# =============================================================================
LEDGER <- if (has("docs/data_objects.tsv")) {
  read.csv(file.path(ROOT, "docs", "data_objects.tsv"), sep = "\t",
           check.names = FALSE, stringsAsFactors = FALSE)
} else {
  data.frame()
}

obj_nodes <- character(0)
if (nrow(LEDGER)) {
  grp <- LEDGER %>%
    dplyr::filter(node != "-") %>%
    dplyr::group_by(node) %>%
    dplyr::summarise(step = dplyr::first(step), generator = dplyr::first(generator),
                     repo = dplyr::first(repo), deps = dplyr::first(deps),
                     regenerated = dplyr::first(regenerated), n_obj = dplyr::n(),
                     .groups = "drop")
  for (i in seq_len(nrow(grp))) {
    g <- grp[i, ]
    id <- paste0("obj:", g$node)
    add_node(id, g$node, "data_object",
             if (g$step == "00_preflight") "stage0" else "stage1",
             origin = g$regenerated,
             detail = sprintf("%d object(s) | generator %s | %s pkg | regenerated: %s",
                              g$n_obj, g$generator, g$repo, g$regenerated))
    obj_nodes <- c(obj_nodes, g$node)
    # tie the object group to the step that builds it
    sid <- paste0("step:", g$step)
    if (!is.null(NODES[[sid]])) add_edge(sid, id, "produces")
    else if (g$step == "00_preflight") {
      pb <- paste0("pf:", g$generator)
      if (!is.null(NODES[[pb]])) add_edge(pb, id, "produces")
      else add_edge(paste0("sh:", BDO), id, "produces")
    }
  }
  # deps: semicolon-separated node names (or an upstream generator filename)
  for (i in seq_len(nrow(grp))) {
    g <- grp[i, ]
    if (is.na(g$deps) || g$deps == "-" || !nzchar(g$deps)) next
    for (d in trimws(unlist(strsplit(g$deps, ";")))) {
      if (d %in% obj_nodes) {
        add_edge(paste0("obj:", d), paste0("obj:", g$node), "feeds")
      } else {
        # a generator that lives upstream of this repo (e.g. clinic_download.R)
        did <- paste0("src:upstream/", d)
        add_node(did, d, "source", "external",
                 detail = "upstream generator, not adapted into this repo")
        add_edge(did, paste0("obj:", g$node), "feeds")
      }
    }
  }
  message("ledger     : ", nrow(LEDGER), " objects in ", nrow(grp), " node groups")
}

# =============================================================================
# 4. Freeze outputs — config/file_versions.json
#
# Every row of FREEZE is one file this pipeline writes, named
#   human-precovid-sed-adu_<tissue_code>_<ome>_<data_category>_<data_details>.txt
# exactly what write_with_path_name() composes before appending _v<version><ext>.
# Those four fields are the four levels of the map's tree, so they arrive already
# split; the rows with no tissue_code are the resources/ tier, whose filenames
# follow no such convention. Files are grouped into one node per (ome, category,
# details); the per-file inventory is written to depgraph_freeze_files.csv.
#
# CURATED: the metabolomics stems build their ome list from OME_TISSUE_CODE at
# run time, so no ome literal appears in their source for the literal scan to
# find. Everything else is auto-attributed from a literal.
# =============================================================================
OME_FAMILY <- c(
  "^metab-[tu]-(?!clinical)" = "metab",     # generate_metab_*  (loops over OME_TISSUE_CODE)
  "clinical$"                = "clinical"   # prot-clinical + metab-t-clinical
)

# Which stems can produce which category. Derived from the step a stem lives in:
# step 06 writes qc-norm + metadata, its *_imputed stems write imputed, step 09
# writes da. A stem's own name settles the rest.
CAT_STEP <- c("qc-norm" = "06_generate_qc_norm", "metadata" = "06_generate_qc_norm",
              "imputed" = "06_generate_qc_norm", "da" = "09_build_da")

stem_for <- function(ome, category, details) {
  if (identical(details, "removed-samples")) return("stem:06_generate_qc_norm/stage_removed_samples")
  step <- CAT_STEP[[category]]
  if (is.null(step) || is.na(step)) return(NA_character_)
  cands <- Filter(function(s) STEPS[[s]]$step == step, names(STEPS))
  if (identical(category, "imputed")) cands <- Filter(function(s) grepl("_imputed$", s), cands)
  else cands <- Filter(function(s) !grepl("_imputed$", s), cands)
  # auto: the stem names this ome as a literal
  hit <- Filter(function(s) ome %in% STEPS[[s]]$omes, cands)
  if (length(hit) == 1) return(hit)
  if (length(hit) > 1) return(hit[1])
  # curated fallback: ome family -> stem family
  for (pat in names(OME_FAMILY)) {
    if (grepl(pat, ome, perl = TRUE)) {
      fam <- OME_FAMILY[[pat]]
      h2 <- Filter(function(s) grepl(paste0("generate_", fam, "_"), s), cands)
      if (length(h2)) return(h2[1])
    }
  }
  NA_character_
}

FREEZE$stem <- NA_character_
FREEZE$node <- NA_character_
for (i in seq_len(nrow(FREEZE))) {
  f <- FREEZE[i, ]
  if (is.na(f$category)) next        # the resources/ tier
  # The seven .html entries are upstream QC reports carried into the release, not
  # built here. Two of them sit at category=metadata/details=qc-report and the rest
  # at category=qc-report, so the carve-out keys on the extension.
  if (identical(f$ext, "html")) next
  FREEZE$stem[i] <- stem_for(f$ome, f$category, f$details)
  FREEZE$node[i] <- paste0("freeze:", f$ome, "|", f$category, "|", f$details)
}

fam <- FREEZE %>%
  dplyr::filter(!is.na(node)) %>%
  dplyr::group_by(node, ome, category, details) %>%
  dplyr::summarise(n_files = dplyr::n(),
                   versions = paste(sort(unique(version)), collapse = ", "),
                   tissues  = paste(sort(unique(tissue_code)), collapse = ", "),
                   stem     = dplyr::first(stats::na.omit(stem)),
                   .groups = "drop")

for (i in seq_len(nrow(fam))) {
  f <- fam[i, ]
  add_node(f$node, paste0(f$category, " ", f$details), "freeze_file", "stage1",
           ome = f$ome,
           origin = if (is.na(f$stem)) "unattributed" else "auto",
           detail = sprintf("%s | %d file(s) | v%s | %s",
                            f$ome, f$n_files, f$versions, f$tissues))
  if (!is.na(f$stem)) add_edge(f$stem, f$node, "produces")
}
n_unattr <- sum(is.na(fam$stem))
message("freeze     : ", nrow(FREEZE), " files -> ", nrow(fam), " grouped nodes",
        if (n_unattr) paste0(" (", n_unattr, " with no producing stem)") else "")
if (n_unattr) message("  unattributed: ", paste(fam$node[is.na(fam$stem)], collapse = ", "))

# The resources/ tier, whose names do not follow the BIC stem convention, plus the
# upstream qc-report HTMLs that are carried through rather than built here.
odd <- FREEZE %>% dplyr::filter(is.na(node))
for (i in seq_len(nrow(odd))) {
  o <- odd[i, ]
  carried <- identical(o$ext, "html")
  id <- paste0("freeze:", o$key)
  add_node(id, sub("[.][a-z]+$", "", o$key), "freeze_file", "stage1",
           ome = if (is.na(o$ome)) "shared" else o$ome,
           origin = if (carried) "carried" else "auto",
           detail = sprintf("v%s%s", o$version,
                            if (carried) " | upstream QC report, carried not built" else ""))
}
# The resource files this repo does build, wired to the stems that write them. Every
# resources/ entry is covered here except motrpac_human-precovid_1kg_pca.csv, which the
# collection carries but this pipeline neither builds nor stages.
RESOURCE_STEMS <- c(
  "motrpac_human-precovid_metabolite-cv.txt"             = "05_stage_metabolomics_cvs/METABOLOMICS_CVS_resource",
  "motrpac-mappings-human-feature-to-gene.txt"           = "07_build_human_feature_to_gene/HUMAN_FEATURE_TO_GENE",
  "txdb_hsapiens_ensembl_v105.sqlite"                    = "06_generate_qc_norm/stage_ensembl_txdb",
  "motrpac_human-precovid_refmet-map.txt"                = "06_generate_qc_norm/stage_refmet_map",
  "motrpac_human-precovid_refmet-map_provenance.json"    = "06_generate_qc_norm/stage_refmet_map")
for (k in names(RESOURCE_STEMS)) {
  fid <- paste0("freeze:", k)
  sid <- paste0("stem:", RESOURCE_STEMS[[k]])
  if (!is.null(NODES[[fid]]) && !is.null(NODES[[sid]])) add_edge(sid, fid, "produces")
}

# ---- who reads the freeze --------------------------------------------------
# CURATED: these reads go through load_qc_local() or a list.files() sweep over
# staging/freeze, so no filename literal exists in the consuming source. Keyed by
# "category|details" first, falling back to the category. The one entry that is
# not a plain category rule is metadata|removed-samples: those files are produced
# upstream and copied verbatim into the freeze by stage_removed_samples.R so the
# release carries them, and nothing in Stage 1 reads them back — they go straight
# to the bucket. (OUTLIERS.R does read them, but from the vendored preflight
# sources, not from the freeze.)
FREEZE_CONSUMERS <- list(
  "metadata|removed-samples" = character(0),
  "metadata" = c("07_build_human_feature_to_gene", "08_build_qc_objects", "09_build_da"),
  "qc-norm"  = c("08_build_qc_objects", "09_build_da"),
  "imputed"  = c("08_build_qc_objects"),
  "da"       = c("11_build_sum_stats")
)
for (i in seq_len(nrow(fam))) {
  key <- paste0(fam$category[i], "|", fam$details[i])
  cons <- FREEZE_CONSUMERS[[key]]
  if (is.null(cons)) cons <- FREEZE_CONSUMERS[[fam$category[i]]]
  for (st in cons) {
    sid <- paste0("step:", st)
    if (!is.null(NODES[[sid]])) add_edge(fam$node[i], sid, "consumed_by")
  }
}

# =============================================================================
# 5. External inputs — external_assets.tsv + gmt_files.tsv
#
# Both manifests already name their consumer (consumed_by / builder_script), so
# the edges come straight out of the columns. Consumers are matched by script
# basename, with glob support so a manifest entry like "generate_*_qc_norm" fans
# out to every matching stem.
# =============================================================================
# name -> node id, for every node a manifest could plausibly name: script
# basenames, and the data-object groups the ledger defines (several manifest rows
# name the object an asset feeds rather than the script that reads it, e.g.
# production_bucket -> "METABOLOMICS_CVS; OUTLIERS").
script_index <- list()
for (id in names(NODES)) {
  n <- NODES[[id]]
  if (n$type %in% c("stem", "preflight_builder", "stage_driver", "test_script")) {
    key <- paste0(n$label, if (n$type %in% c("stem", "test_script")) ".R" else "")
    script_index[[key]] <- c(script_index[[key]], id)
    script_index[[n$label]] <- c(script_index[[n$label]], id)
  }
  if (n$type == "data_object") script_index[[n$label]] <- c(script_index[[n$label]], id)
}
match_consumers <- function(txt) {
  if (is.na(txt) || !nzchar(txt)) return(character(0))
  out <- character(0)
  for (tok in unlist(strsplit(txt, "[;,[:space:]]+"))) {
    tok <- basename(gsub("[()]", "", tok))
    if (nchar(tok) < 3) next                        # "all", "2", ... are not names
    if (grepl("[*]", tok)) {
      pat <- paste0("^", gsub("[*]", ".*", tok), "$")
      out <- c(out, unlist(script_index[grepl(pat, names(script_index))]))
    } else if (!is.null(script_index[[tok]])) {
      out <- c(out, script_index[[tok]])
    }
  }
  unique(out)
}

ext_tsv <- file.path(ROOT, "scripts/00_preflight/config/external_assets.tsv")
if (file.exists(ext_tsv)) {
  EA <- read.csv(ext_tsv, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(EA))) {
    a <- EA[i, ]
    id <- paste0("src:", a$id)
    add_node(id, a$id, "source", "external",
             origin = if (identical(a$gated, "yes")) "gated" else "open",
             detail = sprintf("%s | %s%s", a$kind, a$location,
                              if (identical(a$gated, "yes")) " | consortium-gated" else ""))
    for (cid in match_consumers(a$consumed_by)) add_edge(id, cid, "consumed_by")
    # Stage 0 is exactly "verify every external asset before anything expensive
    # runs", so each asset is checked by the preflight driver whether or not its
    # consumed_by text names something in this repo.
    if (!is.null(NODES[["sh:scripts/00_preflight/00_preflight.sh"]]))
      add_edge(id, "sh:scripts/00_preflight/00_preflight.sh", "checked_by")
  }
  message("assets     : ", nrow(EA), " external assets")
}

gmt_tsv <- file.path(ROOT, "scripts/00_preflight/config/gmt_files.tsv")
if (file.exists(gmt_tsv)) {
  GM <- read.csv(gmt_tsv, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(GM))) {
    g <- GM[i, ]
    id <- paste0("pfout:", g$gmt_file)
    add_node(id, g$gmt_file, "preflight_object", "stage0", origin = "auto",
             detail = sprintf("source: %s | %s", g$source_file, g$license_note))
    for (cid in match_consumers(g$builder_script)) add_edge(cid, id, "produces")
    # every GMT feeds the MOLECULAR_SIGNATURES build (step 03)
    if (!is.null(NODES[["step:03_build_molecular_signatures"]]))
      add_edge(id, "step:03_build_molecular_signatures", "consumed_by")
  }
  message("gmts       : ", nrow(GM), " GMT build assets")
}

# =============================================================================
# 6. check_required_inputs.sh — the Stage 1 hard gate
#
# The gate arrays are `label::absolute/path` pairs plus a few compgen globs. The
# labels alone are the useful node: they name the pre-existing inputs a Stage 1
# run must have before any generator runs.
# =============================================================================
CRI <- "scripts/10_build_data/check_required_inputs.sh"
if (has(CRI)) {
  cl <- rd(CRI)
  gate_id <- "gate:check_required_inputs"
  add_node(gate_id, "check_required_inputs", "gate", "stage1",
           detail = "hard gate: every pre-existing input Stage 1 consumes")
  add_edge(gate_id, paste0("sh:", BUILD_DATA), "requires")
  labels <- grab(cl, '"([A-Za-z0-9_:.-]+)::\\$\\{PIPELINE_ROOT\\}', 2L)
  for (lb in unique(labels)) {
    lid <- paste0("gatein:", lb)
    add_node(lid, lb, "gate_input", "stage0", origin = "required",
             detail = "required before any Stage 1 generator runs")
    add_edge(lid, gate_id, "requires")
    # preflight:X labels point at a preflight object built in stage 0
    if (grepl("^preflight:", lb)) {
      onm <- sub("^preflight:", "", lb)
      if (!is.null(NODES[[paste0("obj:", onm)]])) add_edge(paste0("obj:", onm), lid, "feeds")
    }
  }
  # the compgen-glob checks (methylcap beta-values / DA tables, ATAC passthrough)
  for (v in unique(grab(cl, 'record_check (?:PASS|FAIL) "([a-z_]+):\\$\\{tc\\}"', 2L))) {
    lid <- paste0("gatein:", v)
    add_node(lid, v, "gate_input", "external", origin = "vendored",
             detail = "vendored external table, one per tissue (fetch: make sources)")
    add_edge(lid, gate_id, "requires")
  }
  message("gate       : ", length(unique(labels)), " required-file labels")
}

# =============================================================================
# 7. The intra-repo function call graph (AST, no execution)
#
# Every .R under scripts/ is parsed with R's own parser. Function definitions and
# their line spans are extracted; each SYMBOL_FUNCTION_CALL is attributed to its
# innermost enclosing definition (or to the file, if top-level). The callee
# universe is restricted to functions defined in this repo, so library calls are
# excluded. Shared-library functions (scripts/10_build_data/lib/*.R) are the ones
# that matter here: they are what every stem leans on.
# =============================================================================
r_files <- list.files(file.path(ROOT, "scripts"), pattern = "\\.R$",
                      recursive = TRUE, full.names = TRUE)
r_files <- r_files[!grepl("/(deprecated|sources)/", r_files)]

file_node_id <- function(p) {
  r <- rel(p)
  for (id in names(NODES)) {
    n <- NODES[[id]]
    if (n$type %in% c("stem", "test_script") &&
        grepl(paste0("/", n$label, "\\.R$"), r)) return(id)
    if (n$type == "preflight_builder" && grepl(paste0("/", n$label, "$"), r)) return(id)
  }
  NA_character_
}

defs_all <- list(); calls_all <- list()
for (f in r_files) {
  pd <- tryCatch(utils::getParseData(parse(f, keep.source = TRUE)),
                 error = function(e) NULL)
  if (is.null(pd) || !nrow(pd)) next
  # A definition is `name <- function(...)`. In the parse tree the FUNCTION token's
  # parent is the function's own expr; that expr's parent is the assignment expr,
  # whose LHS SYMBOL sits one level deeper again. Only TOP-LEVEL definitions are
  # kept (assignment expr's parent == 0), so closures local to another function do
  # not become nodes.
  fn_tok <- pd[pd$token == "FUNCTION", ]
  defs <- data.frame(name = character(0), line1 = integer(0), line2 = integer(0))
  for (k in seq_len(nrow(fn_tok))) {
    pid <- fn_tok$parent[k]
    gp <- pd$parent[pd$id == pid]
    if (!length(gp) || gp <= 0) next
    if (!identical(pd$parent[pd$id == gp], 0L)) next          # not top-level
    kids <- pd$id[pd$parent == gp]
    cand <- pd[(pd$parent == gp | pd$parent %in% kids) &
               pd$token %in% c("SYMBOL", "STR_CONST") & pd$id < pid, ]
    if (!nrow(cand)) next
    nm <- gsub('^"|"$', "", cand$text[which.max(cand$id)])    # nearest LHS symbol
    span <- pd[pd$id == gp, ]
    if (!nrow(span)) next
    defs <- rbind(defs, data.frame(name = nm, line1 = span$line1, line2 = span$line2))
  }
  defs <- defs[!duplicated(defs$name), , drop = FALSE]
  if (nrow(defs)) defs_all[[rel(f)]] <- defs
  cl <- pd[pd$token == "SYMBOL_FUNCTION_CALL", c("text", "line1")]
  if (nrow(cl)) calls_all[[rel(f)]] <- cl
}
repo_fns <- unique(unlist(lapply(defs_all, function(d) d$name)))

# A shared library is any .R file that another .R file source()s. Deriving the
# set this way catches both scripts/10_build_data/lib/* and the step-local shared
# files (09_build_da/da_common.R and its modelling engine) without naming either.
SOURCE_LINES <- lapply(stats::setNames(r_files, rel(r_files)), function(f) {
  txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
  grep("source\\(", txt, value = TRUE)
})
sourced_basenames <- unique(unlist(lapply(names(SOURCE_LINES), function(f) {
  hits <- vapply(basename(r_files), function(b)
    any(grepl(b, SOURCE_LINES[[f]], fixed = TRUE)) && b != basename(f), logical(1))
  basename(r_files)[hits]
})))
LIB_FILES <- intersect(names(defs_all), rel(r_files)[basename(r_files) %in% sourced_basenames])

# lib file + lib function nodes, and the defines edges between them. The stage
# comes from the file's own path: 00_preflight/data-raw/lib/helpers.R is Stage 0
# shared code, not Stage 1, and mislabelling it turns every one of its `sourced_by`
# edges into a spurious stage1 -> stage0 crossing.
for (lf in LIB_FILES) {
  lid <- paste0("lib:", lf)
  add_node(lid, basename(lf), "lib_file", STAGE_OF(lf),
           detail = paste("shared source file:", dirname(lf)))
  for (nm in defs_all[[lf]]$name) {
    fid <- paste0("fn:", nm)
    add_node(fid, nm, "lib_fn", STAGE_OF(lf),
             detail = paste("defined in", basename(lf)))
    add_edge(lid, fid, "defines")
  }
}
LIB_FN <- unique(unlist(lapply(defs_all[LIB_FILES], function(d) d$name)))

# Calls into the shared library. A call site inside a top-level definition is
# attributed to that definition (so the lib's own internal wiring shows as fn ->
# fn); everything else is attributed to the file it sits in.
n_call_edges <- 0L
for (f in names(calls_all)) {
  file_src <- file_node_id(file.path(ROOT, f))
  if (is.na(file_src) && f %in% LIB_FILES) file_src <- paste0("lib:", f)
  if (is.na(file_src)) next
  d <- defs_all[[f]]; cc <- calls_all[[f]]
  for (i in seq_len(nrow(cc))) {
    callee <- cc$text[i]
    if (!(callee %in% LIB_FN)) next
    fid <- paste0("fn:", callee)
    if (is.null(NODES[[fid]])) next
    src <- file_src
    if (!is.null(d)) {
      encl <- d[d$line1 <= cc$line1[i] & d$line2 >= cc$line1[i], , drop = FALSE]
      if (nrow(encl)) {
        nm <- encl$name[which.max(encl$line1)]
        if (!is.null(NODES[[paste0("fn:", nm)]])) src <- paste0("fn:", nm)
      }
    }
    if (identical(src, fid)) next
    add_edge(src, fid, "calls"); n_call_edges <- n_call_edges + 1L
  }
}
# source("<lib>.R") edges: which files pull in which shared library
for (f in names(SOURCE_LINES)) {
  src <- file_node_id(file.path(ROOT, f))
  if (is.na(src)) src <- if (f %in% LIB_FILES) paste0("lib:", f) else NA_character_
  if (is.na(src)) next
  for (lf in LIB_FILES) {
    if (identical(lf, f)) next
    if (any(grepl(basename(lf), SOURCE_LINES[[f]], fixed = TRUE)))
      add_edge(paste0("lib:", lf), src, "sourced_by")
  }
}
message("functions  : ", length(repo_fns), " defined in-repo, ",
        length(LIB_FN), " in the shared lib, ", n_call_edges, " call edges")

# =============================================================================
# 8. Stage 2 / Stage 3 sinks — where the freeze ends up
# =============================================================================
penv <- rd("config", "pipeline.env")
bucket_of <- function(v) {
  ln <- grep(paste0("^", v, "="), penv, value = TRUE)
  if (!length(ln)) return(NA_character_)
  val <- gsub('^"|"$', "", sub(paste0("^", v, "="), "", ln[1]))
  val <- gsub("\\$\\{CURRENT_VERSION\\}", gsub('"', "", sub(".*=", "", grep("^CURRENT_VERSION=", penv, value = TRUE)[1])), val)
  val <- gsub("\\$\\{NEW_VERSION\\}", gsub('"', "", sub(".*=", "", grep("^NEW_VERSION=", penv, value = TRUE)[1])), val)
  val
}
for (v in c("STAGING_BUCKET", "NEW_PRODUCTION_BUCKET")) {
  b <- bucket_of(v)
  if (is.na(b)) next
  id <- paste0("sink:", v)
  add_node(id, v, "sink", "stage2", origin = "gated", detail = b)
}
upl <- grep("^sh:scripts/20_upload", names(NODES), value = TRUE)
if (length(upl)) {
  # every freeze file goes to the bucket, including the carried-through QC reports
  for (fid in grep("^freeze:", names(NODES), value = TRUE)) add_edge(fid, upl[1], "consumed_by")
  if (!is.null(NODES[["sink:STAGING_BUCKET"]])) add_edge(upl[1], "sink:STAGING_BUCKET", "produces")
  if (!is.null(NODES[["sink:NEW_PRODUCTION_BUCKET"]]))
    add_edge("sink:STAGING_BUCKET", "sink:NEW_PRODUCTION_BUCKET", "produces")
}
# Stage 3 consumes the object layer, whatever the stage is currently called. Its
# sinks are read off the driver rather than assumed, because the stage's target
# has changed: the figures stage drove ACUTE_REPO, the package-update stage drives
# DATA_PKG_REPO + ANALYSIS_PKG_REPO.
s3 <- grep("^sh:scripts/30_", names(NODES), value = TRUE)
if (length(s3)) {
  for (nm in obj_nodes) add_edge(paste0("obj:", nm), s3[1], "consumed_by")
  s3_txt <- rd(sub("^sh:", "", s3[1]))
  for (v in c("DATA_PKG_REPO", "ANALYSIS_PKG_REPO", "ACUTE_REPO")) {
    if (!any(grepl(v, s3_txt, fixed = TRUE))) next
    id <- paste0("sink:", v)
    add_node(id, sub("_REPO$", "", v), "sink", "stage3",
             detail = "downstream repo this stage writes into")
    add_edge(s3[1], id, "produces")
  }
}

# =============================================================================
# 9. Assemble, score, level
# =============================================================================
nodes <- dplyr::bind_rows(NODES)
edges <- dplyr::bind_rows(EDGES) %>% dplyr::distinct()
edges <- edges %>% dplyr::filter(from %in% nodes$id, to %in% nodes$id, from != to)

# Two graphs, deliberately:
#   flow  — everything except calls/defines. Used for LAYOUT, so the picture reads
#           make -> driver -> step -> stem -> freeze -> bucket left to right.
#   blast — dataflow only: `runs`/`orders`/`tested_by` are control flow, and
#           counting them makes `make precheck` the biggest node in the repo,
#           which is true and useless. n_downstream is the blast radius: change
#           this node, and this many nodes are downstream of the change.
flow <- edges %>% dplyr::filter(!type %in% c("calls", "defines"))
#           `calls` is flipped on the way in: a caller depends on its callee, so
#           the blast radius of a helper is everyone who calls it (and their
#           downstream), not the other way round.
blast <- edges %>%
  dplyr::filter(type %in% c("produces", "consumed_by", "feeds", "requires",
                            "defines", "calls", "sourced_by")) %>%
  dplyr::transmute(from = ifelse(type == "calls", to, from),
                   to   = ifelse(type == "calls", from, to))
gb <- igraph::graph_from_data_frame(blast[, c("from", "to")],
                                    vertices = data.frame(name = nodes$id), directed = TRUE)
reach <- igraph::distances(gb, mode = "out")
nodes$n_downstream <- vapply(nodes$id, function(i)
  sum(is.finite(reach[i, ])) - 1L, integer(1))
# Hotspot = top 5% of reach across ALL nodes. Taking the decile of the positive
# values instead flags a fifth of the graph, because every freeze file has the
# same moderate reach (its consuming steps, then the bucket) and they all pile up
# on the threshold. The 95th percentile isolates the genuinely wide nodes.
thr <- stats::quantile(nodes$n_downstream, 0.95, names = FALSE)
nodes$hotspot <- nodes$n_downstream >= max(thr, 1)

# Levels: a class floor, then relaxed forward so every node sits to the right of
# everything that flows into it. Bounded iteration so a cycle cannot hang it.
FLOOR <- c(make_target = 0, source = 0, preflight_builder = 1, preflight_object = 2,
           gate_input = 2, stage_driver = 3, gate = 3, lib_file = 3, lib_fn = 4,
           step = 5, stem = 6, data_object = 7, freeze_file = 7,
           test_script = 8, sink = 30)
nodes$level <- unname(ifelse(is.na(FLOOR[nodes$type]), 5, FLOOR[nodes$type]))
# steps and their stems fan out left-to-right in build order
si <- STEP_INDEX[sub("^step:", "", nodes$id)]
nodes$level <- ifelse(nodes$type == "step" & !is.na(si), 5 + 2 * si, nodes$level)
stem_step <- vapply(nodes$id, function(i)
  if (!is.null(STEPS[[i]])) STEPS[[i]]$step else NA_character_, "")
ss <- STEP_INDEX[stem_step]
nodes$level <- ifelse(nodes$type == "stem" & !is.na(ss), 6 + 2 * ss, nodes$level)

lv <- stats::setNames(nodes$level, nodes$id)
for (pass in 1:40) {
  moved <- FALSE
  for (k in seq_len(nrow(flow))) {
    need <- lv[[flow$from[k]]] + 1L
    if (lv[[flow$to[k]]] < need && lv[[flow$to[k]]] < 60) {
      lv[[flow$to[k]]] <- need; moved <- TRUE
    }
  }
  if (!moved) break
}
# Dense-rank the result: relaxation leaves long empty stretches (a chain through
# the object layer into Stage 3 reaches the sixties), and empty columns are dead
# horizontal space in a hierarchical layout. Ranking preserves the ordering while
# collapsing the gaps.
nodes$level <- as.integer(factor(as.integer(lv[nodes$id]),
                                 levels = sort(unique(as.integer(lv))))) - 1L

# =============================================================================
# 10. Render
# =============================================================================
STAGE_COLOR <- c(orchestration = "#455A64", stage0 = "#00897B", stage1 = "#7E57C2",
                 stage2 = "#FB8C00", stage3 = "#C2185B", external = "#9E9E9E")
TYPE_SHAPE <- c(make_target = "database", stage_driver = "box", step = "square",
                stem = "square", test_script = "star", lib_file = "square",
                lib_fn = "triangleDown", data_object = "dot", freeze_file = "dot",
                preflight_builder = "square", preflight_object = "dot",
                source = "diamond", gate = "hexagon", gate_input = "diamond",
                sink = "box")
TYPE_LABEL <- c(make_target = "make target", stage_driver = "stage driver",
                step = "build step", stem = "R generator stem",
                test_script = "test script", lib_file = "shared lib file",
                lib_fn = "shared lib function", data_object = "data object group",
                freeze_file = "freeze output", preflight_builder = "preflight builder",
                preflight_object = "preflight asset", source = "external source",
                gate = "required-inputs gate", gate_input = "required input",
                sink = "downstream sink")
EDGE_COLOR <- c(orders = "#455A64", runs = "#5C6BC0", produces = "#2E7D32",
                consumed_by = "#546E7A", feeds = "#EF6C00", requires = "#C62828",
                tested_by = "#00838F", calls = "#B0BEC5", defines = "#CFD8DC",
                sourced_by = "#90A4AE")

vis_nodes <- nodes %>%
  dplyr::mutate(
    shape = unname(ifelse(is.na(TYPE_SHAPE[type]), "dot", TYPE_SHAPE[type])),
    color = unname(ifelse(is.na(STAGE_COLOR[stage]), "#9E9E9E", STAGE_COLOR[stage])),
    borderWidth = ifelse(hotspot, 3, 1),
    size = 12 + pmin(n_downstream, 40) * 0.7,
    title = sprintf("<b>%s</b><br>%s | %s | %s<br>downstream reach: %d%s%s<br><code>%s</code>",
                    label, TYPE_LABEL[type], stage, ome, n_downstream,
                    ifelse(hotspot, " | <b>hotspot</b>", ""),
                    ifelse(nzchar(detail), paste0("<br>", detail), ""), id),
    group = stage
  ) %>%
  dplyr::select(id, label, shape, color, borderWidth, size, title, group, level)
vis_nodes$color <- lapply(seq_len(nrow(vis_nodes)), function(i)
  list(background = vis_nodes$color[i], border = if (nodes$hotspot[i]) "#E91E63" else "#37474F"))

vis_edges <- edges %>%
  dplyr::transmute(from, to, arrows = "to", title = type,
                   color = unname(ifelse(is.na(EDGE_COLOR[type]), "#B0BEC5", EDGE_COLOR[type])),
                   dashes = type %in% c("calls", "defines", "sourced_by"))

legend_nodes <- data.frame(
  label = names(STAGE_COLOR),
  color = unname(STAGE_COLOR),
  shape = "dot", stringsAsFactors = FALSE)

net <- visNetwork(vis_nodes, vis_edges,
                  main = "motrpac-human-presuspension-repro — pipeline dependency graph",
                  submain = sprintf("%d nodes, %d edges | make -> stage -> step -> stem -> freeze -> bucket | generated from source text only",
                                    nrow(nodes), nrow(edges)),
                  width = "100%", height = "1000px") %>%
  visHierarchicalLayout(direction = "LR", sortMethod = "directed",
                        levelSeparation = 260, nodeSpacing = 90) %>%
  visNodes(font = list(size = 15)) %>%
  visEdges(smooth = list(enabled = TRUE, type = "cubicBezier")) %>%
  visOptions(highlightNearest = list(enabled = TRUE, degree = 2, hover = TRUE),
             nodesIdSelection = list(enabled = TRUE, useLabels = TRUE),
             selectedBy = list(variable = "group", multiple = TRUE)) %>%
  visInteraction(hover = TRUE, tooltipDelay = 120, navigationButtons = TRUE) %>%
  visLegend(addNodes = legend_nodes, useGroups = FALSE, width = 0.12,
            main = "stage")
htmlwidgets::saveWidget(net, file.path(SCRIPT_DIR, "depgraph.html"),
                        selfcontained = TRUE, title = "motrpac-human-presuspension-repro dependency graph")
# saveWidget leaves the unbundled dependency dir behind even when it inlines
# everything; drop it so only the single self-contained file is committed.
unlink(file.path(SCRIPT_DIR, "depgraph_files"), recursive = TRUE)

# ---- Layer 2: collapsed to stage super-nodes + the hotspots -----------------
sum_nodes <- nodes %>%
  dplyr::count(stage, name = "n") %>%
  dplyr::transmute(id = stage, label = sprintf("%s\n(%d nodes)", stage, n),
                   color = unname(STAGE_COLOR[stage]), shape = "box",
                   size = 30, level = match(stage, c("external", "orchestration", "stage0",
                                                     "stage1", "stage2", "stage3")))
stage_of <- stats::setNames(nodes$stage, nodes$id)
sum_edges <- edges %>%
  dplyr::filter(!type %in% c("calls", "defines")) %>%
  dplyr::transmute(from = unname(stage_of[from]), to = unname(stage_of[to])) %>%
  dplyr::filter(from != to) %>%
  dplyr::count(from, to, name = "n") %>%
  dplyr::transmute(from, to, label = as.character(n), arrows = "to",
                   width = 1 + log1p(n), color = "#78909C")

# The hotspots hang off their stage one level to the right, on a single-line
# label. Stacked under the stage node with a two-line label they collide with
# each other at any sane node spacing.
top <- nodes %>% dplyr::arrange(dplyr::desc(n_downstream)) %>% utils::head(12)
hot_nodes <- top %>%
  dplyr::transmute(id = id, label = sprintf("%s  (%d)", label, n_downstream),
                   color = unname(STAGE_COLOR[stage]), shape = "dot",
                   size = 12 + pmin(n_downstream, 60) * 0.35,
                   title = sprintf("%s | %s | downstream reach: %d", label, type, n_downstream),
                   level = 6 + match(stage, c("external", "orchestration", "stage0",
                                              "stage1", "stage2", "stage3")))
hot_edges <- data.frame(from = top$stage, to = top$id, arrows = "to",
                        color = "#CFD8DC", dashes = TRUE, stringsAsFactors = FALSE)

net2 <- visNetwork(dplyr::bind_rows(sum_nodes, hot_nodes),
                   dplyr::bind_rows(sum_edges, hot_edges),
                   main = "motrpac-human-presuspension-repro — stage summary + hotspots",
                   submain = "stage edges: label = number of dataflow edges crossing that boundary. Trailing nodes = widest blast radius.",
                   width = "100%", height = "760px") %>%
  visHierarchicalLayout(direction = "LR", sortMethod = "directed",
                        levelSeparation = 320, nodeSpacing = 70) %>%
  visNodes(font = list(size = 16, align = "left")) %>%
  visEdges(smooth = list(enabled = TRUE, type = "cubicBezier"), font = list(size = 18)) %>%
  visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE))
htmlwidgets::saveWidget(net2, file.path(SCRIPT_DIR, "depgraph_summary.html"),
                        selfcontained = TRUE, title = "motrpac-human-presuspension-repro stage summary")
unlink(file.path(SCRIPT_DIR, "depgraph_summary_files"), recursive = TRUE)

# =============================================================================
# 11. Tables + console summary
# =============================================================================
utils::write.table(nodes, file.path(SCRIPT_DIR, "depgraph_nodes.csv"),
                   sep = ",", row.names = FALSE, qmethod = "double")
utils::write.table(edges, file.path(SCRIPT_DIR, "depgraph_edges.csv"),
                   sep = ",", row.names = FALSE, qmethod = "double")
utils::write.table(FREEZE, file.path(SCRIPT_DIR, "depgraph_freeze_files.csv"),
                   sep = ",", row.names = FALSE, qmethod = "double")

message("\n== nodes by type ==")
print(table(nodes$type))
message("\n== edges by type ==")
print(table(edges$type))
message("\n== hotspots (widest blast radius) ==")
print(nodes %>% dplyr::arrange(dplyr::desc(n_downstream)) %>%
        dplyr::select(label, type, stage, n_downstream) %>% utils::head(15),
      row.names = FALSE)
message("\nwrote: depgraph.html, depgraph_summary.html, depgraph_nodes.csv, ",
        "depgraph_edges.csv, depgraph_freeze_files.csv -> ", SCRIPT_DIR)
