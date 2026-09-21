#!/usr/bin/env Rscript
# Stage 1 step 16 — the SCION regulatory networks behind Figure 6 and Extended Data 8.
# (deps: QC objects (08), DA (10), HUMAN_FEATURE_TO_GENE (07), UTORONTO_TFs (14),
#  prot-ph imputed freeze matrix (06))
#
# Adapted from MotrpacHumanPreSuspensionAnalysis's run_SCION() and the driver that
# MotrpacPreSuspensionAcute/figures/landscape/figure_7/scion_figures.Rmd wrapped it in.
# The engine is vendored alongside: run_SCION.R, scion_matrixes.R, scion_run_cmeans.R.
#
# NOT a data object. The output is a set of edge tables under staging/scion/, not an .rda
# routed into either package: 2 networks x 2 exercise groups x (1 + SCION_PERMUTATIONS)
# runs is far past what lazy data can carry, and no accessor reads it. Figure 6B is a
# Cytoscape layout of a hand-merged export of these tables, and ED8A's TFEB target list
# comes out of the same session.
#
# Two networks, both per exercise group:
#   muscle_protph_transcript       prot-ph TF regulators -> DE transcript targets
#   blood_metab_muscle_transcript  blood metabolite regulators -> DE muscle transcript targets
#
# The observed network is inferred first, then SCION_PERMUTATIONS shuffled replicates.
# One permutation is a random-forest fit per cluster per group, so the published 100 is
# hours to days; submit it with EXECUTOR=slurm (config/slurm.json keys this step).
#
# Output: staging/scion/<network>/<group>/ and staging/scion/scion_manifest.tsv

suppressWarnings(suppressMessages({ library(dplyr) }))

for (pkg in c("Mfuzz", "Biobase", "e1071", "parallel", "doParallel", "randomForest"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(pkg, " is required by step 16 — install it before running this step.")
# Mfuzz::mfuzz() resolves e1071::cmeans() and Biobase::exprs() off the search path, and
# Mfuzz Depends on both rather than importing them.
suppressWarnings(suppressMessages({ library(Biobase); library(e1071); library(Mfuzz) }))

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
HERE   <- file.path(ROOT, "scripts", "10_build_data", "16_build_scion")
STEP12 <- file.path(ROOT, "scripts", "10_build_data", "12_build_camera_results")

source(file.path(ROOT, "scripts", "10_build_data", "lib", "qc_helpers.R"))
# .prepare_DA_results() comes from step 12's vendored copy rather than a second one here,
# so the two steps cannot drift on how DA results become z-score matrices.
source(file.path(STEP12, "run_cameraPR.R"))
source(file.path(HERE, "scion_matrixes.R"))
source(file.path(HERE, "scion_run_cmeans.R"))
source(file.path(HERE, "run_SCION.R"))

# CONTRAST_CONVERTER is a global .prepare_DA_results() closes over, as it does in step 13.
# Without it the clustering dies inside run_SCION() with "object 'CONTRAST_CONVERTER' not
# found" - after the matrices have loaded, which is the expensive part.
CONTRAST_CONVERTER <- local({
  f <- file.path(.STAGE1_DATA, "CONTRAST_CONVERTER.rda")
  if (!file.exists(f)) stop("missing CONTRAST_CONVERTER.rda - run step 10 first: ", f)
  e <- new.env(); load(f, envir = e); e[[ls(e)[1]]]
})

out_dir <- file.path(.STAGING, "scion")
manifest_tsv <- file.path(out_dir, "scion_manifest.tsv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- build parameters --------------------------------------------------------------------
# NULL permutations means the observed network and nothing else, which is what run_SCION()
# itself takes for `permute`.
permutations <- local({
  v <- Sys.getenv("SCION_PERMUTATIONS", unset = "")
  if (!nzchar(v) || toupper(v) == "NULL") return(NULL)
  n <- suppressWarnings(as.integer(v))
  if (is.na(n) || n < 0) stop("SCION_PERMUTATIONS must be NULL or a non-negative integer")
  if (n == 0L) NULL else n
})
num_cores <- local({
  n <- suppressWarnings(as.integer(Sys.getenv("SCION_CORES", "1")))
  if (is.na(n) || n < 1) 1L else n
})
force_rerun <- toupper(Sys.getenv("SCION_FORCE", "TRUE")) %in% c("TRUE", "1", "YES")

EXERCISE_GROUPS <- c("ADUResist", "ADUEndur")

# --- the networks ------------------------------------------------------------------------
#
# The two networks and their arguments come from
# MotrpacPreSuspensionAcute/figures/landscape/figure_7/scion_figures.Rmd, which ran the
# prot-ph network unpermuted and the blood-metabolomics network over 100 permutations on 11
# cores. Both are declared here and both run under one permutation plan; the permutation
# count and the core count are build parameters rather than literals.
#
# The prot-ph network is first because it is the one Figure 6 is drawn from: a run stopped
# partway through should leave that network rather than the other.
NETWORKS <- list(
  list(
    name = "muscle_protph_transcript",
    regulators = list(matrixes = "muscle.prot-ph", subset_TFs = TRUE, subset_DE = FALSE),
    targets    = list(matrixes = "muscle.transcript-rna-seq", subset_TFs = FALSE, subset_DE = TRUE)
  ),
  list(
    name = "blood_metab_muscle_transcript",
    regulators = list(matrixes = "blood.metabolomics", subset_TFs = FALSE, subset_DE = FALSE),
    targets    = list(matrixes = "muscle.transcript-rna-seq", subset_TFs = FALSE, subset_DE = TRUE)
  )
)

manifest_rows <- list()

#' Where one run's files land.
#'
#' run_SCION() appends the permutation number to dir.name itself, so the observed run and
#' permutation i share a parent and never overwrite one another.
run_dir <- function(network_name, group, permute) {
  base <- file.path(out_dir, network_name, group)
  if (is.null(permute)) base else file.path(base, as.character(permute))
}

#' The network table run_SCION() wrote for a run, if it wrote one.
#'
#' The prefix is taken from the regulator rownames inside run_SCION() and collapses every
#' metabolomics platform to "metab", so it is not knowable from the arguments here.
run_network_file <- function(dir_path, permute) {
  pattern <- if (is.null(permute)) {
    "-SCION-network-full\\.tsv$"
  } else {
    paste0("-SCION-network-permutation-", permute, "\\.tsv$")
  }
  hits <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (length(hits) == 0) NULL else hits[[1]]
}

#' Run one network x group x permutation, unless it is already on disk.
run_one <- function(network, group, permute) {
  tag <- paste0(network$name, ":", group, ":",
                if (is.null(permute)) "observed" else paste0("perm", permute))
  dir_path <- run_dir(network$name, group, permute)

  existing <- run_network_file(dir_path, permute)
  if (!force_rerun && !is.null(existing)) {
    message("[scion] ", tag, " — already inferred, skipping (SCION_FORCE=FALSE is set)")
    return(invisible(NULL))
  }

  message("[scion] ", tag)
  run_SCION(
    randomGroupCode = group,
    regulators = network$regulators_matrix,
    targets = network$targets_matrix,
    permute = permute,
    dir.name = file.path(out_dir, network$name, group),
    num.cores = num_cores
  )

  written <- run_network_file(dir_path, permute)
  if (is.null(written))
    stop("run_SCION() returned without writing a network table under ", dir_path)

  # colClasses: the hub layer writes an empty Cluster, which a type-guessing read turns into
  # NA — and as.character(NA) is "NA", which is nzchar TRUE, so the hub layer would count
  # below as a fourteenth cluster.
  edges <- read.csv(written, sep = "\t", check.names = FALSE, colClasses = "character")
  manifest_rows[[length(manifest_rows) + 1L]] <<- data.frame(
    network = network$name,
    randomGroupCode = group,
    permutation = if (is.null(permute)) NA_integer_ else permute,
    kind = if (is.null(permute)) "observed" else "null",
    regulators = network$regulators$matrixes,
    targets = network$targets$matrixes,
    file = sub(paste0("^", out_dir, "/"), "", written),
    n_edges = nrow(edges),
    n_regulators = dplyr::n_distinct(edges$Regulator),
    n_targets = dplyr::n_distinct(edges$Target),
    n_clusters = dplyr::n_distinct(edges$Cluster[nzchar(edges$Cluster)]),
    num_cores = num_cores,
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    inferred_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )
  message(sprintf("[scion] %s — %s, %d edges", tag, basename(written), nrow(edges)))
  invisible(NULL)
}

# --- drive -------------------------------------------------------------------------------
for (i in seq_along(NETWORKS)) {
  network <- NETWORKS[[i]]

  network$regulators_matrix <- .load_scion_matrixes(
    desired_matrixes = network$regulators$matrixes,
    subset_TFs = network$regulators$subset_TFs,
    subset_DE = network$regulators$subset_DE
  )
  network$targets_matrix <- .load_scion_matrixes(
    desired_matrixes = network$targets$matrixes,
    subset_TFs = network$targets$subset_TFs,
    subset_DE = network$targets$subset_DE
  )
  message(sprintf("[scion] %s — %d regulators x %d targets", network$name,
                  nrow(network$regulators_matrix), nrow(network$targets_matrix)))

  # The observed network first, then the null. Ordered this way on purpose: a run stopped
  # partway through leaves the network the figure needs, rather than shuffles of it.
  permute_plan <- c(list(NULL), if (!is.null(permutations)) as.list(seq_len(permutations)))

  for (permute in permute_plan) {
    for (group in EXERCISE_GROUPS) {
      run_one(network, group, permute)
    }
  }
}

# --- manifest ----------------------------------------------------------------------------
if (length(manifest_rows) > 0) {
  rows <- do.call(rbind, manifest_rows)
  if (file.exists(manifest_tsv)) {
    prior <- read.csv(manifest_tsv, sep = "\t", check.names = FALSE)
    prior <- prior[!prior$file %in% rows$file, , drop = FALSE]
    rows <- rbind(prior, rows)
  }
  rows <- rows[order(rows$network, rows$randomGroupCode, rows$permutation), , drop = FALSE]
  write.table(rows, file = manifest_tsv, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  message(sprintf("scion_manifest.tsv: %d network(s) recorded -> %s", nrow(rows), out_dir))
} else {
  message("no networks inferred this run — everything was already on disk")
}
