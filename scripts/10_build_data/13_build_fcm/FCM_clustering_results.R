#!/usr/bin/env Rscript
# Stage 1 data objects: FCM_CLUSTERS, FCM_CAMERA, FCM_ORA  (deps: DA_ASSEMBLE,
# MOLECULAR_SIGNATURES, SET_TO_ID, HUMAN_FEATURE_TO_GENE, CONTRAST_CONVERTER)
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/FCM_clustering_results.R.
#
# Fuzzy c-means clustering of the acute-exercise differential-analysis results, then two
# enrichment tests over the resulting clusters:
#
#   FCM_CLUSTERS  run_cmeans()             features clustered by their z-score trajectory
#                                          across contrasts, one `fclust` per tissue
#   FCM_CAMERA    run_cluster_cameraPR()   CAMERA-PR on the cluster membership probabilities —
#                                          which signatures follow each centroid's trajectory
#   FCM_ORA       run_cluster_ORA()        ORA on the disjoint (hard-assigned) clusters
#
# Upstream's data-raw script is three one-line calls, so all the substance is in those three
# functions, vendored alongside this file (see each header for the changes).
#
# Output: scripts/10_build_data/data/{FCM_CLUSTERS,FCM_CAMERA,FCM_ORA}.rda
#
# ---- Adaptations from upstream ----------------------------------------------------------
#
# 1. DA_list is built here from the step-10 *_DA objects and passed into run_cmeans(), which
#    gains a DA_list argument for the purpose. Upstream leaves it NULL, which sends
#    .prepare_DA_results() to load_differential_analysis() and the installed analysis package —
#    the dependency this stage removes. Same change step 12 makes.
#
# 2. The molecular signatures, set metadata, feature->gene map and contrast metadata are the
#    pipeline's own objects (steps 03, 04, 07 and 10), read here as globals that the vendored
#    functions close over. Same pattern step 09 uses for COVARIATES_FILE.
#
# 3. The cluster count per tissue is a build parameter (FCM_K_* below) rather than something
#    typed at a prompt. Upstream's run_cmeans() will, when handed a range, plot Dmin and block
#    on readline(). fcm_diagnostics.R does that sweep non-interactively and writes the curve to
#    disk; the number chosen from it is set here, so the object does not depend on a console.
#
# 4. usethis::use_data() becomes a plain save(), as in every other adapted step.
#
# ---- Scope ------------------------------------------------------------------------------
#
# CLUSTERING uses five omes — run_cmeans()'s own default: transcript-rna-seq, prot-pr, prot-ol,
# prot-ph and metab. Only the `exercise_with_controls` contrasts are clustered, and the two
# `during` contrasts are dropped inside run_cmeans(), so features must be measured at every
# remaining timepoint (complete.cases) to be clustered at all. Adipose prot-pr and prot-ph are
# dropped by run_cmeans() for insufficient timepoints. Epigen and the two clinical DA objects
# are not clustered, and all three are dropped at the file glob below rather than by
# `selected_omes` — epigen because loading ~7 GB only to discard it is wasteful, and
# metab-t-clinical because `selected_omes` cannot tell it apart from the stacked metabolomics
# object for its tissue (see the glob for what that costs).
#
# ENRICHMENT covers four of those five — prot-ol is absent from .prepare_cluster_mem()'s
# choices upstream, so blood prot-ol features are clustered but not enrichment-tested. Both
# tests run over every collection in MOLECULAR_SIGNATURES, which is the upstream default and
# includes PTMSIGDB. PTMSigDB sets are keyed on flanking sequences carrying a ";u"/";d"
# direction suffix that only the PTMsigDB branch of .prepare_sets() strips, and that branch
# fires only when PTMSIGDB is the sole database; mixed in with the others its sets match
# nothing and drop out at the size filter. The released FCM_CAMERA has no PTMSIGDB rows either.
#
# PHOSPHOPROTEOMICS DEPENDS ON A STEP-07 FIX, exactly as step 12 does: cluster membership is
# keyed on the flanking sequence for prot-ph, because PhosphoSitePlus kinase sets are defined
# over flanking sequences. Step 07 retains that column; the check below fails loudly if the
# built object predates it.

suppressWarnings(suppressMessages({library(dplyr); library(data.table)}))

for (pkg in c("Mfuzz", "Biobase", "e1071", "TMSig"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(pkg, " is required by step 13 — install it before running this step.")

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
HERE    <- file.path(ROOT, "scripts", "10_build_data", "13_build_fcm")
STEP12  <- file.path(ROOT, "scripts", "10_build_data", "12_build_camera_results")
out_dir <- file.path(ROOT, "scripts", "10_build_data", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- the objects the vendored functions close over ---------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a
load_one <- function(name) {
  f <- file.path(out_dir, paste0(name, ".rda"))
  if (!file.exists(f)) stop("missing ", name, ".rda — run its build step first: ", f)
  e <- new.env(); load(f, envir = e); e[[name]] %||% e[[ls(e)[1]]]
}

MOLECULAR_SIGNATURES  <- load_one("MOLECULAR_SIGNATURES")    # step 03
SET_TO_ID             <- load_one("SET_TO_ID")               # step 04
HUMAN_FEATURE_TO_GENE <- load_one("HUMAN_FEATURE_TO_GENE")   # step 07
CONTRAST_CONVERTER    <- load_one("CONTRAST_CONVERTER")      # step 10

if (!"flanking_sequence" %in% colnames(HUMAN_FEATURE_TO_GENE))
  stop("HUMAN_FEATURE_TO_GENE has no flanking_sequence column — prot-ph cluster membership ",
       "keys on it (PhosphoSitePlus kinase sets are defined over flanking sequences). ",
       "Rebuild step 07.")

# .prepare_DA_results(), .create_index() and .prepare_sets() come from step 12's vendored
# copy rather than a second one here, so the two steps cannot drift on how DA results become
# z-score matrices or on how sets are filtered to a background.
source(file.path(STEP12, "run_cameraPR.R"))
source(file.path(HERE, "run_cmeans.R"))
source(file.path(HERE, "run_cluster_cameraPR.R"))
source(file.path(HERE, "run_cluster_ORA.R"))

# --- cluster counts ----------------------------------------------------------------------
# The pipeline clusters adipose at 13, blood at 12 and muscle at 12. Override after reading
# the Dmin curves fcm_diagnostics.R writes; a changed k changes FCM_CAMERA and FCM_ORA too.
k_from_env <- function(var, default) {
  v <- Sys.getenv(var)
  if (!nzchar(v)) return(default)
  n <- suppressWarnings(as.integer(v))
  if (is.na(n) || n < 2L) stop(var, " must be an integer >= 2, got ", dQuote(v))
  n
}
K_ADIPOSE <- k_from_env("FCM_K_ADIPOSE", 13L)
K_BLOOD   <- k_from_env("FCM_K_BLOOD",   12L)
K_MUSCLE  <- k_from_env("FCM_K_MUSCLE",  12L)

# --- DA_list: nested tissue -> ome, from the step-10 objects ------------------------------
# .prepare_DA_results() unlists this one level into "tissue.ome" names, so the nesting is the
# interface, not a convenience. Same exclusions as step 12's, plus prot-clinical.
da_files <- list.files(out_dir, pattern = "_DA\\.rda$", full.names = TRUE)
da_files <- da_files[!grepl("SPLICING_DA\\.rda$", da_files)]   # a 00_preflight leaf, not this tier
da_files <- da_files[!grepl("_EPIGEN_.*_DA\\.rda$", da_files)] # 20.4M rows, and not clustered
# The two clinical-chemistry objects are dropped HERE rather than by `selected_omes`.
# metab-t-clinical has to be: it carries assay = "metab", the same ome key as the stacked
# metabolomics object for its tissue, so it is indistinguishable by ome and would land on the
# same DA_list slot — leaving blood metab as 6 analytes, too few for any set to reach min_size.
# prot-clinical would fall out at `selected_omes`, but dropping it alongside keeps one rule for
# clinical chemistry, which is not clustered or enrichment-tested either way.
da_files <- da_files[!grepl("_METAB_T_CLINICAL_DA\\.rda$", da_files)]
da_files <- da_files[!grepl("_PROT_CLINICAL_DA\\.rda$", da_files)]
if (length(da_files) == 0)
  stop("no assembled *_DA objects in ", out_dir, " — run step 10 first")

DA_list <- list()
for (f in da_files) {
  e <- new.env(); load(f, envir = e); o <- as.data.frame(e[[ls(e)[1]]])
  tissue <- unique(as.character(o$tissue)); assay <- unique(as.character(o$assay))
  if (length(tissue) != 1 || length(assay) != 1)
    stop(basename(f), ": expected one tissue and one assay, got [",
         paste(tissue, collapse = ","), "] / [", paste(assay, collapse = ","), "]")
  # Two files claiming one tissue x ome slot would silently leave whichever loaded last, so the
  # collision is an error rather than a quiet substitution.
  if (!is.null(DA_list[[tissue]][[assay]]))
    stop(basename(f), ": ", tissue, " x ", assay, " is already filled by another *_DA object — ",
         "two objects share this ome key; exclude one at the glob above.")
  DA_list[[tissue]][[assay]] <- o
}
message(sprintf("DA_list: %d tissue(s), %d tissue x ome table(s)",
                length(DA_list), sum(lengths(DA_list))))

# --- FCM clustering ----------------------------------------------------------------------
FCM_CLUSTERS <- run_cmeans(DA_list = DA_list,
                           num_clusters_adipose = K_ADIPOSE,
                           num_clusters_blood   = K_BLOOD,
                           num_clusters_muscle  = K_MUSCLE)

save(FCM_CLUSTERS, file = file.path(out_dir, "FCM_CLUSTERS.rda"),
     compress = TRUE, version = 3)
for (t in names(FCM_CLUSTERS)) {
  x <- FCM_CLUSTERS[[t]]
  message(sprintf("FCM_CLUSTERS: %s — %d features x %d clusters | omes %s | sizes %s",
                  t, nrow(x$membership), ncol(x$membership),
                  paste(sort(unique(sub(" .*", "", rownames(x$membership)))), collapse = ","),
                  paste(x$size, collapse = ",")))
}
rm(DA_list); invisible(gc())

# --- enrichment of the clusters ----------------------------------------------------------
FCM_CAMERA <- run_cluster_cameraPR(FCM = FCM_CLUSTERS)
save(FCM_CAMERA, file = file.path(out_dir, "FCM_CAMERA.rda"),
     compress = TRUE, version = 3)
message(sprintf("FCM_CAMERA: %d rows x %d cols | %d tissue x ome | %d database(s)",
                nrow(FCM_CAMERA), ncol(FCM_CAMERA),
                nrow(unique(FCM_CAMERA[, c("tissue", "assay")])),
                dplyr::n_distinct(FCM_CAMERA$database)))

# ORA over the hard-assigned clusters. Not used beyond being compared to FCM_CAMERA, which is
# upstream's own note on this object; it is built because that comparison is the check that
# the fuzzy result is not an artefact of the membership weighting.
FCM_ORA <- run_cluster_ORA(FCM = FCM_CLUSTERS)
save(FCM_ORA, file = file.path(out_dir, "FCM_ORA.rda"),
     compress = TRUE, version = 3)
message(sprintf("FCM_ORA: %d rows x %d cols | %d tissue x ome | %d database(s) -> %s",
                nrow(FCM_ORA), ncol(FCM_ORA),
                nrow(unique(FCM_ORA[, c("tissue", "assay")])),
                dplyr::n_distinct(FCM_ORA$database), out_dir))
