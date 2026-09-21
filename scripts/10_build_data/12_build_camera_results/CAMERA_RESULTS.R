#!/usr/bin/env Rscript
# Stage 1 data objects: CAMERA_RESULTS  (deps: DA_ASSEMBLE, MOLECULAR_SIGNATURES, SET_TO_ID,
# HUMAN_FEATURE_TO_GENE)
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/CAMERA_RESULTS.R.
#
# Pre-ranked CAMERA (Correlation Adjusted MEan RAnk) molecular-signature analysis of the
# differential-analysis results: for each tissue x ome and each contrast, test whether the
# z-statistics of the features in a gene/kinase/metabolite set are shifted relative to the
# rest. Upstream's data-raw script is a single line —
# `CAMERA_RESULTS <- MotrpacHumanPreSuspensionAnalysis::run_cameraPR()` — so the substance is
# in run_cameraPR(), vendored alongside this file (see its header for the changes).
#
# Output: scripts/10_build_data/data/CAMERA_RESULTS.rda
#
# ---- Adaptations from upstream ----------------------------------------------------------
#
# 1. DA_list is built here from the step-10 *_DA objects and passed in explicitly. Upstream
#    leaves it NULL, which makes run_cameraPR() call load_differential_analysis() and read
#    the objects out of the installed analysis package — the dependency this stage removes.
#    run_cameraPR() already takes a DA_list argument for exactly this case, so no behaviour
#    changes; it just gets its input from the local build instead of a library() call.
#
# 2. The molecular signatures, set metadata, feature->gene map and contrast metadata are the
#    pipeline's own objects (steps 03, 04, 07 and 10), read here as globals that the vendored
#    functions close over. Same pattern step 09 uses for COVARIATES_FILE.
#
# 3. usethis::use_data() becomes a plain save(), as in every other adapted step.
#
# ---- Scope ------------------------------------------------------------------------------
#
# Five omes are tested, which is run_cameraPR()'s own default and matches the released
# CAMERA_RESULTS: transcript-rna-seq, prot-pr, prot-ph, prot-ol and metab. Neither clinical
# chemistry object is enrichment-tested; between them they are 9 analytes, too few to reach a
# set of min_size. Both are dropped at the glob below. Epigen is assembled by step 10 and so
# does reach the output dir; it too is dropped at the glob, to avoid loading it first.
#
# Databases are every collection in MOLECULAR_SIGNATURES except PTMSIGDB, again the upstream
# default. PTMSIGDB is excluded there because its sets are keyed on flanking sequences with
# a ";u"/";d" direction suffix that only the .prepare_sets() PTMsigDB branch handles, and the
# released object carries no PTMSIGDB rows either.
#
# PHOSPHOPROTEOMICS DEPENDS ON A STEP-07 FIX. run_cameraPR() keys prot-ph enrichment on the
# flanking sequence rather than the feature_id, because PhosphoSitePlus kinase sets (the PSP
# database) are defined over flanking sequences. Step 07's select() used to drop that column
# from HUMAN_FEATURE_TO_GENE, which would fail this step for prot-ph; it is now retained.

suppressWarnings(suppressMessages({library(dplyr); library(data.table)}))

if (!requireNamespace("TMSig", quietly = TRUE))
  stop("TMSig is required by run_cameraPR() — install it before running step 12.")

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
HERE    <- file.path(ROOT, "scripts", "10_build_data", "12_build_camera_results")
out_dir <- file.path(ROOT, "scripts", "10_build_data", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- the objects the vendored functions close over ---------------------------------------
load_one <- function(name) {
  f <- file.path(out_dir, paste0(name, ".rda"))
  if (!file.exists(f)) stop("missing ", name, ".rda — run its build step first: ", f)
  e <- new.env(); load(f, envir = e); e[[name]] %||% e[[ls(e)[1]]]
}
`%||%` <- function(a, b) if (is.null(a)) b else a

MOLECULAR_SIGNATURES  <- load_one("MOLECULAR_SIGNATURES")    # step 03
SET_TO_ID             <- load_one("SET_TO_ID")               # step 04
HUMAN_FEATURE_TO_GENE <- load_one("HUMAN_FEATURE_TO_GENE")   # step 07
CONTRAST_CONVERTER    <- load_one("CONTRAST_CONVERTER")      # step 10

if (!"flanking_sequence" %in% colnames(HUMAN_FEATURE_TO_GENE))
  stop("HUMAN_FEATURE_TO_GENE has no flanking_sequence column — prot-ph enrichment keys on ",
       "it (PhosphoSitePlus kinase sets are defined over flanking sequences). Rebuild step 07.")

source(file.path(HERE, "run_cameraPR.R"))

# --- DA_list: nested tissue -> ome, from the step-10 objects ------------------------------
# run_cameraPR() unlists this one level into "tissue.ome" names, so the nesting is the
# interface, not a convenience. The ome key is the object's own assay — "metab" for the
# stacked metabolomics objects, whose platform stays in its own column and is not a separate
# dataset here (upstream tests metabolomics as one set of features per tissue).
da_files <- list.files(out_dir, pattern = "_DA\\.rda$", full.names = TRUE)
da_files <- da_files[!grepl("SPLICING_DA\\.rda$", da_files)]   # a 00_preflight leaf, not this tier
# Epigen is dropped HERE rather than by run_cameraPR()'s selected_omes, which is where the
# other non-enrichment objects fall out. It would fall out there too, but only after this
# loop had loaded ~7 GB into DA_list to hand it over — the five epigen objects are 20.4M rows.
da_files <- da_files[!grepl("_EPIGEN_.*_DA\\.rda$", da_files)]
# Clinical chemistry is not enrichment-tested. metab-t-clinical also carries assay = "metab",
# the same ome key as the stacked metabolomics object for its tissue, so it must go before the
# load rather than rely on `selected_omes`, which cannot tell the two apart.
da_files <- da_files[!grepl("_CLINICAL_DA\\.rda$", da_files)]
if (length(da_files) == 0)
  stop("no assembled *_DA objects in ", out_dir, " — run step 10 first")

DA_list <- list()
for (f in da_files) {
  e <- new.env(); load(f, envir = e); o <- as.data.frame(e[[ls(e)[1]]])
  tissue <- unique(as.character(o$tissue)); assay <- unique(as.character(o$assay))
  if (length(tissue) != 1 || length(assay) != 1)
    stop(basename(f), ": expected one tissue and one assay, got [",
         paste(tissue, collapse = ","), "] / [", paste(assay, collapse = ","), "]")
  # tissue x assay is the slot key, so two objects sharing one must not both be loaded: the
  # second would replace the first and the run would continue on the wrong feature set.
  if (!is.null(DA_list[[tissue]][[assay]]))
    stop(basename(f), ": ", tissue, " x ", assay, " is already filled — two DA objects share ",
         "this ome key. Exclude one at the glob above rather than let it overwrite the other.")
  DA_list[[tissue]][[assay]] <- o
}
message(sprintf("DA_list: %d tissue(s), %d tissue x ome table(s) — %s",
                length(DA_list), sum(lengths(DA_list)),
                paste(sprintf("%s[%s]", names(DA_list),
                              vapply(DA_list, function(x) paste(names(x), collapse = ","), character(1))),
                      collapse = " ")))

# --- CAMERA-PR ---------------------------------------------------------------------------
CAMERA_RESULTS <- run_cameraPR(DA_list = DA_list)

save(CAMERA_RESULTS, file = file.path(out_dir, "CAMERA_RESULTS.rda"),
     compress = TRUE, version = 3)
message(sprintf("CAMERA_RESULTS: %d rows x %d cols | %d tissue x ome | %d database(s) -> %s",
                nrow(CAMERA_RESULTS), ncol(CAMERA_RESULTS),
                nrow(unique(CAMERA_RESULTS[, c("tissue", "assay")])),
                dplyr::n_distinct(CAMERA_RESULTS$database), out_dir))
