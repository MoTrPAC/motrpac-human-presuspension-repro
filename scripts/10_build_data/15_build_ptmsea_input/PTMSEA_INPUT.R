#!/usr/bin/env Rscript
# Stage 1 build step 15: PTMSEA_INPUT  (deps: DA_ASSEMBLE (prot-ph), QC_OBJECTS (prot-ph))
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/preprocess_PTMSEA.R.
#
# Formats the phosphoproteomics differential-analysis results as a GCT matrix of
# site-level z-scores (one column per selected contrast, one row per confidently
# localized phosphosite) and writes it as a .gct text file.
#
# ---- WHERE THIS OUTPUT GOES, AND WHAT THIS REPO DOES NOT DO -----------------------------
#
# The .gct written here is an INPUT FILE, not a result. PTM-SEA itself is run from the
# PTM-SEA / ssGSEA2.0 DOCKER IMAGE MAINTAINED BY THE BROAD INSTITUTE, against PTMsigDB.
# THAT RUN IS OUT OF SCOPE FOR THIS REPO: motrpac-human-presuspension-repro regenerates data objects, and this
# step's contract ends at a valid .gct on disk. Nothing downstream in Stage 1, 2 or 3 reads
# the enrichment scores, no step here pulls or runs the image, and there is no
# reproducibility claim over the Broad container's version, its PTMsigDB release, or its
# parameters. Hand the .gct to that container by hand and treat its output as an external
# artifact. See this folder's README for the handoff.
#
# Output (both gitignored, like every other scripts/**/data/ artifact):
#   scripts/10_build_data/data/PTMSEA_INPUT_<tissue>.gct   <- the deliverable
#   scripts/10_build_data/data/PTMSEA_INPUT.rda            <- the same GCT objects, named list
#
# ---- Adaptations from upstream ----------------------------------------------------------
#
# 1. The DA and QC objects are the pipeline's own locally built globals (steps 10 and 08)
#    instead of MotrpacHumanPreSuspensionAnalysis:: / MotrpacHumanPreSuspensionData::
#    lookups. Same pattern step 12 uses. See preprocess_PTMSEA.R for the marked changes.
#
# 2. confident_site is read straight off *_PROT_PH_QC$feature_metadata, as upstream does.
#    That column used to be dropped before the freeze was ever written — step 06's
#    .annotate_prot_ph() select()d down to assay/feature_id/entrez_gene/gene_symbol/
#    ensembl_gene/uniprot/flanking_sequence — which is why upstream's preprocess_PTMSEA()
#    does not run against the shipped Data-package object either. Step 06 now keeps the
#    column (it was always on the raw ratio-results rdesc, so nothing had to be recomputed)
#    and step 08 carries it into the *_QC objects. So this step reads it, and does not
#    reconstruct it. HUMAN_FEATURE_TO_GENE (step 07) does not carry it and is not read here:
#    that map has no tissue column, and confident_site is per tissue.
#
# 3. Both tissues are built by default, not just upstream's `selected_tissues = "muscle"`
#    default. The function still takes one tissue at a time; this stem loops. Narrow it with
#    PTMSEA_TISSUES if you only want one.
#
# 4. usethis::use_data() has no analogue here (upstream never saved this object at all —
#    preprocess_PTMSEA() is a package function whose caller writes the file). The .gct is
#    written with cmapR::write_gct(appenddim = FALSE), which is what upstream's own @examples
#    block does, and a plain save() keeps the GCT objects for the tests.
#
# ---- Build parameters (env) -------------------------------------------------------------
#   PTMSEA_TISSUES             space-separated; default "muscle adipose"
#   PTMSEA_CONTRAST_TYPE       space-separated; default "exercise_with_controls"
#   PTMSEA_CONTRAST_CATEGORY   space-separated; default "EE-CON RE-CON"

suppressWarnings(suppressMessages({ library(dplyr) }))

if (!requireNamespace("cmapR", quietly = TRUE))
  stop("cmapR is required to build a GCT — install it before running step 15.")
if (!requireNamespace("tidyr", quietly = TRUE))
  stop("tidyr is required by preprocess_PTMSEA() — install it before running step 15.")

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
HERE      <- file.path(ROOT, "scripts", "10_build_data", "15_build_ptmsea_input")
out_dir   <- file.path(ROOT, "scripts", "10_build_data", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

.env_words <- function(name, default) {
  v <- Sys.getenv(name)
  if (!nzchar(v)) return(default)
  strsplit(trimws(v), "[[:space:]]+")[[1]]
}
TISSUES            <- .env_words("PTMSEA_TISSUES", c("muscle", "adipose"))
CONTRAST_TYPE      <- .env_words("PTMSEA_CONTRAST_TYPE", "exercise_with_controls")
CONTRAST_CATEGORY  <- .env_words("PTMSEA_CONTRAST_CATEGORY", c("EE-CON", "RE-CON"))

bad <- setdiff(TISSUES, c("muscle", "adipose"))
if (length(bad))
  stop("PTMSEA_TISSUES: prot-ph exists for muscle and adipose only — got ", paste(bad, collapse = ", "))

# --- the objects preprocess_PTMSEA() closes over ------------------------------------------
load_one <- function(name) {
  f <- file.path(out_dir, paste0(name, ".rda"))
  if (!file.exists(f)) stop("missing ", name, ".rda — run its build step first: ", f)
  e <- new.env(); load(f, envir = e); e[[name]] %||% e[[ls(e)[1]]]
}

# --- confident_site gate ------------------------------------------------------------------
# preprocess_PTMSEA() subsets with prot_da_wide[prot_da_wide$confident_site, ]. Both failure
# modes of that expression are silent, so both are checked here instead of downstream:
#   - column absent  -> NULL subscript, which drops EVERY row and yields an empty GCT;
#   - value NA       -> an all-NA ROW that survives into the matrix rather than being
#                       removed, which PTM-SEA then reads as real data.
# A logical column with no NAs is the contract; anything else stops the build.
.require_confident_site <- function(qc, tissue) {
  fm <- qc$feature_metadata
  if (!"confident_site" %in% colnames(fm))
    stop(tissue, ": feature_metadata has no confident_site column — rebuild the prot-ph ",
         "freeze and QC objects (steps 06 then 08). It is written by step 06's ",
         ".annotate_prot_ph(); a freeze built before that column was kept will not have it.")
  cs <- fm$confident_site
  if (!is.logical(cs)) {
    # read_freeze() type-converts TRUE/FALSE to logical on its own; a character or factor
    # here means the freeze file holds something else, and as.logical() of a factor is NA.
    stop(tissue, ": feature_metadata$confident_site is ", class(cs)[1],
         ", expected logical — check what step 06 wrote to the metadata_features freeze file.")
  }
  if (anyNA(cs))
    stop(sprintf("%s: confident_site is NA for %d of %d feature(s) — ", tissue, sum(is.na(cs)), length(cs)),
         "every prot-ph feature must be flagged either way; do not build a GCT from this.")
  message(sprintf("  %s: %d of %d sites confidently localized", tissue, sum(cs), length(cs)))
  qc
}

MUSCLE_PROT_PH_DA  <- load_one("MUSCLE_PROT_PH_DA")     # step 10
ADIPOSE_PROT_PH_DA <- load_one("ADIPOSE_PROT_PH_DA")    # step 10
MUSCLE_PROT_PH_QC  <- .require_confident_site(load_one("MUSCLE_PROT_PH_QC"),  "muscle")    # step 08
ADIPOSE_PROT_PH_QC <- .require_confident_site(load_one("ADIPOSE_PROT_PH_QC"), "adipose")   # step 08

source(file.path(HERE, "preprocess_PTMSEA.R"))

# --- build one GCT per tissue -------------------------------------------------------------
PTMSEA_INPUT <- list()
for (tissue in TISSUES) {
  message(sprintf("PTM-SEA input: %s [%s / %s]", tissue,
                  paste(CONTRAST_TYPE, collapse = ","), paste(CONTRAST_CATEGORY, collapse = ",")))
  gct <- preprocess_PTMSEA(selected_tissues = tissue,
                           selected_contrast_type = CONTRAST_TYPE,
                           selected_contrast_category = CONTRAST_CATEGORY)
  if (nrow(gct@mat) == 0 || ncol(gct@mat) == 0)
    stop(tissue, ": empty GCT (", nrow(gct@mat), " x ", ncol(gct@mat),
         ") — check PTMSEA_CONTRAST_TYPE / PTMSEA_CONTRAST_CATEGORY against the DA object.")

  ofile <- file.path(out_dir, paste0("PTMSEA_INPUT_", tissue, ".gct"))
  # appenddim = FALSE: keep the filename stable across runs. cmapR otherwise stamps
  # _n<cols>x<rows> into it, which would rename the file whenever a contrast is added.
  cmapR::write_gct(gct, ofile = ofile, appenddim = FALSE)
  PTMSEA_INPUT[[tissue]] <- gct
  message(sprintf("  -> %s (%d sites x %d contrasts)",
                  basename(ofile), nrow(gct@mat), ncol(gct@mat)))
}

save(PTMSEA_INPUT, file = file.path(out_dir, "PTMSEA_INPUT.rda"), compress = TRUE, version = 3)
message(sprintf("PTMSEA_INPUT: %d tissue(s) — %s -> %s",
                length(PTMSEA_INPUT),
                paste(sprintf("%s[%dx%d]", names(PTMSEA_INPUT),
                              vapply(PTMSEA_INPUT, function(g) nrow(g@mat), integer(1)),
                              vapply(PTMSEA_INPUT, function(g) ncol(g@mat), integer(1))),
                      collapse = " "),
                out_dir))
message("NOTE: these .gct files are inputs to the Broad Institute PTM-SEA Docker image. ",
        "Running it is NOT part of this repo — see 15_build_ptmsea_input/README.md.")
