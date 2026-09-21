#!/usr/bin/env Rscript
# Step 09 DA stem: epigen-methylcap-seq (blood/t03-edta, muscle/t06-muscle, adipose/t11-adipose).
#
# There is no upstream R build to adapt here, so nothing is vendored commented-out above.
# MethylCap-seq differential analysis is fit with a MALAX GLMM in an external,
# assay-specific pipeline (Yongchao Ge); that modeling code is not part of this repo, and
# the dream/voom engine every other ome shares cannot reproduce it. This stem therefore
# does for the DA tier exactly what generate_methylcap_qc_norm.R does for the qc-norm
# tier: it treats the externally produced tables as a vendored INPUT and copies them into
# the freeze verbatim.
#
# The vendored files already carry the standard BIC
# <tissue_code>_<ome>_da_<model>-<model_type>_<version> token, so once copied they sit
# alongside every other ome's DA output and step 10, the freeze tests and everything
# downstream pick methylcap up without needing a special case for it.
#
# They keep their upstream `malax-glmm-acute_v1.2` name rather than being rewritten to
# this pipeline's `dream-acute` token. Two reasons: the model genuinely is a MALAX GLMM and
# not dream, so a dream- token would misdescribe the contents; and the qc-norm tier
# already vendors its v1.2 beta-values under the same reasoning. The bucket-diff test
# keys on the version-stripped filename, so the pairing still works.
#
# Presence of the sources is enforced up front by
# scripts/10_build_data/check_required_inputs.sh. This stem still checks that every
# methylcap tissue in OME_TISSUE_CODE resolves to exactly one file, because that is what
# catches a PARTIAL download — the gate only proves the folder is non-empty per tissue.
#
# SCHEMA — these tables do NOT have the same columns as the dream-fit omes, because a
# different model produced them. Shared with .convert_dream_output(): assay, feature_id, t,
# AveExpr, p_value, adj_p_value, contrast, full_model. The effect size is
# `methylation_diff`, NOT `logFC`, and the dream-specific z.std /
# degrees_of_freedom / logLik columns are absent. Nothing is renamed here — the tables are
# passed through as produced. Whoever ports step 10 (DA_ASSEMBLE, still a stub) needs to
# decide how to reconcile that, rather than assuming a logFC column exists for every ome.
#
# Output: staging/freeze/epigenomics/da/
#   human-precovid-sed-adu_{t03-edta,t06-muscle,t11-adipose}_epigen-methylcap-seq_da_malax-glmm-acute_v1.2.txt

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
# NOTE: this is the one stem that does NOT source da_common.R, deliberately. That file
# pulls in the dream engine and the whole modelling chain. A verbatim file copy needs none
# of it, and depending on it would gate the methylcap freeze on steps 05-08 for no reason.
# Sourcing qc_helpers.R directly is exactly what the
# qc-norm counterpart generate_methylcap_qc_norm.R does, and it supplies everything used
# here: .STAGING, .freeze_subdir_for_ome and OME_TISSUE_CODE.
source(file.path(ROOT, "scripts", "10_build_data", "lib", "qc_helpers.R"))

generate_methylcap_da <- function(){
  desired_ome <- "epigen-methylcap-seq"
  # Same path da_common.R's .da_freeze_dir() builds; the ome -> ome-group mapping stays
  # shared via .freeze_subdir_for_ome() rather than being hardcoded to "epigenomics".
  da_path <- file.path(.STAGING, "freeze", .freeze_subdir_for_ome(desired_ome), "da")
  dir.create(da_path, recursive = TRUE, showWarnings = FALSE)

  # external artifact (Yongchao Ge; MALAX GLMM modeling not in this repo), vendored
  # locally. See sources/README.md for provenance and how to refresh.
  src_dir <- file.path(ROOT, "scripts", "10_build_data", "09_build_da", "sources", "methylcap_da")
  src_files <- list.files(src_dir, pattern = "methylcap.*_da_.*\\.txt$", full.names = TRUE)

  # tissue (blood/muscle/adipose) <-> tissue_code for methylcap, from the data object
  methyl_tissues <- OME_TISSUE_CODE[OME_TISSUE_CODE$ome == desired_ome, c("tissue", "tissue_code")]

  # Driven by the expected tissue list rather than by whatever happens to be in src_dir,
  # so a missing or duplicated tissue is an error instead of a silently shorter freeze.
  for (i in seq_len(nrow(methyl_tissues))) {
    tissue <- methyl_tissues$tissue[i]
    tissue_code <- methyl_tissues$tissue_code[i]

    hit <- src_files[grepl(tissue_code, basename(src_files), fixed = TRUE)]
    if (length(hit) != 1)
      stop("expected exactly 1 vendored methylcap DA file for ", tissue_code,
           " but found ", length(hit), " in ", src_dir)

    dest <- file.path(da_path, basename(hit))
    if (!file.copy(hit, dest, overwrite = TRUE))
      stop("copy failed: ", hit, " -> ", dest)
    # These tables are ~0.8-1.2 GB each. A truncated copy would still parse as a valid
    # table, just with features missing, so the size is compared rather than trusting
    # file.copy()'s return alone. No read.delim() here by design: the file is passed
    # through untouched, and parsing GBs only to report a row count would cost more
    # than the copy itself.
    if (file.size(dest) != file.size(hit))
      stop("size mismatch after copying ", basename(hit), " — source ", file.size(hit),
           " bytes, freeze copy ", file.size(dest), " bytes")

    message(sprintf("%s / %s: copied %s (%.2f GB) -> %s",
                    tissue, desired_ome, basename(hit), file.size(dest) / 2^30, da_path))
  }
}

generate_methylcap_da()
