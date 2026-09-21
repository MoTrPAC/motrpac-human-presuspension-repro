# Where each built data object belongs, and why.
#
# One rule, applied to the object name, plus one exception list. Stage 3 routes
# from here and `docs/data_objects.tsv` records the same answer, so the two
# cannot drift: `02_route_objects.R` re-derives the map and fails if the
# inventory disagrees.
#
# The rule — the data package carries per-sample measurements and the phenotype
# and clinical tables they are keyed against; the analysis package carries what
# was computed from them.
#   cln_* / pheno / blood_transcript_*  -> MotrpacHumanPreSuspensionData
#   *_QC                                -> MotrpacHumanPreSuspensionData
#   everything else                     -> MotrpacHumanPreSuspensionAnalysis
#
# The exception
#   epigenomics *_QC / *_DA  -> neither; they ship through the GCS bucket.
#
# Why that exception and not a size threshold: the ten epigenomics QC and DA
# objects run 57-315 MB, and GitHub rejects a blob over 100 MB while neither
# package uses git-lfs. Four of the ten would fit under that limit, but
# `load_qc(epigen = TRUE)` and `load_differential_analysis()` already read every
# epigenomics tier from the bucket, so splitting the tier across two transports
# by file size would make the loaders' source depend on how well a table
# happened to compress. The tier moves together.
#
# Epigenomics *_SUM_STATS are a different matter and do ship: they are 0.1-0.9 MB
# because a summary statistic collapses samples, and the package already carries
# them.

EPIGEN_ASSAYS <- c("EPIGEN_ATAC_SEQ", "EPIGEN_METHYLCAP_SEQ")

#' Route one object name to its destination
#'
#' @param object character; the data object's name.
#' @return a one-row data.frame: object, destination, reason. `destination` is
#'   one of "data", "analysis", "bucket".
route_object <- function(object) {
  is_epigen <- any(vapply(EPIGEN_ASSAYS, function(a) grepl(a, object, fixed = TRUE), logical(1)))
  tier <- if (grepl("_QC$", object)) "QC"
          else if (grepl("_DA$", object)) "DA"
          else if (grepl("_SUM_STATS$", object)) "SUM_STATS"
          else "other"

  if (is_epigen && tier %in% c("QC", "DA")) {
    return(data.frame(object = object, destination = "bucket",
                      reason = "epigenomics QC/DA tier; 57-315 MB and loaded from GCS by the package loaders",
                      stringsAsFactors = FALSE))
  }
  if (tier == "QC") {
    return(data.frame(object = object, destination = "data",
                      reason = "QC tier ships in the data package",
                      stringsAsFactors = FALSE))
  }
  # Sample-level tables that are inputs to the analysis rather than results of
  # it: the clinical CRF exports, the phenotype table they key against, and the
  # two blood transcript matrices split only because of their size.
  if (grepl("^cln_", object) || object == "pheno" || grepl("^blood_transcript_", object)) {
    return(data.frame(object = object, destination = "data",
                      reason = "clinical, phenotype or sample-level input table",
                      stringsAsFactors = FALSE))
  }
  data.frame(object = object, destination = "analysis",
             reason = paste0(tolower(tier), " tier ships in the analysis package"),
             stringsAsFactors = FALSE)
}

#' Route a vector of object names
route_objects <- function(objects) {
  do.call(rbind, lapply(objects, route_object))
}
