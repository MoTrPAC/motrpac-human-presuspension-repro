#!/usr/bin/env Rscript
# Per-ome expected column schemas for structural validation (step 05).
#
# Vendored from MotrpacHumanPreSuspensionAnalysis/data-raw/google_cloud_bucket_checks/
# expected_columns.R. The schemas themselves are unchanged; the adaptations are at the
# bottom, in the lookup keys, and are listed there.
#
# Each entry is a named list with:
#   required  — columns that must be present
#   forbidden — columns that must NOT be present (catches stale schema)
#
# DA files all share the same core schema regardless of ome.
# qc-norm and metadata schemas are ome-specific.

# CI.L/CI.R were dropped for v2.0. variancePartition's topTable computes a dream fit's
# interval as se * qt(alpha, df = eb$df.total[top]); df.total is a features x contrasts
# matrix and `top` is a linear index, so every contrast's interval used the first
# contrast's df. Present upstream through 1.41.5, so the columns are not shipped.
DA_COLS_CORE = c(
  "feature_id", "contrast", "full_model",
  "logFC",
  "t", "p_value", "adj_p_value"
)

DA_COLS_PROT = c(DA_COLS_CORE, "AveExpr", "degrees_of_freedom", "logLik")

DA_COLS_METAB = c(DA_COLS_CORE, "AveExpr", "platform")

DA_COLS_METHYL = c(
  "feature_id", "contrast", "full_model",
  "methylation_diff", "p_value", "adj_p_value"
)

# Dropped for v2.0 and forbidden from here on, so a table built by an older pipeline —
# or by a topTable call that reverts to confint = TRUE — fails validation instead of
# shipping intervals computed from the wrong contrast's degrees of freedom.
#
# The dream tables do carry an interval, under CI.L_calculated / CI.R_calculated: the same
# bounds with df.total subset by feature AND contrast, built by run_dream() rather than by
# topTable. Those names are deliberately absent from both lists here. Not forbidden,
# because they are correct; not required, because this schema is also checked against
# bucket tables written before they existed, and a required column would fail every one of
# them. Their presence is asserted where it can be, on the freeze, in Stage 1 da_tests.R.
DA_COLS_DROPPED = c("CI.L", "CI.R")

QC_NORM_COLS_EXPR = c("feature_id") # + sample ID columns (checked by ncol)

# since we combine with pheno upon loading, all we need is the vialLabel for sample mapping.
METADATA_SAMPLE_COLS_CORE = c(
  "vialLabel"
)

# Participant identifiers, demographics and study design. No removed-samples file may ship them.
REMOVED_SAMPLES_COLS_IDENTIFYING = c("pid", "sex", "age", "race", "siteID", "short.id",
                                     "rdgroup", "visit", "timept")

METADATA_FEATURE_COLS_TRANSCRIPT = c(
  "assay", "feature_id", "gene_symbol", "ensembl_gene"
)

METADATA_FEATURE_COLS_PROT = c(
  "assay", "feature_id", "gene_symbol", "uniprot"
)

# prot-ph and prot-pr carry the vendor annotation their own .annotate_*() keeps, on top of the
# shared prot set. prot-ol and prot-clinical have none of it and stay on the shared set.
METADATA_FEATURE_COLS_PROT_PH = c(
  METADATA_FEATURE_COLS_PROT,
  "flanking_sequence", "confident_site", "confident_score", "ptm_score", "redundant_ids"
)

METADATA_FEATURE_COLS_PROT_PR = c(
  METADATA_FEATURE_COLS_PROT,
  "num_peptides", "percent_coverage", "protein_score", "redundant_ids"
)

METADATA_FEATURE_COLS_METAB = c(
  "assay", "feature_id", "refmet_name"
)

METADATA_FEATURE_COLS_EPIGEN = c(
  "assay", "feature_id", "gene_symbol", "ensembl_gene", "entrez_gene", "relationship_to_gene"
)

# Map from (data_category, ome) to the set of required columns.
# Keys are "<data_category>__<ome>" or "<data_category>__*" for shared schemas.
#
# ---- Adapted from upstream --------------------------------------------------------
#
# 1. The three clinical-chemistry keys are replaced by prot-clinical keys. The clinical
#    panels are now ordinary omes: metab-t-clinical is picked up by the ^metab fallback
#    in .get_expected_cols() and needs no key of its own, but prot-clinical does — the
#    fallback only fires for metab. Its DA table carries AveExpr, degrees_of_freedom and
#    logLik, so it takes the proteomics schema, not the core one upstream gave
#    clinical-chemistry.
#
# 2. imputed keys added for prot-pr / prot-ph. The imputed matrices have the same
#    feature_id + sample-column shape as a qc-norm matrix.
#
# 3. metadata__removed-samples forbids REMOVED_SAMPLES_COLS_IDENTIFYING. Upstream forbids
#    nothing there.
#
# 4. metadata__features__prot-ph and __prot-pr take their own required sets rather than the
#    shared proteomics one. Upstream gives all four prot omes METADATA_FEATURE_COLS_PROT, which
#    asserts nothing about the vendor annotation those two files carry from v2.1 on.
EXPECTED_COLUMNS = list(

  # DA schemas
  "da__transcript-rna-seq"  = list(required = DA_COLS_CORE,  forbidden = DA_COLS_DROPPED),
  "da__prot-pr"             = list(required = DA_COLS_PROT,  forbidden = DA_COLS_DROPPED),
  "da__prot-ph"             = list(required = DA_COLS_PROT,  forbidden = DA_COLS_DROPPED),
  "da__prot-ol"             = list(required = DA_COLS_CORE,  forbidden = DA_COLS_DROPPED),
  "da__prot-clinical"       = list(required = DA_COLS_PROT,  forbidden = DA_COLS_DROPPED),
  "da__metab"               = list(required = DA_COLS_METAB, forbidden = DA_COLS_DROPPED),
  "da__epigen-atac-seq"     = list(required = DA_COLS_CORE,  forbidden = DA_COLS_DROPPED),
  "da__epigen-methylcap-seq"= list(required = DA_COLS_METHYL,forbidden = DA_COLS_DROPPED),

  # qc-norm schemas (expression matrix: just need feature_id + at least 1 sample col)
  "qc-norm__transcript-rna-seq"   = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),
  "qc-norm__prot-pr"              = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),
  "qc-norm__prot-ph"              = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),
  "qc-norm__prot-ol"              = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),
  "qc-norm__prot-clinical"        = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),
  "qc-norm__epigen-atac-seq"      = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),
  "qc-norm__epigen-methylcap-seq" = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),
  "qc-norm__metab"                = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),

  # imputed matrices — same shape as qc-norm
  "imputed__prot-pr"              = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),
  "imputed__prot-ph"              = list(required = QC_NORM_COLS_EXPR, forbidden = character(0)),

  # metadata schemas
  "metadata__samples"         = list(required = METADATA_SAMPLE_COLS_CORE, forbidden = character(0)),
  # removed-samples files use either "vialLabel" or "sample" as the ID column
  "metadata__removed-samples" = list(required = character(0), forbidden = REMOVED_SAMPLES_COLS_IDENTIFYING),
  "metadata__features__transcript-rna-seq" = list(required = METADATA_FEATURE_COLS_TRANSCRIPT, forbidden = character(0)),
  "metadata__features__prot-pr"            = list(required = METADATA_FEATURE_COLS_PROT_PR,    forbidden = character(0)),
  "metadata__features__prot-ph"            = list(required = METADATA_FEATURE_COLS_PROT_PH,    forbidden = character(0)),
  "metadata__features__prot-ol"            = list(required = METADATA_FEATURE_COLS_PROT,       forbidden = character(0)),
  "metadata__features__prot-clinical"      = list(required = METADATA_FEATURE_COLS_PROT,       forbidden = character(0)),
  "metadata__features__epigen-atac-seq"      = list(required = METADATA_FEATURE_COLS_EPIGEN, forbidden = character(0)),
  "metadata__features__epigen-methylcap-seq" = list(required = METADATA_FEATURE_COLS_EPIGEN, forbidden = character(0)),
  "metadata__features__metab"              = list(required = METADATA_FEATURE_COLS_METAB,      forbidden = character(0))
)

# Helper: look up the expected columns for a given data_category + ome combination.
# Falls back to the generic metab key for any metab platform.
.get_expected_cols = function(data_category, ome, data_details = NULL) {
  key = paste0(data_category, "__", ome)
  if (key %in% names(EXPECTED_COLUMNS)) {
    return(EXPECTED_COLUMNS[[key]])
  }
  # metab fallback — try with data_details first (e.g. "metadata__features__metab"),
  # then without (e.g. "qc-norm__metab", "da__metab")
  if (grepl("^metab", ome)) {
    if (data_category == "metadata" && !is.null(data_details)) {
      key_metab_detail = paste0("metadata__", data_details, "__metab")
      if (key_metab_detail %in% names(EXPECTED_COLUMNS)) return(EXPECTED_COLUMNS[[key_metab_detail]])
    }
    key_metab = paste0(data_category, "__metab")
    if (key_metab %in% names(EXPECTED_COLUMNS)) {
      return(EXPECTED_COLUMNS[[key_metab]])
    }
  }
  # metadata samples vs features
  if (data_category == "metadata" && !is.null(data_details)) {
    key_meta = paste0("metadata__", data_details, "__", ome)
    if (key_meta %in% names(EXPECTED_COLUMNS)) return(EXPECTED_COLUMNS[[key_meta]])
    key_meta_shared = paste0("metadata__", data_details)
    if (key_meta_shared %in% names(EXPECTED_COLUMNS)) return(EXPECTED_COLUMNS[[key_meta_shared]])
  }
  return(NULL)
}
