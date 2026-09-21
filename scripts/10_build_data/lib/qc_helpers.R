# Vendored shared internals for the Stage 1 QC-norm stems.
# Copied/adapted from the Analysis + Data packages so the stems run without an
# installed package:
#   .../data-raw/gsutil_path_parsing.R
#     write_with_path_name, .match_ome_tissue_code, .find_path_name,
#     .clinical_ome_for_feature, .freeze_subdir_for_ome
#   .../data-raw/generate_differential_analysis/generate_differential_modeling_functions.R
#     process_covariates   (vendored into 09_build_da/da_common.R on this side)
#   MotrpacHumanPreSuspensionData/R/load_pheno.R
#     load_pheno
#
# The in-house data objects these reference (OME_TISSUE_CODE, OUTLIERS,
# COVARIATES_FILE, pheno) are read from local .rda/.rds (preflight outputs +
# this stage's data/) instead of from an installed package. The caller must set
# `.QC_REPO` to the motrpac-human-presuspension-repro root before sourcing this file.
#
# The BIC metabolomics tools are vendored separately, in lib/bic_norm_qc_tools.R.

# Repo root comes from the PRECOVID_ROOT config env var (config/pipeline.env);
# fall back to searching upward for config/pipeline.env so standalone runs work.
precovid_root <- function() {
  r <- Sys.getenv("PRECOVID_ROOT")
  if (nzchar(r)) return(normalizePath(r))
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "config", "pipeline.env")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config", "pipeline.env")))
    stop("cannot find motrpac-human-presuspension-repro root; set PRECOVID_ROOT")
  d
}
suppressWarnings(suppressMessages(library(dplyr)))

.QC_REPO        <- precovid_root()
.PREFLIGHT_DATA <- file.path(.QC_REPO, "scripts", "00_preflight", "data")
.STAGE1_DATA    <- file.path(.QC_REPO, "scripts", "10_build_data", "data")
.STAGING        <- file.path(.QC_REPO, "staging")           # Stage 1 build intermediates
.RAW_FILES      <- file.path(.STAGING, "raw-files")         # scratch for gated dl_read_gcp downloads
dir.create(.RAW_FILES, showWarnings = FALSE, recursive = TRUE)

# --- parallelism switch for every variancePartition path --------------------
# Covers the step-06 ATAC voomWithDreamWeights and the step-09 dream fits alike: both hand
# variancePartition a BiocParallel MulticoreParam (fork) backend. Cores come from
# VARIANCEPARTITION_PARALLEL_CORES in config/pipeline.env, which explains why fork and not the
# SOCK backend it replaced. Deliberately NOT PARALLEL_CORES — that drives mice, unaffected here.
variancepartition_parallel_enabled <- function(){
  n <- suppressWarnings(as.integer(Sys.getenv("VARIANCEPARTITION_PARALLEL_CORES", "")))
  return(!is.na(n) && n > 1)
}

# --- ATAC rerun switch ------------------------------------------------------
# RERUN_ATAC in config/pipeline.env. FALSE makes the step-06 and step-09 ATAC stems
# stage the released qc-norm and DA tables instead of fitting them. Unset means TRUE,
# so a caller that never loaded pipeline.env still gets the full rebuild.
rerun_atac_enabled <- function(){
  raw <- Sys.getenv("RERUN_ATAC", "TRUE")
  v <- toupper(trimws(raw))
  if (!nzchar(v)) return(TRUE)
  if (v %in% c("TRUE", "T", "1", "YES")) return(TRUE)
  if (v %in% c("FALSE", "F", "0", "NO")) return(FALSE)
  stop("RERUN_ATAC must be TRUE or FALSE, got: ", raw)
}

# Per-file freeze versions: provides freeze_version(), freeze_files() and
# FREEZE_STEM_PREFIX, backed by config/file_versions.json — the hand-maintained version
# each file will carry in the next staging-bucket upload, held as a
# [[tissue_code]][[ome]][[data_category]][[data_details]] tree. See that file's header
# for why the versions are a per-file map rather than one NEW_VERSION.
source(file.path(.QC_REPO, "scripts", "10_build_data", "lib", "file_versions.R"))

# --- local data objects (were <pkg>::OBJECT) --------------------------------
OME_TISSUE_CODE <- readRDS(file.path(.PREFLIGHT_DATA, "OME_TISSUE_CODE.rds"))
OUTLIERS        <- readRDS(file.path(.PREFLIGHT_DATA, "OUTLIERS.rds"))
COVARIATES_FILE <- readRDS(file.path(.PREFLIGHT_DATA, "COVARIATES_FILE.rds"))
local({ e <- new.env(); load(file.path(.STAGE1_DATA, "pheno.rda"), envir = e); pheno <<- e$pheno })

# --- resolve raw quant-id input paths from the preflight catalog ------------
# The stems' raw inputs live on the (gated) quant-id bucket; rather than hard-
# coding gs:// paths, look them up in QUANTID_BUCKET_FILES by their tags.
QUANTID_BUCKET_FILES <- readRDS(file.path(.PREFLIGHT_DATA, "QUANTID_BUCKET_FILES.rds"))
quantid_path <- function(assay_class, tissue_code = NULL, data_details = NULL, content_type = NULL){
  hit <- QUANTID_BUCKET_FILES[QUANTID_BUCKET_FILES$assay_class == assay_class, , drop = FALSE]
  if (!is.null(tissue_code))  hit <- hit[hit$tissue_code  == tissue_code,  , drop = FALSE]
  if (!is.null(data_details)) hit <- hit[hit$data_details == data_details, , drop = FALSE]
  if (!is.null(content_type)) hit <- hit[hit$content_type == content_type, , drop = FALSE]
  hit <- hit[!is.na(hit$gcs_path), , drop = FALSE]
  if (nrow(hit) != 1)
    stop(sprintf("quantid_path: expected 1 file, got %d (assay_class=%s tissue_code=%s data_details=%s content_type=%s)",
                 nrow(hit), assay_class, tissue_code, data_details, content_type))
  hit$gcs_path
}

# --- local Ensembl v105 annotation cache ------------------------------------
# The stems annotate features against Ensembl release 105. The live Ensembl
# connection is flaky, so we read a local cache built once by
# sources/ensembl_v105/build_ensembl_v105_cache.R (see that folder's README).
# All Ensembl access at run time goes through the two accessors below — no
# biomaRt::useEnsembl / txdbmaker::makeTxDbFromEnsembl during a pipeline run.
.ENSEMBL_V105_DIR  <- file.path(.QC_REPO, "scripts", "00_preflight", "data-raw", "sources", "ensembl_v105")
.ENSEMBL_V105_TXDB <- file.path(.ENSEMBL_V105_DIR, "txdb_hsapiens_ensembl_v105.sqlite")
.ENSEMBL_V105_GENE <- file.path(.ENSEMBL_V105_DIR, "ensembl_v105_gene_attributes.rds")

# Load the cached ChIPseeker TxDb (was txdbmaker::makeTxDbFromEnsembl(release=105)).
.ensembl_v105_txdb <- function() {
  if (!file.exists(.ENSEMBL_V105_TXDB))
    stop("Ensembl v105 TxDb cache not found: ", .ENSEMBL_V105_TXDB,
         "\n  build it: Rscript ", file.path(.ENSEMBL_V105_DIR, "build_ensembl_v105_cache.R"))
  AnnotationDbi::loadDb(.ENSEMBL_V105_TXDB)
}

# Drop-in for biomaRt::getBM(attributes, filters, values, mart) against the cached
# v105 gene table. `unique()` reproduces the distinct rows a narrower live getBM
# returns (the cached table is the wider uniprot/entrez-expanded union of columns).
ensembl_v105_getBM <- function(attributes, filters, values) {
  if (!file.exists(.ENSEMBL_V105_GENE))
    stop("Ensembl v105 gene cache not found: ", .ENSEMBL_V105_GENE,
         "\n  build it: Rscript ", file.path(.ENSEMBL_V105_DIR, "build_ensembl_v105_cache.R"))
  tab <- readRDS(.ENSEMBL_V105_GENE)
  if (!filters %in% names(tab)) stop("ensembl_v105_getBM: unknown filter column '", filters, "'")
  miss <- setdiff(attributes, names(tab))
  if (length(miss)) stop("ensembl_v105_getBM: attribute(s) not cached: ", paste(miss, collapse = ", "))
  hit <- tab[tab[[filters]] %in% values, attributes, drop = FALSE]
  unique(hit)
}

# --- local RefMet / KEGG annotation snapshot --------------------------------
# The metabolomics stem standardizes metabolite names against RefMet and fills
# residual KEGG IDs from KEGG. Both were live queries, so refmet_name, refmet_id
# and kegg_id reflected whatever those databases held on the day the pipeline ran.
# All access now goes through a snapshot built once by
# sources/refmet/build_refmet_cache.R (see that folder's README). Sourcing
# refmet_lib.R also brings in refmet_fix_names()/refmet_lookup_key(), so the stem
# and the snapshot builder derive lookup keys from one definition.
REFMET_DIR <- file.path(.QC_REPO, "scripts", "00_preflight", "data-raw", "sources", "refmet")
source(file.path(REFMET_DIR, "refmet_lib.R"))

# --- gsutil_path_parsing.R internals ----------------------------------------
.match_ome_tissue_code <- function(desired_ome, input_tissue){
  specific_output <- OME_TISSUE_CODE %>%
    dplyr::filter(ome == desired_ome, tissue == input_tissue)
  if (nrow(specific_output) == 0) stop("No ome matches the desired ome/tissue combo")
  specific_output[["tissue_code"]]
}

# `version = NULL` (the default) resolves the version from config/file_versions.json
# for the stem this call composes, so no caller has to name one. It used to default to
# the literal "1.2" with callers passing "1.4" where that was wrong, which spread the
# per-file versions across ~20 hardcoded call sites. Pass an explicit string only to
# deliberately override the map for one call; to change the version a file ships at,
# edit the map instead.
#
# The four fields this pastes together are the four levels of that map's `files` tree,
# in order, under FREEZE_STEM_PREFIX — which is read from the map rather than written
# here so the composed name and the looked-up name cannot drift apart.
write_with_path_name <- function(actual_data_object = NULL, local_path = NULL, ome = NULL,
                                 tissue = NULL, data_category = NULL, data_details = NULL,
                                 version = NULL, return_name_only = FALSE){
  file_type <- ".txt"
  tissue_code <- .match_ome_tissue_code(desired_ome = ome, input_tissue = tissue)
  file_name <- paste(FREEZE_STEM_PREFIX, tissue_code, ome, data_category, data_details, sep = "_")
  if (is.null(version)) version <- freeze_version(file_name, sub("^\\.", "", file_type))
  file_name <- paste0(local_path, "/", file_name, "_v", version, file_type)
  if (return_name_only) return(file_name)
  utils::write.table(actual_data_object, file = file_name, row.names = FALSE, sep = "\t", quote = FALSE)
}

.clinical_ome_for_feature <- function(feature_id){
  prot_features  <- c("Insulin", "Glucagon", "CK")
  metab_features <- c("Glucose", "Lactate", "Glycerol", "NEFA", "KET", "Cortisol")
  unknown <- setdiff(feature_id, c(prot_features, metab_features))
  if (length(unknown) > 0) stop("Unrecognized clinical analyte(s): ", paste(unknown, collapse = ", "))
  dplyr::case_when(feature_id %in% prot_features  ~ "prot-clinical",
                   feature_id %in% metab_features ~ "metab-t-clinical")
}

.freeze_subdir_for_ome <- function(ome){
  dplyr::case_when(
    grepl("^transcript", ome) ~ "transcriptomics",
    grepl("^prot", ome)       ~ "proteomics",
    grepl("^epigen", ome)     ~ "epigenomics",
    grepl("^metab-u", ome)    ~ "metabolomics-untargeted",
    grepl("^metab-t", ome)    ~ "metabolomics-targeted"
  )
}

# (.find_path_name removed — raw input paths are resolved via quantid_path() from
#  the preflight quant-id catalog instead of a live gsutil ls.)

# --- load_pheno (was MotrpacHumanPreSuspensionData::pheno) -------------------
load_pheno <- function(load_acute_only = TRUE){
  pheno_data <- pheno
  pheno_data$pheno_data <- pheno_data$data
  if (load_acute_only)
    pheno_data$pheno_data <- pheno_data$pheno_data %>% dplyr::filter(visitcode == "ADU_BAS")
  pheno_data
}

# --- GCT normalization helpers (shared by prot-ph / prot-pr) ----------------
# Vendored from .../generate_normalized_expression/generate_prot_pr_qc_norm.R.
.median_mad_norm <- function(x, mad = TRUE){
  if (mad){
    scale.factor <- mean(apply(x@mat, 2, mad, na.rm = TRUE))
    x@mat <- scale(x@mat, center = apply(x@mat, 2, median, na.rm = TRUE),
                   scale = apply(x@mat, 2, mad, na.rm = TRUE))
    x@mat <- x@mat * scale.factor
  } else {
    x@mat <- scale(x@mat, center = apply(x@mat, 2, median, na.rm = TRUE), scale = FALSE)
  }
  x
}

.remove_na <- function(x, pct){
  cmapR::subset_gct(x, rid = which(rowSums(is.na(x@mat)) <= pct * ncol(x@mat)))
}

# Correlation between the two measurements of each muscle replicate pair, across features.
# `pairs` carries one row per pair: `kept` is the column the mean is written into, `dropped`
# is the redundant column. Call this on the matrix as it stands BEFORE any pair is averaged —
# once the intra-site means are written, a column that also belongs to an inter-site pair
# would be correlated against an already-merged measurement.
#
# The correlations are reported, not enforced: they say how well two measurements of one
# sample agree, which is the assumption averaging rests on, but the point below which a pair
# is untrustworthy is a release-specific judgement left to whoever reads the log.
report_replicate_correlations <- function(mat, pairs, pair_type = c("intra-site", "inter-site"),
                                          ome, tissue, method = c("pearson", "spearman")){
  pair_type <- match.arg(pair_type)
  method <- match.arg(method)
  tag <- paste0(tissue, ";", ome, ";replicate-cor;", pair_type)

  if (nrow(pairs) == 0){
    message(tag, ";no pairs")
    return(invisible(pairs))
  }

  missing <- setdiff(c(pairs$kept, pairs$dropped), colnames(mat))
  if (length(missing) > 0){
    stop("replicate pair columns absent from the matrix: ", paste(missing, collapse = ", "))
  }

  pairs$n_features <- NA_integer_
  pairs$correlation <- NA_real_
  for (k in seq_len(nrow(pairs))){
    x <- as.numeric(mat[[pairs$kept[k]]])
    y <- as.numeric(mat[[pairs$dropped[k]]])
    usable <- is.finite(x) & is.finite(y)
    pairs$n_features[k] <- sum(usable)
    # cor() needs two complete observations and non-zero variance in both vectors.
    if (sum(usable) >= 2 && stats::sd(x[usable]) > 0 && stats::sd(y[usable]) > 0){
      pairs$correlation[k] <- stats::cor(x[usable], y[usable], method = method)
    }
    message(sprintf("%s;%s vs %s;n=%d;%s=%s", tag, pairs$kept[k], pairs$dropped[k],
                    pairs$n_features[k], method,
                    ifelse(is.na(pairs$correlation[k]), "NA", sprintf("%.4f", pairs$correlation[k]))))
  }

  scored <- pairs$correlation[!is.na(pairs$correlation)]
  if (length(scored) == 0){
    message(sprintf("%s;%d pair(s);%s=all NA", tag, nrow(pairs), method))
  } else {
    worst <- pairs[which.min(pairs$correlation), ]
    message(sprintf("%s;%d pair(s);%s min=%.4f median=%.4f max=%.4f;lowest=%s vs %s",
                    tag, nrow(pairs), method, min(scored), stats::median(scored), max(scored),
                    worst$kept, worst$dropped))
  }
  invisible(pairs)
}

# --- process_covariates (was ::COVARIATES_FILE) -----------------------------
process_covariates <- function(meta, selected_ome, tissue_input,
                               include_technical = TRUE, custom_covariates = NULL){
  covariates_return <- list()
  covariates_return[["original_meta"]] <- meta
  input_covariates <- if (!is.null(custom_covariates)) custom_covariates else COVARIATES_FILE
  covariates <- input_covariates %>%
    as.data.frame() %>%
    dplyr::filter(ome == selected_ome) %>%
    dplyr::filter(tissue == 'all' | tissue == tissue_input)

  num_cov    <- covariates %>% dplyr::filter(data_type == "numerical")
  factor_cov <- covariates %>% dplyr::filter(data_type == "factor")

  sel_meta <- meta %>%
    dplyr::select(all_of(covariates$covariate)) %>%
    dplyr::mutate(across(all_of(num_cov$covariate), ~ scale(.) %>% as.numeric())) %>%
    dplyr::mutate(across(all_of(factor_cov$covariate), ~ as.factor(.) %>% droplevels())) %>%
    dplyr::mutate(group_timepoint = droplevels(interaction(randomGroupCode, Timepoint))) %>%
    dplyr::mutate(visit_group_timepoint = droplevels(interaction(visitcode, randomGroupCode, Timepoint))) %>%
    dplyr::mutate(sex_group_timepoint = droplevels(interaction(Sex, randomGroupCode, Timepoint)))

  technical_covs <- covariates %>% dplyr::filter(tech_or_design == "Technical")
  full_formula <- names(sel_meta)[!names(sel_meta) %in%
    c("randomGroupCode", "Timepoint", "visitcode", "pid",
      "group_timepoint", "visit_group_timepoint", "sex_group_timepoint")]
  design_covs <- c(full_formula[!full_formula %in% as.character(technical_covs$covariate)], "group_timepoint")
  if (!include_technical) full_formula <- full_formula[!full_formula %in% technical_covs$covariate]

  sex_diff_covs <- full_formula[!full_formula %in% c("Sex", "codedsiteid")]
  sex_formula_string <- paste(sex_diff_covs, collapse = " + ")
  formula_string_sex_differences <- paste("~ 0 + sex_group_timepoint + ", sex_formula_string, "+ (1 | pid)")
  formula_string <- paste(full_formula, collapse = " + ")
  formula_string_full     <- paste("~ 0 + group_timepoint + ", formula_string, "+ (1 | pid)")
  formula_string_training <- paste("~ 0 + visit_group_timepoint + ", formula_string, "+ (visitcode | pid)")
  non_mixed_model         <- paste("~ 0 + group_timepoint + ", formula_string)

  covariates_return[["technical_cov"]] <- technical_covs
  covariates_return[["design_cov"]]    <- design_covs
  covariates_return[["full_formula"]]  <- formula_string_full
  covariates_return[["training_formula"]] <- formula_string_training
  covariates_return[["sex_differences_formula"]] <- formula_string_sex_differences
  covariates_return[["metadata"]] <- sel_meta
  covariates_return[["non_mixed_model"]] <- non_mixed_model
  covariates_return
}

# --- epigen peak annotation (shared by ATAC + methylcap stems) ---------------
# Vendored from .../generate_normalized_expression/generate_atac_qc_norm.R.
# Annotates genomic-range features (chr:start-end) to genes via a TxDb built
# from Ensembl (online) + ChIPseeker::annotatePeak. Both epigen stems call this;
# each supplies its own assay-specific .annotate_*_features wrapper on top.
pre_cawg_get_peak_annotations_hs <- function(counts_dt,
                                             species = "Homo Sapiens",
                                             release = 105,
                                             txdb = NULL) {
  if (!"feature_id" %in% colnames(counts_dt) & !data.table::is.data.table(counts_dt)) {
    genomic_peaks = data.table::data.table(
      feature_id = rownames(counts_dt),
      chrom = gsub(":.*", "", rownames(counts_dt)),
      start = as.numeric(gsub(".*:|-.*", "", rownames(counts_dt))),
      end = as.numeric(gsub(".*-", "", rownames(counts_dt)))
    )
  } else if (!"feature_id" %in% colnames(counts_dt) & data.table::is.data.table(counts_dt)) {
    counts = counts_dt
    counts[, feature_id := paste0(chrom, ':', start, '-', end)]
    genomic_peaks = counts[, .(chrom, start, end, feature_id)]
  } else if ("feature_id" %in% colnames(counts_dt) & data.table::is.data.table(counts_dt)) {
    counts = counts_dt
    genomic_peaks = counts[, .(chrom, start, end, feature_id)]
  } else {
    counts = data.table::as.data.table(counts_dt)
    genomic_peaks = counts[, .(chrom, start, end, feature_id)]
  }

  if (is.null(txdb)) {
    # load the local v105 TxDb cache instead of building from Ensembl online
    txdb = .ensembl_v105_txdb()
  }

  accepted_chrom = GenomeInfoDb::seqlevels(txdb)
  accepted_chrom = accepted_chrom[!grepl("\\.", accepted_chrom)]

  genomic_peaks = genomic_peaks[!grepl("\\.", chrom)]
  genomic_peaks[, chrom := gsub("^chr", "", as.character(chrom))]

  if (!all(unique(genomic_peaks[, chrom]) %in% accepted_chrom)) {
    stop(sprintf(
      "The following chromosomes are found in the input but not in the txdb object: %s",
      paste0(unique(!genomic_peaks[, chrom] %in% accepted_chrom), collapse = ', ')
    ))
  }

  peak = GenomicRanges::GRanges(
    seqnames = genomic_peaks[, chrom],
    ranges = IRanges::IRanges(as.numeric(genomic_peaks[, start]), as.numeric(genomic_peaks[, end]))
  )
  peakAnno = ChIPseeker::annotatePeak(peak,
                                      level = "gene",
                                      tssRegion = c(-2000, 1000),
                                      TxDb = txdb,
                                      overlap = "all")
  pa = data.table::as.data.table(peakAnno@anno)

  if (nrow(pa) == nrow(genomic_peaks)) {
    pa[, feature_id := genomic_peaks[, feature_id]]
  } else {
    cols = c('seqnames', 'start', 'end')
    pa[, (cols) := lapply(.SD, as.character), .SDcols = cols]
    cols = c('chrom', 'start', 'end')
    genomic_peaks[, (cols) := lapply(.SD, as.character), .SDcols = cols]
    pa = merge(pa, genomic_peaks, by.x = c('seqnames', 'start', 'end'), by.y = c('chrom', 'start', 'end'), all.y = TRUE)
  }

  pa[, short_annotation := annotation]
  pa[grepl('Exon', short_annotation), short_annotation := 'Exon']
  pa[grepl('Intron', short_annotation), short_annotation := 'Intron']

  pa[, c('geneChr', 'strand') := NULL]

  cols = c('start', 'end', 'geneStart', 'geneEnd', 'geneStrand')
  pa[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
  pa[, dist_upstream := ifelse(end - geneStart <= 0, end - geneStart, NA_real_)]
  pa[, dist_downstream := ifelse(start - geneEnd >= 0, start - geneEnd, NA_real_)]
  pa[end >= geneStart & start <= geneEnd, dist_downstream := 0]
  pa[end >= geneStart & start <= geneEnd, dist_upstream := 0]
  pa[, relationship_to_gene := ifelse(is.na(dist_downstream), dist_upstream, dist_downstream)]
  pa[, c('dist_upstream', 'dist_downstream') := NULL]

  pa[relationship_to_gene == 0 & grepl("Downstream|Intergenic", short_annotation), short_annotation := "Overlaps Gene"]
  pa[geneStrand == 1 & relationship_to_gene > 0 & relationship_to_gene < 5000, short_annotation := "Downstream (<5kb)"]
  pa[geneStrand == 2 & relationship_to_gene < 0 & relationship_to_gene > -5000, short_annotation := "Downstream (<5kb)"]
  pa[geneStrand == 1 & relationship_to_gene > -5000 & relationship_to_gene < 0 & grepl("Downstream|Intergenic", short_annotation), short_annotation := "Upstream (<5kb)"]
  pa[geneStrand == 2 & relationship_to_gene < 5000 & relationship_to_gene > 0 & grepl("Downstream|Intergenic", short_annotation), short_annotation := "Upstream (<5kb)"]
  pa[abs(relationship_to_gene) >= 5000, short_annotation := "Distal Intergenic"]

  data.table::setnames(pa,
                       c('short_annotation', 'annotation', 'seqnames', 'geneId'),
                       c('custom_annotation', 'chipseeker_annotation', 'chrom', 'ensembl_gene')
  )

  return(pa)
}

# --- MICE imputation (shared by the prot-pr / prot-ph imputed stems) ---------
# Vendored verbatim from .../data-raw/run_mice.R (authors: Natalie M Clark,
# Stephanie Vartany). Sample-wise multiple imputation via mice; returns the
# per-feature mean across `num_imps` imputed datasets. `mice` is only required
# when this is actually called (imputed stems), so it stays a soft dependency.
run_mice <- function(data,
                     na_max = 0.4,
                     num_imps = 15,
                     seed = 2023,
                     num_cores = 1){
  # typically MICE wants the matrix as samples x features (feature-wise imputation),
  # but that is far too slow; sample-wise imputation gives good results here.
  data_filt <- data[rowSums(is.na(data)) / dim(data)[2] <= na_max, ]
  # The imputed tier carries no feature_metadata, so this is the only place the cut is
  # visible; for prot-ph it is about half the features.
  message(sprintf("run_mice: %d of %d features pass na_max = %.2f (%d dropped)",
                  nrow(data_filt), nrow(data), na_max, nrow(data) - nrow(data_filt)))

  print("Imputing using MICE")
  if (num_cores > 1) {
    mice_out <- mice::futuremice(data_filt, m = num_imps, parallelseed = seed, n.core = num_cores, print = TRUE)
  } else {
    mice_out <- mice::mice(data_filt, m = num_imps, seed = seed, print = TRUE)
  }
  print("Imputation done")

  print("Aggregating imputations")
  imp_data <- mice::complete(mice_out, 'long')

  avg_imp_data <- imp_data %>%
    dplyr::group_by(.id) %>%
    dplyr::summarize(dplyr::across(.fns = mean)) %>%
    dplyr::select(-c('.id', '.imp'))
avg_imp_data <- as.data.frame(avg_imp_data)

  # Feature-count guard. futuremice() splits the imputation across workers; a worker that
  # drops a feature would otherwise surface as an opaque "invalid 'row.names' length" on the
  # next line, or not at all. Same failure mode the step-09 dream backend has — see
  # VARIANCEPARTITION_PARALLEL_CORES in config/pipeline.env.
  if (nrow(avg_imp_data) != nrow(data_filt))
    stop(sprintf("run_mice returned %d of %d features (num_cores = %d) — %d lost",
                 nrow(avg_imp_data), nrow(data_filt), num_cores,
                 nrow(data_filt) - nrow(avg_imp_data)))

  rownames(avg_imp_data) <- rownames(data_filt)
  print("Aggregation done")

  return(avg_imp_data)
}

# ── load_qc_local and friends (moved here from lib/da_helpers.R) ──────────────
# Step 08's tests and step 11 need this view of the *_QC objects as much as the step-09
# DA stems do, so it lives with the other QC-object helpers rather than in a DA-only file.
#
# METABOLOMICS_CVS (built in step 05) and HUMAN_FEATURE_TO_GENE (step 07) are read on
# demand rather than at source time: qc_helpers.R is sourced by preflight and by every
# step-06 stem, all of which run before those objects exist.
.metab_cvs <- function() {
  f <- file.path(.STAGE1_DATA, "METABOLOMICS_CVS.rda")
  if (!file.exists(f)) stop("missing ", f, " — build step 05 before remove_redundant_metab")
  e <- new.env(); load(f, envir = e); e$METABOLOMICS_CVS
}
.feature_to_gene <- function() {
  f <- file.path(.STAGE1_DATA, "HUMAN_FEATURE_TO_GENE.rda")
  if (!file.exists(f)) stop("missing ", f, " — build step 07 before remove_unnamed_metab")
  e <- new.env(); load(f, envir = e); e$HUMAN_FEATURE_TO_GENE
}

# --- load_qc_local: assemble the nested [[tissue]][[ome]] list from local .rda ----
# Mirrors the exported load_qc(), but the object source is the local *_QC.rda in
# data/ rather than the lazily-loaded installed package objects. Unlike upstream there is no
# `epigen` argument: step 08 builds the ATAC and methylcap *_QC objects alongside every other
# ome, so they are returned by the ordinary scan below and need no special case.
load_qc_local <- function(selected_tissues = "all",
                          selected_omes = "all",
                          load_acute_only = TRUE,
                          remove_unnamed_metab = FALSE,
                          remove_redundant_metab = TRUE,
                          load_clinical = FALSE) {
  qc_files  <- list.files(.STAGE1_DATA, pattern = "_QC\\.rda$")
  obj_names <- sub("\\.rda$", "", qc_files)

  tissues <- tolower(sub("_.*", "", obj_names))
  omes <- gsub("_", "-", tolower(sub("^[^_]+_(.*)_QC$", "\\1", obj_names)))
  omes[omes == "trnscrpt"] <- "transcript-rna-seq"

  if (identical(selected_tissues, "all")) selected_tissues <- c("adipose", "blood", "muscle")
  avail_omes <- sort(unique(omes))
  if ("metab" %in% selected_omes)
    selected_omes <- unique(c(setdiff(selected_omes, "metab"), grep("^metab-", avail_omes, value = TRUE)))
  if (identical(selected_omes, "all")) selected_omes <- avail_omes
  if (!load_clinical) selected_omes <- setdiff(selected_omes, c("prot-clinical", "metab-t-clinical"))

  keep <- tissues %in% selected_tissues & omes %in% selected_omes
  out <- list()
  for (i in which(keep)) {
    e <- new.env(); load(file.path(.STAGE1_DATA, qc_files[i]), envir = e)
    out[[paste0(tissues[i], ".", omes[i])]] <- get(obj_names[i], envir = e)
  }

  # nest by tissue
  tvec <- sub("\\..*$", "", names(out))
  names(out) <- sub(".*\\.", "", names(out))
  out <- split(out, tvec)

  # recombine the split blood transcriptomics halves
  if ("blood" %in% selected_tissues && "transcript-rna-seq" %in% selected_omes) {
    e <- new.env()
    load(file.path(.STAGE1_DATA, "blood_transcript_1.rda"), envir = e)
    load(file.path(.STAGE1_DATA, "blood_transcript_2.rda"), envir = e)
    bc <- e$blood_transcript_1
    bc$qc_norm <- cbind(bc$qc_norm, e$blood_transcript_2$qc_norm)
    out[["blood"]][["transcript-rna-seq"]] <- bc
  }

  if (load_acute_only)      out <- subset_qc(out, desired_visit = "ADU_BAS")
  if (remove_unnamed_metab) out <- .filter_unnamed_metab(out)
  if (remove_redundant_metab) out <- .eliminate_redundant_metab(out)
  return(out)
}

# --- subset_qc (adapted from the package — divergences noted below) ----------
# Upstream: MotrpacHumanPreSuspensionData/R/load_qc.R:368, where the allowed values were
# literals buried in the body (`if (desired_group == "all") desired_group <- c(...)`).
# They are now the formal defaults, so `args(subset_qc)` documents them, and
# .match_subset_arg() below validates against them. Adapted here:
#   * every filter argument is validated. A typo ("ADUendur") used to match no rows,
#     which upstream treated as "skip this one" and returned that tissue x ome FULLY
#     UNFILTERED — the caller got the whole dataset back as though the filter had
#     applied. It now errors instead.
#   * a vector argument works. `desired_group == "all"` is a length-n condition whenever
#     the caller passes more than one value, which is a hard error on R >= 4.2, so the
#     documented subset_qc(qc, desired_group = c("ADUEndur", "ADUResist")) died.
#   * desired_visit gains "ADU_PAS". Upstream's "all" expanded to the duplicated literal
#     c("ADU_BAS", "ADU_BAS"), silently dropping every post-training sample (3,247 of
#     pheno's 20,180 rows). Unreachable from load_qc_local(), which always names the
#     visit, so no built output moves.
#   * a legitimately empty match writes the empty result back rather than `next`-ing,
#     for the same reason as the typo case above.
#   * sample_metadata is intersected and reordered against the surviving qc_norm columns.
#     The qc-object contract only requires colnames(qc_norm) subset-of vialLabel
#     (test_helpers.R:168), so sample_metadata can carry rows with no matching column,
#     and generate_prot_pr_da.R:86 / generate_metab_da.R:108 hard-error on those via
#     select(metadata$vialLabel). Anchoring on colnames(qc_norm) — the same pattern as
#     all_group_stats.R:96 — leaves the matrix and its column order untouched.
#   * drop = FALSE on every subscript, so a tissue x ome filtered down to one sample or
#     one feature stays a data.frame (test_helpers.R:153 asserts that it is one).
#   * desired_features filters qc_imputed too, and qc_imputed now follows the same
#     samples AND features as qc_norm rather than the samples alone.
#   * upstream's two message() calls, dropped in the original vendoring, are back as
#     warning()s (the convention at generate_metab_da.R:100), so a silently empty subset
#     is visible in a non-interactive run.
# Note the level lists stay closed, as upstream wrote them: a value outside a list is
# rejected rather than passed through, and an NA in the column is dropped under "all"
# because `%in%` reports NA as FALSE.

# Strict match.arg for subset_qc's filters. Reads the allowed values from subset_qc's own
# formals, so the signature is the single source of truth; "all" expands to the full set.
# match.arg(several.ok = TRUE) is deliberately NOT used here: it keeps whichever elements
# matched and silently discards the rest, so c("ADUEndur", "ADUendur") would pass with the
# typo quietly dropped — exactly the failure this is meant to catch.
.match_subset_arg <- function(value, arg_name) {
  choices <- setdiff(eval(formals(subset_qc)[[arg_name]]), "all")
  if (is.factor(value)) value <- as.character(value)
  # identical(value, c("all", choices)) is the argument-not-supplied case: R passes the
  # formal default through untouched.
  if (identical(value, "all") || identical(value, c("all", choices))) return(choices)
  if (!is.character(value) || length(value) == 0L)
    stop(arg_name, ' must be "all" or a non-empty character vector of: ',
         paste(choices, collapse = ", "))
  bad <- setdiff(value, choices)
  if (length(bad))
    stop(arg_name, ": unknown value(s) ", paste0('"', bad, '"', collapse = ", "),
         '; allowed values are "all" or any of: ', paste(choices, collapse = ", "))
  unique(value)
}

subset_qc <- function(load_qc_output,
                      desired_group     = c("all", "ADUControl", "ADUEndur", "ADUResist"),
                      desired_sex       = c("all", "Male", "Female"),
                      desired_timepoint = c("all", "pre_exercise", "during_20_min",
                                            "during_40_min", "post_10_min",
                                            "post_15_30_45_min", "post_3.5_4_hr",
                                            "post_24_hr"),
                      desired_age_group = c("all", "10-30", "30-60", "60-80"),
                      desired_visit     = c("all", "ADU_BAS", "ADU_PAS"),
                      desired_features  = NULL) {
  desired_group     <- .match_subset_arg(desired_group,     "desired_group")
  desired_sex       <- .match_subset_arg(desired_sex,       "desired_sex")
  desired_timepoint <- .match_subset_arg(desired_timepoint, "desired_timepoint")
  desired_age_group <- .match_subset_arg(desired_age_group, "desired_age_group")
  desired_visit     <- .match_subset_arg(desired_visit,     "desired_visit")

  for (tissue in names(load_qc_output)) {
    for (ome in names(load_qc_output[[tissue]])) {
      metadata <- load_qc_output[[tissue]][[ome]][["sample_metadata"]] %>%
        dplyr::filter(randomGroupCode %in% desired_group) %>%
        dplyr::filter(Sex %in% desired_sex) %>%
        dplyr::filter(Timepoint %in% desired_timepoint) %>%
        dplyr::filter(Age_3_groups %in% desired_age_group) %>%
        dplyr::filter(visitcode %in% desired_visit)

      # Anchor on the qc_norm columns so the matrix keeps its incoming order, then pull
      # sample_metadata into that order and drop rows that have no column.
      counts <- load_qc_output[[tissue]][[ome]][["qc_norm"]]
      vials  <- colnames(counts)[colnames(counts) %in% as.character(metadata$vialLabel)]
      counts <- counts[, vials, drop = FALSE]
      metadata <- metadata[match(vials, as.character(metadata$vialLabel)), , drop = FALSE]
      if (length(vials) == 0)
        warning("subset_qc: no samples satisfy the requested characteristics for ",
                tissue, "/", ome, " — qc_norm is now empty")

      if (!is.null(desired_features)) {
        counts <- counts[rownames(counts) %in% desired_features, , drop = FALSE]
        if (nrow(counts) == 0)
          warning("subset_qc: no features match desired_features for ", tissue, "/", ome)
      }

      if ("qc_imputed" %in% names(load_qc_output[[tissue]][[ome]])) {
        ci <- load_qc_output[[tissue]][[ome]][["qc_imputed"]]
        ci <- ci[, colnames(ci)[colnames(ci) %in% vials], drop = FALSE]
        if (!is.null(desired_features))
          ci <- ci[rownames(ci) %in% desired_features, , drop = FALSE]
        load_qc_output[[tissue]][[ome]][["qc_imputed"]] <- ci
      }

      load_qc_output[[tissue]][[ome]][["sample_metadata"]] <- metadata
      load_qc_output[[tissue]][[ome]][["qc_norm"]] <- counts
    }
  }
  return(load_qc_output)
}

# --- metab cleanup (package refs -> local objects) --------------------------
.eliminate_redundant_metab <- function(norm_qc_list) {
  if (!any(grepl("metab", lapply(norm_qc_list, names)))) return(norm_qc_list)
  metab_cvs <- .metab_cvs()   # read once, not per tissue x platform
  for (single_tissue in names(norm_qc_list)) {
    platforms <- names(norm_qc_list[[single_tissue]])
    metab_platforms <- platforms[grepl("metab", platforms) & platforms != "metab-t-clinical"]
    for (mp in metab_platforms) {
      lowest_cv_check <- metab_cvs %>%
        dplyr::filter(tissue == single_tissue, assay == mp) %>%
        dplyr::select(feature_id, refmet_name, lowest_CV)
      norm_qc_list[[single_tissue]][[mp]][["qc_norm"]] <- norm_qc_list[[single_tissue]][[mp]][["qc_norm"]] %>%
        dplyr::mutate(feature_id = rownames(.)) %>%
        dplyr::left_join(lowest_cv_check, by = "feature_id") %>%
        # Drops the redundant platform's copy (lowest_CV == "no") and, via `!NA == "no"`
        # being NA, every feature METABOLOMICS_CVS does not cover. Both intended: the CV
        # table defines the analyzable set, so this count is well below feature_metadata.
        dplyr::filter(!lowest_CV == "no") %>%
        dplyr::mutate(feature_refmet = dplyr::case_when(!is.na(refmet_name) ~ refmet_name, TRUE ~ feature_id)) %>%
        tibble::column_to_rownames("feature_refmet") %>%
        dplyr::select(-c(feature_id, lowest_CV, refmet_name))
    }
  }
  return(norm_qc_list)
}

.filter_unnamed_metab <- function(norm_qc) {
  feature_to_gene <- NULL   # 1.9M rows: read at most once, and only if a metab-u- ome is present
  for (tissue in names(norm_qc)) {
    for (ome in names(norm_qc[[tissue]])) {
      if (grepl("metab-u-", ome)) {
        curr <- norm_qc[[tissue]][[ome]][["qc_norm"]]
        if (is.null(feature_to_gene)) feature_to_gene <- .feature_to_gene()
        # upstream calls this possible_features_metab: the metabolites that HAVE a name
        possible_features_metab <- feature_to_gene %>% dplyr::filter(assay == "metab")
        named <- unique(c(possible_features_metab$feature_id,
                          possible_features_metab$refmet_name))
        norm_qc[[tissue]][[ome]][["qc_norm"]] <- curr[rownames(curr) %in% named, ]
      }
    }
  }
  return(norm_qc)
}

