#!/usr/bin/env Rscript
# Stage 1 imputed stem: prot-pr (protein-ratio proteomics) — multiple imputation
# Adapted from MotrpacHumanPreSuspensionAnalysis/data-raw/generate_normalized_expression/
# generate_prot_pr_imputed.R 
#
# Multiple imputation (mice) of the prot-pr QC-norm matrix. These imputed datasets
# are used ONLY for SCION analyses in the manuscript, never for the primary
# differential analysis. Downstream of the prot-pr qc-norm stem: rather than
# re-download via load_qc(), this reads this stage's OWN prot-pr qc-norm freeze
# file (all samples, not acute-filtered) from staging/freeze. mice imputation is
# via the vendored run_mice() in lib/qc_helpers.R (mice is a soft dependency, only
# needed here). write_with_path_name is vendored too.
#
# Output: BIC-named imputed matrices under staging/freeze/proteomics/qc-norm/
#         (data_category 'imputed', version 1.4).
suppressWarnings(suppressMessages({ library(dplyr); library(tibble) }))
ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
source(file.path(ROOT, "scripts", "10_build_data", "lib", "qc_helpers.R"))

generate_prot_pr_imputed = function(repo_local_dir, num_cores = 1, num_imps = 15){
  if (!requireNamespace("mice", quietly = TRUE))
    stop("package 'mice' is required for the imputed stems")
  desired_ome = 'prot-pr'; tissue_types = c('muscle', 'adipose')
  qc_norm_dir = file.path(repo_local_dir, "freeze", "proteomics", "qc-norm")
  output_path = qc_norm_dir

  for (tissue in tissue_types){
    message(paste("Generating imputed matrixes for", desired_ome, tissue))
    tissue_code = .match_ome_tissue_code(desired_ome, tissue)

    # read this stage's own prot-pr qc-norm (was load_qc(..., load_acute_only=FALSE))
    fp = list.files(qc_norm_dir,
                    pattern = paste0(tissue_code, "_", desired_ome, "_qc-norm_log2-mn"),
                    full.names = TRUE)
    if (length(fp) != 1)
      stop(sprintf("expected 1 prot-pr qc-norm for %s, got %d", tissue_code, length(fp)))
    ome_data = utils::read.csv(fp, sep = "\t", check.names = FALSE)
    rownames(ome_data) = ome_data[[1]]
    ome_data[[1]] = NULL

    # mice requires syntactically valid (non-numeric) column names
    colnames(ome_data) = make.names(colnames(ome_data))
    prot_output_imputed = run_mice(ome_data, num_imps = num_imps, num_cores = num_cores)
    # undo the sanitization (vialLabels are all-numeric -> X-prefixed by make.names)
    colnames(prot_output_imputed) = gsub("X", "", colnames(prot_output_imputed))

    prot_output_imputed = prot_output_imputed %>%
      dplyr::mutate(feature_id = rownames(prot_output_imputed)) %>%
      dplyr::select(feature_id, dplyr::everything())

    write_with_path_name(prot_output_imputed,
                         local_path = output_path,
                         ome = desired_ome,
                         tissue = tissue,
                         data_category = 'imputed',
                         data_details = 'log2-mn')
  }
}

out_base <- .STAGING
# parallel cores from config/pipeline.env (PARALLEL_CORES); >1 -> mice::futuremice
.cores <- suppressWarnings(as.integer(Sys.getenv("PARALLEL_CORES", "1")))
if (is.na(.cores) || .cores < 1) .cores <- 1L
generate_prot_pr_imputed(out_base, num_cores = .cores)
message("prot-pr imputed matrices written under ", file.path(out_base, "freeze", "proteomics", "qc-norm"))
