#!/usr/bin/env Rscript
# Step 09 shared internals — the whole helper surface for the DA stems in one file.
# Every generate_<ome>_da.R stem sources this and nothing else from lib/, except the two
# stems that fit nothing: generate_methylcap_da.R never sources it, and generate_atac_da.R
# sources it only when RERUN_ATAC=TRUE.
#
# Three things live here, in this order:
#   1. the vendored upstream engine — run_dream(), .convert_dream_output(),
#      process_covariates(), .generate_contrasts_acute(), .generate_sex_contrasts()
#   2. the model runner/writer adapted for this pipeline — run_da_models()
#   3. the parallel switch — da_parallel_enabled()
#
# The engine defines process_covariates(), and so does lib/qc_helpers.R, which is sourced
# below it. The qc_helpers definition therefore wins, as it did when the engine was a
# separate file sourced at this same point. That is deliberate and inert: the one
# pipeline-relevant difference — reading the installed package COVARIATES_FILE rather than
# the local object — was changed in place, so the two are equivalent. Verified by running
# both orders and comparing the full output. Keep qc_helpers.R sourced AFTER the engine.

suppressWarnings(suppressMessages({ library(dplyr); library(BiocParallel) }))

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
  Sys.setenv(PRECOVID_ROOT = ROOT)
}
.DA_STEP <- file.path(ROOT, "scripts", "10_build_data", "09_build_da")


# ---------------------------------------------------------------------------
# 1. Vendored upstream engine
# ---------------------------------------------------------------------------

#' Fit linear mixed models using the dream framework
#'
#' This function fits linear mixed models for differential analysis using the
#' \code{dream} framework from the \code{variancePartition} package. It supports
#' both acute and training analyses internally, with contrast definitions and
#' model formulas determined by the supplied metadata and analysis type.
#'
#' For acute analyses, contrasts are generated using a restricted, baseline-anchored
#' contrast set to improve statistical stability and reduce computational burden.
#' Training contrasts are generated separately when requested.
#'
#' While the function is capable of fitting both acute and training models,
#' **only acute exercise results are publicly released** due to limited sample
#' sizes and power constraints in training cohorts.
#'
#' @param expression_object
#' A model-ready expression object (e.g., \code{edgeR::DGEList} or matrix),
#' with columns corresponding to samples.
#'
#' @param model_type
#' Character scalar specifying whether the analysis corresponds to
#' \code{"acute"} or \code{"training"} exercise.
#'
#' @param process_metadata
#' A list returned by \code{process_covariates}, containing processed metadata,
#' model formulas, and design information.
#'
#' @param voom
#' Logical indicating whether voom-based precision weights should be applied
#' prior to model fitting. Typically \code{TRUE} for RNA-seq and ATAC-seq data.
#'
#' @param parallel
#' Logical indicating whether model fitting should be parallelized using
#' \code{BiocParallel}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Reorders metadata to match expression columns
#'   \item Generates contrasts appropriate to the analysis type
#'   \item Constructs contrast matrices using \code{makeContrastsDream}
#'   \item Optionally applies voom-based weighting
#'   \item Fits linear mixed models using \code{dream}
#'   \item Applies empirical Bayes moderation via \code{eBayes}, with the estimator
#'     chosen explicitly from \code{EBAYES_LEGACY} rather than left to limma's
#'     data-dependent default
#' }
#'
#' @return
#' A fitted \code{dream} model object with moderated statistics.
#'
#' @note
#' Training models are fitted for internal analyses only and are not included
#' in public data releases.
#'
#' @keywords internal
#' @author christopher jin

# PARALLEL PATH: this is where it actually happens — `parallel = TRUE` builds the backend
# below and hands it to voomWithDreamWeights() and dream() as BPPARAM.
#
# MulticoreParam, not SnowParam("SOCK"). The one-feature-per-worker loss documented on
# VARIANCEPARTITION_PARALLEL_CORES in config/pipeline.env is a property of SOCK, not of
# parallelism: a SOCK worker is a cold R session that inherits the manager's options(warn = 2)
# via exportglobals, so the first findbars() it evaluates turns lme4's once-per-session
# reformulas deprecation warning into an error and bptry drops that feature. MulticoreParam
# forks, so each worker inherits an already-armed warning cache and never fires it.
#
# Cores come from VARIANCEPARTITION_PARALLEL_CORES so one switch governs both this and the
# step-06 ATAC voom, rather than this site silently taking detectCores() - 2 regardless.
run_dream = function(expression_object = NULL,
                     model_type = NULL,
                     process_metadata,
                     voom = FALSE,
                     parallel = FALSE){
  # check_package_installation("variancePartition")
  if(parallel) {
    num_cores = suppressWarnings(as.integer(Sys.getenv("VARIANCEPARTITION_PARALLEL_CORES", "")))
    if (is.na(num_cores) || num_cores < 1) num_cores = parallel::detectCores() - 2
    param <- MulticoreParam(num_cores, progressbar = TRUE)
  }
  meta_matrix = process_metadata$metadata
  meta_matrix = meta_matrix[match(colnames(expression_object), rownames(meta_matrix)), ] #reorder so they're in the same order as colnames of raw_counts
  #in theory, ^ they should already be the same, because they're matched beforehand...but I'm keeping it just in case.
  if(model_type == "acute") contrast_expressions = .generate_contrasts_acute(meta_matrix) #generate contrasts in a consistent way
  if(model_type == "training") contrast_expressions = .generate_contrasts_training(meta_matrix)
  if(model_type == "sex_differences") contrast_expressions = .generate_sex_contrasts(meta_matrix)

  if(model_type == "acute")  formula = process_metadata$full_formula %>% as.formula()
  if(model_type == "training") formula = process_metadata$training_formula %>% as.formula()
  if(model_type == "sex_differences") formula = process_metadata$sex_differences_formula %>% as.formula()

  L = variancePartition::makeContrastsDream(formula,
                                            meta_matrix,
                                            contrasts = contrast_expressions)
  if (voom){
    if (parallel){
      suppressWarnings({expression_object <- variancePartition::voomWithDreamWeights(expression_object, formula, meta_matrix, BPPARAM = param)}) #only for atac, rna-seq
    }else{
      expression_object <- variancePartition::voomWithDreamWeights(expression_object, formula, meta_matrix) #only for atac, rna-seq
    }
  }

  if (parallel){
    suppressWarnings({fit = variancePartition::dream(expression_object, formula, meta_matrix, L = L, BPPARAM = param)})
  }else{
    fit = variancePartition::dream(expression_object, formula, meta_matrix, L = L)
  }
  fit = variancePartition::eBayes(fit, legacy = .ebayes_legacy())
  return(fit)
}

# Which empirical-Bayes estimator eBayes uses, from EBAYES_LEGACY in config/pipeline.env.
# Always passed explicitly, never left to limma's default: squeezeVar() resolves a NULL
# `legacy` as identical(min(df), max(df)), and dream gives every feature its own fractional
# residual df, so that is always FALSE and limma >= 3.62.0 silently switches estimator. The
# choice moves df.prior/s2.prior, hence every feature's moderated variance, t and p.
# See the EBAYES_LEGACY block in config/pipeline.env for the measured effect.
.ebayes_legacy <- function(){
  if (utils::packageVersion("limma") < "3.62.0")
    stop("limma ", utils::packageVersion("limma"), " has no `legacy` argument to eBayes(); ",
         "this pipeline requires >= 3.62.0 so the estimator is chosen explicitly")
  v <- toupper(Sys.getenv("EBAYES_LEGACY"))
  if (!nzchar(v)) {
    f <- file.path(Sys.getenv("PRECOVID_ROOT"), "config", "pipeline.env")
    ln <- if (file.exists(f)) grep("(^|export )EBAYES_LEGACY=", readLines(f, warn = FALSE), value = TRUE) else character(0)
    if (length(ln)) v <- toupper(sub(".*:-([A-Za-z]+)\\}.*", "\\1", ln[length(ln)]))
  }
  if (!v %in% c("TRUE", "FALSE"))
    stop("EBAYES_LEGACY must be TRUE or FALSE, got '", v, "'")
  as.logical(v)
}


#' Convert dream model fits into a standardized differential analysis table
#'
#' This internal helper function extracts contrast-specific results from a fitted
#' \code{dream} model and converts them into a standardized, long-format table
#' suitable for downstream analysis, visualization, and public release.
#'
#' Only coefficients corresponding to explicit contrasts (i.e., containing
#' subtraction operators) are retained. For each contrast, full summary statistics
#' are extracted and annotated with model, tissue, and omics metadata.
#'
#' @param fit
#' A fitted \code{dream} model object returned by \code{run_dream}.
#'
#' @param metadata
#' Original sample metadata used for model fitting.
#'
#' @param formula
#' The model formula used to fit the mixed model.
#'
#' @param tissue
#' Character scalar specifying the tissue analyzed.
#'
#' @param ome
#' Character scalar specifying the omics layer analyzed.
#'
#' @details
#' For each retained contrast, the function extracts:
#' \itemize{
#'   \item Effect sizes and standard errors
#'   \item Moderated test statistics
#'   \item Confidence intervals
#'   \item Raw and adjusted p-values
#'   \item Model degrees of freedom and log-likelihood
#' }
#'
#' Results are concatenated across contrasts and sorted by adjusted p-value.
#'
#' @return
#' A data frame containing standardized differential analysis results across
#' all contrasts for the given tissue–omics combination.
#'
#' @note
#' This function does not perform filtering beyond contrast selection; all
#' downstream significance thresholds are applied elsewhere.
#'
#' @keywords internal output
#' @author christopher jin

.convert_dream_output = function(fit,
                                 metadata = NULL,
                                 formula = NULL,
                                 tissue = NULL,
                                 ome = NULL){
  full_contrasts = colnames(coef(fit))
  comparison_subset = full_contrasts[grep("-", full_contrasts)]
  #subset to only the models with some type of actual contrast
  res_tissue = data.frame() #final output file
  for(contrast in comparison_subset){
    # confint = FALSE, and the interval is rebuilt below instead. variancePartition's
    # topTable builds CI.L/CI.R as `se * qt(alpha, df = eb$df.total[top])`, and for a dream
    # fit df.total is a features x contrasts matrix while `top` is a linear index — so
    # `df.total[top]` reads down column one whatever `coef` asks for, and every contrast
    # after the first is bounded by the FIRST contrast's df. Still present upstream as of
    # 1.41.5. Asking for those columns and discarding them would only cost time, so the
    # tables carry the corrected pair alone and CI.L/CI.R stay forbidden downstream.
    res = variancePartition::topTable(fit, coef = contrast, number = Inf, p.value = 1, confint = FALSE)
    # CI.L_calculated / CI.R_calculated — the 95% interval with the df bug corrected.
    # df.total is subset by BOTH axes, feature name and contrast name, so each contrast is
    # bounded by its own df instead of column one's. logFC/t is the standard error, which
    # is why the margin needs no separate se: t = logFC/se by construction. Rows follow
    # `res`, which the chain below never reorders, so the two stay aligned.
    df_contrast = fit$df.total[rownames(res), contrast]
    margin = (res$logFC / res$t) * qt(0.975, df_contrast)
    # Keyed by name, NOT positional. topTable defaults to sort.by = "p", so `res` comes back
    # permuted relative to the fit, while fit$rdf and fit$logLik are per-feature NAMED vectors
    # in the fit's original feature order — rdf varies per feature because variancePartition
    # derives an approximate residual df from each feature's variance components. Assigning
    # them positionally gave every row another feature's df and logLik, silently, since
    # feature_id below comes from rownames() and stayed correct.
    res_single = res %>%
      dplyr::mutate(CI.L_calculated = logFC - margin) %>%
      dplyr::mutate(CI.R_calculated = logFC + margin) %>%
      dplyr::mutate(degrees_of_freedom = fit$rdf[rownames(.)]) %>%
      dplyr::mutate(logLik = fit$logLik[rownames(.)]) %>%
      dplyr::mutate(feature_id = rownames(.)) %>%
      dplyr::mutate(full_model = contrast) %>%
      dplyr::mutate(assay = ome) %>%
      dplyr::mutate(contrast = contrast) %>%
      dplyr::mutate(full_model = formula) %>%
      dplyr::mutate(tissue = tissue) %>%
      dplyr::rename(p_value = P.Value) %>%
      dplyr::rename(adj_p_value = adj.P.Val) %>%
      dplyr::select(assay,
                    feature_id,
                    z.std,
                    logFC,
                    CI.L_calculated, CI.R_calculated,
                    degrees_of_freedom,
                    logLik,
                    t,
                    AveExpr,
                    p_value, adj_p_value,
                    contrast,
                    full_model)

    res_tissue = rbind(res_tissue, res_single)
  }
  res_tissue = res_tissue %>% dplyr::arrange(adj_p_value)
  return(res_tissue)
}

#' Construct covariate metadata and model formulas for mixed model analysis
#'
#' This function defines, scales, and encodes covariates for use in linear mixed
#' model differential analysis. Covariate inclusion is controlled by omics layer
#' and tissue, with optional exclusion of technical covariates to support
#' sensitivity analyses.
#'
#' The function returns a structured list containing processed metadata matrices,
#' model formulas for acute and training analyses, and documentation of included
#' technical and design covariates.
#'
#' @param meta
#' A data frame of sample metadata containing participant identifiers, group
#' assignments, timepoints, and candidate covariates.
#'
#' @param selected_ome
#' Character scalar specifying the omics layer being analyzed.
#'
#' @param tissue_input
#' Character scalar specifying the tissue being analyzed.
#'
#' @param include_technical
#' Logical indicating whether technical covariates should be included in the
#' mixed model. Default is \code{TRUE}.
#'
#' @param custom_covariates
#' Optional data frame specifying custom covariate definitions. If \code{NULL},
#' a package-defined covariate configuration file is used. See methods for hose these
#' covariate configurations were originally derived.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Selects covariates relevant to the chosen ome and tissue
#'   \item Scales numerical covariates
#'   \item Encodes categorical covariates as factors
#'   \item Constructs interaction terms for group and timepoint
#'   \item Assembles fixed and random effect model formulas
#' }
#'
#' Acute analyses use participant-level random intercepts, while training
#' analyses include visit-level random slopes.
#'
#' @return
#' A named list containing:
#' \itemize{
#'   \item Processed metadata matrix
#'   \item Acute and training model formulas
#'   \item Lists of technical and design covariates
#'   \item Non–mixed model formula for diagnostic use
#' }
#'
#' @note
#' Although both acute and training formulas are generated, only acute analyses
#' are intended for public release due to training sample size limitations.
#'
#' @keywords modeling covariates
#' @author christopher jin

process_covariates = function(meta,
                              selected_ome,
                              tissue_input,
                              include_technical = TRUE,
                              custom_covariates = NULL){
  covariates_return = list() #output list for the end part
  covariates_return[["original_meta"]] = meta

  if(!is.null(custom_covariates)){
    input_covariates = custom_covariates
  }else{
    # ---- Upstream, verbatim: ----
    # input_covariates = MotrpacHumanPreSuspensionAnalysis::COVARIATES_FILE
    #
    # Changed to the pipeline's own COVARIATES_FILE (the global that lib/qc_helpers.R loads
    # from scripts/00_preflight/data/). Reading the installed package object is exactly the
    # dependency this stage exists to remove, and it also made the result depend on whether
    # this file or lib/qc_helpers.R was sourced last -- a silent hazard, since the wrong order
    # fits the models against a different covariate table without erroring. With both
    # definitions reading the same object the two are equivalent and source order no longer
    # matters. The lookup is resolved at call time, so this file can still be sourced first.
    input_covariates = COVARIATES_FILE
  }
  covariates = input_covariates %>%
    as.data.frame() %>%
    dplyr::filter(ome == selected_ome) %>%
    dplyr::filter(tissue == 'all' | tissue == tissue_input)

  num_cov = covariates %>% dplyr::filter(data_type == "numerical") #numerical covariates
  factor_cov = covariates %>% dplyr::filter(data_type == "factor")

  sel_meta = meta %>%
    dplyr::select(all_of(covariates$covariate)) %>%
    dplyr::mutate(across(all_of(num_cov$covariate), ~ scale(.) %>% as.numeric())) %>%
    dplyr::mutate(across(all_of(factor_cov$covariate), ~ as.factor(.) %>% droplevels())) %>%
    dplyr::mutate(group_timepoint = droplevels(interaction(randomGroupCode, Timepoint))) %>%
    dplyr::mutate(visit_group_timepoint = droplevels(interaction(visitcode, randomGroupCode, Timepoint))) %>%
    dplyr::mutate(sex_group_timepoint = droplevels(interaction(Sex, randomGroupCode, Timepoint)))

  technical_covs = covariates %>% dplyr::filter(tech_or_design == "Technical")
  full_formula = names(sel_meta)[!names(sel_meta) %in% c("randomGroupCode", "Timepoint", "visitcode", "pid", "group_timepoint", "visit_group_timepoint", "sex_group_timepoint")] #remove these from the character vector
  #the purpose of the design covariates section is to make a model.matrix()
  design_covs = c(full_formula[!full_formula %in% as.character(technical_covs$covariate)], "group_timepoint")
  if(!include_technical){ #remove for any modeling where some covariates have been regressed out
    full_formula = full_formula[!full_formula %in% technical_covs$covariate]
  }
  sex_diff_covs = full_formula[!full_formula %in% c("Sex", "codedsiteid")]
  sex_formula_string = paste(sex_diff_covs, collapse = " + ")
  formula_string_sex_differences = paste("~ 0 + sex_group_timepoint + ", sex_formula_string, "+ (1 | pid)")

  formula_string = paste(full_formula, collapse = " + ")

  #we readd group_timepoint first because of the way some contrast matrixes drop values in case of interactions w other levels of factors in the contrast matrixes
  formula_string_full = paste("~ 0 + group_timepoint + ", formula_string, "+ (1 | pid)")
  formula_string_training = paste("~ 0 + visit_group_timepoint + ", formula_string, "+ (visitcode | pid)")
  non_mixed_model = paste("~ 0 + group_timepoint + ", formula_string)

  #so i remove group_timepoint above and then make sure that it comes first because the order of the string can sometimes
  #actually change the contrast matrix formed and which columns are dropped in terms of the contrast comparisons.
  covariates_return[["technical_cov"]] = technical_covs
  covariates_return[["design_cov"]] = design_covs

  covariates_return[["full_formula"]] = formula_string_full
  covariates_return[["training_formula"]] = formula_string_training
  covariates_return[["sex_differences_formula"]] = formula_string_sex_differences

  covariates_return[["metadata"]] = sel_meta #so this is with all the tech/num cov in the correct format
  covariates_return[["non_mixed_model"]] = non_mixed_model

  return(covariates_return)
}


#' Generate a restricted set of acute exercise contrasts for linear mixed models
#'
#' This internal helper function constructs a curated set of contrast expressions
#' for acute exercise analyses, designed to capture biologically interpretable
#' within-group timepoint effects and between-group differences while avoiding
#' the combinatorial explosion associated with fully pairwise contrast generation.
#'
#' Rather than enumerating all possible timepoint × group comparisons, which
#' substantially increases computational burden and multiple-testing penalties,
#' this function defines a semi-manual contrast scheme anchored to the
#' \code{pre_exercise} baseline. This approach balances interpretability,
#' statistical efficiency, and computational feasibility.
#'
#' Only acute exercise contrasts are generated. Training-related contrasts are
#' intentionally excluded, as only acute analyses are publicly released due to
#' sample size and power considerations.
#'
#' @param metadata
#' A data frame of sample metadata containing at least a \code{Timepoint} column.
#' Timepoints must include \code{pre_exercise} and one or more post-baseline
#' acute exercise timepoints.
#'
#' @details
#' For each post-baseline acute timepoint, the following contrasts are generated
#' when applicable:
#' \itemize{
#'   \item Endurance vs Control (change from pre-exercise)
#'   \item Endurance within-group change from pre-exercise
#'   \item Resistance vs Control (change from pre-exercise)
#'   \item Resistance within-group change from pre-exercise
#'   \item Endurance vs Resistance (change from pre-exercise)
#'   \item Control within-group change from pre-exercise
#' }
#'
#' For timepoints lacking resistance blood draws (e.g.,
#' \code{during_20_min}, \code{during_40_min}), resistance-related contrasts
#' are automatically omitted.
#'
#' In addition, baseline (pre-exercise) group contrasts are included to
#' characterize between-group differences prior to exercise onset.
#'
#' @return
#' A character vector of contrast expressions suitable for use with
#' \code{limma}- or \code{dream}-based modeling frameworks.
#'
#' @note
#' This restricted contrast set is intentionally conservative and is used
#' to support stable inference in settings with limited sample sizes.
#' Fully pairwise contrast generation is deliberately avoided.
#'
#' @keywords internal

.generate_contrasts_acute = function(metadata){
  pre_contrast_expressions <- c()
  timepoints = unique(metadata$Timepoint)
  for(tp in timepoints){
    # message(tp)
    if (!tp == 'pre_exercise'){
      # meta_tp = metadata %>% filter(Timepoint == tp)
      # Chris: note - should just make contrast_Endur_Cntrl = paste0(contrast_Endur, " - ", contrast_Cntrls) at some point. This could definitely be refactored...
      contrast_Endur_Cntrl = sprintf("group_timepointADUEndur.%s - group_timepointADUEndur.pre_exercise - group_timepointADUControl.%s + group_timepointADUControl.pre_exercise", tp, tp)
      contrast_Endur = sprintf("group_timepointADUEndur.%s - group_timepointADUEndur.pre_exercise", tp)
      contrast_Resist_Cntrl = sprintf("group_timepointADUResist.%s - group_timepointADUResist.pre_exercise - group_timepointADUControl.%s + group_timepointADUControl.pre_exercise", tp, tp)
      contrast_Resist = sprintf("group_timepointADUResist.%s - group_timepointADUResist.pre_exercise", tp)
      contrast_Endur_Resist = sprintf("group_timepointADUEndur.%s - group_timepointADUEndur.pre_exercise - group_timepointADUResist.%s + group_timepointADUResist.pre_exercise", tp, tp)
      contrast_Cntrls = sprintf("group_timepointADUControl.%s - group_timepointADUControl.pre_exercise", tp)
      if (tp == 'during_20_min' | tp == 'during_40_min') {contrast_Resist_Cntrl = NULL; contrast_Endur_Resist = NULL; contrast_Resist = NULL} #resistance group doesn't get blood draws here
      pre_contrast_expressions = c(pre_contrast_expressions,
                                   contrast_Endur_Cntrl,
                                   contrast_Endur,
                                   contrast_Resist_Cntrl,
                                   contrast_Resist,
                                   contrast_Endur_Resist,
                                   contrast_Cntrls)
    }
  }

  pre_ex_endur_res = "group_timepointADUEndur.pre_exercise - group_timepointADUResist.pre_exercise"
  pre_ex_res_cntrl = "group_timepointADUResist.pre_exercise - group_timepointADUControl.pre_exercise"
  pre_ex_endur_cntrl = "group_timepointADUEndur.pre_exercise - group_timepointADUControl.pre_exercise"
  pre_contrast_expressions = c(pre_contrast_expressions, pre_ex_endur_res, pre_ex_res_cntrl, pre_ex_endur_cntrl)
  return(pre_contrast_expressions)
}

.generate_sex_contrasts = function(metadata){
  contrast_expressions <- c()
  timepoints = unique(metadata$Timepoint)
  for(tp in timepoints){
    # message(tp)
    if (tp == 'pre_exercise') next
    meta_tp = metadata %>% filter(Timepoint == tp)
    for(group in c("ADUEndur", "ADUResist")){
      if ((tp == 'during_20_min' | tp == 'during_40_min') &&  group == "ADUResist") next
      #female tp - female baseline (relative to control)
      female_grp_vs_control = paste0("(sex_group_timepointFemale.", group, ".", tp, " - ", "sex_group_timepointFemale.", group, ".pre_exercise) - ",
                                     "(sex_group_timepointFemale.ADUControl.", tp, " - ", "sex_group_timepointFemale.ADUControl.pre_exercise)"
                                     )
      #male is the same except sub Female for Male
      male_grp_vs_control = gsub("Female", "Male", female_grp_vs_control)

      #female vs male. Positive = female change > male change.
      female_change_vs_male_change = paste0(female_grp_vs_control, " - ", male_grp_vs_control)

      contrast_expressions = c(contrast_expressions, female_grp_vs_control, male_grp_vs_control, female_change_vs_male_change)

    }
  }
  return(contrast_expressions)
}

# ---------------------------------------------------------------------------
# 2. Pipeline-local internals and the adapted runner/writer
# ---------------------------------------------------------------------------

source(file.path(.DA_STEP, "filter_paired_n.R"))
# pipeline-local internals
source(file.path(ROOT, "scripts", "10_build_data", "lib", "qc_helpers.R"))

# Only the acute analysis is publicly released, so the stems never expose model_type.
DA_MODEL_TYPE <- "acute"
# Freeze versions are no longer one constant here. This used to be DA_VERSION <- "1.4"
# on the reasoning that the reference bucket is internally inconsistent (prot-ol and
# transcriptomics DA sit at v1.2, prot-pr/ph at v1.3) so the pipeline should just write
# one current version throughout. Each DA file now carries its own intended version for
# the next staging-bucket upload, resolved from config/file_versions.json by
# write_with_path_name(); passing version = NULL below lets that lookup happen. To move a
# DA table to a new version, edit that map. See lib/file_versions.R.

# ---- Upstream build, verbatim from MotrpacHumanPreSuspensionAnalysis/data-raw/
# ---- generate_differential_analysis/generate_DA_inputs.R (.run_models) ----
#
# .run_models = function(repo_local_dir,
#                        model_type,
#                        expression_object,
#                        process_metadata,
#                        tissue,
#                        ome,
#                        voom = FALSE,
#                        parallel = FALSE){
#
#   da_path = file.path(repo_local_dir, "data", "tmp", "freeze_DA/") #path for output
#   dir.create(da_path, recursive = TRUE, showWarnings = FALSE)
#
#   fit = run_dream(expression_object = expression_object,
#                   model_type = model_type,
#                   process_metadata = process_metadata,
#                   voom = voom,
#                   parallel = parallel)
#
#   # saveRDS(fit, file = file.path(da_path, paste0(tissue, "_", ome, "_", model_type, "_da_fit.rds")))
#
#   if(model_type == "acute") relevant_formula = process_metadata[["full_formula"]]
#   if(model_type == "training") relevant_formula = process_metadata[["training_formula"]]
#   if(model_type == "sex_differences") relevant_formula = process_metadata[["sex_differences_formula"]]
#
#   write_output = .convert_dream_output(fit,
#                                        metadata = process_metadata$original_meta,
#                                        tissue = tissue,
#                                        formula = relevant_formula,
#                                        ome = ome)
#   write_with_path_name(write_output,
#                        local_path = da_path,
#                        ome = ome,
#                        tissue = tissue,
#                        data_category = 'da',
#                        data_details = paste0('dream-', model_type))
#
# }
#
# Adapted below. Two changes: output lands in the pipeline's staging/freeze tree under
# <ome-group>/da/ (mirroring the bucket layout) rather than data/tmp/freeze_DA/, and the
# formula branch collapses to the acute one. The fit object is not persisted, as upstream.

.da_freeze_dir <- function(ome){
  d <- file.path(.STAGING, "freeze", .freeze_subdir_for_ome(ome), "da")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  return(d)
}

# PARALLEL PATH: `parallel` is handed straight to run_dream(). See da_parallel_enabled()
# below — the SOCK backend silently drops NA-carrying features, so every stem passes FALSE.
run_da_models <- function(expression_object,
                          process_metadata,
                          tissue,
                          ome,
                          voom = FALSE,
                          parallel = FALSE,
                          version = NULL){
  da_path <- .da_freeze_dir(ome)

  fit <- run_dream(expression_object = expression_object,
                   model_type = DA_MODEL_TYPE,
                   process_metadata = process_metadata,
                   voom = voom,
                   parallel = parallel)

  write_output <- .convert_dream_output(fit,
                                        metadata = process_metadata$original_meta,
                                        tissue = tissue,
                                        formula = process_metadata[["full_formula"]],
                                        ome = ome)
  write_with_path_name(write_output,
                       local_path = da_path,
                       ome = ome,
                       tissue = tissue,
                       data_category = "da",
                       data_details = paste0("dream-", DA_MODEL_TYPE),
                       version = version)
  message(sprintf("%s / %s: %d rows across %d contrast(s) -> %s",
                  tissue, ome, nrow(write_output), dplyr::n_distinct(write_output$contrast), da_path))
  return(invisible(write_output))
}

# PARALLEL SWITCH FOR EVERY STEP-09 STEM. Thin alias over
# variancepartition_parallel_enabled() in lib/qc_helpers.R, which governs every
# variancePartition path — these dream fits and the step-06 ATAC voom alike. The reasoning
# lives on VARIANCEPARTITION_PARALLEL_CORES in config/pipeline.env.
da_parallel_enabled <- function(){
  return(variancepartition_parallel_enabled())
}
