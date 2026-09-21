#!/usr/bin/env Rscript
# VENDORED from MotrpacHumanPreSuspensionAnalysis/R/run_SCION.R (.load_scion_matrixes,
# .find_shared_scion_matrixes) and MotrpacHumanPreSuspensionData/R/load_qc.R:474
# (combine_qc_matrixes), copied so step 16 runs without either installed package.
#
# CHANGES to .load_scion_matrixes():
#   - load_qc()/subset_qc() come from lib/qc_helpers.R, not from
#     asNamespace("MotrpacHumanPreSuspensionData").
#   - prot-ph and prot-pr read the multiply-imputed matrix from step 06's freeze output.
#     Step 08 skips the imputed freeze files, so the *_QC objects carry no qc_imputed.
#   - the TF pool is step 14's UTORONTO_TFs, which has a feature_id column. The released
#     .rda is the raw 2,765 x 28 download and has none, and subset_qc() skips a NULL
#     desired_features, so subset_TFs = TRUE is a silent no-op against it.
#   - subset_DE reads the step-10 *_DA objects instead of load_differential_analysis().
#   - the tissue/ome vocabularies are inlined constants rather than
#     tissue_available_list() / ome_available_list() / metab_only_list().
#
# CHANGES to .find_shared_scion_matrixes():
#   - pheno comes from load_pheno(). Upstream reads it with
#     get("pheno", envir = asNamespace("MotrpacHumanPreSuspensionData"), inherits = FALSE);
#     pheno is lazy-loaded data with no binding in that namespace, so the call fails and
#     every run_SCION() call dies with "object 'pheno' not found".
#
# CHANGES to combine_qc_matrixes():
#   - HUMAN_FEATURE_TO_GENE is the step-07 object, read through .feature_to_gene().

SCION_TISSUES <- c("adipose", "blood", "muscle")

SCION_METAB_OMES <- c(
  "metab-u-hilicpos", "metab-u-ionpneg", "metab-u-lrpneg",
  "metab-u-lrppos", "metab-u-rpneg", "metab-u-rppos",
  "metab-t-amines", "metab-t-conv", "metab-t-imm-crt",
  "metab-t-oxylipneg", "metab-t-tca", "metab-t-nuc", "metab-t-acoa",
  "metab-t-ka"
)

SCION_OMES <- c(
  "prot-ol", "prot-ph", "prot-pr", "prot-clinical", "transcript-rna-seq",
  "epigen-methylcap-seq", "epigen-atac-seq", SCION_METAB_OMES, "metab-t-clinical"
)

# The omes SCION models on the imputed matrix rather than qc_norm. Imputation is what makes
# a per-feature random forest fittable: qc_norm carries missing values for these two.
SCION_IMPUTED_OMES <- c("prot-ph", "prot-pr")

#' Path to step 06's multiply-imputed freeze matrix for one tissue x ome.
.scion_imputed_path <- function(tissue, ome) {
  freeze_dir <- file.path(.STAGING, "freeze", .freeze_subdir_for_ome(ome), "qc-norm")
  path <- write_with_path_name(local_path = freeze_dir, ome = ome, tissue = tissue,
                               data_category = "imputed", data_details = "log2-mn",
                               return_name_only = TRUE)
  if (!file.exists(path))
    stop("no imputed matrix for ", tissue, " ", ome, " at ", path,
         " - run step 06 (generate_prot_ph_imputed.R / generate_prot_pr_imputed.R)")
  path
}

#' Read an imputed freeze matrix as features x vialLabel.
.scion_read_imputed <- function(tissue, ome) {
  wide <- read.csv(.scion_imputed_path(tissue, ome), sep = "\t", check.names = FALSE)
  matrix_out <- as.matrix(wide[, setdiff(names(wide), "feature_id"), drop = FALSE])
  rownames(matrix_out) <- wide$feature_id
  matrix_out
}

#' The prot-ph feature ids whose gene is a curated transcription factor (step 14).
.scion_tf_features <- function() {
  path <- file.path(.STAGE1_DATA, "UTORONTO_TFs.rda")
  if (!file.exists(path)) stop("UTORONTO_TFs.rda not built - run step 14 first: ", path)
  loaded <- new.env(parent = emptyenv())
  load(path, envir = loaded)
  tf <- get(ls(loaded)[1], envir = loaded)
  if (!"feature_id" %in% names(tf))
    stop("UTORONTO_TFs has no feature_id column - it is the raw download, not the ",
         "regulator pool. Rebuild step 14.")
  unique(as.character(tf$feature_id))
}

#' One assembled *_DA object as a data.frame, keyed by tissue x ome.
#'
#' The step-10 objects are named <TISSUE>_<OME>_DA with the ome upcased and hyphens as
#' underscores; transcript-rna-seq is stored as TRNSCRPT, matching load_qc_local().
.scion_load_da <- function(tissue, ome) {
  ome_stem <- if (ome == "transcript-rna-seq") "TRNSCRPT" else toupper(gsub("-", "_", ome))
  path <- file.path(.STAGE1_DATA, paste0(toupper(tissue), "_", ome_stem, "_DA.rda"))
  if (!file.exists(path))
    stop("no assembled DA object for ", tissue, " ", ome, " at ", path, " - run step 10 first")
  loaded <- new.env(parent = emptyenv())
  load(path, envir = loaded)
  as.data.frame(get(ls(loaded)[1], envir = loaded))
}

#' Differentially abundant feature ids for one tissue x ome, from the step-10 objects.
#'
#' Replaces load_differential_analysis(single_matrix = TRUE) piped into the same filter.
.scion_de_features <- function(tissue, ome) {
  .scion_load_da(tissue, ome) %>%
    dplyr::filter(contrast_type == "exercise_with_controls", adj_p_value < 0.05) %>%
    dplyr::pull(feature_id) %>%
    unique()
}

# --- combine_qc_matrixes (MotrpacHumanPreSuspensionData/R/load_qc.R:474) ----------------
combine_qc_matrixes = function(qc_norm_object,
                               drop_na = TRUE,
                               scale = TRUE,
                               filter_named_metab = TRUE,
                               make_metab_rownames = FALSE){
  full_qc_norm_matrix = data.frame()
  for(tissue in names(qc_norm_object)){
    for(ome in names(qc_norm_object[[tissue]])){
      curr_qc_norm = qc_norm_object[[tissue]][[ome]][["qc_norm"]]
      if (grepl("metab-u-", ome) & filter_named_metab){ #filtering only named metab
        possible_features_metab = .feature_to_gene() %>%
          dplyr::filter(assay == "metab")
        named_features = unique(c(possible_features_metab$feature_id, possible_features_metab$refmet_name))
        curr_qc_norm = curr_qc_norm[rownames(curr_qc_norm) %in% named_features,]
      }
      shortened_ome = ifelse(make_metab_rownames & grepl("metab", ome), "metab", ome)
      rownames(curr_qc_norm) = interaction(rownames(curr_qc_norm), shortened_ome, tissue, sep = "..")
      ncol_before = ncol(curr_qc_norm) #check number of columns before

      curr_metadata = qc_norm_object[[tissue]][[ome]][["sample_metadata"]]
      interaction_vector = with(curr_metadata, paste(pid, Timepoint, sep = "..")) #so here I'm converting viallabels into interactions of timepoint/pid
      name_mapping = setNames(interaction_vector, curr_metadata$vialLabel)
      curr_qc_norm_cols = colnames(curr_qc_norm)
      new_names = name_mapping[curr_qc_norm_cols]
      colnames(curr_qc_norm) = new_names
      if(!ncol(curr_qc_norm) == ncol_before){message(paste(tissue, ome, "not matching"))}  #make sure all the columns have a matching after
      #find shared columns
      if(nrow(full_qc_norm_matrix) == 0){
        full_qc_norm_matrix = curr_qc_norm
      }else{
        shared_cols = intersect(colnames(full_qc_norm_matrix), colnames(curr_qc_norm))
        if(length(shared_cols) == 0 || any(is.na(shared_cols))){
          message(paste("no shared cols between previous omes and", ome))
          next
        }
        full_qc_norm_matrix = full_qc_norm_matrix %>% dplyr::select(dplyr::all_of(shared_cols))
        curr_qc_norm = curr_qc_norm %>% dplyr::select(dplyr::all_of(shared_cols))
        full_qc_norm_matrix = rbind(full_qc_norm_matrix, curr_qc_norm)
      }
    }
  }
  if(drop_na & anyNA(full_qc_norm_matrix)){
    full_qc_norm_matrix = na.omit(full_qc_norm_matrix)
  }
  if(scale){
    full_qc_norm_matrix = t(scale(t(full_qc_norm_matrix)))
  }
  return(as.data.frame(full_qc_norm_matrix))
}

# --- .load_scion_matrixes (run_SCION.R:510) ---------------------------------------------
.load_scion_matrixes = function(desired_matrixes,
                                subset_TFs = TRUE,
                                subset_DE = TRUE){
  final_combined_matrix = data.frame()
  for(desired_input in desired_matrixes){
    desired_tissue = sub("\\..*", "", desired_input)
    if(!desired_tissue %in% SCION_TISSUES)
      stop("The syntax for regulators or targets isn't right. Regulators and targets should be in tissue.ome format (e.g blood.transcript-rna-seq)")
    desired_ome = sub(".*?\\.", "", desired_input)
    if(desired_ome == "metabolomics" | desired_ome == "metab") desired_ome = SCION_METAB_OMES
    if(!all(desired_ome %in% SCION_OMES))
      stop("The syntax for regulators or targets isn't right. Regulators and targets should be in tissue.ome format (e.g blood.transcript-rna-seq)")

    desired_qc_norm_single = load_qc_local(selected_omes = desired_ome,
                                           selected_tissues = desired_tissue,
                                           remove_unnamed_metab = TRUE)
    #-------------------------------------------------------------------------------
    # here we handle parameters and use imputed matrixes for prot-(pr/ph)
    #
    # Step 08 skips the imputed freeze files, so the *_QC objects carry no qc_imputed. The
    # matrix is read from step 06's freeze output and substituted for qc_norm instead. That
    # matrix is built from the all-sample qc-norm files and carries columns the acute object
    # does not, so it is restricted to the object's own samples - which is what keeps
    # qc_norm and sample_metadata aligned for combine_qc_matrixes() below.
    for(imputed_ome in intersect(desired_ome, SCION_IMPUTED_OMES)){
      imputed_matrix = .scion_read_imputed(desired_tissue, imputed_ome)
      sample_metadata = desired_qc_norm_single[[desired_tissue]][[imputed_ome]][["sample_metadata"]]
      shared_samples = intersect(colnames(imputed_matrix), as.character(sample_metadata$vialLabel))
      if(length(shared_samples) == 0)
        stop("the imputed matrix for ", desired_tissue, " ", imputed_ome,
             " shares no samples with the loaded QC object")
      desired_qc_norm_single[[desired_tissue]][[imputed_ome]][["qc_norm"]] =
        as.data.frame(imputed_matrix[, shared_samples, drop = FALSE])
      desired_qc_norm_single[[desired_tissue]][[imputed_ome]][["sample_metadata"]] =
        sample_metadata[as.character(sample_metadata$vialLabel) %in% shared_samples, , drop = FALSE]
      message(paste("Imputed matrix for", desired_tissue, imputed_ome, "-",
                    nrow(imputed_matrix), "features x", length(shared_samples), "samples"))
    }
    #basically just copy over imputed into the name "qc_norm" because then I can just use the combine_qc_matrixes function as usual
    if(any(desired_ome == "prot-ph") & subset_TFs) desired_qc_norm_single = subset_qc(desired_qc_norm_single,
                                                                                      desired_features = .scion_tf_features())

    if(subset_DE){
      DE_features = unique(unlist(lapply(desired_ome, function(one_ome)
        .scion_de_features(desired_tissue, one_ome))))
      desired_qc_norm_single = subset_qc(desired_qc_norm_single, desired_features = DE_features)
    }
    #--------------------------------------------------------------------------------
    participant_format = combine_qc_matrixes(desired_qc_norm_single,
                                             make_metab_rownames = TRUE)
    #just make all the metab stuff just = "metab"

    if(nrow(final_combined_matrix) == 0){
      final_combined_matrix = participant_format
    }else{
      shared_cols = intersect(colnames(final_combined_matrix), colnames(participant_format))
      if(length(shared_cols) == 0 || any(is.na(shared_cols))) stop("desired targets or regulators have no shared participants")
      final_combined_matrix = final_combined_matrix %>% dplyr::select(dplyr::all_of(shared_cols))
      participant_format = participant_format %>% dplyr::select(dplyr::all_of(shared_cols))
      final_combined_matrix = rbind(final_combined_matrix, participant_format)
    }
  }
  unique_parts = sub("\\.\\..*", "", colnames(final_combined_matrix))
  message(paste("Matrix of", paste(desired_matrixes, collapse = ", "),
                "has", length(unique(unique_parts)), "participants,", length(unique_parts), "samples"))
  return(final_combined_matrix)
}

# --- .find_shared_scion_matrixes (run_SCION.R:592) --------------------------------------
.find_shared_scion_matrixes = function(randomGroupCode,
                                       scion_matrix_one,
                                       scion_matrix_two){
  # load_pheno(load_acute_only = TRUE) applies the visitcode == "ADU_BAS" filter itself.
  pheno = load_pheno(load_acute_only = TRUE)[["pheno_data"]] %>%
    dplyr::filter(randomGroupCode == !!randomGroupCode)

  final_output = list()
  shared_cols = intersect(colnames(scion_matrix_one), colnames(scion_matrix_two))
  if(length(shared_cols) == 0 || any(is.na(shared_cols))) stop("desired targets or regulators have no shared participants")

  scion_matrix_one = scion_matrix_one %>%
    dplyr::select(dplyr::all_of(shared_cols)) %>%
    dplyr::select(dplyr::starts_with(as.character(pheno$pid))) #subset to the randomGroupCode

  scion_matrix_two = scion_matrix_two %>%
    dplyr::select(dplyr::all_of(shared_cols)) %>%
    dplyr::select(dplyr::starts_with(as.character(pheno$pid)))  #subset to the randomGroupCode

  if(ncol(scion_matrix_one) == 0 || ncol(scion_matrix_two) == 0)
    stop("desired targets or regulators have no shared participants in this randomGroupCode")

  message(paste("The new matrixes have", ncol(scion_matrix_one),  "shared samples"))
  final_output[["scion_matrix_one"]] = scion_matrix_one
  final_output[["scion_matrix_two"]] = scion_matrix_two
  return(final_output)
}
