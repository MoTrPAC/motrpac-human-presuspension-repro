#!/usr/bin/env Rscript
# VENDORED from MotrpacHumanPreSuspensionAnalysis/R/run_SCION.R:638 (scion_run_cmeans),
# copied so step 16 runs without the installed analysis package.
#
# CHANGES from upstream:
#   - the DA list is nested tissue -> ome -> data.frame before it reaches
#     .prepare_DA_results(). Upstream hands over a FLAT list("tissue.ome" = <data.frame>),
#     and .prepare_DA_results() opens with unlist(DA_list, recursive = FALSE), which
#     flattens each frame into its columns. The tissue/ome name filter then matches
#     nothing, the function returns zero elements, and cmeans dies stacking NULL with
#     "non-numeric matrix extent".
#   - rownames are emitted in the matrix's own feature..ome..tissue form rather than
#     "<ome> <feature_id>", so the cluster lookup in run_SCION() is a direct name match.
#     Upstream reconstructs the sanitised names with a regex; run_SCION() no longer
#     sanitises them (see run_SCION.R in this folder).
#
# .prepare_DA_results() is NOT re-vendored here - step 12 carries it, and SCION_NETWORKS.R
# sources that file.

scion_run_cmeans = function(DA_Input_Cmeans,
                            num_clusters = 13){
  on.exit(gc())

  # Nested tissue -> ome -> data.frame, which is the shape .prepare_DA_results()
  # unlists one level; a flat list of data frames flattens into its columns instead.
  nested_DA_input = list()
  for (flat_name in names(DA_Input_Cmeans)) {
    tissue_key = sub("\\..*", "", flat_name)
    ome_key = sub("^[^.]+\\.", "", flat_name)
    nested_DA_input[[tissue_key]][[ome_key]] = DA_Input_Cmeans[[flat_name]]
  }
  DA_list <- .prepare_DA_results(DA_list = nested_DA_input,
                                 convert_features = FALSE,
                                 .contrast_type = "exercise_with_controls")
  if (length(DA_list) == 0)
    stop("no differential-analysis tables survived preparation for clustering")

  nm <- names(DA_list)
  #this part just formats and removes during
  DA_list <- lapply(names(DA_list), function(name_i) {
    zi <- DA_list[[name_i]]

    zi <- zi[, !grepl("during", colnames(zi))]
    zi <- zi[, !grepl("post_10_min", colnames(zi))]

    ome_i <- sub("^[^.]+\\.", "", name_i)
    # Rownames arrive as feature_id__tissue; emitted as feature..ome..tissue, the form
    # the regulator and target matrices carry.
    feature_i <- sub("__[^_]+$", "", rownames(zi))
    tissue_i <- sub("^.*__", "", rownames(zi))
    rownames(zi) <- paste(feature_i, ome_i, tissue_i, sep = "..")

    return(zi)
  })

  names(DA_list) <- nm
  #so here is one big difference from the "larger" cmeans function
  #we don't split by tissue and just rbind everything as desired
  zmat_stacked = do.call(rbind, DA_list)

  #same eSet concept, just 1 column of zeros because we're just using 1 modality
  zero_mat <- matrix(data = 0, nrow = nrow(zmat_stacked), ncol = 1L)
  colnames(zero_mat) <- "pre_exercise"
  tmp <- cbind(zero_mat, zmat_stacked)
  # Scale features
  sd <- apply(tmp, 1, sd)
  zmat_stacked <- sweep(zmat_stacked, 1, sd, FUN = "/")
  zmat_stacked <- zmat_stacked[order(sd, decreasing = TRUE), ]
  eSet = Biobase::ExpressionSet(assayData = zmat_stacked)

  #This is still the same stuff as the other c-means stuff, we just call
  #the main mfuzz function from the Mfuzz package.
  m_i <- Mfuzz::mestimate(eSet)

  set.seed(0)
  FCM_i <- Mfuzz::mfuzz(
    eset = eSet,
    centers = num_clusters,
    m = m_i
  )
  return(FCM_i)
}
