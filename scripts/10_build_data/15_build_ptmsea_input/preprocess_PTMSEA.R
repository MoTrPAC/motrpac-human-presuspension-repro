#' @title Pre-processing for PTM-SEA
#'
#' VENDORED from MotrpacHumanPreSuspensionAnalysis/data-raw/preprocess_PTMSEA.R
#' (identical to MotrpacHumanPreSuspension/R/preprocess_PTMSEA.R apart from the
#' `MotrpacHumanPreSuspensionAnalysis::` prefixes, which this copy drops anyway).
#'
#' Every change from upstream is marked in place with the upstream original kept
#' commented directly above it. They are all the same kind as step 12's: an
#' installed-package object becomes the pipeline's local global of the same name.
#' PTMSEA_INPUT.R loads MUSCLE_PROT_PH_DA / ADIPOSE_PROT_PH_DA (step 10) and
#' MUSCLE_PROT_PH_QC / ADIPOSE_PROT_PH_QC (step 08) into the calling frame and the
#' functions below close over them.
#'
#' @description This code produces a GCT file of the selected muscle contrasts
#' which is used as input to the PTM-SEA docker.
#'
#' NOTE ON WHAT HAPPENS TO THE OUTPUT: the .gct this returns is an INPUT to
#' PTM-SEA, which is run as a Docker image maintained by the Broad Institute
#' (broadinstitute/ptmsea, the ssGSEA2.0 / PTM-SEA distribution). That run is NOT
#' part of this repo and this pipeline does not attempt it — Stage 1 stops at
#' producing the .gct. See this step's README for where the file goes next.
#'
#' @param selected_tissues character; tissue to use, either adipose or muscle.
#' Default is muscle. Only one tissue may be processed at a time.
#' @param selected_contrast_type character; one or more contrast_types to use. Default is
#' \code{"exercise_with_controls"}.
#' @param selected_contrast_category character; one or more contrast_categories to use.
#' Default is \code{"c("EE-CON","RE-CON")"}
#'
#' @returns A GCT object with the properly formatted DA results. This GCT object
#' can be saved locally for use with PTM-SEA.
#'
#' @author Natalie M Clark
#'
#' @importFrom dplyr %>% select rename full_join bind_rows mutate across
#'   left_join arrange where desc everything

preprocess_PTMSEA <- function(selected_tissues = "muscle",
                              selected_contrast_type = "exercise_with_controls",
                              selected_contrast_category = c("EE-CON","RE-CON"))
{
  # CHANGE 1 — dropped. check_package_installation() is an Analysis-package helper and
  # this stage does not depend on that package being installed. PTMSEA_INPUT.R fails up
  # front if cmapR is missing, the same way CAMERA_RESULTS.R does for TMSig.
  # check cmapR is installed
  # check_package_installation(pkg = "cmapR", fun = "GCT")

  #obtain the correct DA results and filter accordingly
  if(selected_tissues=="muscle"){
    # CHANGE 2 — the installed-package object becomes the local global built by step 10.
    # prot_da <- MotrpacHumanPreSuspensionAnalysis::MUSCLE_PROT_PH_DA %>%
    prot_da <- MUSCLE_PROT_PH_DA %>%
      dplyr::filter(contrast_type %in% selected_contrast_type &
                      contrast_category %in% selected_contrast_category)
  }else{
    # CHANGE 2 (cont.)
    # prot_da <- MotrpacHumanPreSuspensionAnalysis::ADIPOSE_PROT_PH_DA %>%
    prot_da <- ADIPOSE_PROT_PH_DA %>%
      dplyr::filter(contrast_type %in% selected_contrast_type &
                      contrast_category %in% selected_contrast_category)
  }
  #pivot wider
  prot_da_wide <- tidyr::pivot_wider(prot_da,
                              id_cols="feature_id",
                              names_from="contrast_short",
                              values_from=c("logFC","z.std","p_value","adj_p_value"))

  #add feature metadata
  # CHANGE 3 — the `id` -> `feature_id` rename is now conditional. Upstream carries a
  # `TODO: @Chris -> Remove this chunk once id/feature_id is changed`; the locally built
  # *_PROT_PH_QC objects (step 08) already key feature_metadata on `feature_id`, so the
  # rename is a no-op here and errors if applied unconditionally. Kept guarded rather than
  # deleted so this file still works against a package object that has not been changed yet.
  # if(selected_tissues=="muscle"){
  #   prot_meta <- MUSCLE_PROT_PH_QC$feature_metadata %>% dplyr::rename(feature_id=id)
  # }else{
  #   prot_meta <- ADIPOSE_PROT_PH_QC$feature_metadata %>% dplyr::rename(feature_id=id)
  # }
  if(selected_tissues=="muscle"){
    prot_meta <- MUSCLE_PROT_PH_QC$feature_metadata
  }else{
    prot_meta <- ADIPOSE_PROT_PH_QC$feature_metadata
  }
  if("id" %in% colnames(prot_meta) && !"feature_id" %in% colnames(prot_meta))
    prot_meta <- prot_meta %>% dplyr::rename(feature_id=id)
  prot_da_wide <- dplyr::right_join(prot_meta,prot_da_wide,by="feature_id")

  #create GCT object for PTMSEA
  #filter for fully localized sites
  # UNCHANGED. Note that this line did not work against the shipped Data-package object
  # either: step 06's .annotate_prot_ph() used to select() confident_site away before the
  # metadata_features freeze was written, so feature_metadata had no such column and the
  # NULL subscript silently dropped every row. Step 06 now keeps it and step 08 carries it
  # through, so prot_meta really does have it. PTMSEA_INPUT.R gates on that before calling
  # this function rather than letting an absent-or-NA column fail silently here.
  prot_da_wide <- prot_da_wide[prot_da_wide$confident_site,]

  #use z-scores for the values
  rdesc <- prot_da_wide %>% dplyr::select(!dplyr::starts_with('z.std'))
  mat <- prot_da_wide %>% dplyr::select(dplyr::starts_with('z.std'))
  rownames(mat) <- rdesc$feature_id

  #add tissue to the contrast names
  colnames(mat) <- paste(selected_tissues,colnames(mat),sep=".")

  #create and return the GCT object
  # CHANGE 4 — cmapR::GCT rather than a bare GCT(). Upstream relies on the Analysis
  # package's @importFrom; nothing attaches cmapR here.
  # gct <- GCT(mat=as.matrix(mat),
  gct <- cmapR::GCT(mat=as.matrix(mat),
             rdesc=as.data.frame(rdesc),
             rid=rownames(mat),
             cid=colnames(mat))
  return(gct)
}
