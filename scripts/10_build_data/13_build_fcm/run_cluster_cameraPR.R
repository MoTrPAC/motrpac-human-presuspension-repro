#!/usr/bin/env Rscript
# VENDORED from MotrpacHumanPreSuspensionAnalysis/R/run_cluster_cameraPR.R.
#
# run_cluster_cameraPR() plus its internal helpers .prepare_cluster_mem() and
# .split_membership_by_ome(), copied so step 13 can run without the installed analysis
# package. Same convention as 12_build_camera_results/run_cameraPR.R: the upstream original is
# kept commented directly above each changed line.
#
# CHANGES, nearly all of one kind — `MotrpacHumanPreSuspensionAnalysis::<OBJECT>` becomes the
# pipeline's local global of the same name, supplied by FCM_clustering_results.R:
#   MOLECULAR_SIGNATURES  (step 03)   SET_TO_ID  (step 04)   HUMAN_FEATURE_TO_GENE  (step 07)
# plus two structural ones, marked in place:
#   - check_package_installation() is dropped; FCM_clustering_results.R requires TMSig up front.
#   - HUMAN_FEATURE_TO_GENE is wrapped in as.data.frame() before the dplyr chain, because the
#     pipeline's step-07 object is a keyed data.table rather than the package's data.frame.
#     Step 12 makes the same wrap for the same reason.
#
# .create_index() and .prepare_sets() are NOT re-vendored here — step 12 already carries them,
# and FCM_clustering_results.R sources that file so both steps filter sets identically.
#
# The test itself is unchanged: TMSig::cameraPR.matrix() on the membership probability
# matrices, use.ranks = TRUE, upper-tail.

suppressWarnings(suppressMessages({
  library(dplyr); library(tidyr); library(tibble); library(data.table)
}))

run_cluster_cameraPR <- function(FCM,
                                 selected_omes = c("transcript-rna-seq",
                                                   "prot-pr", "prot-ph",
                                                   "metab"),
                                 selected_tissues = "all",
                                 # database = names(MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES),
                                 database = names(MOLECULAR_SIGNATURES),
                                 path_to_gmt = NULL,
                                 min_size = 5L,
                                 overlap_cutoff = 0.7)
{
  on.exit(gc())

  # Allow users with older R versions to still use the package
  # check_package_installation(pkg = "TMSig", fun = "run_cluster_cameraPR")
  # (dropped: FCM_clustering_results.R fails up front if TMSig is missing)

  # List of membership probability matrices by tissue and ome
  mem_list <- .prepare_cluster_mem(FCM = FCM,
                                   selected_omes = selected_omes,
                                   selected_tissues = selected_tissues)

  # Features that contributed to FCM results
  background_list <- lapply(mem_list, rownames)

  ## Prepare molecular signatures ----
  index <- .create_index(database = database,
                         path_to_gmt = path_to_gmt)

  index_list <- .prepare_sets(background_list = background_list,
                              index = index,
                              overlap_cutoff = overlap_cutoff,
                              min_size = min_size)

  res_list <- lapply(names(index_list), function(name_i) {
    TMSig::cameraPR.matrix(
      statistic = mem_list[[name_i]],
      index = index_list[[name_i]],
      use.ranks = TRUE, # modified signed rank test
      inter.gene.cor = 0.01,
      alternative = "greater", # upper-tail test
      adjust.globally = FALSE,
      min.size = min_size
    )
  })
  names(res_list) <- names(index_list)

  ## Reformat results ----
  out <- res_list %>%
    bind_rows(.id = "idcol") %>%
    dplyr::mutate(tissue = sub("\\..*", "", idcol),
           assay = sub(".*\\.", "", idcol),
           idcol = NULL,
           across(.cols = c(tissue, assay),
                  .fns = ~ factor(.x, levels = unique(.x)))) %>%
    dplyr::select(-Direction) %>% # all "Up"
    # Rename columns
    dplyr::rename(cluster = Contrast, set = GeneSet, set_size = NGenes,
           p_value = PValue, adj_p_value = FDR) %>%
    # Include collection, database, set_id, and set_short columns
    # dplyr::left_join(MotrpacHumanPreSuspensionAnalysis::SET_TO_ID, by = "set") %>%
    dplyr::left_join(SET_TO_ID, by = "set") %>%
    droplevels.data.frame() %>%
    dplyr::mutate(set_size_DB = lengths(index)[set],
           size_ratio = round(set_size / set_size_DB, digits = 3L),
           across(.cols = everything(),
                  .fns = ~ structure(.x, names = NULL)),
           across(.cols = c(set_size, set_size_DB),
                  .fns = as.integer)) %>%
    # Convert set columns to factors to reduce the object size
    dplyr::mutate(across(.cols = c(set_id, set, set_short),
                  .fns = ~ factor(.x, levels = sort(unique(.x))))) %>%
    # Adjust p-values separately by tissue, ome, collection, and cluster.
    dplyr::mutate(.by = c(tissue, assay, collection, cluster),
           adj_p_value = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(tissue, assay, cluster, collection, database, p_value) %>%
    # Reorder columns
    dplyr::select(tissue, assay, cluster,
           collection, database, set_id, set, set_short,
           set_size, set_size_DB, size_ratio,
           p_value, adj_p_value) %>%
    # Remove columns with all NA values
    dplyr::select(where(function(x) !all(is.na(x))))

  return(out)
}



## Internal functions ----------------------------------------------------------

.prepare_cluster_mem <- function(FCM,
                                 selected_omes = c("transcript-rna-seq",
                                                   "prot-pr", "prot-ph",
                                                   "metab"),
                                 selected_tissues = "all")
{
  if (!is.vector(FCM, mode = "list") || is.null(names(FCM))) {
    stop("`FCM` must be the output of run_cmeans(). ",
         "See documentation for details.")
  }

  selected_omes <- match.arg(arg = tolower(selected_omes),
                             choices = c("transcript-rna-seq", "prot-pr",
                                         "prot-ph", "metab"),
                             several.ok = TRUE)

  selected_tissues <- match.arg(arg = tolower(selected_tissues),
                                choices = c("all", "adipose",
                                            "blood", "muscle"),
                                several.ok = TRUE)

  if (any(selected_tissues == "all")) {
    selected_tissues <- c("adipose", "blood", "muscle")
  }

  FCM <- FCM[names(FCM) %in% selected_tissues]

  if (!length(FCM)) {
    stop("names(FCM) do not match `selected_tissues`.")
  }

  mem_list <- lapply(FCM, function(fclust) {
    mem <- fclust[["membership"]]

    # Convert features to gene symbols, flanking sequences, or RefMet IDs and
    # split membership matrix rows by ome.
    mem_list <- .split_membership_by_ome(mem = mem,
                                         selected_omes = selected_omes)

    return(mem_list)
  })

  mem_list <- unlist(mem_list, recursive = FALSE)

  keep <- vapply(mem_list, function(mem_i) !is.null(mem_i), logical(1L))

  if (sum(keep) == 0L) {
    stop("Rownames of the membership probability matrices must ",
         "start with one of `selected_omes`.")
  }

  mem_list <- mem_list[keep]

  return(mem_list)
}


.split_membership_by_ome <- function(mem, selected_omes) {
  mem_df <- mem %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature_id") %>%
    dplyr::mutate(assay = sub("(^[^ ]+).*", "\\1", feature_id),
           feature_id = sub("[^ ]+ ", "", feature_id)) %>%
    dplyr::filter(assay %in% selected_omes)

  if (nrow(mem_df) == 0L)
    return(NULL)

  feature_conv <-
    # MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE %>%
    as.data.frame(HUMAN_FEATURE_TO_GENE) %>%
    dplyr::filter(!grepl("^epi", assay)) %>%
    dplyr::select(feature_id, gene_symbol, flanking_sequence) %>%
    dplyr::mutate(across(.cols = everything(),
                  .fns = as.character)) %>%
    dplyr::mutate(
      new_id = case_when(
        !is.na(flanking_sequence) ~ flanking_sequence,
        !is.na(gene_symbol) ~ gene_symbol,
        TRUE ~ feature_id # if N/A or if features are RefMet names
      )
    ) %>%
    dplyr::mutate(new_id = strsplit(new_id, split = "\\|")) %>%
    tidyr::unnest(cols = new_id) %>%
    dplyr::select(feature_id, new_id)

  mem_list <- mem_df %>%
    dplyr::left_join(feature_conv, by = "feature_id") %>%
    # Some features are not in HUMAN_FEATURE_TO_GENE
    dplyr::mutate(new_id = ifelse(is.na(new_id), feature_id, new_id)) %>%
    tidyr::pivot_longer(cols = all_of(colnames(mem)),
                        names_to = "cluster",
                        values_to = "membership") %>%
    dplyr::mutate(cluster = factor(cluster, levels = unique(cluster))) %>%
    # Select highest probability per new_id in each assay/cluster combination
    dplyr::arrange(assay, cluster, desc(membership), new_id) %>%
    dplyr::filter(.by = c(assay, cluster),
           !duplicated(new_id)) %>%
    tidyr::nest(.by = assay) %>%
    dplyr::mutate(data = lapply(data, function(data_i) {
      data_i %>%
        tidyr::pivot_wider(id_cols = new_id,
                           names_from = cluster,
                           values_from = membership) %>%
        tibble::column_to_rownames("new_id") %>%
        as.matrix()
    })) %>%
    tibble::deframe()

  return(mem_list)
}
