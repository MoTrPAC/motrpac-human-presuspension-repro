#!/usr/bin/env Rscript
# VENDORED from MotrpacHumanPreSuspensionAnalysis/R/run_cluster_ORA.R.
#
# run_cluster_ORA() plus its internal helpers .cluster_ORA() and .fast_list_intersect(),
# copied so step 13 can run without the installed analysis package. Same convention as
# 12_build_camera_results/run_cameraPR.R: the upstream original is kept commented directly
# above each changed line.
#
# CHANGES, all of one kind except the last — `MotrpacHumanPreSuspensionAnalysis::<OBJECT>`
# becomes the pipeline's local global of the same name, supplied by FCM_clustering_results.R:
#   MOLECULAR_SIGNATURES  (step 03)   SET_TO_ID  (step 04)
# plus:
#   - check_package_installation() is dropped; FCM_clustering_results.R requires TMSig up front.
#
# .prepare_cluster_mem() comes from run_cluster_cameraPR.R, and .create_index() /
# .prepare_sets() from step 12's run_cameraPR.R; both are sourced by the stem, so the two
# cluster tests see identical membership matrices and identical filtered sets. That shared
# input is what makes the ORA-vs-CAMERA comparison this object exists for meaningful.
#
# The hypergeometric test is unchanged.

suppressWarnings(suppressMessages({
  library(dplyr); library(tidyr); library(tibble); library(data.table)
}))

run_cluster_ORA <- function(FCM,
                            selected_omes = c("transcript-rna-seq",
                                              "prot-pr", "prot-ph",
                                              "metab"),
                            selected_tissues = "all",
                            # database = names(MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES),
                            database = names(MOLECULAR_SIGNATURES),
                            path_to_gmt = NULL,
                            min_size = 5L,
                            overlap_cutoff = 0.7,
                            min_prob = 0.3)
{
  on.exit(gc())

  # Allow users with older R versions to still use the package
  # check_package_installation(pkg = "TMSig", fun = "run_cluster_ORA")
  # (dropped: FCM_clustering_results.R fails up front if TMSig is missing)

  min_prob <- max(0, min(1, min_prob, na.rm = TRUE))

  # List of membership probability matrices by tissue and ome
  mem_list <- .prepare_cluster_mem(FCM = FCM,
                                   selected_omes = selected_omes,
                                   selected_tissues = selected_tissues)

  # Every feature that was clustered, whether or not it reached min_prob below — so the
  # universe here is wider than the clusters built from it. Deliberate: cameraPR scores the
  # same features continuously, and FCM_ORA exists to be read against FCM_CAMERA, which only
  # holds if both start from one universe. min_prob is what ORA needs and cameraPR does not,
  # a discrete set to test.
  background_list <- lapply(mem_list, rownames)

  # Prepare list of clusters. Features are filtered based on membership
  # probability first
  cluster_list <- lapply(mem_list, function(mem_i) {
    cluster_levels <- seq_len(ncol(mem_i))

    # A feature reaching min_prob in no cluster is assigned to none, and stays in the
    # background above. That is 12-57% of features depending on tissue x ome (muscle metab
    # 57%, blood transcript-rna-seq 12%, 35% overall), so cluster sizes here sum to well
    # under background_size in the FCM_ORA table -- the gap is this threshold, not a bug.
    keep <- apply(mem_i >= min_prob, 1, any)
    mem_i <- mem_i[keep, , drop = FALSE]
    cluster_id <- apply(mem_i, 1, which.max)

    # Convert to factor to keep empty clusters
    cluster_id <- factor(cluster_id, levels = cluster_levels)

    split(x = rownames(mem_i), f = cluster_id)
  })

  ## Prepare molecular signatures ----
  index <- .create_index(database = database,
                         path_to_gmt = path_to_gmt)

  index_list <- .prepare_sets(background_list = background_list,
                              index = index,
                              overlap_cutoff = overlap_cutoff,
                              min_size = min_size)

  ## Over-representation analysis ----
  res_list <- lapply(names(index_list), function(name_i) {
    .cluster_ORA(
      cluster_list = cluster_list[[name_i]],
      background = background_list[[name_i]],
      index = index_list[[name_i]],
      min_size = min_size
    )
  })
  names(res_list) <- names(index_list)

  ## Reformat results ----
  out <- res_list %>%
    dplyr::bind_rows(.id = "idcol") %>%
    dplyr::mutate(tissue = sub("\\..*", "", idcol),
           assay = sub(".*\\.", "", idcol),
           idcol = NULL,
           across(.cols = c(tissue, assay),
                  .fns = ~ factor(.x, levels = unique(.x)))) %>%
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
           set_size_in_cluster, cluster_size, background_size,
           p_value, adj_p_value) %>%
    # Remove columns with all NA values
    dplyr::select(where(function(x) !all(is.na(x))))

  return(out)
}



## Internal functions ----------------------------------------------------------

.cluster_ORA <- function(cluster_list,
                         background,
                         index,
                         min_size = 5L)
{
  background <- unique(background)
  background <- background[!is.na(background)]
  size_background <- length(background)

  index <- TMSig::filterSets(
    x = index,
    background = background,
    min_size = min_size,
    max_size = size_background
  )
  set_sizes <- lengths(index)

  dt_list <- lapply(cluster_list, function(cluster_i) {
    cluster_i_elements <- unlist(cluster_i)
    size_cluster_i <- length(cluster_i_elements)

    # Restrict index to cluster elements. Keep empty sets
    index_ci <- .fast_list_intersect(x = index, y = cluster_i_elements)
    set_sizes_ci <- lengths(index_ci)

    dt_i <- data.table(set = names(set_sizes_ci),
                       set_size_in_cluster = set_sizes_ci,
                       stringsAsFactors = FALSE)

    return(dt_i)
  })

  dt <- rbindlist(dt_list, idcol = "cluster")

  dt[, `:=`(set_size = set_sizes[set],
            cluster_size = lengths(cluster_list)[cluster],
            background_size = size_background)]

  dt[, p_value := phyper(q = set_size_in_cluster - 1L,   # successes in cluster
                         m = set_size,                   # total successes
                         n = background_size - set_size, # total failures
                         k = cluster_size,               # sample_size
                         lower.tail = FALSE)]

  # Adjust p-values across clusters and sets
  dt[, adj_p_value := p.adjust(p_value, method = "BH")]
  dt[, cluster := factor(cluster, levels = names(cluster_list))]

  setorderv(dt, cols = c("cluster", "p_value", "set_size", "set"),
            order = c(1, 1, -1, 1))
  setcolorder(dt, neworder = "p_value", before = "adj_p_value")
  setDF(dt)

  return(dt)
}


# Note: does not check that x is a valid named list of character vectors or that
# y is a character vector.
.fast_list_intersect <- function(x, y) {
  # Convert list to data.table for fast filtering
  dt <- data.table(sets = rep(names(x), lengths(x)),
                   elements = unlist(x),
                   stringsAsFactors = FALSE)

  # Convert to factor to preserve order and keep empty sets when splitting
  dt[, sets := factor(sets, levels = unique(names(x)))]

  setorderv(dt, cols = c("sets", "elements"), order = rep(1L, 2L))

  y <- unique(y)
  y <- y[!is.na(y)]
  dt <- subset(dt, subset = elements %in% y)
  dt <- unique(dt)

  out <- split(x = dt[["elements"]], f = dt[["sets"]])

  return(out)
}
