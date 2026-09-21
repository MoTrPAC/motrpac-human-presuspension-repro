#!/usr/bin/env Rscript
# VENDORED from MotrpacHumanPreSuspensionAnalysis/R/run_cmeans.R.
#
# run_cmeans() plus its internal helper .reorder_clusters(), copied so step 13 can run without
# the installed analysis package — the dependency this stage exists to remove. Same convention
# as 12_build_camera_results/run_cameraPR.R: the upstream original is kept commented directly
# above each changed line.
#
# CHANGES:
#   - a `DA_list` argument is added and passed through to .prepare_DA_results(). Upstream
#     hardcodes `DA_list = NULL` there, which makes the helper reach into the installed package
#     via load_differential_analysis(); FCM_clustering_results.R supplies the step-10 objects
#     instead. This mirrors what run_cameraPR() already accepts in step 12.
#   - check_package_installation() is dropped; FCM_clustering_results.R requires Mfuzz and
#     Biobase up front.
#   - the multi-`num_clusters` branch is removed. Upstream plots Mfuzz::Dmin() and then blocks
#     on readline() for the operator to type a cluster count, which a batch build cannot do and
#     which would make the object depend on what was typed. The sweep now lives in
#     fcm_diagnostics.R, which writes the Dmin curve to disk; the chosen k comes back in as a
#     scalar argument. Passing a range here is an error rather than a silent prompt.
#
# .prepare_DA_results() is NOT re-vendored here — step 12 already carries it (with the same
# no-installed-package change), and FCM_clustering_results.R sources that file. One copy, so
# the two steps cannot drift apart on how DA results become z-score matrices.
#
# Clustering itself is unchanged: same scaling, same Mfuzz::mestimate()/mfuzz() calls, same
# set.seed(0), same cluster reordering.

suppressWarnings(suppressMessages({
  library(dplyr); library(data.table)
  # Attached rather than called as Biobase::exprs() below, so the upstream line is untouched.
  library(Biobase)
  # Mfuzz declares e1071 and Biobase in Depends, not Imports, so Mfuzz::mfuzz() resolves
  # cmeans() off the search path. Calling it namespace-qualified — as upstream does, from a
  # package that attaches these itself — fails with "could not find function 'cmeans'" unless
  # they are attached here.
  library(e1071)
}))

run_cmeans <- function(DA_list = NULL,
                       selected_tissues = c("all", "adipose",
                                            "blood", "muscle"),
                       selected_omes = c("transcript-rna-seq",
                                         "prot-pr", "prot-ol",
                                         "prot-ph", "metab"),
                       num_clusters_adipose = 13L,
                       num_clusters_blood = 12L,
                       num_clusters_muscle = 12L,
                       modality = c("both", "Endur", "Resist"))
{
  on.exit(gc())
  set.seed(0)
  # check_package_installation(pkg = "Mfuzz", fun = "run_cmeans")
  # check_package_installation(pkg = "Biobase", fun = "run_cmeans")
  # (dropped: FCM_clustering_results.R fails up front if Mfuzz or Biobase is missing)

  modality <- match.arg(modality,
                        choices = c("both", "Endur", "Resist"))

  selected_omes <- match.arg(arg = tolower(selected_omes),
                             choices = c("transcript-rna-seq",
                                         "prot-pr", "prot-ol",
                                         "prot-ph", "metab"),
                             several.ok = TRUE)

  num_clusters <- list("adipose" = num_clusters_adipose,
                       "blood" = num_clusters_blood,
                       "muscle" = num_clusters_muscle)

  num_clusters <- lapply(num_clusters, function(nc) {
    if (!is.vector(nc, mode = "numeric")) {
      stop("One or more of `num_clusters_adipose`, `num_clusters_blood`, or ",
           "`num_clusters_muscle` are not integer vectors.")
    }

    sort(unique(as.integer(pmax(2L, nc, na.rm = TRUE))))
  })

  # Upstream accepts a range here and resolves it interactively per tissue (see the removed
  # branch below). A batch build has no console, so exactly one k per tissue is required.
  bad_k <- names(num_clusters)[lengths(num_clusters) > 1L]
  if (length(bad_k))
    stop("run_cmeans() needs exactly one cluster count per tissue; got a range for ",
         paste(bad_k, collapse = ", "), ". Run fcm_diagnostics.R to sweep a range and write ",
         "the Dmin curve, then set the chosen k (FCM_K_ADIPOSE / FCM_K_BLOOD / FCM_K_MUSCLE).")

  ## Prepare DA results ----
  # DA_list = NULL upstream, which sends .prepare_DA_results() to load_differential_analysis()
  # and the installed package; the step-10 objects are passed in instead.
  DA_list <- .prepare_DA_results(DA_list = DA_list,
                                 selected_omes = selected_omes,
                                 selected_tissues = selected_tissues,
                                 convert_features = FALSE,
                                 .contrast_type = "exercise_with_controls")

  # Remove adipose global and phosphoproteomics due to insufficient timepoints
  DA_list[c("adipose.prot-pr", "adipose.prot-ph")] <- NULL
  nm <- names(DA_list)
  DA_list <- lapply(names(DA_list), function(name_i) {
    zi <- DA_list[[name_i]]

    zi <- zi[, !grepl("during", colnames(zi))]
    ome_i <- sub("^[^.]+\\.", "", name_i)
    rownames(zi) <- paste(ome_i, rownames(zi))

    return(zi)
  })
  names(DA_list) <- nm

  tissue <- sub("\\..*", "", names(DA_list))

  # Nest by tissue and stack matrices
  DA_list <- split(do.call(list, DA_list), tissue)

  zmat_list <- lapply(DA_list, function(li) {
    zi <- do.call(what = rbind, args = li)
    zi <- zi[complete.cases(zi), ]

    return(zi)
  })

  # Create ExpressionSet objects for mfuzz. Each row is scaled by dividing by
  # the standard deviation calculated after including two columns of zeros to
  # represent the pre-exercise timepoints. If this is not done, features that
  # change similarly in all contrasts will move further away from 0.
  eset_list <- lapply(zmat_list, function(zi) {
    zero_mat <- matrix(data = 0, nrow = nrow(zi), ncol = 2L)
    colnames(zero_mat) <- paste0(c("Endur", "Resist"),
                                 ".pre_exercise")
    tmp <- cbind(zero_mat, zi)

    # Scale features
    sd <- apply(tmp, 1L, sd)

    zi <- sweep(zi, 1L, sd, FUN = "/")

    zi <- zi[order(sd, decreasing = TRUE), ]

    Biobase::ExpressionSet(assayData = zi)
  })

  cmeans_list <- vector(mode = "list", length = length(eset_list))
  names(cmeans_list) <- names(eset_list)

  for (tissue_i in names(eset_list)) {
    eset_i <- eset_list[[tissue_i]]

    # m is determined before subsetting, since it becomes too large with few
    # contrast columns.
    m_i <- Mfuzz::mestimate(eset_i)

    # Optionally subset to specific modality
    if (modality != "both") {
      eset_i <- eset_i[, grepl(modality, colnames(eset_i)), drop = FALSE]
    } else {
      eset_i <- eset_i[, c(grep("Endur", colnames(eset_i)),
                           grep("Resist", colnames(eset_i)))]
    }

    num_clusters_i <- num_clusters[[tissue_i]]

    # The upstream interactive branch, removed above in favour of fcm_diagnostics.R:
    #
    # if (length(num_clusters_i) > 1L) {
    #   message("Determining optimal number of clusters for ", tissue_i, "...")
    #   min_centroid_dist <- Mfuzz::Dmin(eset = eset_i, m = m_i,
    #                                    crange = num_clusters_i, repeats = 1L, visu = FALSE)
    #   plot(x = num_clusters_i, y = min_centroid_dist, xlab = "Number of clusters",
    #        ylab = "Min. centroid dist.", main = tissue_i)
    #   num_clusters_i <- ""
    #   while (num_clusters_i %in% c("", "0", "1")) {
    #     num_clusters_i <- readline(prompt = sprintf(
    #       "Enter the optimal number of clusters (>= 2) for %s: ", tissue_i))
    #     num_clusters_i <- gsub("[^[:digit:].]", "", num_clusters_i)
    #     num_clusters_i <- sub("\\..*", "", num_clusters_i)
    #     num_clusters_i <- sub("^[0]+", "", num_clusters_i)
    #   }
    #   num_clusters_i <- as.integer(num_clusters_i)
    # }

    message(sprintf("FCM: %s — %d features x %d contrasts, k = %d, m = %.4f",
                    tissue_i, nrow(eset_i), ncol(eset_i), num_clusters_i, m_i))

    # FCM clustering
    FCM_i <- Mfuzz::mfuzz(
      eset = eset_i,
      centers = num_clusters_i,
      m = m_i
    )

    # Reorder clusters so those with similar trajectories are consecutive
    FCM_i <- .reorder_clusters(FCM_i)

    # Include input matrix (both modalities) and value of weighting exponent, m
    FCM_i[["input"]] <- exprs(eset_list[[tissue_i]])
    FCM_i[["call"]][["m"]] <- m_i

    cmeans_list[[tissue_i]] <- FCM_i
  }

  return(cmeans_list)
}



## Internal functions ----------------------------------------------------------

.reorder_clusters <- function(fclust) {
  d <- as.dist(1 - cor(t(fclust[["centers"]])))
  hc <- hclust(d)
  neworder <- hc[["order"]]

  fclust[["centers"]] <- fclust[["centers"]][neworder, ]
  rownames(fclust[["centers"]]) <- seq_len(nrow(fclust[["centers"]]))

  fclust[["size"]] <- fclust[["size"]][neworder]

  fclust[["cluster"]][] <- match(fclust[["cluster"]], neworder)
  fclust[["membership"]] <- fclust[["membership"]][, neworder]
  colnames(fclust[["membership"]]) <- seq_len(ncol(fclust[["membership"]]))

  return(fclust)
}
