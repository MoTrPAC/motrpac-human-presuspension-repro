#!/usr/bin/env Rscript
# VENDORED from MotrpacHumanPreSuspensionAnalysis/R/run_cameraPR.R.
#
# run_cameraPR() plus its three internal helpers (.prepare_DA_results, .create_index,
# .prepare_sets, .prepare_DA_results_and_sets), copied so step 12 can run without the
# installed analysis package — the dependency this stage exists to remove. Same convention as
# 09_build_da/da_common.R: the upstream original is kept
# commented directly above each changed line.
#
# CHANGES, all of the same kind — `MotrpacHumanPreSuspensionAnalysis::<OBJECT>` becomes the
# pipeline's local global of the same name, supplied by CAMERA_RESULTS.R:
#   MOLECULAR_SIGNATURES   (step 03)   SET_TO_ID              (step 04)
#   HUMAN_FEATURE_TO_GENE  (step 07)   CONTRAST_CONVERTER     (step 10)
# plus two structural ones, marked in place:
#   - check_package_installation() is dropped; CAMERA_RESULTS.R requires TMSig up front.
#   - .prepare_DA_results() no longer falls back to load_differential_analysis(); a DA_list
#     is always supplied from the step-10 objects, and a NULL is an error rather than a
#     silent reach for an installed package.
# Modelling and reshaping are otherwise unchanged.

suppressWarnings(suppressMessages({
  library(dplyr); library(tidyr); library(tibble); library(data.table)
}))

run_cameraPR <- function(DA_list = NULL,
                         selected_omes = c("transcript-rna-seq",
                                           "prot-pr",
                                           "prot-ph",
                                           "prot-ol",
                                           "metab"),
                         selected_tissues = "all",
                         database = setdiff(names(MOLECULAR_SIGNATURES), "PTMSIGDB"),
                         path_to_gmt = NULL,
                         min_size = 5L,
                         overlap_cutoff = 0.7)
{
  on.exit(gc())

  # Allow users with older R versions to still use the package
  # check_package_installation(pkg = "TMSig", fun = "run_cameraPR")
  # (dropped: CAMERA_RESULTS.R fails up front if TMSig is missing)

  ## Prepare differential analysis results and molecular signatures ----
  ls <- .prepare_DA_results_and_sets(
    DA_list = DA_list,
    selected_omes = selected_omes,
    selected_tissues = selected_tissues,
    database = database,
    path_to_gmt = path_to_gmt,
    min_size = min_size,
    overlap_cutoff = overlap_cutoff
  )

  # Dump contents of ls into function namespace
  for (name_i in names(ls))
    assign(x = name_i, value = ls[[name_i]])

  ## CAMERA-PR ----
  res_list <- lapply(names(index_list), function(name_i) {
    TMSig::cameraPR.matrix(
      statistic = DA_list[[name_i]], # matrix of z-statistics
      index = index_list[[name_i]], # named list of signatures
      use.ranks = FALSE,
      inter.gene.cor = 0.01,
      sort = TRUE,
      alternative = "two.sided",
      adjust.globally = FALSE,
      min.size = min_size
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
                  .fns = ~ factor(.x, levels = sort(unique(.x))))) %>%
    # Rename columns
    dplyr::rename(contrast = Contrast, set = GeneSet, set_size = NGenes,
           direction = Direction, t = TwoSampleT, z.std = ZScore,
           p_value = PValue, adj_p_value = FDR) %>%
    # Include contrast_type and contrast_short columns
    # dplyr::left_join(MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER, by = "contrast") %>%
    dplyr::left_join(CONTRAST_CONVERTER, by = "contrast") %>%
    dplyr::mutate(
      contrast = factor(
        x = contrast,
        levels = intersect(
          # levels(MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER$contrast),
          levels(CONTRAST_CONVERTER$contrast),
          levels(contrast)
        )),
      contrast_short = factor(
        x = contrast_short,
        levels = intersect(
          # levels(MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER$contrast_short),
          levels(CONTRAST_CONVERTER$contrast_short),
          levels(contrast_short)
        ))
    ) %>%
    # Include collection, database, set_id, and set_short columns
    # dplyr::left_join(MotrpacHumanPreSuspensionAnalysis::SET_TO_ID, by = "set") %>%
    dplyr::left_join(SET_TO_ID, by = "set") %>%
    dplyr::mutate(set_size = as.integer(set_size),
           set_size_DB = lengths(index)[set],
           size_ratio = round(set_size / set_size_DB,
                              digits = 3L),
           df = as.integer(df),
           direction = factor(direction, levels = c("Up", "Down")),
           across(.cols = everything(),
                  .fns = ~ structure(.x, names = NULL))) %>%
    # Convert set columns to factors to reduce the object size
    dplyr::mutate(across(.cols = c(set_id, set, set_short),
                  .fns = ~ factor(.x, levels = sort(unique(.x))))) %>%
    droplevels.data.frame() %>%
    # Adjust p-values separately by tissue, ome, contrast, and collection
    dplyr::mutate(.by = c(tissue, assay, contrast, collection),
           adj_p_value = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(contrast_type, tissue, assay, contrast,
            collection, database, p_value) %>%
    # Reorder columns
    dplyr::select(tissue, assay, contrast_type, contrast, contrast_short,
           collection, database, set_id, set, set_short,
           set_size, set_size_DB, size_ratio, direction,
           t, df, z.std, p_value, adj_p_value)

  return(out)
}


## Helper functions ------------------------------------------------------------

.prepare_DA_results <- function(DA_list = NULL,
                                selected_omes = c("transcript-rna-seq",
                                                  "prot-pr",
                                                  "prot-ph",
                                                  "prot-ol",
                                                  "metab"),
                                selected_tissues = "all",
                                convert_features = TRUE,
                                .contrast_type = NULL)
{
  selected_omes <- match.arg(selected_omes,
                             choices = c("transcript-rna-seq", "prot-pr",
                                         "prot-ph", "metab", "prot-ol"),
                             several.ok = TRUE)

  selected_tissues <- match.arg(selected_tissues,
                                choices = c("all", "adipose",
                                            "blood", "muscle"),
                                several.ok = TRUE)

  if ("all" %in% selected_tissues) {
    selected_tissues <- c("adipose", "blood", "muscle")
  }

  # if (is.null(DA_list)) {
  #   DA_list <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(...)
  # }
  # The caller always supplies the step-10 objects; reaching for an installed package here
  # is the dependency this stage removes, so a missing DA_list is an error.
  if (is.null(DA_list))
    stop("DA_list must be supplied — step 12 reads the step-10 *_DA objects, not an ",
         "installed MotrpacHumanPreSuspensionAnalysis.")

  # Unnest list and collapse tissue and ome with a "."
  DA_list <- unlist(DA_list, recursive = FALSE)

  keep_tissues <- sub("\\..*", "", names(DA_list)) %in% selected_tissues
  keep_omes <- sub(".*\\.", "", names(DA_list)) %in% selected_omes

  DA_list <- DA_list[keep_tissues & keep_omes]
  DA_names <- names(DA_list)

  if (is.null(DA_names))
    stop("`DA_list` must be a named list of differential analysis results ",
         "tables with names of the form 'tissue.assay'.")

  # Used by run_cmeans
  if (!is.null(.contrast_type)) {
    DA_list <- lapply(DA_list, function(xi) {
      xi %>%
        dplyr::select(-dplyr::any_of(c("contrast_type", "contrast_short"))) %>%
        dplyr::mutate(
          contrast = factor(contrast, levels = levels(CONTRAST_CONVERTER$contrast))
        ) %>%
        # dplyr::left_join(MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER, by = "contrast") %>%
        dplyr::left_join(CONTRAST_CONVERTER, by = "contrast") %>%
        dplyr::filter(contrast_type %in% .contrast_type) %>%
        dplyr::arrange(contrast) %>%
        droplevels.data.frame()
    })
  }

  # Convert list of data.frames to a list of matrices with features (genes,
  # phosphosites, or metabolites/lipids) as rows and contrasts as columns.
  DA_list <- lapply(DA_names, function(name_i) {
    ome_i <- sub(".*\\.", "", name_i)

    x <- DA_list[[name_i]] %>%
      dplyr::mutate(across(.cols = where(is.factor),
                    .fns = as.character))

    if (convert_features) {
      if (ome_i == "prot-ph") {
        flanking <-
          as.data.frame(HUMAN_FEATURE_TO_GENE)[, c("feature_id", "flanking_sequence")] %>%
          dplyr::mutate(across(.cols = everything(),
                        .fns = as.character))

        # Use flanking sequence as ID
        x <- x %>%
          dplyr::left_join(flanking, by = "feature_id") %>%
          dplyr::mutate(flanking_sequence = strsplit(flanking_sequence,
                                              split = "\\|")) %>%
          # Convert to single-sequence data
          tidyr::unnest(cols = flanking_sequence) %>%
          dplyr::mutate(new_id = flanking_sequence)
      } else if (ome_i == "metab") {
        x <- x %>%
          dplyr::mutate(feature_id = as.character(feature_id),
                 new_id = feature_id) %>%
          dplyr::select(contrast, new_id, z.std)
      } else {
        # Used to convert proteins and transcripts to genes
        feature_to_symbol <-
          as.data.frame(HUMAN_FEATURE_TO_GENE) %>%
          dplyr::select(feature_id, gene_symbol) %>%
          dplyr::mutate(across(.cols = where(is.factor),
                        .fns = as.character)) %>%
          dplyr::distinct()

        ## "prot-pr", "transcript-rna-seq", and default behavior for
        ## user-supplied DA results that are not "prot-ph" or "metab".
        required_cols <- c("feature_id", "contrast", "z.std")

        if (any(!required_cols %in% colnames(x)))
          stop("All `DA_list` tables must contain the following columns: ",
               paste(dQuote(required_cols), collapse = ", "), ".")

        # Include gene_symbol column
        x <- dplyr::left_join(x, feature_to_symbol,
                       by = "feature_id") %>%
          dplyr::mutate(new_id = gene_symbol)

        if (all(is.na(x$gene_symbol)))
          stop("No features in the ", name_i,
               " DA results map to gene symbols.")
      }

      # For each gene (or phosphosite or metabolite/lipid) and contrast pairing,
      # select the most extreme z-score (positive or negative)
      x <- x %>%
        # If the new_id is missing, use the feature ID to avoid unnecessarily
        # collapsing or removing rows.
        dplyr::mutate(new_id = ifelse(!is.na(new_id),
                               new_id,
                               feature_id)) %>%
        dplyr::arrange(contrast, desc(abs(z.std)), z.std) %>%
        dplyr::filter(.by = contrast,
               !duplicated(new_id))
    } else {
      # Do not convert features
      x <- dplyr::mutate(x, new_id = feature_id)
    }

    # Convert to a matrix with contrasts as columns and features as rows. Values
    # of the matrix are z-scores.
    x_wide <- x %>%
      dplyr::select(contrast, new_id, z.std) %>%
      tidyr::pivot_wider(id_cols = new_id,
                         names_from = contrast,
                         values_from = z.std) %>%
      tibble::column_to_rownames("new_id") %>%
      as.matrix()

    return(x_wide)
  })

  names(DA_list) <- DA_names

  return(DA_list)
}


.create_index <- function(
    database = names(MOLECULAR_SIGNATURES),
    path_to_gmt = NULL
) {
  choices <- names(MOLECULAR_SIGNATURES)

  if (is.null(path_to_gmt)) {
    choices <- names(MOLECULAR_SIGNATURES)

    if (!is.character(database))
      stop("`database` must be a character vector of one or more databases ",
           "to test, selected from the following options: ",
           paste(dQuote(choices), collapse = ", "),
           ".")

    database <- match.arg(arg = toupper(database),
                          choices = choices,
                          several.ok = TRUE)

    index <- MOLECULAR_SIGNATURES[database]

  } else {
    if (!is.character(path_to_gmt))
      stop("`database` must be a character vector of one or more databases ",
           "to test, selected from the following options: ",
           paste(dQuote(choices), collapse = ", "),
           ".")

    index <- lapply(path_to_gmt, TMSig::readGMT)
  }

  names(index) <- NULL # index is a nested list
  index <- unlist(index, recursive = FALSE)

  return(index)
}


.prepare_sets <- function(background_list,
                          index,
                          overlap_cutoff = 0.7,
                          min_size = 5L)
{
  overlap_cutoff <- max(0, min(1, overlap_cutoff))

  # Empty lists to store results
  index_list <- similar_sets_list <- list()

  for (name_i in names(background_list)) {
    ome_i <- sub(".*\\.", "", name_i)

    # Background of genes, metabolites, or phosphosites from the DA results
    background_i <- background_list[[name_i]]

    if (all(grepl("^PTMSIGDB", names(index)))) {
      # This code is from TMSig::filterSets. It was repurposed to work with
      # PTMsigDB, where the sites in each set end with ";u" or ";d", but the
      # background vector of IDs do not.
      set_dt <- data.table(
        sets = rep(names(index), lengths(index)),
        elements = unlist(index),
        stringsAsFactors = FALSE
      )

      set_dt[, elements2 := sub(";.*$", "", elements)]
      set_dt <- subset(set_dt, subset = elements2 %in%  background_i)

      index_i <- split(x = set_dt[["elements"]], f = set_dt[["sets"]])
      set_sizes <- lengths(index_i)
      keep_sizes <- (set_sizes >= min_size) & (set_sizes < length(background_i))
      index_i <- index_i[keep_sizes]
    } else {
      # Restrict sets to background, filter by size
      index_i <- TMSig::filterSets(
        x = index,
        background = background_i,
        min_size = min_size,
        max_size = length(background_i) - 1L
      )
    }

    # Overlap filter is not used for phosphoproteomics or
    # metabolomics/lipidomics datasets
    if (!ome_i %in% c("prot-ph", "metab")) {
      # Require a minimum proportion of genes in each set to be in the
      # background
      overlap_prop <- lengths(index_i) / lengths(index)[names(index_i)]

      keep <- which(overlap_prop >= ifelse(ome_i == "prot-ol",
                                           0.1,
                                           overlap_cutoff))

      if (length(keep) == 0L) {
        warning("No sets pass `min_size` and `overlap_cutoff` filters for ",
                ome_i,
                ". Excluding this dataset from the results.")

        next
      }

      index_i <- index_i[keep]
    }

    index_list[[name_i]] <- index_i
  }

  return(index_list)
}


.prepare_DA_results_and_sets <- function(DA_list = NULL,
                                         selected_omes = c("transcript-rna-seq",
                                                           "prot-pr",
                                                           "prot-ph",
                                                           "prot-ol",
                                                           "metab"),
                                         selected_tissues = "all",
                                         database = setdiff(names(MOLECULAR_SIGNATURES), "PTMSIGDB"),
                                         path_to_gmt = NULL,
                                         min_size = 5L,
                                         overlap_cutoff = 0.7)
{
  on.exit(gc())

  index <- .create_index(database = database,
                         path_to_gmt = path_to_gmt)

  ## Prepare DA results and molecular signatures ----
  DA_list <- .prepare_DA_results(DA_list = DA_list,
                                 selected_omes = selected_omes,
                                 selected_tissues = selected_tissues,
                                 convert_features = TRUE)

  index_list <- .prepare_sets(background_list = lapply(DA_list, rownames),
                              index = index,
                              min_size = min_size,
                              overlap_cutoff = overlap_cutoff)

  # Collect results in a list
  out <- list("DA_list" = DA_list,
              "index" = index,
              "index_list" = index_list)

  return(out)
}
