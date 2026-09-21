#!/usr/bin/env Rscript
# Tests for the SCION networks (step 16).
#
# What these are for is the naming contract between the two halves of run_SCION(). The edge
# table is built from the weight matrix's dimnames, and the hub step joins those names back
# against the input matrix rownames, so anything that rewrites one side and not the other
# severs the join silently: hubregdata comes back empty, RS.Get.Weight.Matrix() dies on a
# zero-column matrix, and the run ends BEFORE write.table() — leaving a directory of
# per-cluster CSVs and no network file. Nothing about the numbers looks wrong, because there
# are no numbers to look at.
#
# Observed networks only. The permutations run identical code, so 100 of them would restate
# the same result 100 times.
#
# The input rownames are read from the myregdata.csv / mytargetdata.csv that run_SCION()
# writes in the hub step, so the matrices are not rebuilt here. A run that produced no hub
# layer has neither file, and those checks SKIP.
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
suppressWarnings(suppressMessages({ library(dplyr) }))
cat("== scion network tests ==\n")

SCION_DIR <- file.path(.root, "staging", "scion")
MANIFEST <- file.path(SCION_DIR, "scion_manifest.tsv")

if (!file.exists(MANIFEST)) {
  skip("scion_manifest.tsv present", "no networks inferred — run step 16 without RUN_SCION=FALSE")
  finish()
}

manifest <- read.csv(MANIFEST, sep = "\t", check.names = FALSE)
assert_cols("scion_manifest.tsv", manifest,
            c("network", "randomGroupCode", "permutation", "kind", "regulators", "targets",
              "file", "n_edges", "n_regulators", "n_targets", "n_clusters"))
assert("manifest is non-empty", nrow(manifest) > 0)

observed <- manifest[manifest$kind == "observed", , drop = FALSE]
if (nrow(observed) == 0) {
  skip("observed networks present", "manifest holds permutations only")
  finish()
}
report("observed networks", paste(nrow(observed), "of", nrow(manifest), "runs"))

EXPECTED_COLS <- c("Regulator", "Interaction", "Target", "Weight", "Cluster")

#' Rownames of an input matrix run_SCION() dumped beside the network.
input_rownames <- function(dir_path, which_matrix) {
  path <- file.path(dir_path, paste0(which_matrix, ".csv"))
  if (!file.exists(path)) return(NULL)
  as.character(read.csv(path, check.names = FALSE)[[1]])
}

for (i in seq_len(nrow(observed))) {
  row <- observed[i, ]
  tag <- paste0(row$network, ":", row$randomGroupCode)
  path <- file.path(SCION_DIR, row$file)
  dir_path <- dirname(path)

  if (!file.exists(path)) {
    assert(paste0(tag, ": network file present"), FALSE, paste("missing:", row$file))
    next
  }

  edges <- read.csv(path, sep = "\t", check.names = FALSE, colClasses = "character")
  assert(paste0(tag, ": edge table columns"), identical(colnames(edges), EXPECTED_COLS),
         paste("got:", paste(colnames(edges), collapse = ", ")))
  if (!identical(colnames(edges), EXPECTED_COLS)) next

  assert(paste0(tag, ": edge count matches the manifest"), nrow(edges) == row$n_edges,
         paste("table has", nrow(edges), "rows, manifest says", row$n_edges))

  n_clusters <- dplyr::n_distinct(edges$Cluster[nzchar(edges$Cluster)])
  assert(paste0(tag, ": at least one cluster contributed edges"), n_clusters > 0,
         "no cluster-layer edges")
  report(paste0(tag, ": clusters"), paste(n_clusters, "of the 13 requested"))

  reg_names <- input_rownames(dir_path, "myregdata")
  tgt_names <- input_rownames(dir_path, "mytargetdata")
  if (is.null(reg_names) || is.null(tgt_names)) {
    skip(paste0(tag, ": input-name join"), "no myregdata.csv/mytargetdata.csv — run produced no hub layer")
    next
  }

  # Every name in the table has to be a name that went in. An id rewritten on the way out is
  # unjoinable to the DA tables and to the feature-to-gene map, which is what every
  # downstream use of this network does first.
  in_cluster <- nzchar(edges$Cluster)
  unknown_reg <- setdiff(unique(edges$Regulator), reg_names)
  assert(paste0(tag, ": regulator ids are input rownames"), length(unknown_reg) == 0,
         paste0(length(unknown_reg), " unmatched, e.g. ", paste(head(unknown_reg, 2), collapse = ", ")))

  # Targets are held to the target matrix only in the cluster layer. In the hub layer a
  # target may legitimately be a REGULATOR: when no hub resolves to a cross-ome target —
  # always the case for a metabolomics regulator set — run_SCION() connects the hubs to each
  # other and passes the regulator matrix in as the target matrix.
  unknown_tgt <- setdiff(unique(edges$Target[in_cluster]), tgt_names)
  unknown_hub_tgt <- setdiff(unique(edges$Target[!in_cluster]), c(tgt_names, reg_names))
  assert(paste0(tag, ": target ids are input rownames"),
         length(unknown_tgt) == 0 && length(unknown_hub_tgt) == 0,
         paste0(length(unknown_tgt), " cluster-layer and ", length(unknown_hub_tgt),
                " hub-layer unmatched, e.g. ",
                paste(head(c(unknown_tgt, unknown_hub_tgt), 2), collapse = ", ")))

  # The inputs are meant to show up in the output. Coverage well under half means most of
  # what was handed to the inference never reached an edge.
  reg_seen <- length(intersect(unique(edges$Regulator), reg_names))
  tgt_seen <- length(intersect(unique(edges$Target), tgt_names))
  if (reg_seen >= 0.5 * length(reg_names)) {
    report(paste0(tag, ": input coverage"),
           sprintf("%d/%d regulators (%.0f%%), %d/%d targets (%.0f%%)",
                   reg_seen, length(reg_names), 100 * reg_seen / max(1L, length(reg_names)),
                   tgt_seen, length(tgt_names), 100 * tgt_seen / max(1L, length(tgt_names))))
  } else {
    warn(paste0(tag, ": input coverage"),
         sprintf("only %d of %d input regulators (%.0f%%) reach the network",
                 reg_seen, length(reg_names), 100 * reg_seen / max(1L, length(reg_names))))
  }

  # Hyphenated ids — isoform suffixes, and the "prot-ph" ome tag — are the ones make.names()
  # rewrites, so they are the canary for a sanitising step reappearing anywhere in the path.
  hyphenated_in <- reg_names[grepl("-", sub("\\.\\..*$", "", reg_names))]
  if (length(hyphenated_in) > 0) {
    seen <- intersect(hyphenated_in, unique(edges$Regulator))
    assert(paste0(tag, ": isoform-suffixed regulators reach the network"), length(seen) > 0,
           sprintf("none of %d hyphenated regulator ids appear — make.names() is back somewhere",
                   length(hyphenated_in)))
  }

  # The hub layer: written with an empty Cluster, and joined on the hub names.
  hubs_csv <- file.path(dir_path, "hubs.csv")
  if (file.exists(hubs_csv)) {
    hubs <- tryCatch(as.character(read.csv(hubs_csv, check.names = FALSE)[[2]]),
                     error = function(e) character())
    unmatched <- setdiff(hubs, reg_names)
    assert(paste0(tag, ": hubs join the regulator matrix"),
           length(hubs) > 0 && length(unmatched) == 0,
           paste0(length(unmatched), " of ", length(hubs),
                  " hubs match no regulator rowname — the join that leaves a run with no network file"))
    n_hub_edges <- sum(!in_cluster)
    if (n_hub_edges > 0) {
      report(paste0(tag, ": hub layer"), paste(n_hub_edges, "hub-connection edge(s)"))
    } else {
      warn(paste0(tag, ": hub layer"), "hubs.csv written but no hub-connection edges")
    }
  } else {
    skip(paste0(tag, ": hub layer"), "no hubs.csv")
  }
}

finish()
