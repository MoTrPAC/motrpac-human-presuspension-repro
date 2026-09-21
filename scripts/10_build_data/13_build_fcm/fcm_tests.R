#!/usr/bin/env Rscript
# Tests for FCM_CLUSTERS, FCM_CAMERA and FCM_ORA (step 13) — fuzzy c-means clustering of the
# acute DA results and the two enrichment tests over the clusters.
#
# The expected shape is re-derived from the inputs (FCM_CLUSTERS itself, SET_TO_ID,
# CONTRAST_CONVERTER) rather than imported from the builder, so a mistake in the builder's
# scope or joins shows up as a failure instead of being mirrored.
.root <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(.root)) { .root <- normalizePath(getwd()); while (!file.exists(file.path(.root, "config", "pipeline.env")) && dirname(.root) != .root) .root <- dirname(.root) }
source(file.path(.root, "scripts", "10_build_data", "lib", "test_helpers.R"))
suppressWarnings(suppressMessages({library(dplyr)}))
cat("== FCM tests ==\n")

fcm <- load_built("FCM_CLUSTERS")
assert("FCM_CLUSTERS built", !is.null(fcm), "object missing — run step 13")
if (is.null(fcm)) finish()

TISSUES <- c("adipose", "blood", "muscle")
# Clustering covers five omes; enrichment covers four — prot-ol is absent from
# .prepare_cluster_mem()'s choices upstream, so blood prot-ol is clustered but not tested.
CLUSTER_OMES <- c("transcript-rna-seq", "prot-pr", "prot-ol", "prot-ph", "metab")
ENRICH_OMES  <- c("transcript-rna-seq", "prot-pr", "prot-ph", "metab")

k_expected <- function(tissue) {
  v <- Sys.getenv(paste0("FCM_K_", toupper(tissue)))
  if (nzchar(v)) as.integer(v) else c(adipose = 13L, blood = 12L, muscle = 12L)[[tissue]]
}

# ── FCM_CLUSTERS: structure ───────────────────────────────────────────────────
assert("FCM_CLUSTERS: is a named list", is.list(fcm) && !is.null(names(fcm)))
assert_subset("FCM_CLUSTERS: tissues valid", names(fcm), TISSUES)
assert_completeness("FCM_CLUSTERS: every tissue clustered", names(fcm), TISSUES)

cc <- load_built("CONTRAST_CONVERTER")
ewc <- if (is.null(cc)) NULL else
  as.character(as.data.frame(cc)$contrast[as.character(as.data.frame(cc)$contrast_type) == "exercise_with_controls"])

k_by_tissue <- setNames(integer(0), character(0))
for (tissue in names(fcm)) {
  x <- fcm[[tissue]]
  d <- sprintf("FCM_CLUSTERS[%s]", tissue)
  assert(sprintf("%s: has centers/size/cluster/membership/input", d),
         all(c("centers", "size", "cluster", "membership", "input", "call") %in% names(x)),
         sprintf("has [%s]", paste(names(x), collapse = ", ")))
  if (!all(c("centers", "membership", "input") %in% names(x))) next

  k <- nrow(x$centers); n <- nrow(x$membership)
  k_by_tissue[[tissue]] <- k
  assert(sprintf("%s: k matches the build parameter", d), k == k_expected(tissue),
         sprintf("built %d, expected %d", k, k_expected(tissue)))
  assert(sprintf("%s: membership is %d features x %d clusters", d, n, k),
         ncol(x$membership) == k && ncol(x$centers) == ncol(x$input) &&
           nrow(x$input) == n)

  # .reorder_clusters() relabels centers, membership columns, size and cluster together; if it
  # renumbered any one of them out of step with the others the labels would be meaningless.
  assert(sprintf("%s: centers relabelled 1..k", d),
         identical(rownames(x$centers), as.character(seq_len(k))))
  assert(sprintf("%s: membership columns relabelled 1..k", d),
         identical(colnames(x$membership), as.character(seq_len(k))))
  assert(sprintf("%s: cluster assignment in 1..k", d),
         all(x$cluster >= 1L & x$cluster <= k) && length(x$cluster) == n)
  assert(sprintf("%s: cluster == which.max(membership)", d),
         identical(as.integer(x$cluster), as.integer(apply(x$membership, 1L, which.max))),
         "hard assignment does not match the reordered membership matrix")
  assert(sprintf("%s: size sums to the feature count", d), sum(x$size) == n,
         sprintf("sum(size) = %d, features = %d", sum(x$size), n))

  # Membership probabilities: each feature's row is a distribution over the clusters.
  assert(sprintf("%s: membership in [0,1]", d),
         all(x$membership >= 0 & x$membership <= 1, na.rm = TRUE))
  rs <- rowSums(x$membership)
  assert(sprintf("%s: membership rows sum to 1", d),
         isTRUE(all.equal(rs, setNames(rep(1, n), names(rs)), tolerance = 1e-6)),
         sprintf("range %.6f..%.6f", min(rs), max(rs)))
  assert(sprintf("%s: input is finite", d), all(is.finite(x$input)))
  assert(sprintf("%s: m recorded on the call", d),
         is.numeric(x$call$m) && length(x$call$m) == 1 && x$call$m > 1,
         "fcm_diagnostics.R reads m back out of this")

  # Scope: rownames are "<ome> <feature_id>".
  omes <- sort(unique(sub(" .*", "", rownames(x$membership))))
  assert_subset(sprintf("%s: omes valid", d), omes, CLUSTER_OMES)
  # run_cmeans() drops adipose prot-pr and prot-ph for insufficient timepoints.
  if (tissue == "adipose")
    assert(sprintf("%s: adipose prot-pr/prot-ph excluded", d),
           !any(c("prot-pr", "prot-ph") %in% omes), paste(omes, collapse = ","))

  # Contrasts: exercise_with_controls only, and the two `during` contrasts are dropped.
  assert(sprintf("%s: no `during` contrast clustered", d),
         !any(grepl("during", colnames(x$centers))))
  if (!is.null(ewc))
    assert_subset(sprintf("%s: contrasts are exercise_with_controls", d),
                  colnames(x$centers), ewc)
}

# ── FCM_CAMERA / FCM_ORA: shared contract ─────────────────────────────────────
CAM_COLS <- c("tissue", "assay", "cluster", "collection", "database", "set_id", "set",
              "set_short", "set_size", "set_size_DB", "size_ratio", "p_value", "adj_p_value")
ORA_COLS <- c("tissue", "assay", "cluster", "collection", "database", "set_id", "set",
              "set_short", "set_size", "set_size_DB", "size_ratio",
              "set_size_in_cluster", "cluster_size", "background_size",
              "p_value", "adj_p_value")

s2i <- load_built("SET_TO_ID")
enrich <- list()

for (nm in c("FCM_CAMERA", "FCM_ORA")) {
  obj <- load_built(nm)
  assert(sprintf("%s built", nm), !is.null(obj), "object missing — run step 13")
  if (is.null(obj)) next
  obj <- as.data.frame(obj)
  enrich[[nm]] <- obj
  cols <- if (nm == "FCM_CAMERA") CAM_COLS else ORA_COLS

  assert_cols(nm, obj, cols)
  assert(sprintf("%s: column order", nm), identical(colnames(obj), cols),
         sprintf("got %s", paste(colnames(obj), collapse = " ")))
  assert(sprintf("%s: non-empty", nm), nrow(obj) > 0)

  # Scope.
  assert_subset(sprintf("%s: assays valid", nm), unique(as.character(obj$assay)), ENRICH_OMES)
  assert(sprintf("%s: prot-ol not enrichment-tested", nm),
         !"prot-ol" %in% as.character(obj$assay),
         "prot-ol is clustered but is not in .prepare_cluster_mem()'s choices")
  # prot-ph is the regression guard for the step-07 flanking_sequence fix: cluster membership
  # keys phosphosites on the flanking sequence, so dropping that column would silently reduce
  # this object to the three non-phospho omes rather than failing.
  assert(sprintf("%s: prot-ph was tested (needs HUMAN_FEATURE_TO_GENE$flanking_sequence)", nm),
         "prot-ph" %in% as.character(obj$assay))
  assert_subset(sprintf("%s: tissues valid", nm), unique(as.character(obj$tissue)), names(fcm))

  # Cluster labels must stay inside the k that tissue was clustered at.
  bad_k <- obj %>%
    dplyr::mutate(cl = as.integer(as.character(cluster)),
                  kk = k_by_tissue[as.character(tissue)]) %>%
    dplyr::filter(is.na(cl) | is.na(kk) | cl < 1L | cl > kk)
  assert(sprintf("%s: cluster labels within each tissue's k", nm), nrow(bad_k) == 0,
         sprintf("%d row(s) out of range", nrow(bad_k)))

  # Set metadata came from SET_TO_ID.
  if (is.null(s2i)) skip(sprintf("%s: sets come from SET_TO_ID", nm), "SET_TO_ID missing") else
    assert_subset(sprintf("%s: sets come from SET_TO_ID", nm),
                  unique(as.character(obj$set)), unique(as.character(s2i$set)))
  for (cl in c("collection", "database", "set_id", "set"))
    assert(sprintf("%s: %s fully joined", nm, cl), !anyNA(obj[[cl]]),
           sprintf("%d NA", sum(is.na(obj[[cl]]))))

  # Statistics.
  assert(sprintf("%s: p_value in [0,1]", nm), all(obj$p_value >= 0 & obj$p_value <= 1, na.rm = TRUE))
  assert(sprintf("%s: adj_p_value in [0,1]", nm), all(obj$adj_p_value >= 0 & obj$adj_p_value <= 1, na.rm = TRUE))
  bad <- sum(obj$adj_p_value < obj$p_value - 1e-9, na.rm = TRUE)
  assert(sprintf("%s: adj_p_value >= p_value", nm), bad == 0, sprintf("%d row(s) below", bad))
  assert(sprintf("%s: set_size >= 5", nm), all(obj$set_size >= 5L, na.rm = TRUE),
         sprintf("min = %s", min(obj$set_size, na.rm = TRUE)))
  assert(sprintf("%s: set_size <= set_size_DB", nm), all(obj$set_size <= obj$set_size_DB, na.rm = TRUE),
         sprintf("%d row(s) exceed", sum(obj$set_size > obj$set_size_DB, na.rm = TRUE)))
  sr <- round(obj$set_size / obj$set_size_DB, 3L)
  assert(sprintf("%s: size_ratio = set_size / set_size_DB", nm),
         isTRUE(all.equal(sr, obj$size_ratio, tolerance = 1e-9)))

  # BH is applied within tissue x assay x collection x cluster; re-derive one group and
  # compare, which catches an adjustment done over the wrong grouping (or globally).
  grp <- obj %>% dplyr::count(tissue, assay, collection, cluster) %>%
    dplyr::filter(n >= 20) %>% dplyr::slice(1)
  if (nrow(grp) == 1) {
    g <- obj %>% dplyr::semi_join(grp, by = c("tissue", "assay", "collection", "cluster"))
    assert(sprintf("%s: BH adjusted within tissue x assay x collection x cluster", nm),
           isTRUE(all.equal(p.adjust(g$p_value, method = "BH"), g$adj_p_value, tolerance = 1e-9)),
           sprintf("group %s/%s n=%d", grp$tissue[1], grp$assay[1], grp$n[1]))
  } else skip(sprintf("%s: BH adjusted within tissue x assay x collection x cluster", nm),
              "no group with >= 20 rows")

  assert_unique(sprintf("%s: one row per tissue x assay x cluster x set", nm),
                paste(obj$tissue, obj$assay, obj$cluster, obj$set))
}

# ── FCM_ORA: the counts its p-value is computed from ──────────────────────────
if (!is.null(enrich$FCM_ORA)) {
  o <- enrich$FCM_ORA
  assert("FCM_ORA: set_size_in_cluster <= set_size",
         all(o$set_size_in_cluster <= o$set_size, na.rm = TRUE),
         sprintf("%d row(s) exceed", sum(o$set_size_in_cluster > o$set_size, na.rm = TRUE)))
  assert("FCM_ORA: set_size_in_cluster <= cluster_size",
         all(o$set_size_in_cluster <= o$cluster_size, na.rm = TRUE),
         sprintf("%d row(s) exceed", sum(o$set_size_in_cluster > o$cluster_size, na.rm = TRUE)))
  assert("FCM_ORA: cluster_size <= background_size",
         all(o$cluster_size <= o$background_size, na.rm = TRUE))
  # cluster_size is the count at membership >= min_prob, so summing it over the clusters of one
  # tissue x ome cannot exceed that dataset's background.
  tot <- o %>% dplyr::distinct(tissue, assay, cluster, cluster_size, background_size) %>%
    dplyr::summarise(.by = c(tissue, assay), s = sum(cluster_size), b = dplyr::first(background_size))
  assert("FCM_ORA: cluster sizes sum to at most the background",
         all(tot$s <= tot$b), sprintf("%d dataset(s) exceed", sum(tot$s > tot$b)))
  # Re-derive one hypergeometric p-value from the counts in the row itself.
  r <- o[which.max(o$set_size_in_cluster), ]
  assert("FCM_ORA: p_value is the upper-tail hypergeometric of its own counts",
         isTRUE(all.equal(phyper(q = r$set_size_in_cluster - 1L, m = r$set_size,
                                 n = r$background_size - r$set_size, k = r$cluster_size,
                                 lower.tail = FALSE), r$p_value, tolerance = 1e-9)))
}

# ── the two tests must describe the same grid ─────────────────────────────────
# FCM_CAMERA and FCM_ORA run over the same membership matrices and the same filtered sets, and
# the object exists to be compared against FCM_CAMERA — a differing grid would mean they were
# not testing the same thing.
if (!is.null(enrich$FCM_CAMERA) && !is.null(enrich$FCM_ORA)) {
  key <- function(d) paste(d$tissue, d$assay, d$cluster, d$set)
  kc <- key(enrich$FCM_CAMERA); ko <- key(enrich$FCM_ORA)
  assert("FCM_CAMERA and FCM_ORA cover the same tissue x assay x cluster x set grid",
         setequal(kc, ko),
         sprintf("CAMERA %d rows, ORA %d rows, %d only-CAMERA, %d only-ORA",
                 length(kc), length(ko), length(setdiff(kc, ko)), length(setdiff(ko, kc))))
}

# ── diagnostics artifacts ─────────────────────────────────────────────────────
diag_dir <- Sys.getenv("FCM_DIAG_DIR")
if (!nzchar(diag_dir))
  diag_dir <- file.path(.T_ROOT, "scripts", "10_build_data", "13_build_fcm", "fcm_diagnostics")
sweep_tsv <- file.path(diag_dir, "fcm_cluster_sweep.tsv")
if (!file.exists(sweep_tsv)) {
  skip("FCM diagnostics written", sprintf("no %s (FCM_DIAG=0?)", sweep_tsv))
} else {
  sw <- read.csv(sweep_tsv, sep = "\t", check.names = FALSE)
  assert_cols("fcm_cluster_sweep.tsv", sw,
              c("tissue", "k", "min_centroid_dist", "max_centroid_cor",
                "core_fraction", "min_cluster_size", "k_built"))
  assert_completeness("FCM diagnostics: every tissue swept",
                      unique(as.character(sw$tissue)), names(fcm))
  for (tissue in names(fcm)) {
    for (f in c(sprintf("dmin_%s.png", tissue), sprintf("fcm_diagnostics_%s.pdf", tissue))) {
      p <- file.path(diag_dir, f)
      assert(sprintf("FCM diagnostics: %s written", f),
             file.exists(p) && file.size(p) > 0, "missing or empty")
    }
    # The built k has to land on the swept curves, or the plots cannot speak to it.
    swt <- sw[as.character(sw$tissue) == tissue, ]
    kb <- k_by_tissue[tissue]   # single bracket: NA rather than an error if the tissue bailed
    if (nrow(swt) == 0 || is.na(kb)) {
      skip(sprintf("FCM diagnostics: %s built k is on the swept range", tissue),
           if (is.na(kb)) "tissue has no usable FCM_CLUSTERS entry" else "tissue not in the sweep")
    } else {
      assert(sprintf("FCM diagnostics: %s built k is on the swept range", tissue),
             kb %in% swt$k,
             sprintf("built k = %d, swept %d..%d", kb, min(swt$k), max(swt$k)))
    }
  }
}

# ── diff vs the shipped package objects (INFO only) ───────────────────────────
# Never a FAIL: these are computed from locally refit DA, so they move with the freeze — and
# mfuzz is seeded but its result still depends on the exact input matrix.

# Adjusted Rand index from a contingency table — label-free, so it needs no cluster matching.
.ari <- function(tab) {
  n <- sum(tab); if (n < 2L) return(NA_real_)
  ch2 <- function(x) sum(x * (x - 1) / 2)
  a <- ch2(rowSums(tab)); b <- ch2(colSums(tab)); e <- a * b / ch2(n)
  (ch2(tab) - e) / ((a + b) / 2 - e)
}

# Optimal 1-1 matching of built clusters onto package clusters, maximising total overlap.
# clue::solve_LSAP is exact; the greedy fallback keeps the tests running without it.
.match_clusters <- function(tab) {
  k <- min(dim(tab))
  if (requireNamespace("clue", quietly = TRUE))
    return(as.integer(clue::solve_LSAP(tab[seq_len(k), seq_len(k), drop = FALSE],
                                       maximum = TRUE)))
  out <- rep(NA_integer_, nrow(tab)); tt <- tab
  for (i in seq_len(k)) {
    w <- which(tt == max(tt), arr.ind = TRUE)[1L, ]
    out[w[["row"]]] <- w[["col"]]; tt[w[["row"]], ] <- -1L; tt[, w[["col"]]] <- -1L
  }
  out
}

for (nm in c("FCM_CLUSTERS", "FCM_CAMERA", "FCM_ORA")) {
  rda <- .find_rda(nm)
  if (is.na(rda)) { report(sprintf("diff vs package: %s", nm), "no package .rda"); next }
  p <- .load_rda(rda, nm)
  if (nm == "FCM_CLUSTERS") {
    b <- fcm
    report("diff vs package: FCM_CLUSTERS",
           paste(vapply(intersect(names(b), names(p)), function(t)
             sprintf("%s built %dx%d vs pkg %dx%d", t,
                     nrow(b[[t]]$membership), ncol(b[[t]]$membership),
                     nrow(p[[t]]$membership), ncol(p[[t]]$membership)),
             character(1)), collapse = " | "))

    # A cluster NUMBER is not a cluster IDENTITY. .reorder_clusters() renumbers by an hclust
    # dendrogram over centroid correlations, and that dendrogram's order is only defined up to
    # flipping subtrees — so built cluster 4 and package cluster 4 may describe entirely
    # different trajectories, and a rebuild can permute or even reverse the whole numbering.
    # What the numbering does guarantee is that CONSECUTIVE clusters are adjacent trajectories.
    # Hence: match built to package clusters optimally first, then score
    #   exact   the matched cluster is the package's
    #   +-1     the matched cluster is the package's or either neighbour
    # `+-1 by chance` is what that wider window scores under random assignment, since it spans
    # roughly three of k clusters and would otherwise look impressive on its own.
    report("cluster agreement vs package",
           paste("cluster numbers are a trajectory ORDER, not identities — the same number in",
                 "two builds need not be the same cluster, so agreement is scored after an",
                 "optimal 1-1 match; +-1 also accepts either neighbour, because consecutive",
                 "numbers are adjacent trajectories in the centroid dendrogram"))
    for (t in intersect(names(b), names(p))) {
      kb <- ncol(b[[t]]$membership); kp <- ncol(p[[t]]$membership)
      sh <- intersect(rownames(b[[t]]$membership), rownames(p[[t]]$membership))
      if (!length(sh)) { report(sprintf("cluster agreement: %s", t), "no shared features"); next }
      bc <- factor(apply(b[[t]]$membership[sh, , drop = FALSE], 1L, which.max),
                   levels = seq_len(kb))
      pc <- factor(apply(p[[t]]$membership[sh, , drop = FALSE], 1L, which.max),
                   levels = seq_len(kp))
      tab <- table(bc, pc)
      mapped <- .match_clusters(tab)[as.integer(bc)]
      pcn <- as.integer(pc)
      # Chance baseline: score the package's own cluster distribution against a +-1 window
      # centred on each cluster in turn, so end clusters get their narrower window counted.
      pk <- as.integer(table(pc)) / length(pcn)
      chance <- mean(vapply(seq_len(kp), function(j)
        sum(pk[intersect(seq(j - 1L, j + 1L), seq_len(kp))]), numeric(1)))
      report(sprintf("cluster agreement: %s", t),
             sprintf("shared %d | k %d vs %d | ARI %.3f | exact %.1f%% | +-1 %.1f%% | +-1 by chance %.1f%%",
                     length(sh), kb, kp, .ari(tab),
                     100 * mean(mapped == pcn, na.rm = TRUE),
                     100 * mean(abs(mapped - pcn) <= 1L, na.rm = TRUE),
                     100 * chance))
    }
  } else {
    b <- enrich[[nm]]; if (is.null(b)) next
    p <- as.data.frame(p)
    report(sprintf("diff vs package: %s", nm),
           sprintf("built %d x %d vs pkg %d x %d | cols %s | assays built [%s] pkg [%s]",
                   nrow(b), ncol(b), nrow(p), ncol(p),
                   if (identical(colnames(b), colnames(p))) "identical" else "DIFFER",
                   paste(sort(unique(as.character(b$assay))), collapse = ","),
                   paste(sort(unique(as.character(p$assay))), collapse = ",")))
  }
}

finish()
