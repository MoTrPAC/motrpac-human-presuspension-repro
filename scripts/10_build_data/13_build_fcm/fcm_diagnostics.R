#!/usr/bin/env Rscript
# Diagnostics for the FCM cluster count (step 13).
#
# Upstream's run_cmeans() will, if handed a range of cluster numbers, plot Mfuzz::Dmin() and
# then block on readline() for the operator to type the number to use. A batch build cannot do
# that, and an object whose contents depend on what was typed is not reproducible. This script
# is the other half of that trade: it sweeps the range non-interactively, writes the evidence
# to disk, and leaves the choice to the operator, who feeds it back as FCM_K_<TISSUE> and
# rebuilds.
#
# Input:  scripts/10_build_data/data/FCM_CLUSTERS.rda  (step 13's own first output)
# Output: scripts/10_build_data/13_build_fcm/fcm_diagnostics/
#           fcm_cluster_sweep.tsv        every metric at every k, all tissues
#           dmin_<tissue>.png            the Dmin curve on its own, for a quick look
#           fcm_diagnostics_<tissue>.pdf the full set, one figure per page
#
# The sweep is rerunnable on its own: it reads the scaled input matrix and the weighting
# exponent m back out of FCM_CLUSTERS, so it never has to redo the DA preparation.
#
# ---- What is plotted, and what each plot answers -----------------------------------------
#
# 1. MIN. CENTROID DISTANCE vs k — Mfuzz::Dmin()'s statistic, min(dist(centers)). Falls as
#    clusters are added; the elbow is where extra clusters stop being distinct. This is the
#    curve upstream shows at the prompt.
# 2. MAX. CENTROID CORRELATION vs k — the most similar pair of centroids. Approaching 1 means
#    two clusters have effectively the same trajectory, i.e. k is past useful.
# 3. CORE FRACTION vs k — share of features whose best membership clears `min_prob` (0.3, the
#    threshold run_cluster_ORA() uses for hard assignment). It falls as membership spreads over
#    more clusters; a sharp drop marks the point where clusters stop being well populated.
# 4. SMALLEST CLUSTER vs k — hard-assigned size of the smallest cluster, which catches a k that
#    only adds near-empty clusters.
# 5. CENTROID TRAJECTORIES at the built k — the profiles themselves, for eyeballing duplicates.
# 6. CENTROID CORRELATION HEATMAP at the built k — the same redundancy check as (2), per pair.
# 7. CLUSTER SIZES and MEMBERSHIP DISTRIBUTION at the built k.
#
# Metrics 1-4 come from one fit per k: Mfuzz::Dmin() refits internally and returns only its own
# statistic, so the fit is done here instead and all four are read off it. min_centroid_dist is
# computed exactly as Dmin() computes it, so that curve is the same one.

# Mfuzz declares e1071 and Biobase in Depends, not Imports, so Mfuzz::mfuzz() resolves
# cmeans() off the search path and needs both attached — see run_cmeans.R.
suppressWarnings(suppressMessages({library(Biobase); library(e1071)}))

for (pkg in c("Mfuzz", "Biobase", "e1071"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(pkg, " is required by fcm_diagnostics.R — install it before running step 13.")

ROOT <- Sys.getenv("PRECOVID_ROOT")
if (!nzchar(ROOT)) {
  ROOT <- normalizePath(getwd())
  while (!file.exists(file.path(ROOT, "config", "pipeline.env")) && dirname(ROOT) != ROOT) ROOT <- dirname(ROOT)
}
out_dir  <- file.path(ROOT, "scripts", "10_build_data", "data")
# The sweep is the evidence behind a choice of k, not a data object, so it is written beside
# the step rather than into data/ — which carries only the .rda objects the packages consume.
# check_existing_outputs.sh lists this path as well, so CLEAN=1 clears it with the rest.
diag_dir <- Sys.getenv("FCM_DIAG_DIR")
if (!nzchar(diag_dir))
  diag_dir <- file.path(ROOT, "scripts", "10_build_data", "13_build_fcm", "fcm_diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

# pdf() stamps the wall clock into /CreationDate and /ModDate, so two runs over identical data
# produce files differing in those bytes alone — enough to show up as a change on every
# rebuild. Both are rewritten to a fixed value once the device closes. The replacement is the
# same 14 digits wide, so the xref byte offsets stay valid, and grepRaw works on the raw vector
# rather than a string because compressed PDF streams contain embedded nuls.
PDF_STAMP <- "19700101000000"
freeze_pdf_dates <- function(path) {
  raw <- readBin(path, "raw", file.size(path))
  hit <- FALSE
  for (tag in c("/CreationDate (D:", "/ModDate (D:")) {
    for (at in grepRaw(tag, raw, all = TRUE, fixed = TRUE)) {
      i <- at + nchar(tag)
      raw[i:(i + 13L)] <- charToRaw(PDF_STAMP)
      hit <- TRUE
    }
  }
  if (hit) writeBin(raw, path)
}
load_one <- function(name) {
  f <- file.path(out_dir, paste0(name, ".rda"))
  if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, envir = e); e[[name]] %||% e[[ls(e)[1]]]
}

FCM <- load_one("FCM_CLUSTERS")
if (is.null(FCM))
  stop("missing FCM_CLUSTERS.rda — run FCM_clustering_results.R first: ",
       file.path(out_dir, "FCM_CLUSTERS.rda"))

env_int <- function(var, default, min_value = 1L) {
  v <- Sys.getenv(var)
  if (!nzchar(v)) return(default)
  n <- suppressWarnings(as.integer(v))
  if (is.na(n) || n < min_value) stop(var, " must be an integer >= ", min_value, ", got ", dQuote(v))
  n
}
K_MIN    <- env_int("FCM_DIAG_KMIN", 2L, 2L)
K_MAX    <- env_int("FCM_DIAG_KMAX", 20L, 2L)
REPEATS  <- env_int("FCM_DIAG_REPEATS", 1L)
MIN_PROB <- 0.3   # run_cluster_ORA()'s default; the line drawn on the membership histogram
if (K_MAX < K_MIN) stop("FCM_DIAG_KMAX (", K_MAX, ") is below FCM_DIAG_KMIN (", K_MIN, ")")
KRANGE <- seq.int(K_MIN, K_MAX)

# --- helpers ------------------------------------------------------------------------------

# The matrix run_cmeans() actually clustered: `input` holds both modalities in their original
# column order, and mfuzz saw them reordered Endur-then-Resist. Same expression here, so the
# sweep fits what was built rather than something adjacent to it.
clustered_matrix <- function(fclust) {
  x <- fclust[["input"]]
  x[, c(grep("Endur", colnames(x)), grep("Resist", colnames(x))), drop = FALSE]
}

# Axis labels for the trajectory panels. CONTRAST_CONVERTER$contrast_short is not short enough
# to use here — "Endur.post_15_30_45_min - Control.post_15_30_45_min (delta-delta)" is 66
# characters against a panel a couple of inches wide — so modality and timepoint are taken off
# the contrast itself and abbreviated to "E 15-45m". Anything not matching the expected shape
# falls through unchanged rather than being mangled.
TIMEPOINT_ABBREV <- c(pre_exercise = "pre", during_20_min = "d20m", during_40_min = "d40m",
                      post_10_min = "10m", post_15_30_45_min = "15-45m",
                      post_3.5_4_hr = "3.5-4h", post_24_hr = "24h")
contrast_labels <- function(contrasts) {
  ok <- grepl("^group_timepointADU[A-Za-z]+\\.", contrasts)
  modality  <- sub("^group_timepointADU([A-Za-z]+)\\..*$", "\\1", contrasts)
  timepoint <- sub("^group_timepointADU[A-Za-z]+\\.([^ ]+).*$", "\\1", contrasts)
  tp  <- ifelse(timepoint %in% names(TIMEPOINT_ABBREV),
                unname(TIMEPOINT_ABBREV[timepoint]), timepoint)
  ifelse(ok, paste(substr(modality, 1L, 1L), tp), contrasts)
}

# One mfuzz fit per k, per repeat; four metrics read off each fit, then averaged over repeats.
sweep_k <- function(eset, m, tissue) {
  rows <- list()
  for (rep_i in seq_len(REPEATS)) {
    for (k in KRANGE) {
      set.seed(rep_i)   # fixed per repeat, so the sweep is reproducible
      fit <- Mfuzz::mfuzz(eset = eset, centers = k, m = m)
      centers <- fit[["centers"]]
      cor_mat <- cor(t(centers))
      hard <- table(factor(apply(fit[["membership"]], 1L, which.max), levels = seq_len(k)))
      rows[[length(rows) + 1L]] <- data.frame(
        tissue             = tissue,
        repeat_i           = rep_i,
        k                  = k,
        # Mfuzz::Dmin()'s own statistic: min(dist(cl[[1]])) over the centroids.
        min_centroid_dist  = min(dist(centers)),
        max_centroid_cor   = max(cor_mat[upper.tri(cor_mat)]),
        core_fraction      = mean(apply(fit[["membership"]], 1L, max) >= MIN_PROB),
        min_cluster_size   = min(as.integer(hard)),
        stringsAsFactors   = FALSE
      )
    }
    message(sprintf("  sweep %s: repeat %d/%d done (k = %d..%d)",
                    tissue, rep_i, REPEATS, K_MIN, K_MAX))
  }
  per_rep <- do.call(rbind, rows)
  metrics <- c("min_centroid_dist", "max_centroid_cor", "core_fraction", "min_cluster_size")
  avg <- aggregate(per_rep[, metrics], by = list(k = per_rep$k), FUN = mean)
  avg <- avg[order(avg$k), ]
  data.frame(tissue = tissue, avg, row.names = NULL, stringsAsFactors = FALSE)
}

# A sweep curve with the built k marked. `better` says which end of the y-axis is the good one,
# so the reader is not left to guess which direction on this particular metric is improvement.
plot_sweep <- function(d, metric, ylab, k_built, better) {
  plot(d$k, d[[metric]], type = "b", pch = 19, cex = 0.7,
       xlab = "Number of clusters (k)", ylab = ylab,
       main = sprintf("%s  (%s is better)", ylab, better))
  abline(v = k_built, col = "firebrick", lty = 2, lwd = 2)
  y <- d[[metric]][d$k == k_built]
  if (length(y) == 1) {
    points(k_built, y, col = "firebrick", pch = 19, cex = 1.3)
    mtext(sprintf("built k = %d", k_built), side = 3, col = "firebrick", cex = 0.8, adj = 1)
  }
  grid(col = "grey88")
}

plot_trajectories <- function(centers, tissue) {
  k <- nrow(centers)
  labs <- contrast_labels(colnames(centers))
  # Endur columns come first, Resist second (run_cmeans orders them that way); the divider
  # says where one ends and the other begins without repeating the modality on every tick.
  n_endur <- sum(grepl("Endur", colnames(centers)))
  ncol_panel <- min(k, 4L)
  op <- par(mfrow = c(ceiling(k / ncol_panel), ncol_panel),
            mar = c(3.8, 3.4, 2, 0.6), mgp = c(2.1, 0.6, 0), oma = c(0, 0, 2.4, 0))
  on.exit(par(op))
  ylim <- range(centers)
  for (i in seq_len(k)) {
    plot(seq_len(ncol(centers)), centers[i, ], type = "b", pch = 19, ylim = ylim,
         xaxt = "n", xlab = "", ylab = "centroid", main = sprintf("cluster %d", i),
         cex.main = 0.95, cex.lab = 0.85)
    axis(1, at = seq_len(ncol(centers)), labels = labs, las = 2, cex.axis = 0.62)
    abline(h = 0, col = "grey60", lty = 3)
    if (n_endur > 0 && n_endur < ncol(centers))
      abline(v = n_endur + 0.5, col = "grey60", lty = 2)
  }
  mtext(sprintf("%s - cluster centroids (k = %d); E = Endur, R = Resist", tissue, k),
        outer = TRUE, cex = 0.95, font = 2)
}

plot_centroid_cor <- function(centers, tissue) {
  k <- nrow(centers)
  cm <- cor(t(centers))
  op <- par(mar = c(4, 4, 4, 5))
  on.exit(par(op))
  cols <- hcl.colors(101, "RdBu", rev = TRUE)
  image(seq_len(k), seq_len(k), cm[, rev(seq_len(k)), drop = FALSE], zlim = c(-1, 1),
        col = cols, xlab = "cluster", ylab = "cluster", axes = FALSE,
        main = sprintf("%s - centroid correlation (k = %d)", tissue, k))
  axis(1, at = seq_len(k), cex.axis = 0.7); axis(2, at = seq_len(k), labels = rev(seq_len(k)), las = 1, cex.axis = 0.7)
  box()
  # The diagonal is 1 by construction, so the largest off-diagonal entry is the redundancy
  # summary for this k — a pair above ~0.9 is effectively one trajectory split in two.
  offdiag <- cm; diag(offdiag) <- NA
  mtext(sprintf("max off-diagonal r = %.3f", max(offdiag, na.rm = TRUE)), side = 3, cex = 0.8)
}

plot_sizes_and_membership <- function(fclust, tissue) {
  op <- par(mfrow = c(2, 1), mar = c(4.2, 4.2, 3, 1))
  on.exit(par(op))
  barplot(fclust[["size"]], names.arg = seq_along(fclust[["size"]]),
          xlab = "cluster", ylab = "features (hard-assigned)",
          main = sprintf("%s - cluster sizes", tissue), col = "grey75", border = NA)
  best <- apply(fclust[["membership"]], 1L, max)
  hist(best, breaks = 50, col = "grey75", border = NA,
       xlab = "max membership probability per feature", ylab = "features",
       main = sprintf("%s - membership (%.1f%% >= %.2f)",
                      tissue, 100 * mean(best >= MIN_PROB), MIN_PROB))
  abline(v = MIN_PROB, col = "firebrick", lty = 2, lwd = 2)
}

# --- run ----------------------------------------------------------------------------------

message(sprintf("FCM diagnostics: k = %d..%d, %d repeat(s), %d tissue(s) -> %s",
                K_MIN, K_MAX, REPEATS, length(FCM), diag_dir))

sweeps <- list()
for (tissue in names(FCM)) {
  fclust <- FCM[[tissue]]
  mat <- clustered_matrix(fclust)
  m   <- fclust[["call"]][["m"]]
  k_built <- nrow(fclust[["centers"]])
  if (is.null(m))
    stop("FCM_CLUSTERS[['", tissue, "']]$call$m is missing — run_cmeans() records it; ",
         "rebuild FCM_CLUSTERS with the current step 13.")

  message(sprintf("%s: %d features x %d contrasts, m = %.4f, built k = %d",
                  tissue, nrow(mat), ncol(mat), m, k_built))

  eset <- Biobase::ExpressionSet(assayData = mat)
  d <- sweep_k(eset = eset, m = m, tissue = tissue)
  d$k_built <- k_built
  sweeps[[tissue]] <- d

  # Dmin curve alone, as a PNG — the one plot that answers "why this k".
  png(file.path(diag_dir, sprintf("dmin_%s.png", tissue)),
      width = 1200, height = 900, res = 150)
  plot_sweep(d, "min_centroid_dist", "Min. centroid distance", k_built, better = "elbow")
  title(sub = sprintf("%s - Mfuzz::Dmin statistic, %d repeat(s)", tissue, REPEATS))
  dev.off()

  pdf_path <- file.path(diag_dir, sprintf("fcm_diagnostics_%s.pdf", tissue))
  pdf(pdf_path, width = 9, height = 7, onefile = TRUE)
  plot_sweep(d, "min_centroid_dist", "Min. centroid distance", k_built, better = "elbow")
  title(sub = sprintf("%s - %d features, m = %.4f", tissue, nrow(mat), m))
  plot_sweep(d, "max_centroid_cor",  "Max. centroid correlation", k_built, better = "lower")
  plot_sweep(d, "core_fraction",     sprintf("Fraction with membership >= %.2f", MIN_PROB),
             k_built, better = "higher")
  plot_sweep(d, "min_cluster_size",  "Smallest cluster (features)", k_built, better = "higher")
  plot_trajectories(fclust[["centers"]], tissue)
  plot_centroid_cor(fclust[["centers"]], tissue)
  plot_sizes_and_membership(fclust, tissue)
  dev.off()
  freeze_pdf_dates(pdf_path)
}

sweep_tsv <- do.call(rbind, sweeps)
write.table(sweep_tsv, file = file.path(diag_dir, "fcm_cluster_sweep.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

for (tissue in names(sweeps)) {
  d <- sweeps[[tissue]]
  k_built <- d$k_built[1]
  at <- d[d$k == k_built, ]
  if (nrow(at) == 1)
    message(sprintf("%s @ built k = %d: min centroid dist %.3f | max centroid cor %.3f | core %.1f%% | smallest cluster %d",
                    tissue, k_built, at$min_centroid_dist, at$max_centroid_cor,
                    100 * at$core_fraction, round(at$min_cluster_size)))
  else
    message(sprintf("%s: built k = %d is outside the swept range %d..%d — widen FCM_DIAG_KMIN/KMAX to place it on the curves",
                    tissue, k_built, K_MIN, K_MAX))
}
message("FCM diagnostics written to ", diag_dir)
