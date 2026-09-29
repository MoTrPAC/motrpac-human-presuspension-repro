# Stage 3 carry, step 2 — decide where every built object goes, and account for
# everything that is left over on both sides.
#
# Two ledgers have to agree before anything is copied: the routing rule in
# lib/routing.R, and docs/data_objects.tsv. They are maintained separately, so
# this step re-derives the map and fails if the inventory disagrees — that
# disagreement is exactly the bug that made the inventory assign 42 QC objects to
# the wrong package.
#
# Usage:
#   Rscript 02_route_objects.R --data-pkg <dir> --analysis-pkg <dir> \
#     --build-data <dir> --preflight-data <dir> --inventory <tsv> --out <dir>
#
# Writes:
#   <out>/routing.tsv    one row per object, on either side, with an action
#   <out>/route_report.tsv
# Exits non-zero if the inventory and the rule disagree, or if the row counts
# do not balance.

.here <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
source(file.path(.here, "lib", "carry_helpers.R"))
source(file.path(.here, "lib", "routing.R"))

opt <- parse_opts(commandArgs(trailingOnly = TRUE),
                  c("data-pkg", "analysis-pkg", "build-data", "preflight-data",
                    "inventory", "out"))
for (need in c("data-pkg", "analysis-pkg", "build-data", "out")) {
  if (is.null(opt[[need]])) stop("missing required option --", need)
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
rpt <- new_report()

pkg_dir <- c(data = file.path(opt$`data-pkg`, "data"),
             analysis = file.path(opt$`analysis-pkg`, "data"))

# ---- Sources ---------------------------------------------------------------
# Stage 1 writes .rda; Stage 0 writes .rds for the objects it builds from
# vendored sources. Both are package data once they land, but an .rds dropped
# into data/ would load under a name taken from the file rather than from the
# object, so the format difference has to survive to step 3.

sources <- list()
rda <- list.files(opt$`build-data`, pattern = "\\.rda$")
if (length(rda)) {
  sources[[length(sources) + 1L]] <-
    data.frame(object = sub("\\.rda$", "", rda),
               source_path = file.path(opt$`build-data`, rda),
               source_format = "rda", source_stage = "10_build_data",
               stringsAsFactors = FALSE)
}
if (!is.null(opt$`preflight-data`) && dir.exists(opt$`preflight-data`)) {
  rds <- list.files(opt$`preflight-data`, pattern = "\\.rds$")
  # QUANTID_BUCKET_FILES is pipeline configuration that step 06 of Stage 1 reads;
  # it has never been package data and must not become it.
  rds <- rds[sub("\\.rds$", "", rds) != "QUANTID_BUCKET_FILES"]
  if (length(rds)) {
    sources[[length(sources) + 1L]] <-
      data.frame(object = sub("\\.rds$", "", rds),
                 source_path = file.path(opt$`preflight-data`, rds),
                 source_format = "rds", source_stage = "00_preflight",
                 stringsAsFactors = FALSE)
  }
}
src <- do.call(rbind, sources)
if (is.null(src) || !nrow(src)) stop("no source objects found")

routed <- route_objects(src$object)
src$destination <- routed$destination
src$reason <- routed$reason

record(rpt, "PASS", "sources",
       sprintf("%d object(s): %d .rda from Stage 1, %d .rds from Stage 0",
               nrow(src), sum(src$source_format == "rda"), sum(src$source_format == "rds")))

# ---- Cross-check against docs/data_objects.tsv ------------------------------

if (!is.null(opt$inventory) && file.exists(opt$inventory)) {
  inv <- read_tsv(opt$inventory)
  # Build assets that are not package data objects; the inventory says so in its
  # own notes and the routing rule has no opinion about them.
  inv <- inv[!inv$object %in% c("gmt_files", "DA"), , drop = FALSE]
  m <- merge(src[, c("object", "destination")], inv[, c("object", "repo")],
             by = "object", all.x = TRUE)
  missing_rows <- m$object[is.na(m$repo)]
  disagree <- m[!is.na(m$repo) & m$repo != m$destination, ]

  if (length(missing_rows)) {
    record(rpt, "FAIL", "inventory:coverage",
           sprintf("%d built object(s) absent from the inventory: %s",
                   length(missing_rows), paste(missing_rows, collapse = ", ")))
  } else {
    record(rpt, "PASS", "inventory:coverage", "every built object has an inventory row")
  }

  if (nrow(disagree)) {
    record(rpt, "FAIL", "inventory:agreement",
           sprintf("%d object(s) routed against the inventory, e.g. %s (rule=%s, inventory=%s)",
                   nrow(disagree), disagree$object[1], disagree$destination[1], disagree$repo[1]))
  } else {
    record(rpt, "PASS", "inventory:agreement",
           "the routing rule and docs/data_objects.tsv agree on every object")
  }
} else {
  record(rpt, "WARN", "inventory", "not supplied — routing is unchecked against the inventory")
}

# ---- Actions ----------------------------------------------------------------

src$action <- NA_character_
src$action[src$destination == "bucket"] <- "skip-bucket-only"

for (i in which(src$destination %in% names(pkg_dir))) {
  target <- file.path(pkg_dir[[src$destination[i]]], paste0(src$object[i], ".rda"))
  src$action[i] <- if (file.exists(target)) "replace" else "add"
}

# ---- Objects the packages ship that nothing produces -------------------------
# Named individually rather than counted. An orphan is not necessarily wrong —
# some are deliberately kept — but an unnamed one is indistinguishable from an
# object the build silently dropped.
#
# An orphan goes one of two ways. Kept, and the package's existing bytes are
# carried into the release unchanged; or dropped, and the release ships without
# it. Dropping is the answer when the object would otherwise go out at a version
# older than everything around it.

# None at present. Analysis 2.0.7 removed the assay_codes re-export; callers read
# MotrpacBicQC::assay_codes.
ORPHAN_REASON <- character(0)

# Withdrawn from the release. Three different reasons, and the distinction is
# worth keeping visible in the NEWS entry.
#
# The clinical-chemistry pair is a RENAME, not a loss. v1.3 carried clinical
# chemistry as a single assay; v2.0 splits it into a metabolomics and a
# proteomics assay, and all four split objects are built this cycle. Keeping the
# combined pair beside them would ship the same measurements twice under two
# schemas, with only the v1.3 pair carrying stale values — so the split replaces
# them rather than joining them. Every shipped call site that read the combined
# objects is rewired onto the split in step 4.
#
# Blood ATAC summary statistics are a genuine absence: epigenomics summary
# statistics are trimmed to significant features only, for file size, and blood
# ATAC has none this cycle — 0 of 5,279,211 features reach adj_p_value < 0.05,
# the smallest being 0.0501 — so step 11 produces no object rather than an empty
# one.
#
# The per-platform metabolomics summary statistics are the third case, and also a rename.
# Step 11 now stacks the research platforms into one {TISSUE}_METAB_SUM_STATS per tissue,
# which is the key the *_DA objects have always used, so every row the per-platform objects
# held is in the stack with its platform in the platform column. Read off what the packages
# actually ship rather than hand-kept: a platform that was never carried needs no entry, and
# a platform added later needs no edit here. BLOOD_METAB_T_CLINICAL_SUM_STATS is excluded —
# clinical chemistry stays its own object on both tiers.
metab_per_platform <- unique(unlist(unname(lapply(pkg_dir, function(d) {
  o <- sub("\\.rda$", "", list.files(d, pattern = "^[A-Z]+_METAB_[A-Z_]+_SUM_STATS\\.rda$"))
  setdiff(o, "BLOOD_METAB_T_CLINICAL_SUM_STATS")
}))))

ORPHAN_DROP <- c(
  BLOOD_EPIGEN_ATAC_SEQ_SUM_STATS = "no blood ATAC feature is significant at adj_p_value < 0.05 this cycle, and epigenomics summary statistics carry significant features only, so no object is built",
  BLOOD_CLINICAL_CHEMISTRY_SUM_STATS = "replaced by the v2.0 split into BLOOD_METAB_T_CLINICAL_SUM_STATS and BLOOD_PROT_CLINICAL_SUM_STATS",
  CLIN_CHEMISTRY_DA                 = "replaced by the v2.0 split into BLOOD_METAB_T_CLINICAL_DA and BLOOD_PROT_CLINICAL_DA",
  stats::setNames(
    rep("replaced by the per-tissue {TISSUE}_METAB_SUM_STATS stack, which carries the same rows with the platform in its own column, matching {TISSUE}_METAB_DA",
        length(metab_per_platform)),
    metab_per_platform)
)

orphan_reasons <- c(ORPHAN_REASON, ORPHAN_DROP)

orphans <- list()
for (dest in names(pkg_dir)) {
  have <- sub("\\.rda$", "", list.files(pkg_dir[[dest]], pattern = "\\.rda$"))
  left <- setdiff(have, src$object[src$destination == dest])
  for (o in left) {
    orphans[[length(orphans) + 1L]] <- data.frame(
      object = o, source_path = NA_character_, source_format = NA_character_,
      source_stage = NA_character_, destination = dest,
      reason = if (o %in% names(orphan_reasons)) orphan_reasons[[o]] else "no source produces this object",
      action = if (o %in% names(ORPHAN_DROP)) "orphan-dropped"
               else if (o %in% names(ORPHAN_REASON)) "orphan-intended"
               else "orphan-unexplained",
      stringsAsFactors = FALSE)
  }
}
if (length(orphans)) src <- rbind(src, do.call(rbind, orphans))

unexplained <- src$object[src$action == "orphan-unexplained"]
if (length(unexplained)) {
  record(rpt, "FAIL", "orphans",
         sprintf("%d packaged object(s) with no source and no reason: %s",
                 length(unexplained), paste(unexplained, collapse = ", ")))
} else {
  record(rpt, "PASS", "orphans",
         sprintf("%d packaged object(s) without a source, each with a recorded reason: %d kept, %d dropped",
                 sum(src$action %in% c("orphan-intended", "orphan-dropped")),
                 sum(src$action == "orphan-intended"), sum(src$action == "orphan-dropped")))
}

dropped <- src$object[src$action == "orphan-dropped"]
if (length(dropped)) {
  record(rpt, "WARN", "orphans:dropped",
         sprintf("%d object(s) withdrawn from the release: %s",
                 length(dropped), paste(dropped, collapse = ", ")))
}

# An ORPHAN_DROP entry only fires for an object still sitting in a package's data/, because
# that is the loop above's input. Once a withdrawal has been promoted the object is gone
# from the checkout, the entry stops matching, and every step that keys on it -- the man
# page removal, doc:withdrawn, the NEWS "Removed data" section -- does nothing at all. That
# is correct, there being nothing left to remove, but it is indistinguishable from a
# withdrawal that never ran, so name the entries that matched nothing rather than leaving
# the difference to silence.
unmatched <- setdiff(names(ORPHAN_DROP), dropped)
if (length(unmatched)) {
  record(rpt, "PASS", "orphans:dropped-already",
         sprintf("%d withdrawal(s) with nothing left to withdraw, already absent from both packages: %s",
                 length(unmatched), paste(unmatched, collapse = ", ")))
}

# ---- Balance ----------------------------------------------------------------
# The counts have to add up, or something was dropped silently.

tab <- table(src$action)
for (a in names(tab)) record(rpt, "PASS", paste0("action:", a), sprintf("%d object(s)", tab[[a]]))

n_sources <- sum(!is.na(src$source_path))
n_accounted <- sum(src$action %in% c("replace", "add", "skip-bucket-only"))
if (n_sources != n_accounted) {
  record(rpt, "FAIL", "balance",
         sprintf("%d source object(s) but %d routed", n_sources, n_accounted))
} else {
  record(rpt, "PASS", "balance", sprintf("all %d source object(s) routed", n_sources))
}

write_tsv(src[order(src$destination, src$object), ], file.path(opt$out, "routing.tsv"))
write_tsv(report_frame(rpt), file.path(opt$out, "route_report.tsv"))

message(sprintf("\n%d check(s): %d FAIL, %d WARN", length(rpt$rows), rpt$fails, rpt$warns))
if (rpt$fails > 0L) quit(status = 1L)
