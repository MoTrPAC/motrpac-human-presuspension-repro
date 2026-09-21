#!/usr/bin/env Rscript
# Tests for the Stage 2 validator: that the file name the manifest composes is the file
# name Stage 1 actually wrote, and that every manifest row has a column schema to be
# checked against.
#
# The composition and lookup tiers are pure. The tier that matches composed prefixes
# against real paths reads the freeze and SKIPs when there is none.
#
# Run: RUN_TESTS=1 is the default in build.sh; RUN_TESTS=0 skips.

source(file.path(Sys.getenv("STAGE_LIB"), "bucket_helpers.R"))
source(file.path(Sys.getenv("STAGE_LIB"), "required_structure.R"))
source(file.path(Sys.getenv("STAGE_LIB"), "expected_columns.R"))
source(file.path(.BH_ROOT, "scripts", "10_build_data", "lib", "test_helpers.R"))

# Same two helpers the step defines. Kept in step order so a change there shows up here.
.category_dir <- function(data_category) {
  if (data_category %in% c("qc-norm", "imputed")) "qc-norm" else data_category
}
.expected_prefix <- function(row) {
  paste0(row$gcs_subdir, "/", .category_dir(row$data_category), "/",
         FILE_HEADER, "_", row$tissue_code, "_", row$ome, "_",
         row$data_category, "_", row$data_details, "_v")
}

cat("\n== name composition ==\n")

assert("imputed files are looked for in qc-norm/",
       identical(.category_dir("imputed"), "qc-norm"),
       "there is no imputed/ directory in the bucket")
assert("other categories map to their own dir",
       identical(.category_dir("da"), "da") && identical(.category_dir("metadata"), "metadata"))

prefixes <- vapply(seq_len(nrow(required_structure)),
                   function(i) .expected_prefix(required_structure[i, ]), character(1))

assert("no composed prefix contains NA",
       !any(grepl("NA", prefixes, fixed = TRUE)),
       paste("offending:", paste(utils::head(prefixes[grepl("NA", prefixes, fixed = TRUE)], 3), collapse = ", ")))
assert("composed prefixes are unique",
       !any(duplicated(prefixes)),
       paste("duplicated:", paste(utils::head(prefixes[duplicated(prefixes)], 3), collapse = ", ")))

# startsWith() matching means one prefix being a prefix of another would let a row claim
# another row's file. Nothing in the current grammar does this, and it must stay that way.
sorted <- sort(prefixes)
collisions <- sorted[-1][startsWith(sorted[-1], utils::head(sorted, -1))]
assert("no prefix is a prefix of another",
       length(collisions) == 0,
       paste("collides:", paste(utils::head(collisions, 3), collapse = ", ")))

cat("\n== column schemas ==\n")

# A manifest row with no schema is validated for presence only. That is legitimate for
# removed-samples, and a silent gap for anything else.
no_schema <- vapply(seq_len(nrow(required_structure)), function(i) {
  r <- required_structure[i, ]
  is.null(.get_expected_cols(r$data_category, r$ome, r$data_details))
}, logical(1))
unschemad <- unique(paste(required_structure$data_category, required_structure$ome)[no_schema])
assert("every manifest row resolves to a column schema",
       length(unschemad) == 0,
       paste("no schema for:", paste(unschemad, collapse = "; ")))

# Which DA schema an ome resolves to is Stage 1's da_tests.R to catch: it checks every DA
# file's columns by exact vector, prot-clinical and metab-t-clinical included, so a
# routing mistake shows up there against the real table. The two below have no Stage 1
# counterpart — Stage 1 does not test imputed, and never sees the bucket-resident table
# the forbidden columns guard.
assert("imputed matrices take the qc-norm shape",
       identical(.get_expected_cols("imputed", "prot-pr")$required, QC_NORM_COLS_EXPR))
assert("CI.L/CI.R are forbidden on every DA schema",
       all(vapply(grep("^da__", names(EXPECTED_COLUMNS), value = TRUE),
                  function(k) all(DA_COLS_DROPPED %in% EXPECTED_COLUMNS[[k]]$forbidden), logical(1))),
       "dropped for v2.0 — an older table must fail rather than ship")
# Every removed-samples manifest row, metab platforms included (they resolve through the
# ^metab fallback), must reach the identifier list. validate_structure.R checks forbidden
# columns on any file present, so a row that resolves elsewhere would go unchecked.
rs_rows <- required_structure[required_structure$data_details == "removed-samples", , drop = FALSE]
rs_unguarded <- unique(rs_rows$ome[!vapply(seq_len(nrow(rs_rows)), function(i)
  all(REMOVED_SAMPLES_COLS_IDENTIFYING %in%
      .get_expected_cols("metadata", rs_rows$ome[i], "removed-samples")$forbidden), logical(1))])
assert("participant identifiers are forbidden on every removed-samples row",
       nrow(rs_rows) > 0 && length(rs_unguarded) == 0,
       paste("unguarded omes:", paste(rs_unguarded, collapse = ", ")))

cat("\n== QC report scope ==\n")
assert("QC-report check covers only subdirs that ship one",
       all(QC_REPORT_SUBDIRS %in% required_structure$gcs_subdir) &&
       !any(grepl("^metabolomics", QC_REPORT_SUBDIRS)),
       "no metabolomics QC report exists anywhere; checking for one would warn forever")

# ---- Against the freeze Stage 1 wrote ----------------------------------------------
# The manifest and the freeze must describe the same set of files. Asserted directly,
# against real files, this catches an ome whose DA rows were never generated, a transform
# that composes the wrong name, and a lab-* leak.
#
# Against the freeze, not the bucket: the only snapshot available here is step 01's, taken
# before step 04 uploads, so every genuinely new file would fail a bucket version of this.
# The bucket is step 05's own job, against a listing it takes fresh. Duplicate claims need
# no test — a path starts with two prefixes only if one is a prefix of the other, which
# the pure test above already forbids.

cat("\n== composed names vs the freeze ==\n")
# freeze_inventory() stops on a missing dir and on a dir holding no upload artifacts, and
# those are different problems. Report the message it raised rather than a fixed reason,
# so a real read failure is not filed under "run make data".
fz <- tryCatch(freeze_inventory(), error = function(e) conditionMessage(e))
if (!is.data.frame(fz)) {
  skip("prefix match against freeze", fz)
} else {
  hits  <- lapply(prefixes, function(p) fz$rel_path[startsWith(fz$rel_path, p)])
  n_hit <- lengths(hits)

  # Every REQUIRED row must resolve. Optional rows legitimately match nothing —
  # removed-samples is only written where a sample was actually removed.
  req_misses <- prefixes[required_structure$required & n_hit == 0]
  assert("every required manifest row matches a freeze file",
         length(req_misses) == 0,
         paste(length(req_misses), "unmatched, e.g.", paste(utils::head(req_misses, 3), collapse = ", ")))

  res_hits <- unlist(lapply(seq_len(nrow(REQUIRED_RESOURCES)), function(i) {
    r <- REQUIRED_RESOURCES[i, ]
    fz$rel_path[startsWith(fz$rel_path, paste0(r$rel_stem, if (r$versioned) "_v" else "."))]
  }))
  # Carry-forward files are bucket residents Stage 1 never regenerates, so a freeze that
  # happens to hold one is not an unclaimed file.
  unclaimed <- setdiff(fz$rel_path, c(unique(unlist(hits)), res_hits))
  unclaimed <- unclaimed[!is_carried_forward(unclaimed)]
  assert("every freeze file is claimed by a manifest row",
         length(unclaimed) == 0,
         paste(length(unclaimed), "unclaimed, e.g.", paste(utils::head(unclaimed, 3), collapse = ", ")))

  report("optional rows matching nothing",
         sprintf("%d (removed-samples not written for every ome x tissue)",
                 sum(!required_structure$required & n_hit == 0)))
}

finish()
