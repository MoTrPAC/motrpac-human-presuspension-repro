#!/usr/bin/env Rscript
# Tests for the Stage 2 diff: the name/version parsing that decides which staging object a
# freeze file matches, and the classification that decides what step 04 uploads and deletes.
#
# classify_freeze_file() is the one decision no Stage 1 test stands in for — see this
# step's README on CONFLICT — so its truth table and the parsing feeding it are asserted
# here in full.
#
# All of it is pure — no bucket is contacted. The tiers that want real file names read
# step 01's snapshot and SKIP when it is absent, so this runs on a fresh clone.
#
# Run: RUN_TESTS=1 is the default in build.sh; RUN_TESTS=0 skips.

source(file.path(Sys.getenv("STAGE_LIB"), "bucket_helpers.R"))
source(file.path(Sys.getenv("STAGE_LIB"), "required_structure.R"))
source(file.path(.BH_ROOT, "scripts", "10_build_data", "lib", "test_helpers.R"))

cat("\n== version parsing ==\n")

assert("strip_to_version: standard BIC name",
       identical(strip_to_version("proteomics/da/human-precovid-sed-adu_t10-muscle_prot-pr_da_dream-acute_v2.0.txt"), "2.0"))
assert("strip_to_version: .txt.gz",
       identical(strip_to_version("x/y_v1.4.txt.gz"), "1.4"))
assert("strip_to_version: .csv",
       identical(strip_to_version("resources/motrpac_human-precovid_1kg_pca.csv"), NA_character_),
       "an unversioned name must yield NA, not a spurious version")
assert("strip_to_version: unversioned .sqlite",
       identical(strip_to_version("resources/txdb_hsapiens_ensembl_v105.sqlite"), NA_character_),
       "v105 is part of the stem, not a _v<major>.<minor> suffix")

assert("join_key: strips the version, keeps the path",
       identical(join_key("epigenomics/qc-norm/a_v1.2.txt"), "epigenomics/qc-norm/a"))
assert("join_key: leaves an unversioned path alone, extension and all",
       identical(join_key("resources/txdb_hsapiens_ensembl_v105.sqlite"),
                 "resources/txdb_hsapiens_ensembl_v105.sqlite"),
       "both sides of the diff get the same treatment, so an untouched path still matches")
assert("join_key: a moved file does not match its old location",
       !identical(join_key("clinical_chemistry/da/x_v1.3.txt"), join_key("proteomics/da/x_v1.3.txt")),
       "keying on basename alone would hide the clinical split")

# The reason .ver_rank exists at all: a string sort puts v1.10 below v1.9, and the higher
# version is the one every downstream step treats as the incumbent.
assert("ver_rank: 1.10 outranks 1.9",  .ver_rank("1.10") > .ver_rank("1.9"))
assert("ver_rank: 2.0 outranks 1.4",   .ver_rank("2.0")  > .ver_rank("1.4"))
assert("ver_rank: NA ranks lowest",    .ver_rank(NA_character_) < .ver_rank("1.2"))
assert("ver_rank: vectorised",         identical(.ver_rank(c("1.2", "2.0")), c(1002, 2000)))

cat("\n== file-name parsing ==\n")
# parse_bic_name() fills the diff report's descriptive columns and nothing on the
# upload/delete path, so it is covered at the grammar level.

fields <- function(path) {
  p <- parse_bic_name(path)
  c(p$tissue_code, p$ome, p$data_category, p$data_details, p$version)
}
bic <- "metabolomics-targeted/da/human-precovid-sed-adu_t02-plasma_metab-t-amines_da_dream-acute_v1.3.txt"
assert("parse_bic_name: every field of a standard BIC name",
       identical(fields(bic), c("t02-plasma", "metab-t-amines", "da", "dream-acute", "1.3")),
       paste("got:", paste(fields(bic), collapse = " | ")))

# data_details carries embedded underscores in some names and must not be truncated.
p2 <- parse_bic_name("epigenomics/metadata/human-precovid-sed-adu_t03-edta_epigen-methylcap-seq_metadata_removed-samples_v2.0.txt")
assert("parse_bic_name: hyphenated details survive", identical(p2$data_details, "removed-samples"))

p3 <- parse_bic_name("resources/motrpac_human-precovid_1kg_pca.csv")
assert("parse_bic_name: non-BIC name yields NA fields", is.na(p3$ome) && is.na(p3$tissue_code),
       "resources/ files do not follow the BIC grammar and must not be forced into it")

cat("\n== md5 conversion ==\n")
# gsutil reports base64; tools::md5sum returns hex. Everything downstream compares hex.
assert("md5_base64_to_hex: known value",
       identical(md5_base64_to_hex("mNtROaDrXEslTHK8AgqmCw=="), "98db5139a0eb5c4b254c72bc020aa60b"))
assert("md5_base64_to_hex: NA in, NA out", is.na(md5_base64_to_hex(NA_character_)))
assert("md5_base64_to_hex: garbage in, NA out", is.na(md5_base64_to_hex("not base64 !!")))

cat("\n== change classification ==\n")
# The truth table step 04 acts on. MODIFIED and REVERSIONED are the rows that delete
# anything — both moved the version, so both leave an object at the old path behind.
#
# The table is the two comparisons crossed. Version same/moved x bytes same/differ:
#
#            bytes same    bytes differ   md5 unknown
#   same     UNCHANGED     CONFLICT       CONFLICT
#   moved    REVERSIONED   MODIFIED       MODIFIED
#
# The two off-diagonal cells are the mistakes config/file_versions.json can make, and they
# are not symmetric in cost: CONFLICT would publish new bytes under a number already
# meaning something else, so it stops the run; REVERSIONED only churns a version for no
# content change, so it proceeds and warns.
cf <- function(...) classify_freeze_file(...)
assert("ADDED: stem absent from the bucket",
       identical(cf("2.0", NA_character_, "aa", NA_character_), "ADDED"))
assert("MODIFIED: version moved, bytes differ",
       identical(cf("2.0", "1.2", "aa", "bb"), "MODIFIED"))
assert("REVERSIONED: version moved, bytes identical",
       identical(cf("2.0", "1.2", "aa", "aa"), "REVERSIONED"),
       "a version bump with no content change is worth seeing, not filing under MODIFIED")
assert("MODIFIED: version moved, md5 unknown",
       identical(cf("2.0", "1.2", "aa", NA_character_), "MODIFIED"),
       "cannot tell is not the same as identical — a composite object must not become REVERSIONED")
assert("UNCHANGED: same version, same bytes",
       identical(cf("1.2", "1.2", "aa", "aa"), "UNCHANGED"))
assert("CONFLICT: same version, different bytes",
       identical(cf("1.2", "1.2", "aa", "bb"), "CONFLICT"))
assert("CONFLICT: same version, md5 unknown",
       identical(cf("1.2", "1.2", "aa", NA_character_), "CONFLICT"),
       "cannot tell is not the same as unchanged")
assert("REPLACE: CONFLICT under ALLOW_CONFLICT",
       identical(cf("1.2", "1.2", "aa", "bb", allow_conflict = TRUE), "REPLACE"))
assert("ALLOW_CONFLICT does not disturb the other verdicts",
       identical(cf("1.2", "1.2", "aa", "aa", allow_conflict = TRUE), "UNCHANGED") &&
       identical(cf("2.0", "1.2", "aa", "bb", allow_conflict = TRUE), "MODIFIED") &&
       identical(cf("2.0", "1.2", "aa", "aa", allow_conflict = TRUE), "REVERSIONED"),
       "ALLOW_CONFLICT governs the same-version case only")
# A file with no version token (a pre-v2.0 incumbent) has NA on both sides, so it must land in the
# same-version row rather than reading as a moved version against itself.
assert("unversioned: NA == NA is not a moved version",
       identical(cf(NA_character_, NA_character_, "aa", "aa"), "UNCHANGED") &&
       identical(cf(NA_character_, NA_character_, "aa", "bb"), "CONFLICT"))

cat("\n== overwrite record ==\n")
# write_replace_readme() is the ONLY trace an in-place overwrite leaves — the bucket keeps
# none — so what is asserted is that it records enough to identify the object afterwards:
# the path, and BOTH md5s. Losing either would make the file unable to answer the one
# question it exists for, which of the two versions you are holding.
.repl <- data.frame(
  rel_path = "proteomics/qc-norm/x_v1.2.txt",
  gcs_path = "gs://b/proteomics/qc-norm/x_v1.2.txt",
  new_version = "1.2", old_md5 = "aaa111", new_md5 = "bbb222",
  stringsAsFactors = FALSE
)
.tmp <- tempfile(fileext = ".md")
write_replace_readme(.repl, .tmp, applied = TRUE, bucket = "gs://b", verified = TRUE)
.txt <- paste(readLines(.tmp), collapse = "\n")
assert("replace README: records both md5s",
       grepl("aaa111", .txt, fixed = TRUE) && grepl("bbb222", .txt, fixed = TRUE),
       "without the pair you cannot tell which version an object is")
assert("replace README: records the object path and that the version did not move",
       grepl("gs://b/proteomics/qc-norm/x_v1.2.txt", .txt, fixed = TRUE) &&
       grepl("unchanged", .txt, fixed = TRUE))
assert("replace README: says plainly it was a dry run when it was",
       {
         t2 <- tempfile(fileext = ".md")
         write_replace_readme(.repl, t2, applied = FALSE, bucket = "gs://b")
         grepl("DRY RUN", paste(readLines(t2), collapse = "\n"), fixed = TRUE)
       },
       "a planned overwrite must not read as one that happened")
.none <- tempfile(fileext = ".md")
write_replace_readme(.repl[0, ], .none, applied = TRUE, bucket = "gs://b")
assert("replace README: no REPLACE rows writes no file",
       !file.exists(.none),
       "a stale record from a previous cycle must not be mistaken for this one")
unlink(c(.tmp, .none))

cat("\n== manifest ==\n")

# The manifest agreeing with the OME_TISSUE_CODE catalog is Stage 1's to catch, against
# real files: qc_norm_freeze_tests.R covers every row's qc-norm, features and samples
# files under that ome's transform, da_tests.R every row's DA file. Asserted here: the
# `required` flag, which decides what this stage FAILs on when a file is absent.
#
# The exception is DA coverage. Upstream's da_details lookup named only the five non-metab
# omes, so every metab pair fell through to NA and got no DA row at all — 34 of the 48 DA
# files in the bucket went unvalidated. Stated as "every ome carries a DA row" rather than
# as a count, it needs no total to keep current and holds on a fresh clone, where the tier
# that compares the manifest to the freeze has no freeze to read.
assert("manifest: every ome carries a DA row",
       setequal(unique(required_structure$ome),
                unique(required_structure$ome[required_structure$data_category == "da"])),
       paste("no DA row for:",
             paste(setdiff(unique(required_structure$ome),
                           required_structure$ome[required_structure$data_category == "da"]),
                   collapse = ", ")))
assert("manifest: imputed rows are proteomics-only and required",
       all(required_structure$ome[required_structure$data_category == "imputed"] %in% c("prot-pr", "prot-ph")) &&
       all(required_structure$required[required_structure$data_category == "imputed"]),
       "the row is only generated for omes that produce one, so it must be satisfiable")
assert("manifest: removed-samples is optional",
       !any(required_structure$required[required_structure$data_details == "removed-samples"]),
       "not every ome x tissue removed a sample")
assert("manifest: the clinical assays are ordinary omes",
       all(c("prot-clinical", "metab-t-clinical") %in% required_structure$ome) &&
       !any(required_structure$gcs_subdir == "clinical_chemistry"),
       "clinical_chemistry/ is retired")

cat("\n== carry-forward list ==\n")
assert("carried: qc-report HTML",
       is_carried_forward("transcriptomics/metadata/human-precovid-sed-adu_all_transcript-rna-seq_qc-report_v1.2.html"))
assert("carried: 1kg PCA, unversioned and versioned",
       is_carried_forward("resources/motrpac_human-precovid_1kg_pca.csv") &&
       is_carried_forward("resources/motrpac_human-precovid_1kg_pca_v2.0.csv"),
       "config/file_versions.json ships it at v2.0; the pre-v2.0 name must still match")
assert("carried: an ordinary data file is NOT carried",
       !is_carried_forward("proteomics/da/human-precovid-sed-adu_t10-muscle_prot-pr_da_dream-acute_v2.0.txt"),
       "a regenerated file must never be excused as carry-forward")

# ---- Against the real bucket listing, when one has been taken ----------------------

snap <- file.path(SNAPSHOT_DIR, "latest_staging.tsv")
cat("\n== real bucket names ==\n")
if (!file.exists(snap)) {
  skip("bucket name round-trip", "no staging snapshot — run step 01")
} else {
  st <- read.csv(snap, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
                 colClasses = "character")
  vers <- vapply(st$rel_path, strip_to_version, character(1), USE.NAMES = FALSE)
  keys <- join_key(st$rel_path)

  assert("every bucket path parses to a key shorter than itself, or is unversioned",
         all(nchar(keys) < nchar(st$rel_path) | is.na(vers)),
         paste("unparsed:", paste(utils::head(st$rel_path[nchar(keys) == nchar(st$rel_path) & !is.na(vers)], 3), collapse = ", ")))

  # A duplicated key means the bucket carries two versions of one file. That is a real
  # state the diff reports as DUPLICATE, so this is INFO rather than an assertion.
  dup <- unique(keys[duplicated(keys)])
  if (length(dup)) report("bucket holds multiple versions of a stem",
                          sprintf("%d stem(s), e.g. %s", length(dup), dup[1]))
  else report("bucket holds one version per stem", sprintf("%d files", nrow(st)))

  unversioned <- st$rel_path[is.na(vers)]
  report("unversioned bucket files",
         if (length(unversioned)) paste(unversioned, collapse = ", ") else "none")
}

finish()
