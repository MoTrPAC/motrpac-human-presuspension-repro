# Step 09 — build_da

Fits the acute differential-analysis models and writes one BIC-named DA table per tissue ×
ome under `staging/freeze/<ome-group>/da/`. `build.sh` runs each `generate_<ome>_da.R` stem in
turn; `STEMS="a b"` restricts the run.

## What is in this folder

| File | Kind | Notes |
|------|------|-------|
| `build.sh` | driver | stem loop; per-stem logs + `logs/da_report.tsv` |
| `da_common.R` | shared — the only helper the stems source | the vendored engine (`run_dream()`, `.convert_dream_output()`, `.generate_contrasts_acute()`, `.generate_sex_contrasts()`, `process_covariates()`), plus `run_da_models()`, `DA_MODEL_TYPE` and `da_parallel_enabled()`. The freeze version is not a constant here — see `config/file_versions.json` |
| `filter_paired_n.R` | **vendored, one line changed** | paired-sample-count feature filter (prot-pr / prot-ph only) |
| `generate_*_da.R` | adapted stems | one per ome, each carrying its upstream body commented out above the adaptation |
| `sources/` | **raw inputs** | vendored external DA tables — see `sources/README.md` |

Vendored **code** sits alongside the stems; vendored **raw input files** go under `sources/`,
the same convention as `06_generate_qc_norm/sources/`. `sources/` holds `methylcap_da/` (the
MALAX GLMM tables described under Scope below) and, only when `RERUN_ATAC=FALSE`, `atac_da/`.

## Vendored code — provenance

Copied 2026-07-29 from
`MotrpacHumanPreSuspensionAnalysis/data-raw/generate_differential_analysis/`. To refresh,
re-copy from that directory and re-apply the two notes below.

Every change is marked in place with the upstream original kept commented directly above it.
Most are one line, but two are not, and both move results rather than plumbing:
`da_common.R` adds `.ebayes_legacy()`, which selects limma's empirical-Bayes estimator and so
moves `df.prior`/`s2.prior` and every feature's moderated variance, t and p; and
`filter_paired_n.R` reimplements the cell grid in about 25 lines, because upstream's
`pivot_wider`/`c_across` minimum was non-monotonic.

The vendored engine in `da_common.R` (419 lines upstream) — its `process_covariates()`
read `MotrpacHumanPreSuspensionAnalysis::COVARIATES_FILE`, the installed package object, which
is the dependency this stage exists to remove. It now reads the pipeline's local
`COVARIATES_FILE` global.

`filter_paired_n.R` (122 lines upstream) — its final statement filtered `feature_metadata` on a
column named `id`, which exists neither on this pipeline's `*_QC` objects nor on the installed
package's; both carry `feature_id`, which is also the name the function uses everywhere else.
The `qc_norm` subsetting on the preceding line is keyed on rownames and was always correct.

## `BiocParallel` must be attached

The vendored `run_dream()` calls `SnowParam()` unqualified, so `BiocParallel` is attached in
`da_common.R` rather than merely installed. `da_common.R` then runs `MulticoreParam` instead —
see its header for why SOCK loses a feature per worker — and takes its core count from
`VARIANCEPARTITION_PARALLEL_CORES`, not `PARALLEL_CORES`, which drives mice and does nothing
here. Either way a serial run will not catch the attach.

## Scope

Adapted: prot-ol, prot-pr, prot-ph, metabolomics (all platforms, per tissue), transcriptomics,
ATAC. Transcriptomics and ATAC are fit on raw counts with voom weighting rather than on the
qc-norm matrix.

Not modelled here: **methylcap**, which is fit with a MALAX GLMM in Yongchao Ge's separate
pipeline that is not part of this package at all. Rather than leave the ome absent from the
freeze, `generate_methylcap_da.R` copies the externally produced tables in verbatim from
`sources/methylcap_da/` — the same treatment `generate_methylcap_qc_norm.R` gives the
beta-value matrices one tier up. The copies keep their upstream
`malax-glmm-acute_v1.2` filenames rather than this pipeline's `dream-acute` token, because the
model is not dream; the bucket-diff test keys on the version-stripped name, so they still pair.
Their presence is enforced by `../check_required_inputs.sh`.

The ATAC stem reads its qc-norm and sample metadata directly from `staging/freeze` instead of
through `load_qc_local()`, because step 08 does not package epigen into `*_QC.rda` and
`load_qc_local()` has no `epigen` argument. See the header of `generate_atac_da.R`.
