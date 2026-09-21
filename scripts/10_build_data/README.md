# Stage 1 — build data (step by step)

Regenerates the MoTrPAC Pre-Suspension data objects as a sequence of numbered step
subfolders. **The folder sequence is the build order** — `10_build_data.sh` runs each
`NN_*/build.sh` in order. `docs/data_objects.tsv` is now the object inventory /
regeneration ledger (it has a `step` column pointing back here); it no longer drives the
build.

Each subfolder exposes one `build.sh` entry point that sources `scripts/lib/common.sh`
and does one step. Its exit code is the step status: **0 = PASS**, **77 = SKIP**
(placeholder / not-yet-adapted), anything else = FAIL.

## Steps

| Step | Builds | Kind |
|------|--------|------|
| `01_stage_clinical` | `cln_*` (76) | staged from Data pkg (TODO: local build) |
| `02_stage_pheno` | `pheno` | staged from Data pkg (TODO: local build) |
| `03_build_molecular_signatures` | `MOLECULAR_SIGNATURES` | adapted local build (reads preflight GMTs) |
| `04_build_set_to_id` | `SET_TO_ID` | adapted local build |
| `05_stage_metabolomics_cvs` | `METABOLOMICS_CVS` | staged from Analysis pkg (TODO: local build) |
| `06_generate_qc_norm` | qc-norm **freeze** (all omes) → `staging/freeze/` | adapted; 10 stems, looped in its own `build.sh` |
| `07_build_human_feature_to_gene` | `HUMAN_FEATURE_TO_GENE` | adapted (reads freeze `metadata_features`) |
| `08_build_qc_objects` | the 42 `*_QC` + `blood_transcript_1/2` | adapted (packages freeze into `.rda`) |
| `09_build_da` … `12_build_camera_results` | DA, DA_ASSEMBLE, SUM_STATS, CAMERA_RESULTS | adapted (09 DA chain still in migration) |
| `13_build_fcm` | `FCM_CLUSTERS`, `FCM_CAMERA`, `FCM_ORA` | adapted local build + cluster-count diagnostics |
| `14_build_utoronto_tfs` | `UTORONTO_TFs` | adapted local build (vendored Human TFs extract → prot-ph sites) |
| `15_build_ptmsea_input` | `PTMSEA_INPUT_<tissue>.gct` (+ `PTMSEA_INPUT`) | adapted local build; **the file is an input, the PTM-SEA run is out of scope** |
| `16_build_scion` | SCION edge tables → `staging/scion/` | adapted local build; **not a data object**, on by default, `RUN_SCION=FALSE` opts out |

Step 06 produces the qc-norm **freeze** files (not package `.rda`); steps 07 and 08
consume that freeze, so 06 must run before them. (Stage 0 / `00_preflight` builds the leaf
objects — colors, GMTs, `OME_TISSUE_CODE`, `OUTLIERS`, `COVARIATES_FILE`, `SPLICING_DA` —
before this stage.)

Step 15 is the one step whose output is not consumed anywhere else in this repo. It writes a
`.gct` of prot-ph DA z-scores for **PTM-SEA**, which is run from the **Broad Institute's
PTM-SEA / ssGSEA2.0 Docker image** against PTMsigDB. **That run is not part of this
pipeline**: nothing here pulls or runs the container, no stage reads its enrichment scores,
and there is no reproducibility claim over the container version, the PTMsigDB release, or
its parameters. Stage 1's contract ends at a valid `.gct` on disk. See
`15_build_ptmsea_input/README.md` for the handoff.

Step 16 is the other output that is not a package object: it writes SCION network edge
tables under `staging/scion/`, which Figure 6B is drawn from after a manual Cytoscape merge.
It runs by default (`RUN_SCION=FALSE` opts out) and is measured in hours to days, so submit
it (`EXECUTOR=slurm`) rather than running it inline. See `16_build_scion/README.md`.

## Run

```bash
make data                      # run all steps in order (implies preflight)
STEPS="03 04" make data        # run only steps whose folder name matches a token
STOP_ON_FAIL=false make data   # keep going past failures
```

Per-step logs: `logs/data_<step>.log`. Summary: `logs/build_data_report.tsv`
(`status` ∈ PASS/SKIP/FAIL). Step 06 additionally writes per-stem logs
`logs/qcnorm_<stem>.log` and `logs/qc_norm_report.tsv`.

## Shared resources (stage level)

- `lib/qc_helpers.R` — vendored internals (BIC naming, `load_pheno`, `process_covariates`,
  offline Ensembl v105 accessors, `run_mice`, …) sourced by the step-06 stems and
  step-08 `qc_norm_results.R`. Reads `OME_TISSUE_CODE`/`OUTLIERS`/`COVARIATES_FILE`/
  `QUANTID_BUCKET_FILES` from preflight and `pheno` from `data/`.
- `data/` — staging output; every step reads/writes its `.rda` here. These will eventually
  be organized as the canonical objects in the `MotrpacHumanPreSuspensionData` /
  `-Analysis` packages.
- `check_required_inputs.sh` — the pre-build gate (preflight objects, vendored methylcap
  sources under `06_generate_qc_norm/sources/`, Ensembl v105 cache).
- `06_generate_qc_norm/sources/methylcap_qc_norm/` — vendored MethylCap beta-value inputs
  (only the methylcap stem consumes them).

## Still to do

- Turn the staged steps (01, 02, 05) into real local builds.
- Adapt the placeholder steps (09–14) — port the corresponding
  `MotrpacHumanPreSuspensionAnalysis/data-raw/` generators to local
  `staging/freeze`-reading stems (see each stub's `build.sh` for the upstream source).
## `confident_site` (prot-ph)

Step 06's `.annotate_prot_ph()` used to `select()` the PTM site-localization flag away before
writing `metadata_features`, so it reached neither the freeze nor any `.rda`. The shipped
Data-package object has the same gap, which is why upstream's `preprocess_PTMSEA()` returns
zero rows against it. It is now kept in step 06 and carried into `*_PROT_PH_QC` (step 08), and is what step 15
filters PTM-SEA input on.

It is **not** carried into `HUMAN_FEATURE_TO_GENE` or its `resources` `.txt` (step 07). The
flag is **per tissue** and that map is keyed on `(assay, feature_id)` with no tissue column,
so a value there could only be a collapse across tissues rather than the measurement — 859 of
the 7,865 sites shared by muscle and adipose disagree. **Read
`*_PROT_PH_QC$feature_metadata` for it.** Full reasoning in
`15_build_ptmsea_input/README.md`.

## `custom_annotation` / `relationship_to_gene` (epigen)

Step 06's `.annotate_atac_features()` and `.annotate_methylcap_features()` write both peak
annotation columns onto the epigen `metadata_features`, but step 07's `select()` dropped them
before they reached `HUMAN_FEATURE_TO_GENE` — even though the object's own documentation
(`carry/04_document.R`) describes both. The map therefore carried an ATAC or MethylCap feature
to a gene without saying what the relationship to that gene *is*: a promoter peak and a peak
40kb into an intron mapped to the same gene, indistinguishable. Both are now retained.

- `custom_annotation` — factor; the ChIPseeker region class, renamed by
  `pre_cawg_get_peak_annotations_hs()` (`lib/qc_helpers.R`): `Promoter (<=1kb)`,
  `Promoter (1-2kb)`, `5' UTR`, `3' UTR`, `Exon`, `Intron`, `Overlaps Gene`,
  `Upstream (<5kb)`, `Downstream (<5kb)`, `Distal Intergenic`.
- `relationship_to_gene` — numeric; the signed base-pair distance to that gene, `0` where the
  peak overlaps it. Cast to numeric in step 07: it is read as character and is in the blanket
  `as.factor()` exclusion list, so without the cast every comparison it exists for would
  compare lexicographically (`"40000" < "5000"` is `TRUE` as a string).

Both are `NA` for every non-epigen assay. **They need no cross-tissue collapse** — they are
derived from the peak *coordinates* (the `feature_id` itself), not
measured per tissue, so they agree on every `feature_id` shared between tissues (0 of 306,788
ATAC and 0 of 1,545,930 MethylCap shared keys disagree). Step 07 has a duplicate-key guard
that fails the build rather than trusting that, and the step-07 tests pin it per assay.
