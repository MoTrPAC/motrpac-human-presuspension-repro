# Step 14 — `UTORONTO_TFs`

The prot-ph phosphosites whose gene is a curated human transcription factor. Two columns:

| column | |
|---|---|
| `feature_id` | prot-ph site id, e.g. `O00409_S352s` |
| `gene_symbol` | the site's gene, always a curated TF |

Built from the vendored Human TFs database extract (`sources/DatabaseExtract_v_1.01.txt`,
Lambert et al. 2018) restricted to `Is TF? == "Yes"`, mapped through step 07's
`HUMAN_FEATURE_TO_GENE` to every prot-ph feature sharing the gene symbol.

**1,381 sites across 383 TF genes**, from 1,639 curated TFs — most TFs have no confidently
quantified phosphosite in this study.

```bash
bash scripts/10_build_data/14_build_utoronto_tfs/build.sh    # or: STEPS=14 bash scripts/10_build_data/10_build_data.sh
```

Requires step 07. Output: `scripts/10_build_data/data/UTORONTO_TFs.rda`.

## Why the column set is the whole contract

This object is a **regulator pool**, not an annotation table. Step 16's
`.load_scion_matrixes()` reads it to apply `subset_TFs = TRUE`, where upstream read

```r
private_subset_qc(desired_qc_norm_single,
                  desired_features = MotrpacHumanPreSuspensionAnalysis::UTORONTO_TFs$feature_id)
```

and `subset_qc()` guards its feature filter with
`if (!is.null(desired_features))`. So an object **without** a `feature_id` column does not
error — it silently turns the TF subsetting into a no-op and hands SCION every phosphosite
in the matrix.

That is the state the object shipped by `MotrpacHumanPreSuspensionAnalysis` (through
v2.0.0) is in: it is the raw 2,765 x 28 download, not the two-column mapping its own
`data-raw/UTORONTO_TFs.R` describes. Against it, a `subset_TFs = TRUE` SCION run models on
all ~9,300 imputed muscle phosphosites rather than the ~230 TF sites that survive into that
matrix — a 40x regulator pool with no error and no warning. `utoronto_tests.R` asserts the
column set for exactly this reason.

## Adaptations from upstream

`MotrpacHumanPreSuspensionAnalysis/data-raw/UTORONTO_TFs.R`, with four changes:

1. **Reads `sources/`, not `~/Downloads/`** — so the build is offline and reproducible.
2. **`assay`, not `platform`.** `HUMAN_FEATURE_TO_GENE` has no `platform` column, so
   upstream's `filter(platform == "prot-ph")` errors against it. That is why the upstream
   script cannot regenerate its own object as written, and why the shipped `.rda` is the raw
   download.
3. **Filter before dedup.** Upstream runs `distinct(feature_id, .keep_all = TRUE)` before
   `filter(platform == ...)` / `filter(Is.TF. == "Yes")`. A feature that maps to two genes
   gets a row per gene, so deduping first keeps whichever row sorts first and can drop a
   real prot-ph TF site whose surviving row is a non-TF gene or a different assay.
   Defensive rather than corrective: on the current inputs both orders give the same 1,381
   sites.
4. **`save()`**, not `usethis::use_data()`, into `scripts/10_build_data/data/`.

## Tests

`utoronto_tests.R` re-derives the expected feature set straight from the vendored extract
and `HUMAN_FEATURE_TO_GENE` rather than importing it from the builder, so a filtering or
dedup mistake fails instead of being mirrored. It also asserts every `(feature_id,
gene_symbol)` pair is real, every id is a prot-ph feature, and the object is not the raw
download shape.

`diff_vs_package` is expected to report a difference — the shipped object is the un-mapped
raw extract — so it is recorded as INFO, not a failure.
