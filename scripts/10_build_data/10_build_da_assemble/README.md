# Step 10 — build_da_assemble

Turns the step-09 DA freeze into the package-facing differential-analysis objects. Step 09
**fits**, one BIC-named table per tissue × ome under `staging/freeze/<ome-group>/da/`; this
step **assembles**, reading those tables back, attaching the contrast metadata, and reshaping
them into `CONTRAST_CONVERTER` plus one `{TISSUE}_{OME}_DA` object per fitted table.

Nothing is refit. No statistic changes value here — only which object a row lands in, what
the columns are called, and what order they sit in. `da_assemble_tests.R` asserts exactly
that: every statistic in every object is compared cell-for-cell against the freeze table it
came from.

Output: `scripts/10_build_data/data/*.rda` (19 objects). Requires step 09.

## What is in this folder

| File | Kind | Notes |
|------|------|-------|
| `build.sh` | driver | runs the builder, then the tests unless `RUN_TESTS=0` |
| `differential_analysis_results.R` | adapted | the build; upstream's `data-raw/differential_analysis_results.R` with its loader replaced |
| `da_assemble_tests.R` | tests | structure, conservation, value fidelity, diff-vs-package |
| `sources/contrast_converter.txt` | **raw input** | the 33 contrasts, authored upstream — see `sources/README.md` |

## What it builds

`CONTRAST_CONVERTER` — the 33 acute contrasts and everything derived from them
(`contrast_type`, `contrast_category`, `randomGroupCode`, `Timepoint`). Built from the
vendored `sources/contrast_converter.txt`, whose **row order is significant**: it sets
`contrast_order` and the factor level order every `*_DA` object inherits.

18 `*_DA` objects, one per tissue × ome:

| tissue | objects |
|---|---|
| adipose | `METAB` (12 platforms), `PROT_PH`, `PROT_PR`, `TRNSCRPT`, `EPIGEN_METHYLCAP_SEQ` |
| blood | `METAB` (10 platforms), `METAB_T_CLINICAL`, `PROT_CLINICAL`, `PROT_OL`, `TRNSCRPT`, `EPIGEN_ATAC_SEQ`, `EPIGEN_METHYLCAP_SEQ` |
| muscle | `METAB` (10 platforms), `PROT_PH`, `PROT_PR`, `TRNSCRPT`, `EPIGEN_ATAC_SEQ`, `EPIGEN_METHYLCAP_SEQ` |

Metabolomics platforms are stacked per tissue into one `METAB` object carrying `assay =
"metab"` and the platform in its own `platform` column — upstream's `bind_rows(.id =
"platform")`. Every other ome is one object per fitted table. `transcript-rna-seq` is the one
ome whose object name is not its assay code spelled out: upstream abbreviates it `TRNSCRPT`
and downstream code uses that name.

Two BIC tissue codes map to muscle (`t06-muscle`, `t10-muscle`) and four to blood
(`t02-plasma`, `t03-edta`, `t04-blood-rna`, `t05-pbmc`), so tissue is resolved through
`OME_TISSUE_CODE` rather than read off the filename.

## Clinical chemistry is split — deliberate divergence from the package

This is the one place the local objects differ from the released package by design.

The nine clinical analytes are fit as **two separate tables** in the freeze:

- `prot-clinical` → `BLOOD_PROT_CLINICAL_DA` — CK, Glucagon, Insulin
- `metab-t-clinical` → `BLOOD_METAB_T_CLINICAL_DA` — Cortisol, Glucose, Glycerol, KET,
  Lactate, NEFA

The released package concatenates them into a single `CLIN_CHEMISTRY_DA` under the synthetic
assay name `"clinical-chemistry"` (9 features × 33 contrasts = 297 rows). They are kept apart
here, which is what `docs/data_objects.tsv` lists, so that **every object maps back to exactly
one fitted freeze table** (or, for metab, one stack of platform files) and keeps the assay
name it was actually fit under. The tests assert the split is clean — disjoint features,
distinct assays — and report the union against the package object (99 + 198 = 297 rows,
9 features) as the check that splitting lost nothing.

## Epigen

Assembled here too, unlike upstream's `epigen = FALSE`, one object per freeze table:

| object | tissue code | model | rows |
|---|---|---|---|
| `BLOOD_EPIGEN_ATAC_SEQ_DA` | `t05-pbmc` | `dream-acute` | 5 272 134 |
| `MUSCLE_EPIGEN_ATAC_SEQ_DA` | `t06-muscle` | `dream-acute` | 5 361 699 |
| `BLOOD_EPIGEN_METHYLCAP_SEQ_DA` | `t03-edta` | `malax-glmm-acute` | 2 779 517 |
| `MUSCLE_EPIGEN_METHYLCAP_SEQ_DA` | `t06-muscle` | `malax-glmm-acute` | 4 373 304 |
| `ADIPOSE_EPIGEN_METHYLCAP_SEQ_DA` | `t11-adipose` | `malax-glmm-acute` | 2 637 529 |

Getting these to AWS, where `load_differential_analysis(epigen = TRUE)` fetches from, is
handled separately.

This is where the methylcap schema warning in `../09_build_da/sources/README.md` — *"whoever
ports it must reconcile the two schemas rather than assume every ome has a `logFC`"* — is
exercised. ATAC is dream-fit and carries the standard statistics; methylcap carries
`methylation_diff` in place of `logFC` and none of dream's `z.std` /
`degrees_of_freedom` / `logLik`. The builder orders and keys columns through `intersect()` and
never names a statistic, so both reshape without a branch. The tests assert each family's
columns separately rather than trusting that.

Contrast counts differ per tissue and ome (ATAC 21, methylcap 6 blood / 4 muscle / 2 adipose),
and the methylcap feature × contrast grid is incomplete — both properties of an external
artifact step 09 passes through, which is why step 09 reports that grid rather than asserting
it.

**Two consumers glob `_DA\.rda$` over this output dir**, and both had to be told about epigen:

- step 11 reads it from these objects instead of re-reading the freeze, which it used to do
  while epigen was unassembled — doing both would count every epigen feature twice
- step 12 drops it at the file glob rather than letting `run_cameraPR()`'s `selected_omes`
  drop it, which would mean loading 20.4M rows into `DA_list` first and then discarding them

## Platforms the CV dedup empties

**`metab-t-imm-crt`** — assembled like every other platform, and left with no rows. Its
single feature is Cortisol, the same analyte `metab-t-clinical` carries, measured by
immunoassay instead; blood `metab-u-hilicpos` measures it at a lower CV, so `METABOLOMICS_CVS`
flags no `metab-t-imm-crt` row `lowest_CV == "yes"` and nothing survives the collapse. The
same dedup that reduces every other duplicated analyte to one copy is what keeps a second
Cortisol row out of `BLOOD_METAB_DA` — the outcome upstream reaches by excluding the platform
outright, since it has no room for it in either `BLOOD_METAB_DA` (10 platforms) or
`CLIN_CHEMISTRY_DA` (9 features).

`collapse_metab_cv()` reports an emptied platform on stderr rather than erroring, and the
tests report the same set back from `METABOLOMICS_CVS`, so the gap stays visible without any
ome being named in the code.

## Adaptations from upstream

`data-raw/differential_analysis_results.R` is ported with its loader replaced; everything
downstream of the load — the metabolomics stacking, `.process_raw_DA()`, the naming — is
ported as-is.

1. **Input.** Upstream calls `:::.load_differential_analysis()`, which `gsutil`-copies the
   released DA repo into a tempdir. The freeze *is* that repo's content, already local, so
   the download stem is dropped.
2. **Tissue.** The freeze carries `assay` but not `tissue` (`.convert_dream_output()` sets a
   tissue column and drops it again in its closing `select()`). Upstream recovers it from the
   `"<tissue>.<ome>"` list names its loader builds; here it comes from the BIC tissue code in
   each filename via `OME_TISSUE_CODE`, the same way steps 08 and 11 read the freeze.
3. **`OME_TISSUE_CODE` is read straight from the preflight `.rds`**, not via
   `lib/qc_helpers.R`. Sourcing `qc_helpers` loads `pheno.rda` and the QC objects, none of
   which this step needs — it never touches sample-level data. Reading only the lookup it
   does use keeps the dependency honest.
4. **Output.** `usethis::use_data()` is replaced by a plain `save()` into
   `scripts/10_build_data/data/`, matching every other adapted step. `use_data()` writes into
   an installed package source tree; this pipeline stages objects and promotes them
   separately.
5. **Reader.** `data.table::fread()` rather than `read.csv()`. This step keeps every column,
   so no subset is possible, and the non-epigen DA tier is ~890 MB across 43 files — the
   blood transcriptomics table alone is 281 MB. That is the "several very large matrices at
   once" exception in the workspace `CLAUDE.md`.

## Comparison against the released package

`da_assemble_tests.R` reports a diff against each shipped object. It is **INFO, never FAIL**,
because step 09 refits every dream ome locally, so freeze content legitimately moves ahead of
the release.

As built against freeze `v2.0`:

- `CONTRAST_CONVERTER` is `all.equal`-identical to the package object.
- Column names, column order, `data.table` key, class and `contrast` factor levels are
  identical for every object.
- Row and feature counts match exactly for the proteomics and transcriptomics objects (e.g.
  `BLOOD_PROT_OL_DA` 46 761 rows / 1 417 features).
- The three `METAB` objects carry more features than the release (blood 1 528 vs 1 140), and
  `MUSCLE_TRNSCRPT_DA` seven fewer. These are step-09 / freeze-content differences, not
  assembly differences — the value-fidelity test confirms every statistic present is
  bit-identical to its freeze cell.
