# Local RefMet / KEGG annotation snapshot

Local, offline snapshot of the **Metabolomics Workbench RefMet** batch endpoint and
the **KEGG compound list**, used by the Stage 1 metabolomics QC-norm stem to
standardize metabolite names and attach `refmet_name`, `refmet_id` and `kegg_id`.

Both were live calls inside `.build_metab_refmet_map()`: a POST to the RefMet batch
endpoint and a `KEGGREST::keggList("compound")` fetch, on every run. RefMet is a
rolling database with no releases, so those three columns silently tracked whatever
it held on the day the pipeline happened to run — the header of
`generate_metab_qc_norm.R` used to carry a hand-written "RefMet map last run for
v1.4: 2026-07-24" precisely because nothing else recorded it. Freezing the responses
makes the metabolomics annotation reproducible and auditable; the only script that
contacts either service is the builder below.

## Artifacts

| File | Provenance | Consumed by |
|------|-----------|-------------|
| `refmet_name_map.txt` | one POST of every queried name to `name_to_refmet_new_minID.php`, all 14 returned columns kept | `.annotate_refmet()` in `generate_metab_qc_norm.R`, via `refmet_map()` in `refmet_lib.R` |
| `kegg_compound_list.txt` | `KEGGREST::keggList("compound")` verbatim (entry + `;`-delimited synonyms) | `.get_kegg_ids_from_snapshot()`, via `kegg_compound_list()` in `refmet_lib.R` |
| `refmet_snapshot.json` | written by the builder: build date, endpoint, row counts | provenance; staged to `resources/` beside the map |
| `refmet_lib.R` | — | sourced by both the builder and `lib/qc_helpers.R` |

`refmet_name_map.txt` is keyed by **`Input name`** — the exact string the pipeline
asks RefMet about, which is the submitted `refmet_name` after
`refmet_fix_names()` (the curated corrections for typos, capitalization, outdated
LIPID MAPS names and lab-internal aliases) and after `refmet_lookup_key()` strips
the untargeted platforms' LC suffixes (`_hp_a`, `_rp_b`, …). Both functions live in
`refmet_lib.R`, which the builder and the stem share, so the keys the snapshot is
built with and the keys looked up at run time cannot drift apart.

## Why the response and not the whole database

RefMet publishes a full-table REST dump, and snapshotting that would be the closer
analogue of `ensembl_v105_gene_attributes.rds`. It is not usable here: the batch
endpoint resolves lab-supplied **synonyms and aliases** to standardized names
server-side, and that synonym resolution is not in the dump. Reproducing it locally
would change results. So the snapshot is the endpoint's own answer for the union of
every name the study asks about — the same thing the Ensembl cache does (the union
of every attribute set the stems request), just keyed by query rather than by gene.

The cost of that choice: coverage is bounded by the name universe at build time. A
name absent from the snapshot is therefore a hard error in `.annotate_refmet()`, not
a silent `NA` — RefMet returns a row for every name submitted, including ones it
cannot resolve, so a missing key means the snapshot is stale, not that RefMet has no
entry. Refresh it when the metabolomics feature set changes.

## Build / refresh

Both artifacts are small and committed to git, so a fresh clone already has them —
unlike `txdb_hsapiens_ensembl_v105.sqlite`, nothing needs fetching from the bucket.
Rebuild only to pick up RefMet changes or to cover new metabolite names:

```bash
Rscript build_refmet_cache.R            # skips artifacts that already exist
FORCE=1 Rscript build_refmet_cache.R    # re-query and overwrite both
```

The query universe is read from the raw metabolomics metabolite-metadata files
under `staging/raw-files/metabolomics_qc_norm` (override with `METAB_RAW_DIR`);
those are the gated quant-id inputs, so a refresh needs them downloaded. A refresh
takes the **union** of the names found there and the keys already in the snapshot,
so rebuilding with only some platforms downloaded never shrinks coverage.

Requires network access to `metabolomicsworkbench.org` and `rest.kegg.jp`, and the
`curl`, `KEGGREST`, `jsonlite` packages. The builder errors loudly if RefMet fails
to echo back every submitted name, rather than writing a partial map.

**Built: 2026-08-04** — 1,700 distinct RefMet queries (1,655 resolved to a RefMet
name), 19,620 KEGG compounds.

## Published copy

`stage_refmet_map.R` (a step 06 stem) copies the map and its provenance JSON into
the freeze as `resources/motrpac_human-precovid_refmet-map_v<version>.txt` and
`resources/motrpac_human-precovid_refmet-map_provenance_v<version>.json`, so Stage 2
uploads them with everything else. Versioned from v2.0 on, same treatment as
`resources/txdb_hsapiens_ensembl_v105.sqlite` — the version is a release stamp, and
the freeze version tracks the study release, not this annotation snapshot, whose date
travels in the JSON. Both versions are set under `resources` in
`config/file_versions.json`.

The per-platform `metadata_features` files carry only `refmet_name`, `refmet_id` and
`kegg_id` for features that survived QC. The published map is the whole table —
every name asked about, plus the class hierarchy (super/main/sub class), formula,
exact mass, and the HMDB, ChEBI, LIPID MAPS, PubChem and InChIKey cross-references.
