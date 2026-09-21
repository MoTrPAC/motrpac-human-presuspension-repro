# motrpac-human-presuspension-repro dependency graph

A "targets-style" dependency graph of this repo: which Make target runs which
stage driver, which driver runs which step, which step runs which R stem, what
each stem reads and writes, and which nodes have the widest blast radius if
changed.

This is the pipeline-shaped sibling of the package-shaped graph in
`notes/PreCovid Comp/dependency_graph/`, which maps the three MoTrPAC R packages
(`data-raw/*.R` → `.rda` objects → NAMESPACE exports → figure scripts). That
topology does not exist here. This repo is a Make/shell pipeline whose currency
is BIC-named freeze files, so the node and edge taxonomies are rebuilt from the
ground up; the rendering layer (visNetwork, the two CSVs, `app.R`) is the same.

The shape it draws:

```
external sources → preflight objects → make target → stage driver → build step
    → R generator stem → freeze files → staging bucket / downstream packages
```

## Regenerate

```bash
make depgraph                                  # from the repo root
Rscript docs/dependency_graph/build_depgraph.R  # equivalent
Rscript docs/dependency_graph/build_depgraph.R --root /path/to/another/checkout
```

Needs R (≥4.4) with `visNetwork`, `htmlwidgets`, `igraph`, `dplyr`, `jsonlite`.

The script parses **source text and config only**. It does not run the pipeline,
does not read `staging/freeze`, and does not need the gitignored build outputs
(`scripts/00_preflight/data/*.rds`, `staging/**`) or consortium data access — so
it works on a fresh clone. Re-run it whenever the pipeline changes; the committed
HTML and CSVs are a snapshot, not a source of truth.

## Outputs

| File | What |
|------|------|
| `depgraph.html` | **Layer 1** — the full interactive graph. Hierarchical left→right: external sources → preflight → make → drivers → steps → stems → freeze files → bucket/packages. Includes the shared-library call graph. |
| `depgraph_summary.html` | **Layer 2** — collapsed to six stage super-nodes plus the top hotspots; edge labels = number of dataflow edges crossing each stage boundary. |
| `depgraph_nodes.csv` | Node table (`id, label, type, stage, ome, origin, detail, n_downstream, hotspot, level`). |
| `depgraph_edges.csv` | Edge table (`from, to, type`). |
| `depgraph_freeze_files.csv` | The per-file freeze inventory behind the grouped freeze nodes — one row per `config/file_versions.json` file, with its tissue_code / ome / category / details, its version, its `why_bumped` note, and the stem attributed to it. |

### Reading `depgraph.html`

- **Color = stage.** slate = Makefile orchestration, teal = Stage 0 preflight,
  purple = Stage 1 build-data, orange = Stage 2 upload, pink = Stage 3
  downstream packages, grey = outside this repo.
- **Shape = node type:**
  - database = make target
  - box = stage driver / downstream sink
  - square = build step, R generator stem, shared lib file, preflight builder
  - dot = freeze output, data-object group, preflight asset
  - triangle-down = shared library function
  - diamond = external source / required input
  - hexagon = the Stage 1 required-inputs gate
  - star = test script
- **Size** = downstream blast radius. **Pink border** = hotspot (top decile).
- **Filter** with the stage dropdown; **click** a node to highlight its 2-hop
  neighbourhood; use the id dropdown to jump to one.

### Node ids

| Kind | Id form | Example |
|------|---------|---------|
| make target | `make:<target>` | `make:data` |
| stage driver | `sh:<repo-relative path>` | `sh:scripts/10_build_data/10_build_data.sh` |
| build step | `step:<folder>` | `step:09_build_da` |
| R generator stem | `stem:<folder>/<stem>` | `stem:06_generate_qc_norm/generate_prot_pr_qc_norm` |
| freeze output | `freeze:<ome>\|<category>\|<details>` | `freeze:prot-pr\|qc-norm\|log2-mn` |
| data-object group | `obj:<node>` | `obj:QC_NORM` |
| shared lib file / function | `lib:<path>` / `fn:<name>` | `fn:write_with_path_name` |
| external source | `src:<id>` | `src:quant_id_bucket` |
| required input | `gatein:<label>` | `gatein:preflight:OME_TISSUE_CODE` |

Freeze nodes are **grouped**: one node per (ome, category, details) family across
all tissues, because 226 individual file nodes is a wiring diagram, not a
schematic. The node tooltip carries the file count, versions and tissue codes;
`depgraph_freeze_files.csv` has the per-file rows.

## How edges are derived (hybrid)

**Auto (reproducible):**

- `orders` — the Makefile target/prereq graph, resolving the `$(STAMPS)/x → x`
  stamp indirection so `data` depends on `preflight` rather than on a file.
- `runs` — target → the script its recipe invokes; `10_build_data.sh` → each of
  the 14 steps in its `STEP_DIRS` array; each `build.sh` → its stems, from the
  `ALL_STEMS` array (steps 06/09) or the `"${HERE}/<name>.R"` single-stem form.
  Deriving drivers from the recipes rather than hardcoding them is why renaming
  a stage (`30_figures` → `30_update_relevant_packages`) needs no edit here.
- `produces` — stem → freeze family. Every file in `config/file_versions.json`
  arrives already split into `tissue_code / ome / category / details` — the four
  levels of that map's tree, which are exactly what `write_with_path_name()`
  composes before appending `_v<version><ext>` — and is attributed to the stem
  that owns that ome. This is the one place the script sources repo code
  (`lib/file_versions.R`, for `freeze_files()`) rather than parsing it, so the
  filename rule has a single definition. Attribution is by **string-literal
  scan**: the stem's own source is parsed and its ome literals matched against
  the ome vocabulary. Also preflight builder → its GMT, and step → its ledger
  objects.
- `consumed_by` — external asset → the script its manifest names; freeze family →
  the later steps that read it; every freeze file → the Stage 2 upload driver.
- `feeds` — the object-level DAG, read straight out of the `deps` column of
  `docs/data_objects.tsv` (semicolon-separated node names). That ledger is
  hand-maintained upstream truth, so it is ingested rather than re-derived.
- `requires` — the `label::path` arrays and `compgen` glob checks in
  `check_required_inputs.sh`, the hard gate every Stage 1 run passes first.
- `checked_by` — every row of `external_assets.tsv` → the preflight driver, which
  is precisely the stage that verifies them.
- `tested_by` — each step → the `*_tests.R` beside it.
- `defines` / `calls` / `sourced_by` — the intra-repo function call graph,
  recovered by statically parsing every `.R` under `scripts/` with R's own parser
  (`getParseData` — AST tokens, **no code execution**). Only **top-level**
  function definitions become nodes, so closures local to another function are
  skipped. A call site inside a definition is attributed to that definition (so
  the library's internal wiring shows as fn → fn); everything else is attributed
  to its file. The callee universe is restricted to functions defined in a
  **shared** file, and "shared" is itself derived — any `.R` that another `.R`
  `source()`s — which picks up both `10_build_data/lib/*` and the step-local
  `09_build_da/da_common.R` without naming either.

**Curated (two tables, both in `build_depgraph.R`):**

- `OME_FAMILY` — which stem family owns an ome whose stem never names it as a
  string literal. Only the metabolomics stems need this: they build their ome
  list out of `OME_TISSUE_CODE` at run time, so there is no literal to grep.
  Every other ome is auto-attributed.
- `FREEZE_CONSUMERS` — which later steps read which freeze category. These reads
  go through `load_qc_local()` or a `list.files()` sweep over `staging/freeze`,
  so no filename literal exists in the consuming source either. The one entry
  that is not a plain category rule is `metadata|removed-samples`: those files are
  produced upstream and copied verbatim into the freeze by
  `stage_removed_samples.R` so the release carries them, and nothing in Stage 1
  reads them back. (`OUTLIERS.R` does read them — from the vendored preflight
  sources, not from the freeze.)

## Current hotspots (re-derived each run; see the console summary)

`n_downstream` is the blast radius: how many nodes are downstream of a change to
this one. It is computed on **dataflow edges only** — `runs`/`orders`/`tested_by`
are control flow, and counting them makes `make precheck` the biggest node in the
repo, which is true and useless. `calls` is flipped on the way in, because the
blast radius of a helper is everyone who calls it.

Highest reach → change these with the most care (counts from the run that
produced the committed CSVs; 20 nodes clear the 95th-percentile bar):

| Node | Reach | Why |
|------|-------|-----|
| `00_preflight/data-raw/lib/helpers.R` | 193 | the Stage 0 shared library — every preflight builder sources it |
| `10_build_data/lib/file_versions.R` | 181 | every freeze filename's version resolves through `freeze_version()` |
| `10_build_data/lib/qc_helpers.R` | 174 | `write_with_path_name`, `load_qc_local`, `quantid_path`, … |
| `src:quant_id_bucket` | 92 | the gated raw tier every qc-norm stem reads |
| `src:refmet_api`, `src:keggrest` | 61 | live metabolite annotation, and non-deterministic across runs |
| `generate_metab_qc_norm` | 60 | the widest stem — 14 of the 23 omes |
| `generate_differential_modeling_functions.R` | 45 | the `dream` modelling engine behind every DA stem |
| `filter_paired_n.R`, `da_common.R` | 39–41 | the rest of the DA engine |

## Interactive explorer (`app.R`)

`depgraph.html` is a static snapshot of the whole graph. `app.R` is an
interactive explorer on the same two CSVs: pick a node and it traces that node's
**upstream dependencies** — every source, driver, step, stem and shared helper
that must be correct for it to reproduce — or its **downstream blast radius**.

```r
shiny::runApp("docs/dependency_graph")
# or, from inside this folder:  shiny::runApp()
```

Needs `shiny`, `visNetwork`, `igraph`, `dplyr`, `DT`. Reads only the two CSVs.
Re-run `build_depgraph.R` first if the pipeline changed, then relaunch.

- **Direction** — *Upstream* (what it depends on) or *Downstream* (what depends
  on it). Upstream reverses every edge type except `calls`, so the view flows
  external source → driver → step → stem → freeze → bucket.
- **Focal node** — defaults to the freeze-output list (switch the type dropdown
  to root on any node); narrow it with the ome filter.
- **Max hops** — bounds how far the lineage is traced.
- **Include shared-lib calls / defines** — adds the `calls`/`defines` wiring.
- **Include make/driver control flow** — adds `runs`/`orders`; turn it off for a
  pure dataflow view.
- **Click a node** to re-root on it. The side panel counts reachable nodes, the
  freeze outputs downstream, and the consortium-gated inputs upstream; the
  *Dependency table* tab lists them with hop distance.

## Caveats

- Freeze attribution is **name-based**, not execution-based: a stem is credited
  with a file because it owns that file's ome, not because the run was observed
  writing it. It is calibrated for relative blast radius, not for proving a stem
  emitted a given byte.
- `FREEZE_CONSUMERS` is at **category** granularity. Step 09 is shown consuming
  every `metadata` and `qc-norm` family, which is what `load_qc_local()` does in
  aggregate, but an individual DA stem only touches its own ome.
- The `calls` graph is **static**: it sees `f()` syntactic calls, not functions
  reached via `do.call`/`match.fun`/`get()` or S3/S4 dispatch. Helper names that
  repeat across shared files merge into one node.
- Stub steps (10, 12, 13, 14 — the ones that `exit 77`) appear with
  `origin = "stub"` and produce nothing. Their ledger objects still carry `feeds`
  edges, because the ledger records the intended DAG, not the adapted subset.
- The seven `.html` QC reports in `file_versions.json` are produced upstream and
  carried into the release; they are marked `origin = "carried"` and have no
  producing stem. The two keyed as `metadata / qc-report` are the reason the
  carve-out keys on the file extension rather than the category field.
- The `resources/` tier has no tissue_code, so it sits outside the tree in the map
  and outside the grouped freeze nodes here: one node per file, wired to its
  staging stem. `motrpac_human-precovid_1kg_pca.csv` is the exception with no
  producing stem — the collection carries it, this pipeline never writes it.
- `make clean` and the `docker-*` targets have no dataflow edges at all. That is
  correct — they are out-of-band, and no pipeline stage depends on them.
- The committed `.html`/`.csv` outputs are a snapshot of the commit that
  generated them. Regenerate with `make depgraph` rather than trusting them after
  a pipeline change.
