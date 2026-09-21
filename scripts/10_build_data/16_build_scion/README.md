# Step 16 — SCION networks

The regulatory networks behind Figure 6 and Extended Data 8. Two networks, each inferred
per exercise group (`ADUResist`, `ADUEndur`):

| network | regulators | targets |
|---|---|---|
| `muscle_protph_transcript` | muscle prot-ph, TF sites only (step 14) | muscle transcript-rna-seq, DE only |
| `blood_metab_muscle_transcript` | blood metabolomics, all features | muscle transcript-rna-seq, DE only |

Features are clustered by the shape of their z-score trajectory across the
`exercise_with_controls` contrasts, then a random forest per cluster scores every
regulator against every target in that cluster. The cluster hubs are connected to each
other in a second pass.

```bash
bash scripts/10_build_data/16_build_scion/build.sh
# or: STEPS=16 bash scripts/10_build_data/10_build_data.sh
```

Requires steps 06, 07, 08, 10 and 14.
Output: `staging/scion/<network>/<group>/` and `staging/scion/scion_manifest.tsv`.

## Not a data object

The output is edge tables under `staging/`, not an `.rda` routed into either package.
2 networks x 2 groups x (1 + `SCION_PERMUTATIONS`) runs is far past what lazy data can
carry, and no accessor in either package reads a network. Figure 6B is a Cytoscape layout
of a hand-merged export of these tables, and ED8A's TFEB target list is a spreadsheet taken
out of that same session. This step regenerates the input those manual steps started from.

It is therefore absent from `docs/data_objects.tsv`. It runs by default like every other
step; `RUN_SCION=FALSE` opts out, and `10_build_data.sh` then records it SKIP (exit 77).

## Cost

One run is a random-forest fit per cluster per group, which is hours to days across two
exercise groups and two networks. `SCION_PERMUTATIONS` stays 0: the networks are too dense
for a permutation null to be the boundary, so edges are cut by a stricter weight cutoff
instead, and every permutation multiplies that cost. `config/slurm.json` keys this step, so
`EXECUTOR=slurm` submits it as its own job rather than blocking the stage:

```bash
STEPS=16 EXECUTOR=slurm bash scripts/10_build_data/10_build_data.sh
```

The observed network for a group is always inferred before its permutations, so a run
stopped partway through leaves the network the figure needs rather than shuffles of it.
Re-entering the step re-infers every run; `SCION_FORCE=FALSE` skips any already on disk.

## `SCION_CORES` changes the numbers, not just the runtime

Above 2 cores, `RS.Get.Weight.Matrix()` computes the weight matrix through a
`parallel::makeCluster(num.cores - 1, type = "FORK")` cluster. Each worker inherits the
parent's RNG state at the fork point and `clusterSetRNGStream()` is never called, so which
target genes land on which worker is what fixes each forest's seed. **The same inputs at a
different `SCION_CORES` produce different edge weights** — and 1 or 2 cores take the serial
branch, which is a third result again.

This is upstream's behaviour, carried over unchanged rather than fixed, because changing it
would move every published weight. It means `SCION_CORES` is part of the result and not a
performance knob: to compare two runs, set it the same. `scion_manifest.tsv` records
`num_cores` per run for that reason.

Reproduced on 2026-08-13 against acute-repro's 2026-08-12 Stage 07 output, the last run
before inference moved here — `muscle_protph_transcript` / `ADUResist`, `num.cores = 10`,
R 4.4.1, randomForest 4.7.1.2. The regulator and target matrices matched to float
round-trip (233 x 156 and 8,897 x 156, same row and column order), and
`network_cluster_1.csv` came out byte-identical: sha256 `a20023a8…`, 9,798 edges, zero
weight difference.

## Where this came from

`MotrpacHumanPreSuspensionAnalysis::run_SCION()` and its four helpers, which this step
replaces. That file (`R/run_SCION.R`) is removed from the package: it was the only caller of
all four helpers, and it held most of the package's undeclared Analysis-to-Data runtime calls.

The driver is the loop from
`MotrpacPreSuspensionAcute/figures/landscape/figure_7/scion_figures.Rmd`, quoted verbatim
in `SCION_NETWORKS.R` above the `NETWORKS` list.

| file | vendored from |
|---|---|
| `run_SCION.R` | `run_SCION()`, `scion_data_processing()`, `RS.Get.Weight.Matrix()`, `RSGWM2()` |
| `scion_matrixes.R` | `.load_scion_matrixes()`, `.find_shared_scion_matrixes()`, and `combine_qc_matrixes()` from the Data package |
| `scion_run_cmeans.R` | `scion_run_cmeans()` |

`.prepare_DA_results()` is not re-vendored — step 12 carries it and `SCION_NETWORKS.R`
sources that copy, so the two steps cannot drift on how DA results become z-score matrices.

## Three upstream bugs fixed on the way

Each file's header lists what it diverges from upstream on; `git log` is the record of what
the original looked like.

1. **`make.names()` on the matrix rownames** (`run_SCION.R`). It rewrote every rowname
   carrying a hyphen, a space or a leading digit — which is both isoform suffixes and the
   `prot-ph` ome tag — while the DA table the clustering keys on kept the original ids. Those
   features matched nothing, were never clustered, and never reached inference. Removing it
   also removes the regex that reconstructed the sanitised names, and the two `data.frame()`
   calls building the edge tables now pass `check.names = FALSE` so the `Regulator` column
   stays joinable against the matrix rownames.

2. **`pheno` was unreachable** (`scion_matrixes.R`). `.find_shared_scion_matrixes()` fetched
   it with `get("pheno", envir = asNamespace("MotrpacHumanPreSuspensionData"), inherits = FALSE)`.
   `pheno` is lazy-loaded data with no binding in that namespace, so the call could not find
   it and every `run_SCION()` call died with *object 'pheno' not found*, regardless of what
   was attached. It reads `load_pheno()` here.

3. **A flat list handed to a function expecting a nested one** (`scion_run_cmeans.R`).
   `run_SCION()` built `list("tissue.ome" = <data.frame>)` and passed it to
   `.prepare_DA_results()`, which opens with `unlist(DA_list, recursive = FALSE)` because it
   expects `tissue -> ome -> data.frame`. Unlisting a flat list of data frames flattens each
   frame into its columns, the name filter then matches nothing, and cmeans dies stacking
   NULL with *non-numeric matrix extent*.

One upstream decision worth knowing, noted at `RSGWM2()`: target genes are deliberately
**not** removed from the input matrix. Removing them drops autoregulation from the network,
which is what it was written for, but it breaks a network with a single regulator.

## Inputs that are no longer vendored

Both files this step needs were previously vendored by hand into `acute-repro`:

- the multiply-imputed muscle prot-ph matrix is **step 06's own freeze output**
  (`generate_prot_ph_imputed.R`), read from `staging/freeze/proteomics/qc-norm/`. Step 08
  skips the imputed freeze files, so the `*_QC` objects carry no `qc_imputed` component and
  upstream's `qc_imputed` substitution has nothing to substitute.
- the TF regulator pool is **step 14's `UTORONTO_TFs`**. The object shipped by the Analysis
  package is the raw 2,765 x 28 download with no `feature_id` column, and `subset_qc()`
  skips a NULL `desired_features` — so `subset_TFs = TRUE` is a silent no-op against it and
  SCION gets every phosphosite rather than the TF pool.

## Tests

`scion_tests.R` checks the naming contract between the two halves of `run_SCION()`: every
`Regulator` and `Target` id in the edge table is a rowname that went in, the hubs join the
regulator matrix, and hyphenated (isoform-suffixed) regulator ids reach the network — the
canary for a sanitising step reappearing anywhere in the path. That failure mode is silent:
a severed join leaves `hubregdata` empty, `RS.Get.Weight.Matrix()` dies on a zero-column
matrix, and the run ends *before* `write.table()`, leaving a directory of per-cluster CSVs
and no network file at all.

Observed networks only — the permutations run identical code. Input rownames are read from
the `myregdata.csv` / `mytargetdata.csv` that `run_SCION()` dumps in the hub step, so the
matrices are not rebuilt; a run with no hub layer SKIPs those checks.
