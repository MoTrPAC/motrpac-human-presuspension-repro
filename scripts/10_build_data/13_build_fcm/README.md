# Step 13 — build_fcm

Fuzzy c-means (FCM) clustering of the acute differential-analysis results, and enrichment of
the clusters that come out of it. Features are clustered by the *shape* of their z-score
trajectory across the `exercise_with_controls` contrasts, so a cluster is a response pattern —
"up early, back to baseline by 24 h" — rather than a set of features that merely moved.

Three objects, all in `scripts/10_build_data/data/`:

| Object | Function | What it is |
|---|---|---|
| `FCM_CLUSTERS` | `run_cmeans()` | one `fclust` per tissue: centroids, membership probabilities, hard assignment, and the scaled input matrix |
| `FCM_CAMERA` | `run_cluster_cameraPR()` | CAMERA-PR on the membership probabilities — which signatures follow each centroid's trajectory |
| `FCM_ORA` | `run_cluster_ORA()` | ORA on the disjoint (hard-assigned) clusters |

Requires steps 03, 04, 07 and 10.

`FCM_ORA` is not used downstream. Upstream's own note is that it exists to be compared against
`FCM_CAMERA`; that comparison is the check that the fuzzy result is not an artefact of the
membership weighting, which is why it is built rather than dropped.

That comparison is also why the two objects share one universe. `run_cluster_ORA()` takes its
background from every clustered feature, then hard-assigns only those reaching `min_prob = 0.3`
in some cluster — so in the `FCM_ORA` table the cluster sizes sum to well under
`background_size`, by 12–57% depending on tissue × ome (muscle metab 57%, blood
transcript-rna-seq 12%, 35% across the eight). The gap is that threshold. Narrowing the
background to the assigned features instead would change every p-value and would stop
`FCM_CAMERA`, which scores all of them continuously, from being the comparison it is built for.

## What is in this folder

| File | Kind | Notes |
|------|------|-------|
| `build.sh` | driver | builder, then diagnostics unless `FCM_DIAG=0`, then tests unless `RUN_TESTS=0` |
| `FCM_clustering_results.R` | adapted | builds `DA_list`, calls the three functions, saves |
| `run_cmeans.R` | **vendored** | `run_cmeans()` + `.reorder_clusters()` |
| `run_cluster_cameraPR.R` | **vendored** | `run_cluster_cameraPR()` + `.prepare_cluster_mem()`, `.split_membership_by_ome()` |
| `run_cluster_ORA.R` | **vendored** | `run_cluster_ORA()` + `.cluster_ORA()`, `.fast_list_intersect()` |
| `fcm_diagnostics.R` | new | sweeps the cluster count and writes the plots behind it |
| `fcm_tests.R` | tests | structure, scope, statistics, BH grouping, diagnostics, diff-vs-package |

Vendored **code** sits alongside the stem, the same convention as step 12 and
`09_build_da/da_common.R`. Upstream's
`data-raw/FCM_clustering_results.R` is three one-line calls, so all the substance is in the
functions, which is why they are vendored rather than reimplemented.

`.prepare_DA_results()`, `.create_index()` and `.prepare_sets()` are **not** re-vendored here.
Step 12 already carries them with the same no-installed-package change, and
`FCM_clustering_results.R` sources `../12_build_camera_results/run_cameraPR.R`. One copy, so
the two steps cannot drift on how DA results become z-score matrices or on how sets are
filtered to a background.

## Vendored code — provenance

Copied from `MotrpacHumanPreSuspensionAnalysis/R/{run_cmeans,run_cluster_cameraPR,run_cluster_ORA}.R`.
Every change is marked in place with the upstream original kept commented directly above it,
and they are nearly all the same kind: `MotrpacHumanPreSuspensionAnalysis::<OBJECT>` becomes
the pipeline's local global of the same name — `MOLECULAR_SIGNATURES` (03), `SET_TO_ID` (04),
`HUMAN_FEATURE_TO_GENE` (07), `CONTRAST_CONVERTER` (10).

The structural changes:

- `run_cmeans()` gains a `DA_list` argument, passed through to `.prepare_DA_results()`.
  Upstream hardcodes `DA_list = NULL` there, which sends the helper to
  `load_differential_analysis()` and the installed analysis package — the dependency this stage
  removes. `run_cameraPR()` already took such an argument in step 12.
- **The interactive cluster-number branch is removed.** See below.
- `check_package_installation()` is dropped; `FCM_clustering_results.R` fails up front if
  Mfuzz, Biobase, e1071 or TMSig is missing.
- `HUMAN_FEATURE_TO_GENE` is wrapped in `as.data.frame()` before the dplyr chain in
  `.split_membership_by_ome()`, because the pipeline's step-07 object is a keyed `data.table`.
  Step 12 makes the same wrap.
- `e1071` and `Biobase` are attached. Mfuzz declares both in `Depends` rather than `Imports`,
  so `Mfuzz::mfuzz()` resolves `cmeans()` off the search path; called namespace-qualified from
  a plain script — as upstream does from a package that attaches them itself — it fails with
  `could not find function "cmeans"`.

Clustering and both tests are otherwise unchanged: same scaling, same `set.seed(0)`, same
`mestimate()`/`mfuzz()` calls, same cluster reordering, same `cameraPR.matrix()` and
hypergeometric calls.

## Choosing the number of clusters

This is the one place where upstream cannot be reproduced as written. `run_cmeans()` accepts a
*range* of cluster numbers, and when given one it plots `Mfuzz::Dmin()` and then blocks on
`readline()` for the operator to type the number to use. A batch build has no console, and an
object whose contents depend on what someone typed is not reproducible.

The two halves are separated instead:

- **`fcm_diagnostics.R`** does the sweep non-interactively and writes the evidence to
  `scripts/10_build_data/13_build_fcm/fcm_diagnostics/` — beside the step, not in `data/`,
  which carries only the `.rda` objects the packages consume. `CLEAN=1` clears it.
- **`FCM_clustering_results.R`** takes the chosen number as a build parameter —
  `FCM_K_ADIPOSE`, `FCM_K_BLOOD`, `FCM_K_MUSCLE`, defaulting to **13 / 12 / 12**.

Passing a range to `run_cmeans()` is an error rather than a silent prompt. The loop is:
run the build, read the plots, set `FCM_K_*` if the default is wrong, rebuild. A changed `k`
changes `FCM_CAMERA` and `FCM_ORA` too.

### What the diagnostics produce

Per tissue: `dmin_<tissue>.png` (the headline curve on its own) and
`fcm_diagnostics_<tissue>.pdf` (all of it, one figure per page). Across tissues:
`fcm_cluster_sweep.tsv`, every metric at every k.

Sweep metrics, all read off **one fit per k** — `Mfuzz::Dmin()` refits internally and returns
only its own statistic, so the fit is done in the script and all four metrics come off it.
`min_centroid_dist` is computed exactly as `Dmin()` computes it, so that curve is the same one
upstream shows at the prompt.

| Plot | Question it answers |
|---|---|
| Min. centroid distance vs k | `Dmin()`'s statistic. Falls as clusters are added; the **elbow** is where extra clusters stop being distinct. |
| Max. centroid correlation vs k | The most similar pair of centroids. Approaching 1 means two clusters have effectively the same trajectory — k is past useful. |
| Core fraction vs k | Share of features whose best membership clears 0.3, the threshold `run_cluster_ORA()` hard-assigns at. A sharp drop marks where clusters stop being well populated. |
| Smallest cluster vs k | Catches a k that only adds near-empty clusters. |
| Centroid trajectories at the built k | The profiles themselves — duplicates are visible directly. |
| Centroid correlation heatmap at the built k | The same redundancy check as the max-correlation curve, per pair. |
| Cluster sizes + membership histogram at the built k | How evenly features distribute, and how much of the data clears 0.3. |

Knobs: `FCM_DIAG=0` skips the sweep entirely, `FCM_DIAG_KMIN` / `FCM_DIAG_KMAX` (default
2..20) set the range, `FCM_DIAG_REPEATS` (default 1, upstream's own choice) averages over
random initialisations, and `FCM_DIAG_DIR` relocates the output directory that
`fcm_diagnostics.R` writes and `fcm_tests.R` reads. The sweep reads the scaled input matrix and the weighting exponent `m`
back out of `FCM_CLUSTERS`, so it reruns on its own without redoing the DA preparation.

The sweep is the expensive part of this step — one `mfuzz()` fit per k per tissue per repeat,
against the same 20k–42k feature matrices the build clusters once.

## Scope

**Clustering** uses five omes, `run_cmeans()`'s own default: `transcript-rna-seq`, `prot-pr`,
`prot-ol`, `prot-ph`, `metab`. Three further restrictions are upstream's and are kept:

- only `exercise_with_controls` contrasts, with the two `during` contrasts dropped inside
  `run_cmeans()`;
- features must be measured at **every** remaining timepoint (`complete.cases`) — clustering a
  trajectory requires the whole trajectory;
- adipose `prot-pr` and `prot-ph` are dropped for insufficient timepoints.

Epigen and the two clinical DA objects are not clustered. Epigen is dropped at the file glob
rather than by `selected_omes` — it would fall out there too, but only after loading ~7 GB into
`DA_list` to hand it over.

**Enrichment** covers four of those five. `prot-ol` is absent from
`.prepare_cluster_mem()`'s choices upstream, so blood `prot-ol` features are clustered — they
contribute to the blood centroids — but are not enrichment-tested. The released `FCM_CAMERA`
has no `prot-ol` rows either.

Both tests run over every collection in `MOLECULAR_SIGNATURES`, the upstream default, which
*includes* `PTMSIGDB` — unlike step 12, which excludes it explicitly. It makes no difference:
PTMSigDB sets are keyed on flanking sequences carrying a `;u`/`;d` direction suffix that only
the PTMsigDB branch of `.prepare_sets()` strips, and that branch fires only when `PTMSIGDB` is
the *sole* database. Mixed in with the others its sets match nothing and drop out at the size
filter, leaving 12 databases. The released object has no `PTMSIGDB` rows either.

## Phosphoproteomics depends on the step-07 fix

Exactly as step 12 does. Cluster membership is keyed on the **flanking sequence** for
`prot-ph`, because PhosphoSitePlus kinase sets (the `PSP` database) are defined over flanking
sequences — `.split_membership_by_ome()` selects the column unconditionally, so without it the
step fails for *every* ome, not just prot-ph. `FCM_clustering_results.R` checks for it up front
and says to rebuild step 07; `fcm_tests.R` asserts `prot-ph` is present in both enrichment
objects, which is the regression guard.

## Comparison against the released package

`fcm_tests.R` reports a diff against the shipped objects. It is **INFO, never FAIL** — these
are computed from locally refit DA, so they move with the freeze, and `mfuzz()` is seeded but
its result still depends on the exact input matrix.

As built against freeze `v2.0`:

| | built | package |
|---|---|---|
| adipose | 20,757 features × 13 clusters | 20,757 × 13 |
| blood | 23,086 × 12 | 23,086 × 12 |
| muscle | 41,775 × 12 | 41,775 × 12 |
| `FCM_CAMERA` | 434,337 × 13 | 434,337 × 13 |
| `FCM_ORA` | 434,337 × 16 | 434,337 × 16 |
