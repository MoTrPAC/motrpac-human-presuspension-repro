# Stage 3 — update the relevant packages

Carries the regenerated data into the two packages that distribute it, and gates that on the packages
being in a state worth releasing into.

```
make update-relevant-packages       # audit, then carry into test packages
AUDIT_ONLY=1 bash scripts/30_update_relevant_packages/30_update_relevant_packages.sh
STRICT=0 make update-relevant-packages    # report findings without failing the stage
RUN_CHECK=1 make update-relevant-packages # also R CMD check each test package (slow)
```

## Two halves

**Audit — read-only.** Validates both package trees, checks the dependency graph for
cycles, and inventories code no longer needed. Reads the two package repos and writes only under
`logs/`. A hard cycle or a structural ERROR fails the stage, but at the verdict at the end, not
before the carry: the carry runs either way and its findings are reported alongside the audit's.
Nothing is at stake in that ordering, because the carry only ever writes under `staging/` — a
package that cannot install yields a test package that cannot install, never a broken release. The
cost of a failed audit is the carry's runtime, so use `AUDIT_ONLY=1` to stop at the audit.

**Carry — writes only under `staging/`.** Routes every built object to its package, assembles a
**test package** per repo, re-documents it, bumps its version and NEWS, and runs its own test suite
against the new payload. Neither package checkout is modified. Promoting a test package into its
repo — copying `data/`, committing, pushing, opening the PR — stays manual, the same posture
`make promote` takes toward bucket promotion.

### Why test packages rather than the repos

A carry that writes into `MotrpacHumanPreSuspensionData/data/` is 584 MB of binary churn in a
checkout that also holds hand-written code, and it cannot be reviewed before it exists. Building the
same thing under `staging/test-packages/` makes the result inspectable, installable and testable
while leaving both repos exactly as they were, and it means a failed carry needs no `git checkout --`
to recover from. The test package is a copy of its source tree with a `data/` assembled from the
routing table, so it inherits the branch, version and uncommitted edits of whatever it was built
from — step 1 records all three.

### The carry steps

| Step | Does |
|---|---|
| `01_validate_targets.R` | both repos are git checkouts; records branch, SHA, version and uncommitted paths; fails if a version is behind the v1.3 comparison baseline or the build is missing |
| `02_route_objects.R` | routes every object by `carry/lib/routing.R`, cross-checks against `docs/data_objects.tsv`, and names every packaged object no source produces — as kept or as withdrawn |
| `03_build_test_packages.R` | assembles each test package; writes every object through `save(compress = "bzip2")` to match how the packages store `data/`; records a content verdict per object |
| `04_document.R` | regenerates `@format` where the payload changed shape, documents objects new this cycle, removes the documentation of withdrawn ones, removes any `inst/PROVENANCE.tsv`, re-roxygenises, and verifies every data object has a man page and no withdrawn object still has one |
| `05_version_and_news.R` | bumps `Version` and `Date`, prepends a NEWS entry naming the breaking column removals and the withdrawn object |
| `06_check_and_test.R` | re-pins the expectations the new payload legitimately breaks, adds a regression test for the removed columns, runs both suites |

### What the routing rule is

`carry/lib/routing.R` is the single source of truth, and `docs/data_objects.tsv` records the same
answer so the two can be checked against each other:

- `cln_*`, `pheno`, `blood_transcript_*`, and every `*_QC` → `MotrpacHumanPreSuspensionData`
- everything else → `MotrpacHumanPreSuspensionAnalysis`
- the ten epigenomics `*_QC` / `*_DA` objects → **neither**; they ship through the GCS bucket

That last exception is not a size threshold. Those objects run 57–315 MB against GitHub's 100 MB
blob limit with no git-lfs in either package, and four of the ten would fit under it — but
`load_qc(epigen = TRUE)` and `load_differential_analysis()` already read the whole epigenomics tier
from the bucket, so splitting it by how well a table happened to compress would make the loaders'
source depend on compression luck. The tier moves together. Epigenomics `*_SUM_STATS` are 0.1–0.9 MB
and do ship.

The exception is only safe while those loaders resolve, and each pins the release directory as a
literal that nothing in this pipeline reads or rewrites:

| Constant | File |
|----------|------|
| `gsutil_base` | `MotrpacHumanPreSuspensionData/R/load_qc.R:266` |
| `AWS_header` | `MotrpacHumanPreSuspensionAnalysis/R/load_DA_from_AWS.R:103` |

Both name `v1.3`. The data hub's release directory is `c1.3` (`config/pipeline.env`), so
`load_qc(epigen = TRUE)` already reads a path that is not there — excluding the ten objects from
both packages leaves them shipping nowhere and loading from nowhere. Repointing is manual and
ordered: `make promote` publishes the `c2.0` tier, the two constants move to it, and only then can
either package release as 2.0.0. Whether the CloudFront mirror followed the `v`→`c` rename is not
recorded here; confirm it against the CDN rather than assuming it matches the bucket.

The rule is checked rather than asserted: it reproduces where all 183 currently-packaged objects
actually live.

### Objects the packages ship that no source produces

Step 2 names every one of them, and each is either kept or withdrawn. Kept means the package's
existing bytes go into the release untouched; none are kept at present. Withdrawn means the
release ships without it, and step 4 removes its roxygen so no man
page is left aliasing an object that will not load. An object with no source and no reason is a
FAIL, because that is indistinguishable from one the build silently dropped.

Thirty-five are withdrawn, for three different reasons, and the distinction is worth keeping visible
in the NEWS entry.

**A rename.** `CLIN_CHEMISTRY_DA` and `BLOOD_CLINICAL_CHEMISTRY_SUM_STATS` are replaced by the v2.0
split of clinical chemistry into a metabolomics and a proteomics assay. All four split objects are
built this cycle, so keeping the combined pair beside them would ship the same measurements twice
under two schemas, with only the v1.3 pair carrying stale values. The split is an exact partition of
what they held — 198 + 99 = 297 DA rows, 114 + 57 = 171 summary-stat rows.

Withdrawing them is not just a deletion: `plot_single_feature()` read both at three shipped call
sites, and `Pkg::MISSING` errors at call time, so step 4 moves those onto the split. The halves also
relabel themselves — `metab` (with `platform`) and `prot-clinical`, where the combined objects said
`clinical-chemistry` — so the DA half needs the same metab-platform normalisation the function
already applies to its main DA frame, or `assay` disagrees between the DA and the summary statistics
and the later `full_join` matches nothing; and the `"Clin. Chem."` facet label has to recognise the
two new assay strings, or clinical panels render labelled `NA`.

**A second rename.** The thirty-two `{TISSUE}_METAB_{PLATFORM}_SUM_STATS` objects are replaced by
one `{TISSUE}_METAB_SUM_STATS` per tissue. Step 11 now stacks the research metabolomics platforms
the way step 10 has always keyed the differential analysis — `assay = "metab"` with the platform in
its own column — so the two tiers nest the same way and a caller joining them no longer needs to
know that one said `metab-u-rppos` where the other said `metab`. It is an exact partition: 8,004 +
21,717 + 10,692 rows in, the same rows out, and no `(randomGroupCode, feature_id, Timepoint)` key
appears on two platforms, because step 10's CV collapse keeps each RefMet name on one platform per
tissue. `BLOOD_METAB_T_CLINICAL_SUM_STATS` is not part of the stack, on either tier.

**A genuine absence.** `BLOOD_EPIGEN_ATAC_SEQ_SUM_STATS`. Epigenomics summary statistics keep only
features at `adj_p_value < 0.05`, for file size; blood ATAC has none this cycle — 0 of 5,279,211,
the smallest being 0.0501 — so step 11 builds no object. Muscle ATAC, with 6,584 significant
features, is built as usual.

### Versions

Step 5 releases both packages at the collection's version — `NEW_VERSION` from
`config/pipeline.env`, so v2.0 of the data ships as 2.0.0 of each package. Continuing each
package's own numbering (0.0.1.101 and 0.2.4) would say nothing about which release of the
collection an install holds, which is the only question anyone asks of these two.

### Run it against a pre-promotion tree

Verdicts and NEWS are computed against whatever the package checkouts currently hold. Run the carry
twice into the same release and the second NEWS entry describes the delta against the first rather
than against the last release, and the two entries stack. If a promotion has already happened, point
`GITHUB_ROOT` at worktrees of the pre-promotion commits rather than re-running over the result.

## Outputs

| File | What |
|---|---|
| `logs/update_relevant_packages_report.tsv` | stage-level PASS/WARN/FAIL rows |
| `logs/package_audit/structure_<pkg>.tsv` | structural findings, one row each, severity-ranked |
| `logs/package_audit/cycles_verdict.tsv` | the cycle verdict — read this one first |
| `logs/package_audit/cycles_package_edges.tsv` | declared DESCRIPTION edges between the audited packages |
| `logs/package_audit/cycles_cross_references.tsv` | every cross-package mention, classified by capacity |
| `logs/package_audit/cycles_build_order.tsv` | `data-raw` generator edges |
| `logs/package_audit/cycles_bootstrap_edges.tsv` | install and build edges as one graph |
| `logs/package_audit/cycles_function_level.tsv` | call-graph cycles and self-recursion |
| `logs/package_audit/unused_code_inventory.tsv` | every function and data object, classified |
| `logs/package_audit/unused_code_commented_blocks.tsv` | runs of ≥5 commented-out code lines |
| `logs/package_carry/routing.tsv` | one row per object: destination, action, and why |
| `logs/package_carry/carry_manifest.tsv` | one row per carried object, with a content verdict |
| `logs/package_carry/test_results.tsv` | one row per test, both packages |
| `staging/test-packages/<Package>/` | the assembled test package |

## The checks

`audit/01_validate_structure.R` — DESCRIPTION placeholders and roxygen-version conflicts; `data/` ↔
man-page pairing in both directions; declared-but-unused dependencies; `pkg::` calls and NAMESPACE
imports that DESCRIPTION never declares; unguarded `Suggests`; exports resolving to nothing;
`.Rbuildignore` coverage. Exits non-zero on any ERROR.

`audit/02_check_cycles.R` — three graphs, because a cycle in each breaks something different. The
**package graph** (DESCRIPTION `Depends`/`Imports`/`LinkingTo`) is fatal: a cycle means neither package
installs into a clean library. The **build graph** (`data-raw` generators) is a bootstrap deadlock
rather than an install failure. The **call graph** is reported, not failed on — mutual recursion is
legal R. `Suggests` is deliberately not an edge; R permits a cycle through it precisely because the
dependency is optional.

Two checks matter more than the graphs themselves. A **de-facto cycle** is shipped `R/` code calling a
package DESCRIPTION never declares — including via `asNamespace("Pkg")` string lookups, which no
resolver and no `R CMD check` can see. If the declared graph already runs the other way, the metadata
reports acyclic while the code is not. And the **bootstrap** check reads the install and build graphs
as one relation, because each is separately acyclic here while their union is not.

`audit/03_report_unused_code.R` — classifies every function and data object as LIVE, DEAD, API-ONLY,
DEPRECATED, BUILD-ONLY, CALLED-UNEXPORTED or DYNAMIC-ACCESS, counting call sites in the package, its
tests and vignettes, and the consumer repos. Never a gate: dead code is a backlog, not a release
blocker, and a false DEAD would otherwise stop a good release.

## Why the checks parse rather than grep

Every usage signal comes from R's own parser (`utils::getParseData`), and each substitution fixed a
real false verdict:

- **`pkg::` detection** uses `SYMBOL_PACKAGE` tokens. Grepping `R/` for `pkg::` also matches comments
  and roxygen prose, which invented five dependencies that did not exist.
- **Usage counting** counts `SYMBOL` tokens and excludes only the definition *line*. Excluding the whole
  defining file hides the common case of a helper called by the exported function directly above it.
  The call graph alone is not enough either: `codetools` classifies a name as local once it is assigned,
  so a helper that is read then rebound (`environment(.f) <- ns`) vanishes from `findGlobals`.
- **Function-valued arguments** count as uses. `heatmap_color_fun = .feature_color_function` lands in
  `findGlobals`' `$variables`, not `$functions`; reading only `$functions` reports every callback as dead.
- **Roxygen inline code** counts as uses. `` @format `r .describe_data_format()` `` is executable code
  inside a comment, and 70+ doc stubs call it.
- **S3 methods** are entry points. `print.motrdat` is reached by dispatch and named only in
  `S3method(print,motrdat)`.
- **Boundaries are explicit lookaround**, not `\b`. In `pkg:::.find_ome` both `:` and `.` are non-word
  characters, so `\b` never matches and every dot-prefixed internal reads as unreferenced.
- **Vendored trees are excluded** — `staging/pkgs-*` and `.claude/worktrees/` hold verbatim copies of
  the packages under audit, and counting them makes every symbol look live.
