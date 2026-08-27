# Troubleshooting on SOL

Work from the first complete error message. Later “object not found” or recursive wrap-up errors are often consequences of the first failure rather than independent problems.

## A package is listed as installed but will not load

`installed.packages()` confirms only that a package directory is registered. It does not prove that every imported namespace can load. Use:

```r
pkgs <- c("Seurat", "Signac", "scDblFinder")

for (p in pkgs) {
  tryCatch(
    {
      loadNamespace(p)
      cat(p, as.character(packageVersion(p)), "OK\n")
    },
    error = function(e) {
      cat(p, "FAILED:", conditionMessage(e), "\n")
    }
  )
}
```

Run this in the same R version and `R_LIBS_USER` as the RStudio analysis session. Step 0 performs a broader version of this check.

## RStudio works but an Rscript/Slurm job does not

The two sessions are probably using different R versions or library search paths. Compare:

```r
R.version.string
Sys.getenv("R_LIBS_USER")
.libPaths()
Sys.which("R")
```

For the current project, select R 4.4.2 in Open OnDemand. Packages installed under `/home/dsaiz/R/x86_64-pc-linux-gnu-library/4.4` will not automatically appear under R 4.6.x.

## Missing dependencies appear one at a time

This usually means the user library contains a partially completed dependency tree. A top-level package directory can exist even though an imported package such as `RANN`, `spam`, `cluster`, `igraph`, `fastmatch`, or `biovizBase` is missing.

Use `setup/install_packages.R` from a clean R 4.4.2 session. It tests namespace loading, repairs known top-level and transitive failures, and then repeats the full preflight.

## A `00LOCK-*` directory prevents installation

R creates a lock directory while installing a package. A killed or interrupted installation can leave it behind.

1. Confirm no R process is currently installing into the same user library.
2. Identify the exact stale lock directory reported in the error.
3. Remove only that exact `00LOCK-<package>` directory.
4. Install the failed package again in a clean R session.

Do not delete the whole R library.

## The notebook stops because `approved=FALSE`

That is an analysis gate, not a code error. Step 2 refuses to filter cells until every sample in the active cohort has complete numeric thresholds and `approved=TRUE` in `config/datasets/<dataset_version>/cohorts/<cohort_id>/qc_thresholds.csv`. Review the plots, document the decision, save the CSV, and rerun the approval chunk.

## An RStudio session runs out of memory

1. Check whether the current sample checkpoint was written completely.
2. Start a new RStudio session with the memory recommendation in `docs/SOL_RSTUDIO_SETUP.md`.
3. Open the RStudio Project and rerun the notebook from the top.
4. Leave checkpoint reuse enabled so completed samples are loaded rather than recalculated.

Avoid holding two large Seurat objects in the Environment pane. Restart R between steps.

## `object 'x' not found` appears during error wrap-up

Find the first earlier error. For example, if common-peak parsing fails before `peak_sets` is created, a later cleanup or summary line can also report `object 'peak_sets' not found`. Fixing the first error prevents the downstream object-not-found message.

## A checkpoint exists but lacks required metadata

The code checks expected assays and metadata before reuse. If it reports that an old checkpoint lacks a newer field, set the notebook's reuse flag to `FALSE` for that stage and rebuild from the previous documented checkpoint. Do not manually add blank metadata merely to bypass the check.
