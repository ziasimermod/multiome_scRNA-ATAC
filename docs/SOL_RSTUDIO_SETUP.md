# Running the workflow in ASU SOL RStudio Server

The analysis code is interactive and teaching-oriented, but the computation still occurs on a SOL compute node. Resource requests are made before launching RStudio from Open OnDemand.

## Session fields

For the current Mana Lab account and public partition:

| Open OnDemand field | Project setting |
|---|---|
| Account | `grp_mmana` |
| Partition | `public` |
| QOS | `public` |
| R version | `4.4.2` |
| GPU resources | leave blank |
| Additional modules | leave blank unless a later step explicitly requires one |
| Additional sbatch options | leave blank for this workflow |

R package libraries are version-specific. Packages installed under R 4.4.2 are not automatically available to R 4.6.x. Always confirm `R.version.string`, `Sys.getenv("R_LIBS_USER")`, and `.libPaths()` in Step 0.

## Starting resource requests

These are conservative starting points for the current two libraries, not guarantees for future datasets. Resource use scales with cells, features, peaks, and parallel workers.

| Step | CPU cores | Memory | Wall time | Why |
|---:|---:|---:|---:|---|
| 0: environment and files | 2 | 16 GiB | 2 h | Package/file checks and mouse annotation |
| 1: import and QC metrics | 8 | 96 GiB | 12 h | H5 import, ATAC QC, and doublet detection |
| 2: review and filter | 4 | 64 GiB | 8 h | One large Seurat checkpoint loaded at a time |
| 3: common peak quantification | 8 | 128 GiB | 12–16 h | `FeatureMatrix()` over a shared peak set |
| 4: SCTransform, LSI, WNN | 8 | 192 GiB | 16–24 h | Current expected peak-memory step |
| 5–7: review/annotation/summaries | 4 | 64–96 GiB | 8–12 h | Depends on marker calculations and object size |

Steps 3–7 are roadmap estimates and will be finalized with their notebooks.

## Launch procedure

1. Log in to the SOL Open OnDemand portal.
2. Choose **Interactive Apps → RStudio Server**.
3. Enter the settings for the step you plan to run.
4. Launch and wait for the session to start.
5. In RStudio, use **File → Open Project** and select `multiome_scRNA-ATAC.Rproj` from the cloned repository.
6. Open the next numbered notebook under `analysis/`.
7. Run chunks in order. Read the text above a chunk before running it.
8. When the notebook is complete, use **Session → Restart R** before opening the next step.

## Where to clone the repository

A practical current location is:

```text
/scratch/dsaiz/Yesenia_scData2026/multiome_scRNA-ATAC
```

The analysis does not depend on that exact repository location. The data paths are currently explicit in `config/project_config.R` and `config/datasets/<dataset_version>/samples.csv`.

## Select the dataset and cohort

The project configuration reads two environment variables. Set them in the
RStudio Console before running a notebook's setup chunk:

```r
Sys.setenv(
  MULTIOME_DATASET_VERSION = "v2_resequenced",
  MULTIOME_COHORT_ID = "ppar"
)
```

Valid v2 cohort IDs are `ppar`, `wt`, and `il17`. The v1 initial-depth dataset
contains only `ppar`. If the variables are unset, the workflow defaults to the
v2 PPAR cohort. Restarting R clears console-only selections, so confirm the
active configuration table in Step 0 before creating outputs.

Each independent run writes beneath:

```text
/scratch/dsaiz/Yesenia_scData2026/results/<dataset_version>/independent/<cohort_id>/
```

## RStudio is not the storage location

RStudio holds objects in node memory only for the life of the session. Durable inputs and checkpoints remain on the SOL filesystem. `readRDS()` reads a checkpoint into memory; it does not modify the file. `saveRDS()` is the operation that creates or replaces a checkpoint.

## Session safety

- Do not open multiple copies of a heavy notebook in different RStudio sessions against the same output directory.
- Do not install packages while a large analysis object is loaded.
- Do not run package installation concurrently from two sessions; this can leave a `00LOCK-*` directory.
- Avoid `as.matrix()` on full sparse count matrices.
- Confirm that each numbered completion marker exists before proceeding.
- If a session ends, start a new one and reuse completed checkpoints rather than rerunning earlier samples.

## If the RStudio session disappears

First determine whether the computation created its expected checkpoint. A crashed browser tab does not necessarily mean the compute process failed, but an out-of-memory Slurm event normally ends the R process.

For Steps 1–2:

1. start a new R 4.4.2 RStudio session;
2. open the RStudio Project;
3. confirm completed checkpoint files under the configured `OUTPUT_DIR`;
4. leave checkpoint reuse enabled;
5. rerun the notebook from the top so it reconstructs the small in-memory variables safely.
