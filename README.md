# Mana Lab 10x Multiome analysis

This repository is a stepwise, teaching-oriented workflow for paired single-cell RNA and ATAC data generated with 10x Genomics Multiome and `cellranger-arc count`.

The project includes three mouse-colon cohorts, each with one pooled CON
library and one pooled HFD library:

| Cohort ID | Group | CON sample | HFD sample | Pool metadata |
|---|---|---|---|---|
| `ppar` | PPAR | `Ppar-CON-1` | `Ppar-HFD-1` | CON: 6 mice (2F/4M); HFD: 10 mice (2F/8M) |
| `wt` | WT | `VilB-CON-1` | `VilB-HFD-1` | CON: 7 mice (4F/3M); HFD: 6 mice (3F/3M) |
| `il17` | IL17 | `IL17-CON-1` | `IL17-HFD-1` | CON: 4 mice (2F/2M); HFD: 4 mice (3F/1M) |

Each Cell Ranger library contains colon epithelial cells pooled from multiple mice. Individual mouse identities are not recoverable after pooling. Therefore, pooled-library identity and diet condition are completely confounded.

Cell-level and cluster-level differences can be explored descriptively, but nuclei are **not** independent biological replicates for diet-level inference. This workflow does not treat cell-level P values as replicated evidence for HFD-versus-CON effects.

## Teaching workflow and reference provenance

Start with [Yesenia's handoff guide](docs/YESENIA_HANDOFF.md).
The main independent workflow is analysis/00-07, followed by the whole-WNN-cluster
07D companion. 05B and the Gourab/refinement/trajectory/08-series notebooks are
optional analyses with their own prerequisites. Presentation figures are support
outputs and do not define an analysis step 09.

Code/configuration presence is verified by cleanup. Analysis execution status must
be checked from the saved run records, not inferred from a filename or this README.

## The short answer to "Is Multiome QC separate or together?"

Both:

1. Review each Cell Ranger library separately because capture quality and sequencing depth can differ between libraries.
2. Evaluate RNA-specific and ATAC-specific metrics separately because they measure different failure modes.
3. Make one final joint decision for each barcode because RNA and ATAC originate from the same nucleus and downstream multimodal analysis requires both measurements to be usable.

This is encoded in Step 2 as separate `qc_pass_rna`, `qc_pass_atac`, and `qc_pass_doublet` decisions followed by one final `qc_pass` decision.

## Current analysis strategy

After QC, the two libraries within the selected cohort are placed into a shared
feature space without batch integration because library identity and diet are
completely confounded.

RNA is analyzed using log normalization, highly variable genes, PCA, and an RNA-only exploratory UMAP.

ATAC is analyzed using a common peak set, TF-IDF normalization, SVD/LSI, and an ATAC-only exploratory UMAP. LSI dimension 1 is excluded from downstream neighborhood construction in this dataset because it is almost perfectly associated with ATAC sequencing depth.

RNA PCs and ATAC LSI dimensions are then combined using Seurat weighted nearest neighbors (WNN). A clustering-resolution grid is reviewed manually before selecting the working clustering resolution.

## Start here on ASU SOL

1. Clone this repository into your working directory on SOL.
2. Start an Open OnDemand **RStudio Server** session using **R 4.4.2**. Do not request a GPU.
3. Open `multiome_scRNA-ATAC.Rproj` in RStudio.
4. Read `docs/SOL_RSTUDIO_SETUP.md` before selecting CPU, memory, and wall time.
5. If this is the first setup of the R 4.4 library, run `setup/install_packages.R` once.
6. Set `MULTIOME_DATASET_VERSION` and `MULTIOME_COHORT_ID` before sourcing the project configuration. The defaults are `v2_resequenced` and `ppar`.
7. Open the numbered notebooks under `analysis/` and run them in order, one chunk at a time.
8. Restart R between major notebooks when appropriate to release memory.

For example, select the WT resequenced cohort in a fresh R session with:

```r
Sys.setenv(
  MULTIOME_DATASET_VERSION = "v2_resequenced",
  MULTIOME_COHORT_ID = "wt"
)

Sys.getenv(c(
  "MULTIOME_DATASET_VERSION",
  "MULTIOME_COHORT_ID"
))
```

Valid v2 cohort IDs are `ppar`, `wt`, and `il17`; v1 contains only `ppar`.
Environment selections are session-specific, so set and confirm them again
after restarting R.

The current project-specific data location is configured in `config/project_config.R`:

```text
/scratch/dsaiz/Yesenia_scData2026
```

Input sample metadata live in `config/datasets/<dataset_version>/samples.csv`. Sample-specific paths should not be hard-coded inside individual analysis notebooks.

Reviewed QC decisions are cohort-specific:

```text
config/datasets/<dataset_version>/cohorts/<cohort_id>/qc_thresholds.csv
```

Reviewed clustering and annotation decisions use the same cohort directory:

```text
config/datasets/<dataset_version>/cohorts/<cohort_id>/clustering_decision.csv
config/datasets/<dataset_version>/cohorts/<cohort_id>/cell_state_annotations.csv
```

## What belongs in GitHub

Commit code, documentation, the sample sheet, reviewed QC decisions, and future annotation records.

Do **not** commit FASTQ files, Cell Ranger outputs, fragment files, BAM files, H5 matrices, large RDS checkpoints, generated result directories, credentials, access tokens, or other sensitive information.

The `.gitignore` file provides guardrails, but always inspect `git status` before staging or committing changes.

## Checkpoints and results

Large generated files are written outside this repository to:

```text
/scratch/dsaiz/Yesenia_scData2026/results/<dataset_version>/independent/<cohort_id>
```

Current major checkpoint directories include:

```text
00_run_info
01_qc
02_objects
03_common_peaks
03_objects
04_reduction
05_wnn
06_annotation
07_diet_response
07b_trajectories
07c_candidate_prioritization
```

Annotation and downstream results are written under `06_annotation` and
`07_*` directories.

See `docs/WORKFLOW_AND_CHECKPOINTS.md` for the exact input/output contract for each step.

Known SOL package-library and checkpoint failures are collected in `docs/TROUBLESHOOTING.md`.

## Primary method references

The analysis design follows official 10x Genomics, Seurat, Signac, and Bioconductor documentation. See `docs/REFERENCES.md` for the specific sources and which analysis decisions they support.
