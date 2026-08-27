# Mana Lab 10x Multiome analysis

This repository is a stepwise, teaching-oriented workflow for paired single-cell RNA and ATAC data generated with 10x Genomics Multiome and `cellranger-arc count`.

The project includes three mouse-colon cohorts, each with one pooled CON
library and one pooled HFD library:

| Cohort ID | Group | CON sample | HFD sample | Pool metadata |
|---|---|---|---|---|
| `ppar` | PPAR | `Ppar-CON-1` | `Ppar-HFD-1` | CON: 6 mice (2F/4M); HFD: 10 mice (2F/8M) |
| `wt` | WT | `VilB-CON-1` | `VilB-HFD-1` | pending wet-lab confirmation |
| `il17` | IL17 | `IL17-CON-1` | `IL17-HFD-1` | pending wet-lab confirmation |

Each Cell Ranger library contains colon epithelial cells pooled from multiple mice. Individual mouse identities are not recoverable after pooling. Therefore, pooled-library identity and diet condition are completely confounded.

Cell-level and cluster-level differences can be explored descriptively, but nuclei are **not** independent biological replicates for diet-level inference. This workflow does not treat cell-level P values as replicated evidence for HFD-versus-CON effects.

## Why the workflow is split into steps

Each numbered notebook has one job, reads documented inputs or completed checkpoints, and writes a new checkpoint. A later step can therefore be repeated without rerunning every expensive calculation.

| Step | Notebook | Main purpose | Current status |
|---:|---|---|---|
| 0 | `analysis/00_environment_and_inputs.Rmd` | Validate R, packages, sample metadata, reference annotation, and Cell Ranger ARC inputs | Ready |
| 1 | `analysis/01_import_and_calculate_qc.Rmd` | Import each library separately and calculate RNA/ATAC QC plus doublet scores | Ready |
| 2 | `analysis/02_review_qc_and_filter.Rmd` | Review QC distributions, record thresholds, and make one joint keep/remove decision per paired barcode | Ready |
| 3 | `analysis/03_build_common_atac_peak_set.Rmd` | Build a common GRCm39 ATAC peak space and requantify post-QC nuclei | Ready |
| 4 | `analysis/04_merge_and_reduce_dimensions.Rmd` | Merge libraries, perform RNA log-normalization/PCA and ATAC TF-IDF/LSI, and inspect modality-specific structure | Ready |
| 5 | `analysis/05_build_wnn_and_evaluate_clusters.Rmd` | Construct the WNN graph, compare clustering resolutions, and select a reviewed working resolution | Ready |
| 6 | `analysis/06_annotate_mouse_colon_cell_states.Rmd` | Assign and document mouse-colon epithelial cell states using RNA and ATAC evidence | PPAR complete; WT/IL17 require cohort-specific review |
| 7 | composition/pseudobulk workflow | Descriptive composition summaries and preparation of pseudobulk-ready downstream matrices | Planned |

The PPAR development run is complete through Step 6. Steps 00-05 are shared
across cohorts; Step 06 remains guarded against reusing PPAR cluster numbers
for WT or IL17 before their annotations are independently reviewed.

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

The current project-specific data location is configured in `config/project_config.R`:

```text
/scratch/dsaiz/Yesenia_scData2026
```

Input sample metadata live in `config/datasets/<dataset_version>/samples.csv`. Sample-specific paths should not be hard-coded inside individual analysis notebooks.

Reviewed QC decisions are cohort-specific:

```text
config/datasets/<dataset_version>/cohorts/<cohort_id>/qc_thresholds.csv
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
```

Future annotation and downstream results will be written under `06_annotation` and `07_*` directories.

See `docs/WORKFLOW_AND_CHECKPOINTS.md` for the exact input/output contract for each step.

Known SOL package-library and checkpoint failures are collected in `docs/TROUBLESHOOTING.md`.

## Primary method references

The analysis design follows official 10x Genomics, Seurat, Signac, and Bioconductor documentation. See `docs/REFERENCES.md` for the specific sources and which analysis decisions they support.
