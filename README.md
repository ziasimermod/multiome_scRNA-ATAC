# Mana Lab 10x Multiome analysis

This repository is a stepwise, teaching-oriented workflow for paired single-cell RNA and ATAC data generated with 10x Genomics Multiome and `cellranger-arc count`.

The first project-specific implementation analyzes two mouse-colon PPARalpha/PPARdelta iKO libraries:

| Sample | Diet | Current biological replicate count |
|---|---|---:|
| `Ppar-CON-1` | control diet (CON) | 1 |
| `Ppar-HFD-1` | high-fat diet (HFD) | 1 |

Because there is currently one biological sample per diet, sample and diet are completely confounded. Cell-level differences can be explored descriptively, but cells are **not** independent diet replicates and this repository does not perform a formal HFD-versus-CON significance test at this stage.

## Why the workflow is split into steps

Each numbered notebook has one job, reads documented inputs, and writes a checkpoint. A later step can therefore be repeated without rerunning every expensive calculation.

| Step | Notebook | Main purpose | Current status |
|---:|---|---|---|
| 0 | `analysis/00_environment_and_inputs.Rmd` | Check R, packages, sample metadata, and Cell Ranger ARC files | Ready |
| 1 | `analysis/01_import_and_calculate_qc.Rmd` | Import each library separately and calculate RNA/ATAC QC plus doublet scores | Ready |
| 2 | `analysis/02_review_qc_and_filter.Rmd` | Review distributions, record thresholds, and make one joint keep/remove decision per paired barcode | Ready |
| 3 | common-peak notebook | Build and quantify a common ATAC peak set | Planned |
| 4 | RNA/ATAC processing notebook | SCTransform/PCA, TF-IDF/LSI, and WNN | Planned |
| 5 | cluster-review notebook | Compare clustering resolutions and calculate markers | Planned |
| 6 | annotation notebook | Assign and document mouse-colon cell populations | Planned |
| 7 | composition notebook | Descriptive composition and pseudobulk preparation | Planned |

The repository intentionally begins with Steps 0–2. Those steps establish the file contract and the most important manual decision point before downstream analyses are added.

## The short answer to “Is Multiome QC separate or together?”

Both:

1. Review each Cell Ranger library separately because capture quality and sequencing depth can differ between libraries.
2. Evaluate RNA-specific and ATAC-specific metrics separately because they measure different failure modes.
3. Make one final joint decision for each barcode because RNA and ATAC came from the same nucleus and WNN requires both measurements to be useful.

This is encoded in Step 2 as separate `qc_pass_rna`, `qc_pass_atac`, and `qc_pass_doublet` columns followed by one `qc_pass` column.

## Start here on ASU SOL

1. Clone this repository into the project directory on SOL.
2. Start an Open OnDemand **RStudio Server** session using **R 4.4.2**. Do not request a GPU.
3. Open `multiome_scRNA-ATAC.Rproj` in RStudio.
4. Read `docs/SOL_RSTUDIO_SETUP.md` before selecting CPU, memory, and wall time.
5. If this is the first setup of the R 4.4 library, run `setup/install_packages.R` once.
6. Open the numbered notebooks under `analysis/` and run them in order, one chunk at a time.
7. Restart R between notebooks to release memory.

The current data location is configured explicitly in `config/project_config.R`:

```text
/scratch/dsaiz/Yesenia_scData2026
```

Input sample metadata live in `config/samples.csv`. Do not hard-code additional sample paths inside analysis notebooks.

## What belongs in GitHub

Commit code, documentation, the sample sheet, and reviewed analysis decisions. Do **not** commit FASTQ files, Cell Ranger outputs, fragment files, BAM files, large RDS checkpoints, or generated result directories. The `.gitignore` file provides guardrails, but always inspect `git status` before committing.

## Checkpoints and results

Generated files are written outside this repository to:

```text
/scratch/dsaiz/Yesenia_scData2026/results/Ppar_iKO_multiome_stepwise
```

See `docs/WORKFLOW_AND_CHECKPOINTS.md` for the exact input/output contract for each step.

Known SOL package-library and checkpoint failures are collected in `docs/TROUBLESHOOTING.md`.

## Primary method references

The analysis design follows official 10x Genomics, Seurat, Signac, and Bioconductor documentation. See `docs/REFERENCES.md` for the specific sources and which decisions they support.
