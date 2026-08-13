# Primary references

These are the main official or primary method sources used to structure the workflow. Links were reviewed on 2026-08-13.

## 10x Genomics Cell Ranger ARC

- [Understanding Cell Ranger ARC outputs](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/outputs/understanding-output) — identifies the feature-barcode matrices, ATAC peaks, fragment files, and per-barcode outputs produced by `cellranger-arc count`.
- [Feature-barcode matrices](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/outputs/feature-barcode-matrices) — explains that the H5 contains both `Gene Expression` and `Peaks` feature types.
- [Per-barcode QC metrics](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/outputs/per-barcode-qc-metrics) — documents the paired RNA/ATAC barcode metrics used for the Cell Ranger FRiP calculation and library review.
- [Analyzing GEX and ATAC libraries with Cell Ranger ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/running-pipelines/single-library-analysis) — describes paired GEX/ATAC processing for one Multiome library.

## Seurat and Signac

- [Joint RNA and ATAC analysis: 10x Multiome](https://stuartlab.org/signac/articles/pbmc_multiomic) — primary implementation model for storing paired RNA and ATAC assays and performing joint downstream analysis.
- [Analyzing scATAC-seq data with Signac](https://stuartlab.org/signac/articles/pbmc_vignette.html) — source for TSS enrichment, nucleosome signal, and common scATAC quality-control concepts.
- [Merging Signac objects](https://stuartlab.org/signac/articles/merging) — supports creating a common peak set and quantifying it in each experiment before merging.
- [Parallel and distributed processing in Signac](https://stuartlab.org/signac/articles/future) — documents parallelization behavior for expensive operations such as `FeatureMatrix()`.
- [Weighted nearest-neighbor analysis](https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis) — basis for combining modality-specific representations into a joint neighbor graph, UMAP, and clustering.
- [Using Seurat with multimodal data](https://satijalab.org/seurat/articles/multimodal_vignette) — overview of multimodal object structure and WNN concepts.

## Doublet detection

- [scDblFinder Bioconductor vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html) — method and recommendations for doublet detection, including expected rates and per-capture handling.

## Statistical unit and pseudobulk analysis

- Crowell HL et al. [muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data](https://doi.org/10.1038/s41467-020-19894-4). *Nature Communications* (2020) — supports replicate-aware differential-state reasoning rather than treating cells as independent experimental replicates.
- Squair JW et al. [Confronting false discoveries in single-cell differential expression](https://doi.org/10.1038/s41467-021-25960-2). *Nature Communications* (2021) — demonstrates inflated false discoveries from cell-level methods that ignore biological replication and supports pseudobulk approaches.

References do not provide universal cutoff values. QC cutoffs remain experiment-specific decisions documented in `config/qc_thresholds.csv`.

