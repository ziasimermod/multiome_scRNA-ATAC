# Analysis design decisions

This file records choices that are easy to lose when a pipeline is represented only as code.

## 1. Work from Cell Ranger ARC outputs, not FASTQ files

The current R workflow begins after `cellranger-arc count`. R imports the paired filtered feature-barcode H5, ATAC fragment file, per-barcode metrics, and sample-specific peaks. FASTQ processing and alignment are upstream responsibilities.

## 2. Start with Cell Ranger's filtered matrix

The filtered matrix is a practical first-pass set of cell-associated barcodes. We retain ordinary low-quality cells for QC modeling and remove only near-empty barcodes before doublet detection. A later extension may evaluate raw-matrix ambient RNA/cell-calling methods, but it is not mixed into the initial teaching workflow.

## 3. QC is per sample, per modality, then joint

- **Per sample:** depth and quality distributions can differ among GEM wells/libraries.
- **Per modality:** RNA and ATAC metrics measure different failure modes.
- **Joint final decision:** both assays originated from the same nucleus, and downstream WNN relies on both.

This means “separate or together?” is not an either/or choice.

## 4. Thresholds are reviewed, not universal

Automated robust suggestions identify extreme distribution tails, but rare colon populations can have unusual RNA content or chromatin complexity. Each threshold is therefore a recorded, version-controlled decision. Thresholds should not be copied to a future experiment without reviewing that experiment's distributions.

## 5. Doublets are called within library before final QC

Doublets can have unusually high RNA and ATAC counts. `scDblFinder` is run on each library separately using RNA counts, before the final multi-metric filter. The doublet call is one component of the joint barcode decision.

## 6. Samples will share an ATAC feature space before merging

Cell Ranger calls peaks independently for each sample, so peak coordinates do not match exactly. The next milestone will reduce sample peak calls to a common peak set and re-quantify that same set in every sample before merging.

## 7. Do not integrate away the experimental contrast by default

The current two libraries are also the two diet groups. An aggressive sample correction could erase a real diet-associated shift because sample and diet cannot be separated statistically. Initial processing will merge shared features without Seurat integration by sample. Future replicated designs can evaluate batch correction using known technical covariates while preserving diet.

## 8. WNN is a joint representation, not a third assay

RNA will be represented by PCA and ATAC by LSI. WNN learns cell-specific contributions from the two representations and builds a joint neighbor graph. RNA UMAP, ATAC UMAP, and WNN UMAP will all be retained for diagnostic comparison.

## 9. Clusters are hypotheses, not cell-type labels

Colon annotation will combine cluster markers, canonical molecular signatures, RNA expression patterns, ATAC support, and manual review. The repository will keep the annotation worksheet and uncertainty notes rather than silently replacing cluster numbers with labels.

## 10. Diet comparisons are descriptive at the current sample size

There is one biological sample per diet. Cell proportions, pseudobulk counts, and fold-change summaries can be generated as exploratory descriptions. P-values treating nuclei as independent diet replicates would be pseudoreplication and are excluded.

