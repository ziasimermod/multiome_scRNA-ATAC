------------------------------------------------------------------------

editor_options: markdown: wrap: 72 ---

# Reviewing a Candidate Biological Offshoot in WNN Space

## Purpose

This guide describes how to investigate a cluster that appears to be a real, stable biological offshoot during weighted nearest-neighbor (WNN) analysis but is strongly enriched for one sample or experimental condition.

The goal is to distinguish among four possibilities:

1.  a coherent cell type or cell state;
2.  a genuine marker-positive subpopulation embedded within a broader cluster;
3.  a low-complexity or mixed cluster that should receive a conservative label;
4.  a technical structure caused by poor-quality nuclei, doublets, or one library.

A WNN cluster is evidence of neighborhood structure. It is not, by itself, evidence of a distinct cell identity or a treatment effect.

## When this review is warranted

Initiate a focused offshoot review when one or more of the following are true:

- a cluster remains spatially distinct across several clustering resolutions;
- a cluster is almost entirely populated by one library or condition;
- increasing resolution divides a broad parent population into plausible biological children;
- one modality supports the cluster more strongly than the other;
- the marker list is unusually sparse, conflicting, or dominated by mitochondrial, ribosomal, or stress-associated genes;
- a rare-cell marker is strongly enriched but detected in only a minority of the cluster.

## Evidence hierarchy

Review the candidate in the following order. This prevents marker enthusiasm from overriding sample structure or quality problems.

### 1. Establish cross-resolution persistence

Determine whether the population exists at several resolutions and identify its parent and child populations.

``` r
resolution_crosswalk <- merged[[]] |>
  dplyr::count(
    parent_cluster = .data[["wnn_res_0.2"]],
    child_cluster = .data[["wnn_res_0.8"]],
    name = "nuclei"
  ) |>
  dplyr::group_by(parent_cluster) |>
  dplyr::mutate(
    parent_total = sum(nuclei),
    percent_of_parent = 100 * nuclei / parent_total
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(parent_cluster, dplyr::desc(nuclei))
```

Interpretation:

- A stable child that receives nearly all its nuclei from one parent is unlikely to be a numbering artifact.
- A parent that separates into two marker-coherent children supports using a resolution high enough to preserve that distinction.
- A tiny child with no new marker program may represent overclustering.

### 2. Quantify library and condition composition

Report both the percentage of the cluster contributed by each sample and the percentage of each sample assigned to the cluster.

``` r
cluster_composition <- merged[[]] |>
  dplyr::count(
    sample_id,
    cluster = .data[["wnn_res_0.8"]],
    name = "nuclei"
  ) |>
  dplyr::group_by(cluster) |>
  dplyr::mutate(percent_of_cluster = 100 * nuclei / sum(nuclei)) |>
  dplyr::ungroup() |>
  dplyr::group_by(sample_id) |>
  dplyr::mutate(percent_of_sample = 100 * nuclei / sum(nuclei)) |>
  dplyr::ungroup()
```

With one pooled library per diet, use terms such as **CON-associated**, **HFD-associated**, or **library-skewed**. Do not call the observation a diet effect, depletion, or expansion without biological replication.

### 3. Compare quality within each sample

Compare the candidate with the remaining nuclei from the same library. Pooling samples at this stage can make a sample-quality difference look like a cell-state difference.

At minimum, summarize:

- RNA counts and features;
- ATAC counts;
- mitochondrial percentage;
- TSS enrichment;
- nucleosome signal;
- FRiP;
- continuous doublet score and categorical doublet call;
- learned RNA and ATAC weights.

Passing the original QC gates does not guarantee that a cluster has the same information content as the rest of the dataset. A valid rare population can still have lower complexity, but this should reduce annotation confidence.

### 4. Require modality-aware support

Inspect the candidate on the RNA, ATAC, and WNN embeddings. In this workflow, the reduction names are:

``` r
c("umap.rna", "umap.atac", "wnn.umap")
```

Interpretation:

- RNA, ATAC, and WNN separation supports a multimodally distinct state.
- RNA and WNN separation with ATAC dispersion supports an RNA-defined state, but not yet a distinct chromatin state.
- WNN-only separation warrants examination of modality weights and neighbor construction.
- Broad dispersion in all embeddings argues against a coherent cluster.

### 5. Use two marker contrasts

First calculate conventional one-versus-rest markers. Then compare the candidate directly with its closest cross-resolution sibling within the same sample.

``` r
sample_review <- subset(
  merged,
  subset = sample_id == "SAMPLE-ID"
)

Idents(sample_review) <- sample_review[[]][["wnn_res_0.8"]]

sibling_markers <- FindMarkers(
  sample_review,
  ident.1 = "CANDIDATE",
  ident.2 = "SIBLING",
  assay = "RNA",
  only.pos = FALSE,
  min.pct = 0.05,
  logfc.threshold = 0.10
) |>
  tibble::rownames_to_column("gene") |>
  dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))
```

The sibling contrast is often more informative than one-versus-rest because it asks what caused a specific parent population to divide.

### 6. Test complete marker programs and co-detection

Do not assign an identity from one highly significant gene. For a proposed cell type, audit a biologically coherent program and determine whether its genes occur in the same nuclei.

For a tuft-cell hypothesis, an initial review panel can include:

``` r
tuft_genes_requested <- c(
  "Dclk1", "St18", "Pou2f3", "Trpm5", "Gnat3", "Plcb2",
  "Rgs13", "Sh2d6", "Fyb", "Ltc4s", "Avil", "Hck"
)

tuft_genes <- intersect(
  tuft_genes_requested,
  rownames(merged[["RNA"]])
)
```

Extract normalized expression and retain cell barcodes explicitly:

``` r
program_expression <- SeuratObject::FetchData(
  merged,
  vars = tuft_genes,
  layer = "data"
)

program_expression$barcode <- rownames(program_expression)

program_metadata <- merged[[]] |>
  tibble::rownames_to_column("barcode") |>
  dplyr::transmute(
    barcode,
    sample_id,
    diet,
    cluster = as.character(.data[["wnn_res_0.8"]])
  )

program_review <- dplyr::left_join(
  program_metadata,
  program_expression,
  by = "barcode"
)

program_review$program_genes_detected <- rowSums(
  program_review[, tuft_genes, drop = FALSE] > 0
)
```

Summarize detection by sample and cluster:

``` r
program_detection_summary <- program_review |>
  dplyr::group_by(sample_id, cluster) |>
  dplyr::summarise(
    nuclei = dplyr::n(),
    median_program_genes_detected = median(program_genes_detected),
    percent_two_or_more = 100 * mean(program_genes_detected >= 2),
    dplyr::across(
      dplyr::all_of(tuft_genes),
      ~ 100 * mean(.x > 0),
      .names = "{.col}_percent"
    ),
    .groups = "drop"
  )
```

High fold change in a small minority of nuclei can mean that a real rare-cell core is embedded within a larger low-complexity cluster. Co-detection separates that possibility from uniform but sparse dropout across a coherent population.

### 7. Compare marker-positive and marker-negative cells within the candidate

If most cells do not express the proposed program, compare the program-positive core with the remainder of the candidate cluster.

``` r
candidate_review <- program_review |>
  dplyr::filter(
    sample_id == "SAMPLE-ID",
    cluster == "CANDIDATE"
  ) |>
  dplyr::mutate(
    program_group = dplyr::if_else(
      program_genes_detected >= 2,
      "program_positive",
      "program_negative"
    )
  )

table(candidate_review$program_group)
```

Map these groups back onto all three embeddings. If the positive nuclei form a localized core while the negative nuclei form a diffuse halo, use a state-enriched or mixed-cluster label and consider focused subclustering. If positive and negative cells are spatially intermixed, dropout and low RNA complexity may be the main explanation.

### 8. Seek independent ATAC support

For an RNA-defined candidate, inspect:

- gene activity at defining loci;
- accessibility around defining promoters and enhancers;
- relevant transcription-factor motifs;
- aggregate coverage after stratifying by sample and candidate status.

Sparse RNA expression plus coherent accessibility can rescue confidence. Sparse RNA expression plus dispersed ATAC structure should remain a cautious annotation.

## Annotation outcomes

Use one of the following outcomes rather than forcing every cluster into a fully resolved identity.

| Outcome | Evidence | Recommended annotation |
|----|----|----|
| Coherent cell state | Stable structure, broad program, acceptable QC, preferably multimodal support | Specific state; moderate or high confidence |
| State-enriched cluster | Real marker-positive core, but most cluster members lack the full program | `State-enriched epithelial` or `State-like epithelial`; low/moderate confidence |
| Low-complexity population | Stable cluster but weak program and systematically low information content | `Low-complexity epithelial`; low confidence |
| Mixed population | Incompatible programs or distinct internal marker-positive groups | Conservative mixed label; subcluster before final interpretation |
| Technical structure | Poor QC, strong doublet enrichment, or library-specific technical separation | Exclude or flag; do not assign biological identity |

## IL17 cluster 11 case study

### Structural evidence

At resolution 0.2, cluster 7 contained 544 nuclei. At resolution 0.8 it divided almost completely into:

| Resolution 0.8 child | Nuclei inherited from cluster 7 |
|----------------------|--------------------------------:|
| Cluster 11           |                             397 |
| Cluster 17           |                             146 |
| Cluster 22           |                               1 |

Only three cluster-11 nuclei originated outside the resolution-0.2 parent. Cluster 11 is therefore a persistent biological or technical structure rather than a high-resolution numbering artifact.

### Sample association

At resolution 0.8, cluster 11 contains 390 IL17-CON-1 nuclei and 10 IL17-HFD-1 nuclei.

- 5.56% of IL17-CON-1 nuclei belong to cluster 11.
- 0.21% of IL17-HFD-1 nuclei belong to cluster 11.
- The observed prevalence is approximately 26-fold higher in the CON library.

This is a strong **CON-associated/library-associated** observation. It is not a replicated diet effect because each diet is represented by one pooled library.

### Quality and modality evidence

The 390 CON cluster-11 nuclei pass the established QC gates but have lower information content than other CON nuclei:

| Metric                          | CON cluster 11 | Other CON clusters |
|---------------------------------|---------------:|-------------------:|
| Median RNA counts               |            452 |              1,221 |
| Median RNA features             |            400 |                875 |
| Median ATAC counts              |          1,406 |              3,398 |
| Median mitochondrial percentage |           2.50 |              0.238 |
| Median TSS enrichment           |           7.80 |               8.26 |
| Median FRiP                     |          0.623 |              0.742 |
| Median scDblFinder score        |       0.000559 |            0.00324 |
| Median RNA weight               |          0.614 |              0.548 |

This population is lower-complexity but is not explained by elevated doublet scores in the CON library. Its higher RNA weight and its compact RNA/WNN position, together with dispersion in ATAC space, indicate that the separation is primarily RNA-supported.

The ten HFD cluster-11 nuclei have a median scDblFinder score of 0.240. They are too few and too doublet-score-enriched to define the population, even if their categorical QC call was singlet.

### Marker-program evidence

The initial one-versus-rest analysis found two strong positive markers:

| Gene | Average log2 fold change | Detection in cluster 11 | Detection outside cluster 11 | Adjusted P value |
|----|---:|---:|---:|---:|
| `Dclk1` | 7.33 | 12.6% | 0.5% | 1.91e-96 |
| `St18` | 5.83 | 13.1% | 1.7% | 2.35e-43 |

These genes support a tuft-like hypothesis, but co-detection shows that the program is confined to a minority of the CON cluster:

| IL17-CON-1 cluster 11 pattern | Nuclei | Percent |
|-------------------------------|-------:|--------:|
| Neither `Dclk1` nor `St18`    |    325 |   83.3% |
| `St18` only                   |     16 |   4.10% |
| `Dclk1` only                  |     14 |   3.59% |
| Both genes                    |     35 |   8.97% |

Only 11.8% of the CON cluster detected at least two genes from the reviewed tuft panel, and the median number of detected tuft genes was zero. The current evidence therefore does not support labeling all 390 CON nuclei as definitive tuft cells.

In contrast, nine of the ten HFD cluster-11 nuclei co-detected `Dclk1` and `St18`, and those cells detected a broader tuft program. This intriguing result must be interpreted alongside their small number and elevated continuous doublet scores.

### Sibling-cluster interpretation

Cluster 17 is not simply a second tuft-like cluster. Within IL17-CON-1, cluster 17 strongly expresses a neuroendocrine/enterochromaffin program including `Chga`, `Chgb`, `Tph1`, `Ddc`, `Slc18a1`, `Rfx6`, `Lmx1a`, `Snap25`, and multiple neuronal secretion genes.

Although `St18` is detected in 51.3% of CON cluster 17, `Dclk1` and the broader reviewed tuft program are absent or nearly absent. In this dataset, `St18` alone is therefore not tuft-specific. The cross-resolution split separates an enteroendocrine/enterochromaffin child from the candidate cluster-11 population.

### Current conclusion

The most defensible working interpretation is:

> Cluster 11 is a stable, CON-associated, RNA-defined, low-complexity epithelial offshoot containing a `Dclk1`/`St18` tuft-like core. The tuft program is not sufficiently pervasive to assign a definitive tuft-cell label to the entire cluster.

Recommended provisional annotation:

``` text
cell_state: Tuft-like enriched
cell_state_full: Dclk1+/St18+ tuft-like-enriched low-complexity epithelial state
annotation_confidence: low
```

This label should remain provisional until the marker-positive core is mapped, compared with the marker-negative remainder, and tested for independent ATAC support.

## Next decision gate for IL17 cluster 11

Before approving the IL17 clustering resolution or completing annotation:

1.  Map nuclei with at least two tuft-panel genes onto RNA, ATAC, and WNN UMAPs.
2.  Compare program-positive and program-negative cluster-11 nuclei within IL17-CON-1.
3.  Compare the program-positive core directly with enteroendocrine cluster 17.
4.  Inspect accessibility or gene activity at `Pou2f3`, `Trpm5`, `Dclk1`, `Rgs13`, and `St18`.
5.  Determine whether focused subclustering cleanly recovers the approximately 46 program-positive CON nuclei without creating unstable tiny groups.
6.  Retain sample-associated language until biological replicates are added.

If the tuft-positive cells form a localized and multimodally supported core, annotate that core separately through focused subclustering. If they remain intermixed with marker-negative low-complexity nuclei, retain the conservative cluster-level annotation and treat tuft abundance as a marker-program result rather than a cluster-count result.

## Suggested review outputs

Write focused diagnostics under a cluster-specific directory rather than mixing them with the primary Step 5 plots:

``` text
05_wnn/offshoot_review/cluster_11/
├── cross_resolution_membership.csv
├── sample_composition.csv
├── sample_stratified_qc_summary.csv
├── candidate_vs_sibling_RNA_markers.csv
├── marker_program_detection_summary.csv
├── marker_co_detection.csv
├── RNA_ATAC_WNN_support.pdf
└── focused_marker_programs.pdf
```

These files are descriptive review artifacts. They should not be treated as replicated differential-abundance or diet-response results.

## Reproducibility checklist

Before finalizing the cluster decision, record:

- dataset version and cohort ID;
- WNN dimensions and clustering resolution;
- cluster column;
- parent and sibling clusters;
- sample-stratified counts and percentages;
- QC and doublet-score summaries;
- RNA, ATAC, and WNN support;
- marker-test contrast and sample restriction;
- requested and detected marker-program genes;
- co-detection rule used to define program-positive nuclei;
- final annotation, confidence, reviewer, date, and decision notes;
- explicit statement that causal diet inference was not performed without replication.

## Selected biological references

- Westphalen CB et al. Long-lived intestinal tuft cells serve as colon cancer-initiating cells. *Journal of Clinical Investigation* (2014). <https://doi.org/10.1172/JCI73434>
- Intestinal Tuft-2 antimicrobial-immunity study. *Immunity* (2022). <https://doi.org/10.1016/j.immuni.2022.03.001>
