## Current interpretation status

Steps 8, 8B, and 8C preserve complementary analyses of the WT
persistent HFD-associated gene sets. Step 8 uses annotation unit 8 as a
strict ISC-like sensitivity analysis. Steps 8B and 8C use units 8 and 6
as a source-matched crypt composite and decompose the contribution of
each unit.

The `broad_stem` label used in these notebooks should be interpreted as
a source-matched crypt composite rather than a purified ISC population.
Unit 8 provides an ISC-like component, whereas unit 6 shows a
noncycling, absorptive-biased crypt-progenitor program. Unit 0 remains
under focused review because it may contain both cycling ISCs and early
transit-amplifying cells.

Primary biological interpretation of the WT stem compartment is paused
pending focused RNA, ATAC, and WNN review of global units 0, 6, and 8.
The existing analyses are retained as prespecified sensitivity analyses
and provenance records.

Step 8D is implemented but its motif results should be interpreted as
unit-8- and unit-6-specific regulatory analyses, not yet as a complete
analysis of the WT stem compartment.

# Step 8 persistent HFD analyses

This analysis section deliberately preserves a sequence of related but distinct WT audits.

## Why Step 8 and Step 8B are separate

`08_validate_persistent_hfd_programs.Rmd` is the strict provenance audit. It tests Gourab's persistent day-21/day-50 HFD-associated genes in the single WT annotation unit 8 (`ISC-like crypt`) and integrates those RNA effects with exact plus-or-minus 2-kb TSS accessibility.

`08b_refine_persistent_hfd_programs_in_matched_WT_populations.Rmd` does not replace Step 8. It recomputes the benchmark after aligning the WT annotation units more closely to the source annotation breadth. Its strict unit-8 result must numerically reproduce Step 8 before regrouped results are interpreted.

This separation prevents a scientifically motivated regrouping from erasing the original result or looking like a post hoc replacement.

`08c_decompose_wt_broad_stem_and_prioritize_regulatory_loci.Rmd` addresses a limitation exposed by Step 8B: the broad-stem composite contains very different proportions of unit 8 and unit 6 in CON and HFD. It partitions each linear RNA or ATAC difference into an average within-state term and a composition term, creates component-resolved regulatory peak sets, and reprioritizes exact-promoter candidates. It does not test motifs.

`08d_enrich_wt_broad_stem_regulatory_motifs.Rmd` is a manually gated sequence analysis. It tests only reviewed Step 8C peak sets against backgrounds matched for baseline accessibility, GC content, peak width, and genomic context. The motif decision is stored in `config/datasets/v2_resequenced/cohorts/wt/motif_peak_set_decision.csv`.

## Source-matched population map

| Analysis | WT units | Role | Rationale |
|-----------------|--------------------:|-----------------|-----------------|
| Broad stem | 8 + 6 | Primary | Matches a source stem group broad enough to include ISC-like and crypt-progenitor states. |
| Strict ISC-like | 8 | Sensitivity | Preserves the original Step 8 anchor. |
| Crypt progenitor | 6 | Component | Shows the noncycling progenitor contribution to broad stem. |
| Reg4+ DCS | 9 | Primary | Direct approved DCS-like comparison. |
| All goblet | 4 + 5 + 7 + 13 | Primary | Includes every approved immature and mature goblet unit. |
| Core goblet | 4 + 5 + 7 | Sensitivity | Tests dependence on CON-dominant unit 13. |
| Goblet components | 4, 5, 7, 13 separately | Component | Locates state-specific or composition-driven effects. |

Units 0 (`Cycling ISC/TA`) and 18 (`Cycling secretory progenitor`) remain important biological context, but they are not part of the source-matched composites. Including them would make proliferation a major part of the population definition and would answer a different question.

## Three concepts that must remain separate

- `annotation_confidence` records the consistency and specificity of marker evidence supporting a cell-state label. It is a reviewed categorical judgment (`high`, `moderate`, or `low`), not a percent match.
- `cell_support_tier` records the minimum nuclei available in either diet: Tier A has at least 200 per diet, Tier B at least 100, and Tier C fewer than
  100. It measures sampling support, not biological truth.
- `analysis_role` records why a population is present: `source_matched_primary`, `sensitivity`, or `component`.

The Step 7 labels `secondary_compartment` and `exclude_primary` are also scope rules rather than confidence percentages. `secondary_compartment` keeps non-epithelial populations available for their own analysis while excluding them from the epithelial benchmark. `exclude_primary` retains low-confidence populations for sensitivity review without allowing them to drive primary RNA/ATAC conclusions.

## Run order

1.  Complete and approve Steps 6 and 7.
2.  Run Step 8 and confirm `STEP_08_COMPLETE.txt`.
3.  Run Step 8B in a fresh WT v2 R session.
4.  Review the Step 8B completion marker, primary summary, candidate tables, sensitivity tables, and PDFs.
5.  Run Step 8C and review the observed-versus-standardized RNA effects, component-resolved ATAC classes, motif-ready peak-set sizes, and exact promoter priorities.
6.  Complete `motif_peak_set_decision.csv` only for peak sets that pass that review; every approved row requires reviewer and review date.
7.  Run Step 8D, inspect background matching before motif enrichment, and only then nominate TF families for follow-up.
8.  Carry the finalized observed and composition-standardized WT benchmarks into descriptive PPAR and IL17 genotype-by-diet comparisons.

Step 8B requires two analysis packages beyond the shared workflow preflight. The current `msigdbr` release also uses the R package `curl` for a one-time, versioned MSigDB database download:

``` r
install.packages(c("curl", "msigdbr"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("fgsea")
```

Restart R after installing or updating packages so the notebook does not mix old and new package namespaces.

Before running the MSigDB chunk, this should return a curl handle without an error:

``` r
curl::new_handle(timeout = 600)
```

An error containing `length(keys) == length(values)` indicates a broken or mixed `curl` R package/DLL installation. Reinstall `curl` into the active R library and restart R; changing the requested MSigDB collection cannot repair that download-layer failure.

Step 8D additionally requires:

``` r
BiocManager::install(
  c(
    "JASPAR2024",
    "TFBSTools",
    "motifmatchr",
    "Biostrings",
    "Rsamtools",
    "DBI",
    "RSQLite"
  ),
  ask = FALSE,
  update = FALSE
)
```

The exact indexed GRCm39 FASTA used by Cell Ranger ARC is expected at `PROJECT_DIR/reference/GRCm39_2024-A/fasta/genome.fa`. A different location can be supplied with `MULTIOME_GRCM39_FASTA`; Step 8D will not silently use a different genome build.

## Main Step 8B outputs

The output directory is `08b_source_matched_persistent_hfd/` within the configured WT results tree.

Core provenance and support files:

- `matched_population_definitions.csv`
- `tables/matched_population_counts_by_diet.csv`
- `tables/matched_population_support.csv`
- `tables/strict_Step8_unit8_reproduction_audit.csv`

RNA benchmark files:

- `tables/RNA_pooled_library_effects_by_matched_population.csv.gz`
- `tables/persistent_gene_RNA_evidence_by_matched_population.csv.gz`
- `tables/persistent_program_ranked_RNA_enrichment.csv`
- `tables/persistent_gene_concordant_RNA_candidates.csv`
- `tables/persistent_gene_discordant_RNA_candidates.csv`
- `tables/population_sensitivity_summary.csv`

Targeted interpretation files:

- `tables/targeted_ORA_input_audit.csv`
- `tables/targeted_MSigDB_ORA_results.csv`

ATAC and integrated files:

- `tables/ATAC_pooled_library_effects_by_matched_population.csv.gz`
- `tables/ATAC_top_diet_candidates_by_matched_population.csv`
- `tables/ATAC_effect_summary_by_population_direction_context.csv`
- `tables/persistent_gene_exact_promoter_evidence_by_population_peak.csv.gz`
- `tables/persistent_gene_integrated_primary_candidates.csv`

The d21-versus-d21 source scatter is the primary reference plot. The d50 scatter is retained with an `S` prefix as supplementary context.

## Main Step 8C outputs

The output directory is `08c_broad_stem_decomposition/` within the configured WT results tree.

- `tables/broad_stem_component_counts_and_proportions.csv`
- `tables/RNA_broad_stem_composition_decomposition.csv.gz`
- `tables/persistent_gene_RNA_composition_decomposition.csv`
- `tables/ATAC_broad_stem_composition_decomposition.csv.gz`
- `tables/ATAC_binary_detection_downsampling_stability.csv.gz`
- `tables/motif_ready_peak_sets.csv.gz`
- `tables/motif_peak_set_summary.csv`
- `tables/motif_background_universe.csv.gz`
- `tables/persistent_gene_exact_promoter_decomposition.csv.gz`
- `tables/persistent_gene_prioritized_promoter_candidates.csv`

The midpoint decomposition is exact on its stated linear scale:

``` text
observed HFD-minus-CON difference = within-state term + composition term
```

RNA uses mean per-cell counts per 10,000. ATAC uses the percentage of nuclei with at least one fragment in a peak. Existing pooled pseudobulk log2FCs are retained alongside those quantities; log2FCs themselves are not averaged or additively decomposed.

## Main Step 8D outputs

The output directory is `08d_broad_stem_motif_enrichment/`.

- `approved_motif_peak_sets_applied.csv`
- `tables/motif_background_matching_audit.csv`
- `tables/motif_peak_membership.csv.gz`
- `tables/JASPAR2024_CORE_motif_metadata.csv`
- `tables/motif_enrichment_results.csv.gz`
- `tables/motif_top_candidates.csv`

Step 8D stops if no reviewed peak set is approved, an approved set has fewer than 25 peaks, the exact FASTA or its index is unavailable, or matched background coverage is inadequate.

Review archives can be generated from the repository root with:

``` bash
bash make_wt_step08c_review_pack.sh "$wt_output"
bash make_wt_step08d_review_pack.sh "$wt_output"
```

## Interpretation boundaries

WT has one pooled library per diet. Therefore, HFD-versus-CON changes, rank-based enrichment, ORA, and RNA/ATAC concordance are descriptive. Cell downsampling evaluates sensitivity to nuclei composition but cannot recover mouse-level variance.

Composite pseudobulk effects can combine two phenomena: altered expression within a component state and altered representation of component states between diets. The strict/core and per-unit outputs must be reviewed before a broad-stem or all-goblet effect is described as cell intrinsic.

Step 8C retains both views intentionally. The observed broad-stem effect describes what changed in the pooled compartment, including redistribution. The midpoint-standardized and component-resolved effects are the appropriate inputs for a narrower regulatory-mechanism question.

An exact promoter overlap is stronger genomic context than a nearest-gene assignment, but it does not prove regulatory causality. Step 8B explicitly defers distal peak-to-gene linkage. HyperMet and HypoImmune are also deferred until after the persistent-gene benchmark so predefined program overlays do not alter primary candidate selection.

Likewise, an enriched motif identifies a compatible DNA-binding preference, not a uniquely identified TF or a binding event. Closely related TFs can share motifs. TF RNA expression, chromVAR-like activity, footprinting where coverage permits, distal peak-to-gene evidence, and orthogonal experiments remain follow-up evidence rather than prerequisites silently folded into the motif P value.
