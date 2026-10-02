#!/usr/bin/env Rscript
# Summarize completed 08G outputs for review; no biological replicate inference.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L)
  stop("Usage: Rscript 08g_summarize_C0_C2_C0W0_review.R STEP08G_DIR STEP08F_DIR REVIEW_TABLE_DIR",
       call. = FALSE)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})
src <- normalizePath(args[[1]], mustWork = TRUE)
prior <- normalizePath(args[[2]], mustWork = TRUE)
dest <- args[[3]]
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
read_input <- function(root, file) {
  path <- file.path(root, file)
  if (!file.exists(path) || file.info(path)$size <= 0L)
    stop("Missing or empty review input: ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}
write_review <- function(x, filename) {
  readr::write_csv(x, file.path(dest, filename), na = "NA")
}
require_fields <- function(x, fields, label) {
  missing <- setdiff(fields, names(x))
  if (length(missing))
    stop(label, " missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
}

contrast <- read_input(src, "03_contrast_definitions_and_diet_balance.csv")
effects <- read_input(src, "04_descriptive_effect_counts.csv")
gates <- read_input(src, "02_C0W0_gate_size_and_score_sensitivity.csv")
cells <- read_input(src, "02_all_nuclei_fixed_labels_and_C0W0_gates.csv.gz")
marker <- read_input(src, "01_core_marker_detection_by_state_diet.csv")
candidate <- read_input(src, "06_C0W0_top15_candidates_with_depth_check.csv.gz")
examples <- read_input(src, "07_example_loci_selected_without_Gourab_overlay.csv")
coverage <- read_input(src, "08_example_coverage_status.csv")
source <- read_input(prior, "08_Gourab_gene_assay_audit.csv")
rna <- read_input(src, "04_RNA_evaluable_pooled_effects.csv.gz")
require_fields(contrast, c("contrast", "role", "n_CON", "n_HFD", "eligible"), "contrasts")
require_fields(effects, c("contrast", "modality", "consensus_features",
  "evaluable_features", "clear_features", "HFD_higher_clear",
  "CON_higher_clear"), "effect counts")
require_fields(gates, c("gate", "target_n", "selected_n", "CON_n", "HFD_n",
  "jaccard_detected4", "jaccard_rank600"), "score gate audit")
require_fields(cells, c("cell_barcode", "diet", "focused_cluster",
  "parent_wnn_cluster", "ISC_core6_UCell", "ISC_core_genes_detected",
  "nCount_ATAC", "C0W0_top15", "C0W0_top10", "C0W0_top05",
  "C0W0_top20", "C0W0_top25"), "nucleus table")
require_fields(marker, c("gene", "diet", "nuclei", "detected_n"), "marker detection")
require_fields(candidate, c("feature", "tier", "direction", "anchor_log2FC_HFD_vs_CON",
  "C0_same", "C2_same", "top10_same", "proxy_same",
  "depth_matched_valid_draws", "depth_matched_sign_agreement_percent"),
  "ATAC candidates")
require_fields(examples, c("feature", "gene", "tier", "tier_priority",
  "RNA_evaluable", "RNA_clear", "RNA_ATAC_same_direction",
  "anchor_log2FC_HFD_vs_CON", "RNA_log2FC"), "four example loci")
require_fields(coverage, c("contrast", "feature", "status", "detail"), "coverage")
require_fields(source, c("gene", "program", "expected_HFD_sign",
  "in_RNA_assay"), "Gourab gene audit")
require_fields(rna, c("contrast", "feature", "log2FC_HFD_vs_CON",
  "clear_effect", "CON_percent", "HFD_percent"), "RNA effects")
if (anyDuplicated(contrast$contrast) || anyDuplicated(effects[c("contrast", "modality")]) ||
    anyDuplicated(cells$cell_barcode) || anyDuplicated(candidate$feature) ||
    anyDuplicated(source$gene) || anyDuplicated(rna[c("contrast", "feature")]))
  stop("Unexpected duplicate keys in 08G/08F outputs.", call. = FALSE)
if (!all(examples$feature %in% candidate$feature) || nrow(examples) > 4L)
  stop("Example loci are not a subset of the fixed 08G ATAC candidates.", call. = FALSE)
if (nrow(coverage) != 2L * nrow(examples) ||
    !setequal(coverage$feature, examples$feature) ||
    !setequal(coverage$contrast, c("C0_all", "C0W0_top15")))
  stop("Coverage status does not describe both comparison groups for every example.",
       call. = FALSE)
saved <- dplyr::filter(coverage, status == "saved")
for (j in seq_len(nrow(saved))) {
  locus_index <- match(saved$feature[[j]], examples$feature)
  pdf <- file.path(src, "plots", paste0("08_", saved$contrast[[j]],
    "_example_", locus_index, ".pdf"))
  if (!file.exists(pdf) || file.info(pdf)$size <= 0L)
    stop("Coverage status says saved but PDF is missing/empty: ", pdf,
         call. = FALSE)
}
if (!setequal(unique(effects$contrast), contrast$contrast[contrast$eligible]) ||
    !setequal(unique(rna$contrast), contrast$contrast[contrast$eligible]))
  stop("Effect availability differs from contrast eligibility.", call. = FALSE)
if (any(effects$clear_features != effects$HFD_higher_clear +
        effects$CON_higher_clear))
  stop("Descriptive clear-effect directions do not sum to totals.", call. = FALSE)

# An ineligible group has missing effect counts, not zero changed peaks.
viability <- tidyr::expand_grid(contrast = contrast$contrast,
  modality = c("ATAC", "RNA")) |>
  dplyr::left_join(dplyr::select(contrast, contrast, role, n_CON, n_HFD, eligible),
    by = "contrast", relationship = "many-to-one") |>
  dplyr::left_join(dplyr::select(effects, contrast, modality,
    consensus_features, evaluable_features, clear_features,
    HFD_higher_clear, CON_higher_clear),
    by = c("contrast", "modality"), relationship = "one-to-one") |>
  dplyr::mutate(clear_fraction_of_evaluable = dplyr::if_else(
    !is.na(evaluable_features) & evaluable_features > 0,
    clear_features / evaluable_features, NA_real_))
if (any(viability$eligible & is.na(viability$consensus_features)) ||
    any(!viability$eligible & !is.na(viability$consensus_features)))
  stop("Eligibility and effect counts disagree.", call. = FALSE)
write_review(viability, "01_effect_evaluability_and_diet_balance.csv")
write_review(dplyr::filter(viability, modality == "ATAC"),
  "01_ATAC_effect_counts_first.csv")

# Rank sensitivity is a membership comparison; do not compare four- and
# six-gene score magnitudes as though they share a scale.
gate_diet <- cells |>
  dplyr::filter(parent_wnn_cluster == "0", focused_cluster == "0") |>
  dplyr::group_by(diet) |>
  dplyr::summarise(base_n = dplyr::n(),
    base_median_six_score = median(ISC_core6_UCell),
    base_median_six_genes_detected = median(ISC_core_genes_detected),
    base_median_ATAC_counts = median(nCount_ATAC), .groups = "drop")
write_review(gate_diet, "02_C0W0_source_group_by_diet.csv")
write_review(gates, "02_nested_gate_balance_and_rank_sensitivity.csv")
write_review(marker, "02_six_marker_detection_by_state_diet.csv")

tier_counts <- candidate |>
  dplyr::group_by(tier, direction) |>
  dplyr::summarise(peaks = dplyr::n(),
    also_C0_same = sum(C0_same), also_C2_same = sum(C2_same),
    also_top10_same = sum(top10_same),
    also_proxy_same = sum(proxy_same),
    peaks_with_depth_draws = sum(!is.na(depth_matched_valid_draws)),
    depth_sign_ge80_and_15_draws = sum(
      depth_matched_valid_draws >= 15L &
      depth_matched_sign_agreement_percent >= 80, na.rm = TRUE),
    .groups = "drop")
write_review(tier_counts, "03_ATAC_candidate_support_tiers.csv")
resampling <- read_input(src, "06_top200_ATAC_equal_cell_and_depth_sensitivity.csv.gz")
require_fields(resampling, c("contrast", "feature", "method",
  "valid_iterations", "sign_agreement_percent"), "resampling")
stability <- resampling |>
  dplyr::filter(method == "depth_matched") |>
  dplyr::group_by(contrast) |>
  dplyr::summarise(top_evaluable_peaks_sampled = dplyr::n(),
    with_15_valid_draws = sum(valid_iterations >= 15L),
    with_15_draws_and_80pct_same_sign = sum(valid_iterations >= 15L &
      sign_agreement_percent >= 80, na.rm = TRUE), .groups = "drop")
write_review(stability, "03_depth_sensitivity_top200_evaluable_only.csv")

# The overlay is generated after fixed population/peak selection. Missing
# RNA effect rows are labeled not evaluable (or absent assay gene), never 0.
program_detail <- tidyr::expand_grid(contrast = contrast$contrast,
  gene = source$gene) |>
  dplyr::left_join(dplyr::select(contrast, contrast, eligible),
    by = "contrast", relationship = "many-to-one") |>
  dplyr::left_join(dplyr::select(source, gene, program, original_symbol,
    expected_HFD_sign, in_RNA_assay), by = "gene",
    relationship = "many-to-one") |>
  dplyr::left_join(rna |>
    dplyr::transmute(contrast, gene = feature,
      RNA_log2FC_HFD_vs_CON = log2FC_HFD_vs_CON,
      RNA_clear_effect = clear_effect,
      RNA_CON_percent = CON_percent, RNA_HFD_percent = HFD_percent),
    by = c("contrast", "gene"), relationship = "many-to-one") |>
  dplyr::mutate(RNA_evaluable = !is.na(RNA_log2FC_HFD_vs_CON) & eligible,
    source_direction_match = RNA_evaluable &
      sign(RNA_log2FC_HFD_vs_CON) == expected_HFD_sign,
    source_clear_same_direction = source_direction_match &
      dplyr::coalesce(RNA_clear_effect, FALSE))
write_review(program_detail, "04_posthoc_Gourab_RNA_by_contrast.csv")
program_counts <- program_detail |>
  dplyr::group_by(contrast, program) |>
  dplyr::summarise(gene_list_n = dplyr::n(),
    assay_present = sum(in_RNA_assay),
    RNA_evaluable = sum(RNA_evaluable),
    expected_direction = sum(source_direction_match),
    clear_same_direction = sum(source_clear_same_direction),
    .groups = "drop") |>
  dplyr::left_join(dplyr::select(contrast, contrast, n_CON, n_HFD, eligible),
    by = "contrast", relationship = "many-to-one") |>
  dplyr::mutate(dplyr::across(c(RNA_evaluable, expected_direction,
    clear_same_direction), ~ dplyr::if_else(eligible, .x, NA_integer_)))
write_review(program_counts, "04_posthoc_Gourab_program_counts.csv")

example_review <- examples |>
  dplyr::select(feature, gene, tier, tier_priority,
    RNA_evaluable, RNA_clear, RNA_ATAC_same_direction,
    anchor_log2FC_HFD_vs_CON, RNA_log2FC, dplyr::any_of("program")) |>
  dplyr::left_join(coverage, by = "feature", relationship = "one-to-many") |>
  dplyr::arrange(dplyr::desc(tier_priority),
    dplyr::desc(RNA_ATAC_same_direction),
    dplyr::desc(abs(anchor_log2FC_HFD_vs_CON)), contrast)
write_review(example_review, "05_four_example_loci_selection_and_coverage.csv")

review <- c(
  "# WT Step 08G review: C0, C2, and C0/W0 ISC score gates", "",
  "Read `01_ATAC_effect_counts_first.csv`, `02_nested_gate_balance_and_rank_sensitivity.csv`,",
  "`03_ATAC_candidate_support_tiers.csv`, and `05_four_example_loci_selection_and_coverage.csv`.",
  "Then inspect `08g/plots/` and detailed CSV files. C0 and C2 remain separate.", "",
  "- One pooled CON library and one pooled HFD library; effect thresholds are descriptive,",
  "  with no mouse-level P values or significant DAR claims.",
  "- Ineligible contrasts (<30 nuclei in either diet) have NA effect counts.",
  "- A clear ATAC effect is evaluable with |pooled log2FC| >= 0.5;",
  "  evaluability requires >=30 total counts and >=1% detection in either diet.",
  "- Four plotted loci: C0/W0 top15 clear ATAC peaks overlapping annotated",
  "  gene TSS +/-2 kb, with evaluable RNA; ranked by C0/top10 support tier,",
  "  same RNA/ATAC direction, then absolute ATAC effect; one row per peak.",
  "  RNA clear-effect status and Gourab-program membership were not required.",
  "  The four examples are illustrations, not a motif/DAR shortlist.",
  "- Depth resampling covers only up to 200 highest absolute-effect evaluable",
  "  peaks per eligible contrast. Missing depth values elsewhere are untested.",
  "- BED sets are exploratory region inputs; the provided universe is unmatched",
  "  and no motif enrichment was performed in 08G.", "",
  "## Counts at packaging time", "",
  capture.output(print(as.data.frame(dplyr::filter(viability,
    modality == "ATAC") |>
    dplyr::select(contrast, n_CON, n_HFD, eligible, evaluable_features,
      clear_features, HFD_higher_clear, CON_higher_clear)), row.names = FALSE)),
  "")
writeLines(review, file.path(dest, "README_review_first.md"))
message("Wrote 08G review summaries to ", normalizePath(dest))
