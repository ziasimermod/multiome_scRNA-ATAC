#!/usr/bin/env Rscript
# Summarize completed 08F outputs for a compact, auditable review pack.
# One pooled WT CON and one pooled WT HFD library: descriptive effects only.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 08f_summarize_isc_ucell_review.R STEP08F_DIR REVIEW_TABLE_DIR",
       call. = FALSE)
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})
src <- normalizePath(args[[1]], mustWork = TRUE)
dest <- args[[2]]
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
read_input <- function(file) {
  path <- file.path(src, file)
  if (!file.exists(path) || file.info(path)$size == 0)
    stop("Missing/nonempty 08F input required: ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}
save_table <- function(x, name) {
  readr::write_csv(x, file.path(dest, name), na = "NA")
}
require_columns <- function(x, names_required, label) {
  missing <- setdiff(names_required, names(x))
  if (length(missing)) stop(label, " missing: ", paste(missing, collapse = ", "),
                            call. = FALSE)
}

contrast <- read_input("02_contrast_definitions_and_eligibility.csv")
peaks <- read_input("05_ATAC_all_eligible_pooled_effects.csv.gz")
cells <- read_input("02_all_gate_memberships_with_05f_labels.csv.gz")
comparison <- read_input("06_proxy_vs_UCell_union_all_peaks.csv.gz")
require_columns(contrast,
  c("contrast", "role", "gate", "parent", "n_CON", "n_HFD", "eligible", "status"),
  "contrast definitions")
require_columns(peaks, c("contrast", "feature", "evaluable", "clear_effect",
  "role",
  "log2FC_HFD_vs_CON", "CON_percent", "HFD_percent", "CON_counts", "HFD_counts"),
  "ATAC effects")
require_columns(cells, c("gate", "diet", "selected", "parent_wnn_cluster",
  "focused_cluster", "ISC_core6_UCell", "ISC_core_genes_detected",
  "nCount_ATAC", "Cell_cycle_UCell", "Secretory_DSC_UCell"), "gate memberships")
require_columns(comparison, c("gate", "feature", "proxy_union_class",
  "parent_pattern", "union_evaluable", "union_clear", "proxy_clear",
  "union_log2FC", "proxy_log2FC"), "proxy comparison")

focus <- c("proxy_C0_C2", "ISC_union_top15", "ISC_union_top10",
           "ISC_union_top05")
if (!all(focus %in% contrast$contrast))
  stop("Expected the proxy and three union contrasts in definitions.", call. = FALSE)
if (anyDuplicated(peaks[c("contrast", "feature")]))
  stop("ATAC effects have duplicate contrast-feature pairs.", call. = FALSE)
eligible_keys <- contrast$contrast[contrast$eligible]
if (!setequal(unique(peaks$contrast), eligible_keys))
  stop("ATAC effect contrasts do not match eligible contrast definitions.", call. = FALSE)
feature_count <- dplyr::count(peaks, contrast, name = "features")
if (length(unique(feature_count$features)) != 1L)
  stop("Eligible contrasts do not share the same consensus-peak universe.",
       call. = FALSE)

# "Evaluable" and "clear_effect" are calculated in 08F's pooled_effect().
# Clear ATAC effect: evaluable and absolute pooled log2FC >= 0.5.
count_effects <- peaks |>
  dplyr::group_by(contrast) |>
  dplyr::summarise(
    consensus_peaks = dplyr::n(),
    evaluable_peaks = sum(evaluable, na.rm = TRUE),
    clear_peaks = sum(clear_effect, na.rm = TRUE),
    clear_HFD_open = sum(clear_effect & log2FC_HFD_vs_CON > 0, na.rm = TRUE),
    clear_CON_open = sum(clear_effect & log2FC_HFD_vs_CON < 0, na.rm = TRUE),
    evaluable_abs_log2FC_ge_1 = sum(evaluable &
      abs(log2FC_HFD_vs_CON) >= 1, na.rm = TRUE),
    any_detection_ge_1pct = sum(pmax(CON_percent, HFD_percent) >= 1,
                                 na.rm = TRUE),
    .groups = "drop")
summary <- contrast |>
  dplyr::select(contrast, role, gate, parent, n_CON, n_HFD, eligible, status) |>
  dplyr::left_join(count_effects, by = "contrast", relationship = "one-to-one") |>
  dplyr::mutate(clear_fraction_of_evaluable = dplyr::if_else(
    !is.na(evaluable_peaks) & evaluable_peaks > 0,
    clear_peaks / evaluable_peaks, NA_real_),
    clear_per_1000_evaluable = 1000 * clear_fraction_of_evaluable)
if (any(summary$eligible & is.na(summary$consensus_peaks)) ||
    any(!summary$eligible & !is.na(summary$consensus_peaks)))
  stop("Effect availability and eligibility are inconsistent.", call. = FALSE)
if (any(summary$clear_peaks != summary$clear_HFD_open + summary$clear_CON_open,
        na.rm = TRUE))
  stop("Clear peaks do not resolve into the two diet directions.", call. = FALSE)
save_table(summary, "01_viable_ATAC_peaks_all_13_populations.csv")
save_table(summary |> dplyr::filter(contrast %in% focus) |>
  dplyr::mutate(contrast = factor(contrast, levels = focus)) |>
  dplyr::arrange(contrast) |>
  dplyr::mutate(contrast = as.character(contrast)),
  "01_proxy_and_union_15_10_05_viable_peak_counts.csv")

# A diet-blind gate can still have different selection rates by diet.
gate_balance <- cells |>
  dplyr::group_by(gate, diet, parent_wnn_cluster) |>
  dplyr::summarise(parent_n = dplyr::n(), selected_n = sum(selected),
    selected_fraction = mean(selected),
    selected_median_ISC_score = if (any(selected)) median(ISC_core6_UCell[selected])
      else NA_real_,
    selected_median_ISC_genes_detected = if (any(selected))
      median(ISC_core_genes_detected[selected]) else NA_real_,
    selected_median_ATAC_counts = if (any(selected))
      median(nCount_ATAC[selected]) else NA_real_,
    .groups = "drop")
save_table(gate_balance, "02_gate_balance_parent_by_diet.csv")
subcluster <- cells |>
  dplyr::filter(selected) |>
  dplyr::count(gate, parent_wnn_cluster, diet, focused_cluster,
               name = "selected_n") |>
  dplyr::group_by(gate, parent_wnn_cluster, diet) |>
  dplyr::mutate(fraction_within_parent_diet_gate = selected_n / sum(selected_n)) |>
  dplyr::ungroup()
save_table(subcluster, "02_selected_C0_C5_composition.csv")
overlap <- cells |>
  dplyr::mutate(is_proxy = as.character(focused_cluster) %in% c("0", "2")) |>
  dplyr::group_by(gate, diet) |>
  dplyr::summarise(proxy_n = sum(is_proxy),
    union_gate_n = sum(selected), shared_n = sum(is_proxy & selected),
    union_only_n = sum(selected & !is_proxy),
    proxy_only_n = sum(is_proxy & !selected),
    jaccard_cells = shared_n / (proxy_n + union_gate_n - shared_n),
    .groups = "drop")
save_table(overlap, "02_proxy_union_nucleus_overlap_by_gate_diet.csv")

# Evaluable peaks in four requested comparisons preserve the exact 08F
# detection, count, effect and direction fields; the complete denominator
# including nonevaluable peaks is retained in the summary table above.
readr::write_csv(peaks |>
  dplyr::filter(contrast %in% focus, evaluable) |>
  dplyr::select(contrast, feature, CON_n, HFD_n, CON_counts, HFD_counts,
    CON_CPM, HFD_CPM, CON_percent, HFD_percent, detection_delta_pp,
    log2FC_HFD_vs_CON, evaluable, clear_effect, direction),
  file.path(dest, "03_ATAC_evaluable_proxy_and_union_peaks.csv.gz"), na = "NA")
readr::write_csv(peaks |>
  dplyr::filter(role == "parent_gate", clear_effect),
  file.path(dest, "03_parent_clear_peaks.csv.gz"), na = "NA")
overlap_counts <- comparison |>
  dplyr::count(gate, proxy_union_class, parent_pattern, name = "peaks")
save_table(overlap_counts, "04_proxy_union_peak_overlap_by_gate.csv")
readr::write_csv(comparison |>
  dplyr::filter(union_clear | proxy_clear) |>
  dplyr::select(gate, feature, proxy_union_class, parent_pattern,
    union_evaluable, union_clear, union_log2FC, union_CON_percent,
    union_HFD_percent, proxy_evaluable, proxy_clear, proxy_log2FC,
    proxy_CON_percent, proxy_HFD_percent, n_parent_evaluable,
    n_HFD_open, n_CON_open, HFD_open_parents, CON_open_parents),
  file.path(dest, "04_proxy_union_clear_in_either_arm.csv.gz"), na = "NA")

# The 08F depth resampling is only for the 200 highest |log2FC| evaluable
# peaks per selected contrast. Do not extrapolate its stability to all peaks.
stability <- read_input("07_ATAC_equal_cell_and_depth_sensitivity.csv.gz")
if (nrow(stability)) {
  require_columns(stability, c("contrast", "feature", "method",
    "valid_iterations", "sign_agreement_percent"), "depth sensitivity")
  sampled <- stability |>
    dplyr::filter(method == "depth_matched") |>
    dplyr::group_by(contrast) |>
    dplyr::summarise(peaks_sampled = dplyr::n(),
      peaks_with_at_least_15_valid_draws = sum(valid_iterations >= 15L),
      peaks_with_at_least_80pct_sign_agreement = sum(valid_iterations >= 15L &
        sign_agreement_percent >= 80, na.rm = TRUE), .groups = "drop")
  save_table(sampled, "05_depth_sensitivity_top_200_only.csv")
}

plot_data <- summary |>
  dplyr::filter(contrast %in% focus, eligible) |>
  dplyr::mutate(contrast = factor(contrast, levels = focus)) |>
  tidyr::pivot_longer(c(clear_HFD_open, clear_CON_open),
    names_to = "diet_direction", values_to = "direction_peak_count") |>
  dplyr::mutate(diet_direction = dplyr::recode(diet_direction,
    clear_HFD_open = "HFD higher", clear_CON_open = "CON higher"))

if (nrow(plot_data)) {
  p <- ggplot2::ggplot(plot_data,
    ggplot2::aes(x = contrast, y = direction_peak_count, fill = diet_direction)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("CON higher" = "#E69F00",
      "HFD higher" = "#56B4E9")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)) +
    ggplot2::labs(title = "WT crypt: descriptive accessibility candidates",
      subtitle = "Identical consensus peaks; evaluable and |pooled log2FC| ≥ 0.5; one pooled library per diet",
      x = NULL, y = "Clear-effect consensus peaks", fill = NULL)
  ggplot2::ggsave(file.path(dest, "01_proxy_and_union_clear_peak_counts.pdf"),
    p, width = 11, height = 6, useDingbats = FALSE)
}

review <- c(
  "# WT Step 08F: C0+C2 proxy and UCell ISC gates",
  "",
  "Open `01_proxy_and_union_15_10_05_viable_peak_counts.csv` first, then",
  "`02_gate_balance_parent_by_diet.csv` and the 05F composition and overlap tables.",
  "The four requested comparisons use the same fixed consensus-peak universe.",
  "A missing value means the contrast had fewer than 30 nuclei in one diet and",
  "was not evaluated; it does not mean zero affected peaks.",
  "",
  "- `evaluable_peaks`: total peak counts >=30 across selected CON+HFD nuclei",
  "  and >=1% detection in at least one diet (the original 08F rule).",
  "- `clear_peaks`: evaluable and absolute pooled ATAC log2FC >=0.5.",
  "  CON/HFD direction splits sum to `clear_peaks`.",
  "- `clear_per_1000_evaluable` compares candidate yield despite changing",
  "  gate sizes; interpret with cell balance, depth and marker detection.",
  "- Depth resampling covers up to 200 peaks per specified contrast, chosen",
  "  by absolute effect, not every clear peak. It is a sensitivity audit.",
  "- Gourab gene overlays were calculated after gates and peaks were defined.",
  "",
  "The one pooled CON and one pooled HFD library do not provide biological",
  "replication. No nucleus-level P value or peak FDR is inferred. A clear",
  "effect is a descriptive threshold, not a statistically significant DAR.",
  "The UCell union is an exploratory combination until W0/W6/W8 identity,",
  "diet balance and effect directions have been checked separately.",
  "",
  "## Summary at packaging time",
  "",
  capture.output(print(as.data.frame(summary |>
    dplyr::filter(contrast %in% focus) |>
    dplyr::select(contrast, n_CON, n_HFD, eligible, evaluable_peaks,
                  clear_peaks, clear_HFD_open, clear_CON_open,
                  clear_per_1000_evaluable)), row.names = FALSE)),
  "")
writeLines(review, file.path(dest, "README_review_first.md"))
message("Wrote 08F review summaries to ", normalizePath(dest))
