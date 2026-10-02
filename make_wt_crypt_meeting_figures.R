#!/usr/bin/env Rscript
# Four presentation-sized WT crypt figures from the completed Step 05D/05E CSVs.
# Standalone provenance script for Yesenia; no Seurat object or GitHub checkout needed.
#
# Usage on SOL:
#   Rscript make_wt_crypt_meeting_figures.R [WT_OUTPUT_DIR] [FIGURE_DIR]
# Defaults below use Dom's v2 resequenced WT outputs. Pass Yesenia's mirror
# as WT_OUTPUT_DIR if those result files have been copied there.
# Requires ggplot2; grid and all other functions are in base/recommended R.
# An extracted review pack can also be supplied as WT_OUTPUT_DIR.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Usage: Rscript make_wt_crypt_meeting_figures.R [WT_OUTPUT_DIR] [FIGURE_DIR]", call. = FALSE)
wt_output <- if (length(args) >= 1L) args[[1L]] else
  "/scratch/dsaiz/Yesenia_scData2026/results/v2_resequenced/independent/wt"
out_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path(wt_output, "review_packs", "meeting_figures")
if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("ggplot2 is required in this R environment.", call. = FALSE)

step5d <- file.path(wt_output, "05_wnn", "gourab_d21_rna_reference_mapping",
  "05d_symmetric_pseudobulk_validation", "tables")
step5e <- file.path(wt_output, "05_wnn", "focused_wt_crypt_stem_progenitor", "tables")
inputs <- c(
  similarity = file.path(step5d, "primary_similarity_with_bootstrap_support.csv"),
  summary = file.path(step5e, "primary_r0.3_cluster_summary.csv"),
  parents = file.path(step5e, "primary_r0.3_parent_composition.csv"),
  markers = file.path(step5e, "primary_r0.3_targeted_RNA_promoter_summary_by_cluster.csv")
)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) stop("Missing input files:\n", paste(missing_inputs, collapse = "\n"), call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_input <- function(path) utils::read.csv(path, check.names = FALSE,
  stringsAsFactors = FALSE)
similarity <- read_input(inputs[["similarity"]])
summary <- read_input(inputs[["summary"]])
parents <- read_input(inputs[["parents"]])
markers <- read_input(inputs[["markers"]])
clusters <- paste0("C", 0:5)
ink <- "#20354A"
muted <- "#526577"
parent_colors <- c(W0 = "#2B6F97", W6 = "#80B1D3", W8 = "#B3B3B3")
cluster_colors <- c(C0 = "#1B658A", C1 = "#72AFC9", C2 = "#52A68B",
  C3 = "#D69E5C", C4 = "#BE84A7", C5 = "#868E96")

theme_slide <- function() {
  ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(color = ink, size = 12),
      axis.text = ggplot2::element_text(color = ink, size = 10),
      legend.title = ggplot2::element_text(size = 10, color = ink),
      legend.text = ggplot2::element_text(size = 10, color = ink),
      plot.margin = ggplot2::margin(8, 12, 8, 12)
    )
}

# Slide 1: four clusters sharing Gourab's Stem label do not have one WT match.
source_ids <- c("G6", "G7", "G8", "G9")
wt_ids <- c("W0", "W6", "W8", "W9", "W18", "W19")
source_map <- similarity[similarity$source_cluster_id %in% source_ids &
  similarity$wt_cluster_id %in% wt_ids, c("source_cluster_id", "wt_cluster_id", "spearman_rho")]
if (nrow(source_map) != length(source_ids) * length(wt_ids) ||
    anyDuplicated(source_map[c("source_cluster_id", "wt_cluster_id")]))
  stop("Unexpected Step 05D similarity matrix; expected one row per selected G/W pair.", call. = FALSE)
source_map$source_cluster_id <- factor(source_map$source_cluster_id,
  levels = rev(source_ids))
source_map$wt_cluster_id <- factor(source_map$wt_cluster_id, levels = wt_ids,
  labels = c("W0\nCycling\nISC/TA", "W6\nCrypt\nprogenitor",
    "W8\nISC-like\ncrypt", "W9\nReg4+\nDCS",
    "W18\nCycling sec.\nprogenitor", "W19\nSurface\ncolonocyte"))
source_map$label_color <- ifelse(abs(source_map$spearman_rho) >= 0.48,
  "white", ink)
source_map$label <- sprintf("%.2f", source_map$spearman_rho)
p_source <- ggplot2::ggplot(source_map,
    ggplot2::aes(x = wt_cluster_id, y = source_cluster_id, fill = spearman_rho)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = label, color = label_color),
    size = 4, show.legend = FALSE) +
  ggplot2::scale_color_identity() +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
    high = "#B2182B", midpoint = 0, limits = c(-0.8, 0.8),
    name = "Spearman rho") +
  ggplot2::scale_y_discrete(labels = function(x) paste(x, "·  Stem")) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::coord_equal() + theme_slide() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9),
    axis.text.y = ggplot2::element_text(size = 12),
    legend.position = "right")

# Slide 2 left: focused WNN parent composition.
parent_data <- parents[parents$focused_subcluster %in% 0:5 &
  parents$parent_wnn_cluster %in% c(0, 6, 8), , drop = FALSE]
if (nrow(parent_data) != 18L) stop("Expected 18 focused cluster / parent rows.", call. = FALSE)
parent_data$cluster <- factor(paste0("C", parent_data$focused_subcluster),
  levels = rev(clusters))
parent_data$parent <- factor(paste0("W", parent_data$parent_wnn_cluster),
  levels = rev(names(parent_colors)))
parent_data$share <- parent_data$percent_of_focused_subcluster
parent_data$label <- ifelse(parent_data$share >= 12,
  sprintf("%.0f%%", parent_data$share), "")
parent_data$label_color <- ifelse(parent_data$parent == "W0", "white", ink)
n_by_cluster <- setNames(summary$nuclei, paste0("C", summary$focused_subcluster))
p_parent <- ggplot2::ggplot(parent_data,
    ggplot2::aes(x = cluster, y = share, fill = parent)) +
  ggplot2::geom_col(width = .7, color = "white", linewidth = .2) +
  ggplot2::geom_text(ggplot2::aes(label = label, color = label_color),
    position = ggplot2::position_stack(vjust = .5), size = 3,
    show.legend = FALSE) +
  ggplot2::scale_color_identity() +
  ggplot2::scale_fill_manual(values = parent_colors,
    breaks = names(parent_colors), drop = FALSE,
    name = "Parent") +
  ggplot2::scale_x_discrete(labels = function(x)
    sprintf("%s  (n=%s)", x, format(n_by_cluster[x], big.mark = ",", trim = TRUE))) +
  ggplot2::scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Parent WT global cluster (%)") + theme_slide() +
  ggplot2::theme(legend.position = "top")

# Slide 2 right: explicitly computed column z-scores over C0-C5.
program_columns <- c(
  "curated_canonical_ISC_score_median",
  "source_G6_canonical_ISC_score_median",
  "source_G9_crypt_progenitor_score_median",
  "source_G8_absorptive_like_legacy_stem_score_median",
  "curated_secretory_commitment_score_median"
)
program_names <- c("Canonical ISC", "Gourab G6", "Gourab G9",
  "Gourab G8", "Secretory")
summary <- summary[match(0:5, summary$focused_subcluster), , drop = FALSE]
if (nrow(summary) != 6L || anyNA(summary$focused_subcluster) ||
    !identical(as.integer(summary$focused_subcluster), 0:5))
  stop("Expected one Step 05E summary row for each C0-C5 cluster.", call. = FALSE)
program_data <- do.call(rbind, lapply(seq_along(program_columns), function(i) {
  raw <- as.numeric(summary[[program_columns[[i]]]])
  sigma <- sqrt(mean((raw - mean(raw))^2))
  data.frame(cluster = factor(clusters, levels = rev(clusters)),
    program = program_names[[i]],
    z = if (sigma == 0) 0 else (raw - mean(raw)) / sigma)
}))
program_data$program <- factor(program_data$program, levels = program_names)
p_program <- ggplot2::ggplot(program_data,
    ggplot2::aes(x = program, y = cluster, fill = z)) +
  ggplot2::geom_tile(color = "white", linewidth = .25) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
    high = "#B2182B", midpoint = 0, limits = c(-1.75, 1.75),
    name = "Relative score\n(column z)") +
  ggplot2::labs(x = NULL, y = NULL) + theme_slide() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 32,
    hjust = 1, size = 9), axis.text.y = ggplot2::element_text(size = 11))

# Slide 3: per-cluster detection for selected biologically interpretable RNA markers.
genes <- c("Lgr5", "Smoc2", "Lrig1", "Rnf43", "Dmbt1", "Hmgcs2", "Muc2")
rna <- markers[markers$focused_subcluster %in% 0:5 &
  markers$gene %in% genes, c("focused_subcluster", "gene", "RNA_percent_detected")]
if (nrow(rna) != 42L || anyDuplicated(rna[c("focused_subcluster", "gene")]))
  stop("Expected 42 unique focused cluster / marker rows.", call. = FALSE)
rna$cluster <- factor(paste0("C", rna$focused_subcluster),
  levels = rev(clusters))
rna$gene <- factor(rna$gene, levels = genes)
rna$label <- sprintf("%.0f%%", rna$RNA_percent_detected)
rna$label_color <- ifelse(rna$RNA_percent_detected >= 32, "white", ink)
p_rna <- ggplot2::ggplot(rna,
    ggplot2::aes(x = gene, y = cluster, fill = RNA_percent_detected)) +
  ggplot2::geom_tile(color = "white", linewidth = .25) +
  ggplot2::geom_text(ggplot2::aes(label = label, color = label_color),
    size = 3.8, show.legend = FALSE) +
  ggplot2::scale_color_identity() +
  ggplot2::scale_fill_gradient(low = "#FFFFCC", high = "#08306B",
    limits = c(0, max(50, max(rna$RNA_percent_detected))),
    name = "RNA detected (%)") +
  ggplot2::labs(x = NULL, y = NULL) + theme_slide() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12),
    axis.text.y = ggplot2::element_text(size = 12))

# Slide 4: descriptive proportions in the two pooled libraries.
total_con <- sum(summary$CON_nuclei)
total_hfd <- sum(summary$HFD_nuclei)
composition <- rbind(
  data.frame(diet = "CON", cluster = paste0("C", summary$focused_subcluster),
    nuclei = summary$CON_nuclei, percent = 100 * summary$CON_nuclei / total_con),
  data.frame(diet = "HFD", cluster = paste0("C", summary$focused_subcluster),
    nuclei = summary$HFD_nuclei, percent = 100 * summary$HFD_nuclei / total_hfd)
)
composition$diet <- factor(composition$diet, levels = c("CON", "HFD"),
  labels = c(sprintf("CON  (n=%s)", format(total_con, big.mark = ",")),
    sprintf("HFD  (n=%s)", format(total_hfd, big.mark = ","))))
composition$cluster <- factor(composition$cluster, levels = rev(clusters))
composition$label <- ifelse(composition$percent >= 5,
  sprintf("%s  %.1f%%", composition$cluster, composition$percent), "")
composition$label_color <- ifelse(composition$cluster %in% c("C0", "C3", "C5"),
  "white", ink)
p_composition <- ggplot2::ggplot(composition,
    ggplot2::aes(x = diet, y = percent, fill = cluster)) +
  ggplot2::geom_col(width = .62, color = "white", linewidth = .4) +
  ggplot2::geom_text(ggplot2::aes(label = label, color = label_color),
    position = ggplot2::position_stack(vjust = .5), size = 4,
    fontface = "bold", show.legend = FALSE) +
  ggplot2::scale_color_identity() +
  ggplot2::scale_fill_manual(values = cluster_colors,
    breaks = rev(clusters), guide = "none") +
  ggplot2::scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  ggplot2::labs(x = NULL, y = "Share of focused crypt nuclei (%)") +
  theme_slide() +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "#E5EAF0"),
    axis.text.x = ggplot2::element_text(size = 13))

slides <- list(
  list(filename = "01_source_stem_heterogeneity.pdf",
    title = "One source label spans several WT programs",
    subtitle = "Diet-balanced RNA pseudobulk similarity, 370 shared identity genes, cell-cycle genes excluded",
    plots = list(p_source),
    callout = "G6 -> W0   |   G9 -> W6   |   G8 -> W19   |   G7 -> W18/W9",
    source = "Step 05D primary similarity. Spearman rho; profile similarity is not cell transfer."),
  list(filename = "02_refined_crypt_population_map.pdf",
    title = "The WT crypt core resolves distinct states",
    subtitle = "Focused WNN r0.3 re-clustering of global W0, W6 and W8; 4,468 nuclei",
    plots = list(p_parent, p_program),
    callout = "C0: primary ISC-enriched     C2: progenitor-like     C5: stem-associated sensitivity",
    source = "Step 05E parent composition and program scores. Working labels; no Step 06 override."),
  list(filename = "03_direct_RNA_marker_evidence.pdf",
    title = "Detected transcripts anchor the focused labels",
    subtitle = "Fraction of nuclei with detected RNA; markers do not identify every nucleus in a cluster",
    plots = list(p_rna),
    callout = "C0: Lgr5 / Lrig1 / Rnf43    |    C2: Dmbt1 / Hmgcs2    |    C4: Muc2",
    source = "Step 05E targeted RNA audit. Sparse nuclear RNA and regional expression affect detection."),
  list(filename = "04_diet_composition_context.pdf",
    title = "CON and HFD libraries differ in crypt-state mixture",
    subtitle = "The same 4,468 focused WT nuclei; one pooled library per diet",
    plots = list(p_composition),
    callout = "C1: 10.9% in CON and 35.8% in HFD; composition does not establish a within-state response.",
    source = "Step 05E cluster summary. Percentages are descriptive, not biological-replicate estimates.")
)

draw_slide <- function(slide) {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  grid::grid.text(slide$title, x = .055, y = .955, just = c("left", "top"),
    gp = grid::gpar(fontsize = 22, fontface = "bold", col = ink))
  grid::grid.text(slide$subtitle, x = .055, y = .875, just = c("left", "top"),
    gp = grid::gpar(fontsize = 12, col = muted))
  if (length(slide$plots) == 1L) {
    print(slide$plots[[1L]], newpage = FALSE,
      vp = grid::viewport(x = .51, y = .525, width = .87, height = .65))
  } else {
    print(slide$plots[[1L]], newpage = FALSE,
      vp = grid::viewport(x = .29, y = .54, width = .47, height = .61))
    print(slide$plots[[2L]], newpage = FALSE,
      vp = grid::viewport(x = .755, y = .55, width = .45, height = .59))
  }
  grid::grid.text(slide$callout, x = .055, y = .105,
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 12, fontface = "bold", col = ink))
  grid::grid.text(slide$source, x = .055, y = .035,
    just = c("left", "center"),
    gp = grid::gpar(fontsize = 8.5, col = muted))
}

open_pdf <- function(path, onefile = FALSE) grDevices::pdf(path,
  width = 13.333, height = 7.5, onefile = onefile, useDingbats = FALSE,
  family = "Helvetica")
for (slide in slides) {
  open_pdf(file.path(out_dir, slide$filename))
  tryCatch(draw_slide(slide), finally = grDevices::dev.off())
}
open_pdf(file.path(out_dir, "WT_crypt_working_hypothesis_figures_R.pdf"),
  onefile = TRUE)
tryCatch(for (slide in slides) draw_slide(slide),
  finally = grDevices::dev.off())

input_manifest <- data.frame(
  input = names(inputs), path = normalizePath(inputs),
  bytes = file.info(inputs)$size,
  md5 = unname(tools::md5sum(inputs)), stringsAsFactors = FALSE)
utils::write.csv(input_manifest, file.path(out_dir, "input_manifest.csv"),
  row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Wrote four slide-sized PDFs and a combined PDF to: ", out_dir)
