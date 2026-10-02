#!/usr/bin/env Rscript
# Presentation-only VilB WT d21 walkthrough. Run on SOL with R 4.4.2:
#   Rscript plot_vilb_pseudostem_walkthrough.R
# Optional: VILB_PROJECT_ROOT=/scratch/dsaiz/Yesenia_scData2026 Rscript ...
# Requires completed WT Step 05D, 05E, 05F, 06, and 08E outputs.
# Nothing in this script changes the saved Seurat objects or GitHub workflow.
# C0+C2 is a crypt-core proxy for the G6/G9 portion of Gourab's broad Stem label.
# It is not a cell transfer or a reconstruction of all G6/G7/G8/G9 cells.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Signac)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
})

PROJECT <- Sys.getenv("VILB_PROJECT_ROOT", "/scratch/dsaiz/Yesenia_scData2026")
WT <- file.path(PROJECT, "results", "v2_resequenced", "independent", "wt")
OUT <- file.path(WT, "review_packs", "vilb_pseudostem_walkthrough")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
RUN_COVERAGE <- TRUE             # Set FALSE if fragments cannot be accessed on this node.
COVERAGE_GENES_OVERRIDE <- character() # Example: c("Hmgcs2", "Mgst1")
MAX_COVERAGE_GENES <- 4L

step6_path <- file.path(WT, "06_annotation", "multiome_annotated.rds")
focused_path <- file.path(WT, "05_wnn", "focused_wt_crypt_stem_progenitor",
                         "objects", "WT_crypt_core_0_6_8_focused.rds")
labels_path <- file.path(WT, "05_wnn", "wt_crypt_followup",
                         "05f_review", "refined_working_labels.csv")
mapping_path <- file.path(WT, "05_wnn", "gourab_d21_rna_reference_mapping",
                          "05d_symmetric_pseudobulk_validation", "tables",
                          "primary_similarity_with_bootstrap_support.csv")
promoter_path <- file.path(WT, "08e_refined_wt_crypt_diet",
                           "target_gene_promoter_overlapping_peak_map.csv.gz")
inputs <- c(step6_path, focused_path, labels_path, mapping_path, promoter_path)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) stop("Missing required input(s):\n", paste(missing_inputs, collapse = "\n"))

# Exact user-supplied Figure 2 program lists; Atpif1 is looked up as Atp5if1
# in this GRCm39 RNA assay. Keep absent genes in the audit and plot index.
programs <- list(
  HyperMet = strsplit(paste(
    "Mgst1 Gstm1 Idh1 Gss Gsta4 Gclm Sod1 Sult1a1 Hmgcs2 Acsm3 Csad",
    "Acss2 Shmt1 Impdh2 Akr7a5 Adk Aprt Hao2 Hadh Mecr Scd2 Mif Fads2",
    "Acads Sox9 Aldob Atpif1 Uqcrc1 Slc4a4 Hspd1 Eef1a1 Nox1 Cyba",
    "Tkt Gpx2 Gsto1 Adh1"), " ")[[1]],
  HypoImmune = strsplit(paste(
    "Nupr1 Crip1 Atp2a3 Cdkn1a Cd74 Herpud1 Dnajc10 Ern2 Itpr1 Xbp1",
    "P4hb Selenok Selenos Creb3l1 Serinc3 Sh3glb1 H2-K1 H2-D1 B2m",
    "H2-T23 Psmb8 Calr Ceacam1 Klf4 Mia3 Cd24a Rap1gap Pla2g2f",
    "Spint2 Trp53inp1 Gabarap Gabarapl2 Camk2n1 S100a13 Oit1 Syt7",
    "Serp1 Rab3d Rab15 Agr2 Hspa5 Manf Dnajc3 Creb3l4 Pdia3 Pdia6",
    "Sdf2l1 Sec61b Ramp1 Galnt7 C1galt1 Gfpt1 Gcnt3 St3gal6 Tmem59",
    "St6galnac6 Galnt10 C1galt1c1 Fut2 Galnt3 Rpn1 Krtcap2 Pmm2",
    "Prdx6 Qsox1 Sh3bgrl3 Txndc5 Pdia5 Dap Foxa3 Clca1 Tspan13",
    "Atp2c2 P2rx4 Stim2 Asph Fxyd3 Lrrc26 Kit Plpp1 Pdxdc1 Rrbp1",
    "Tram1 Sec62 Ssr1 Bag1 Lman1 Actb Fos Il13ra1 Ifngr1"), " ")[[1]]
)
stopifnot(length(programs$HyperMet) == 37L, length(programs$HypoImmune) == 91L)
gene_list <- bind_rows(lapply(names(programs), function(p) {
  tibble(program = p, gene = programs[[p]],
         feature = ifelse(programs[[p]] == "Atpif1", "Atp5if1", programs[[p]]),
         expected = if (p == "HyperMet") "up" else "down")
}))
stopifnot(!anyDuplicated(gene_list$feature))

COL <- c(CON = "#E69F00", HFD = "#56B4E9")
FCOL <- c(`0` = "#007F86", `1` = "#8B5CF6", `2` = "#B34A6A",
          `3` = "#75807B", `4` = "#D38C36", `5` = "#3454A5")
theme_slide <- function() theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 17),
        plot.subtitle = element_text(size = 11),
        plot.caption = element_text(size = 9, color = "#464646"),
        legend.position = "right", legend.key.height = grid::unit(0.4, "cm"))
save_pdf <- function(plot, filename, width = 13, height = 8) {
  ggsave(file.path(OUT, filename), plot = plot, width = width, height = height,
         device = grDevices::pdf, useDingbats = FALSE, limitsize = FALSE)
}
require_cols <- function(x, cols, label) {
  m <- setdiff(cols, colnames(x))
  if (length(m)) stop(label, " is missing: ", paste(m, collapse = ", "))
}

message("Loading annotated global WNN and focused Step 05E objects")
global <- readRDS(step6_path)
focused <- readRDS(focused_path)
gm <- global[[]]
fm <- focused[[]]
require_cols(gm, c("wnn_cluster", "cell_state_full", "diet"), "Step 06 metadata")
require_cols(fm, c("crypt_core_wnn_res_0.3", "parent_wnn_cluster", "diet"),
             "Step 05E metadata")
if (!"wnn.umap" %in% names(global@reductions)) stop("Missing global wnn.umap")
if (!all(rownames(fm) %in% rownames(gm))) stop("Focused barcodes absent from global object")
if (!identical(as.character(fm$diet), as.character(gm[rownames(fm), "diet"])))
  stop("Focused and global diet labels do not agree for aligned cells")

labels <- read_csv(labels_path, show_col_types = FALSE)
require_cols(labels, c("focused_cluster", "working_label", "label_status"), "05F map")
labels$focused_cluster <- as.character(labels$focused_cluster)
if (!setequal(labels$focused_cluster, as.character(0:5)) ||
    anyDuplicated(labels$focused_cluster)) stop("Expected reviewed C0-C5 map")
if (!all(labels$label_status == "provisional_working"))
  stop("05F annotation status changed; inspect labels before presenting")
focus_names <- setNames(paste0("C", labels$focused_cluster, "  ", labels$working_label),
                        labels$focused_cluster)

source_map <- read_csv(mapping_path, show_col_types = FALSE) |>
  filter(source_cluster_id %in% paste0("G", 6:9),
         source_to_WT_rank == 1 | (source_cluster_id == "G7" & wt_cluster_id == "W9")) |>
  select(source_cluster_id, wt_cluster_id, wt_cluster_label, spearman_rho,
         source_to_WT_rank, source_to_WT_top1_frequency)
if (!all(paste0("G", 6:9) %in% source_map$source_cluster_id))
  stop("05D source map does not contain all G6-G9 top matches")
expected_top_matches <- c(G6 = "W0", G7 = "W18", G8 = "W19", G9 = "W6")
top_matches <- source_map |> filter(source_to_WT_rank == 1) |>
  distinct(source_cluster_id, wt_cluster_id)
if (!identical(unname(setNames(top_matches$wt_cluster_id,
                             top_matches$source_cluster_id)[names(expected_top_matches)]),
               unname(expected_top_matches)))
  stop("05D source top matches changed; revise proxy rationale and UMAP caption")
write_csv(source_map, file.path(OUT, "00_source_mapping_review.csv"))

emb <- as.data.frame(SeuratObject::Embeddings(global, "wnn.umap"))
if (ncol(emb) < 2L) stop("WNN UMAP must have two dimensions")
names(emb)[1:2] <- c("x", "y")
emb$cell <- rownames(emb)
umap <- emb |>
  left_join(gm |> tibble::rownames_to_column("cell") |>
              select(cell, wnn_cluster, cell_state_full, diet), by = "cell") |>
  mutate(wnn_cluster = as.character(wnn_cluster),
         wnn_label = paste0("W", wnn_cluster),
         annotation = as.character(cell_state_full))
if (anyNA(umap$wnn_cluster) || anyNA(umap$annotation)) stop("Missing WNN labels")
umap$focused_cluster <- as.character(fm[umap$cell, "crypt_core_wnn_res_0.3"])
centers <- umap |>
  group_by(wnn_cluster, wnn_label) |>
  summarise(x = stats::median(x), y = stats::median(y), .groups = "drop")
wnn_levels <- as.character(sort(as.integer(unique(umap$wnn_cluster))))
umap$wnn_label <- factor(umap$wnn_label, levels = paste0("W", wnn_levels))
wnn_colors <- setNames(grDevices::hcl.colors(length(wnn_levels), "Dynamic"),
                       paste0("W", wnn_levels))
base_umap <- function(title, subtitle, caption = NULL) {
  ggplot() + coord_equal() + labs(title = title, subtitle = subtitle,
    caption = caption, x = "WNN UMAP 1", y = "WNN UMAP 2") + theme_slide()
}
wnn_text <- geom_label(data = centers, aes(x, y, label = wnn_label), size = 3.2,
                       inherit.aes = FALSE, alpha = .86, label.size = 0,
                       fill = "white", show.legend = FALSE)

# 01: exact final WNN clustering; 02: reviewed global Step 06 annotations.
p01 <- base_umap("VilB: WNN clusters", "Every nucleus on the saved joint RNA + ATAC embedding") +
  geom_point(data = umap, aes(x, y, color = wnn_label), size = .17, alpha = .7) +
  scale_color_manual(values = wnn_colors, drop = FALSE, name = "WNN cluster") + wnn_text
save_pdf(p01, "01_all_WNN_clusters.pdf", 13, 8)
ann_levels <- sort(unique(umap$annotation))
ann_colors <- setNames(grDevices::hcl.colors(length(ann_levels), "Dark 3"), ann_levels)
p02 <- base_umap("VilB: reviewed global annotations",
  "Full Step 06 identity labels; WNN cluster numbers remain on the map") +
  geom_point(data = umap, aes(x, y, color = annotation), size = .17, alpha = .7) +
  scale_color_manual(values = ann_colors, name = "Step 06 identity") + wnn_text +
  theme(legend.text = element_text(size = 8))
save_pdf(p02, "02_global_annotation.pdf", 16, 9)
marker_panel <- c("Lgr5", "Smoc2", "Mki67", "Top2a", "Dmbt1", "Hmgcs2",
                  "Reg4", "Muc2", "Krt20")
marker_panel <- marker_panel[marker_panel %in% rownames(global[["RNA"]])]
if (length(marker_panel) >= 3L) {
  SeuratObject::DefaultAssay(global) <- "RNA"
  p02b <- Seurat::DotPlot(global, features = marker_panel,
                         group.by = "wnn_cluster") +
    scale_color_gradient(low = "#D7E2E3", high = "#A62B5D") +
    labs(title = "Marker context for global WNN identities",
         subtitle = "Stem, cycle, crypt progenitor, secretory and surface RNA markers",
         x = "Gene", y = "WNN cluster", color = "Mean expression",
         size = "Detected (%)") +
    theme_slide() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_pdf(p02b, "02b_marker_context.pdf", 13, 7)
}

# 03 and 04: overlay focused C0-C5 assignments on the UNCHANGED global UMAP.
# Nonselected cells stay pale grey; WNN number labels stay visible.
focus <- umap |> filter(!is.na(focused_cluster)) |>
  mutate(focus_label = factor(focus_names[focused_cluster],
                              levels = focus_names[as.character(0:5)]))
if (!setequal(unique(focus$focused_cluster), as.character(0:5)))
  stop("Focused C0-C5 assignments do not match expected object")
focus_colors <- setNames(unname(FCOL[as.character(0:5)]),
                         focus_names[as.character(0:5)])
grey_back <- geom_point(data = umap, aes(x, y), color = "#D9D9D9",
                        size = .15, alpha = .4)
p03 <- base_umap("VilB: focused crypt subclustering",
  "C0-C5 projected onto the original global WNN UMAP; provisional 05F identities",
  "Focused WNN fit uses parent W0/W6/W8. Grey nuclei are outside that fit.") +
  grey_back + geom_point(data = focus, aes(x, y, color = focus_label),
                         size = .28, alpha = .85) +
  scale_color_manual(values = focus_colors, name = "Focused subcluster") + wnn_text
save_pdf(p03, "03_focused_C0_to_C5_on_global_WNN.pdf", 15, 9)

proxy <- focus |> filter(focused_cluster %in% c("0", "2"))
p04 <- base_umap("Crypt-core proxy for part of Gourab's 'Stem' label",
  "G6 best matches W0 (C0); G9 best matches W6 (C2). Other WNN states remain grey.",
  "G7 best matches W18/W9 and G8 best matches W19. This proxy does not reconstruct the full source label.") +
  grey_back + geom_point(data = proxy,
    aes(x, y, color = factor(focused_cluster, levels = c("0", "2"))),
    size = .35, alpha = .9) +
  scale_color_manual(values = FCOL[c("0", "2")],
    labels = unname(focus_names[c("0", "2")]), name = "Included nuclei") + wnn_text
save_pdf(p04, "04_C0_C2_crypt_core_proxy_on_global_WNN.pdf", 15, 9)

# Extract only the feature rows needed for plotting from the joined RNA.crypt
# layer. Use the saved normalized data layer for violins, counts for pooled CPM.
if (!all(c("RNA.crypt", "ATAC.crypt") %in% names(focused@assays)))
  stop("Focused object must have RNA.crypt and ATAC.crypt assays")
rna_assay <- focused[["RNA.crypt"]]
atac_assay <- focused[["ATAC.crypt"]]
if (!all(c("counts", "data") %in% SeuratObject::Layers(rna_assay)))
  stop("RNA.crypt needs joined counts and data layers")
rna_counts <- SeuratObject::LayerData(rna_assay, layer = "counts")
rna_data <- SeuratObject::LayerData(rna_assay, layer = "data")
atac_counts <- SeuratObject::LayerData(atac_assay, layer = "counts")
cells <- rownames(fm)[as.character(fm$crypt_core_wnn_res_0.3) %in% c("0", "2")]
if (!all(cells %in% colnames(rna_counts)) || !all(cells %in% colnames(rna_data)) ||
    !all(cells %in% colnames(atac_counts))) stop("Focused matrix/barcode mismatch")
md <- fm[cells, , drop = FALSE]
md$diet <- as.character(md$diet)
md$focused_cluster <- as.character(md$crypt_core_wnn_res_0.3)
if (!setequal(unique(md$diet), c("CON", "HFD"))) stop("Expected both pooled diets")

gene_list$in_RNA <- gene_list$feature %in% rownames(rna_counts)
gene_list$CON_detected_percent <- NA_real_
gene_list$HFD_detected_percent <- NA_real_
for (i in which(gene_list$in_RNA)) {
  g <- gene_list$feature[i]
  gene_list$CON_detected_percent[i] <- 100 * mean(rna_counts[g, cells[md$diet == "CON"]] > 0)
  gene_list$HFD_detected_percent[i] <- 100 * mean(rna_counts[g, cells[md$diet == "HFD"]] > 0)
}
write_csv(gene_list, file.path(OUT, "05_gene_assay_and_detection_audit.csv"))

# Per-cell module score = mean gene-wise z score of saved log-normalized RNA,
# standardized across C0+C2 combined. Skip zero-variance genes and report them.
# These scores describe RNA distributions; they are not tests of diet effects.
score_rows <- list()
score_audit <- list()
for (program in names(programs)) {
  sel <- gene_list |> filter(.data$program == .env$program, .data$in_RNA)
  expr <- as.matrix(rna_data[sel$feature, cells, drop = FALSE])
  gene_sd <- apply(expr, 1L, stats::sd)
  keep <- is.finite(gene_sd) & gene_sd > 0
  score_audit[[program]] <- tibble(program = program, gene = sel$gene,
    feature = sel$feature, used_for_module_score = keep)
  if (sum(keep) < 3L) stop("Fewer than 3 variable genes in ", program)
  z <- scale(t(expr[keep, , drop = FALSE])) # cells x genes; column means 0, SD 1
  scores <- rowMeans(z)
  score_rows[[program]] <- tibble(cell = cells, diet = md$diet,
    component = paste0("C", md$focused_cluster), program = program,
    module_score = as.numeric(scores))
}
write_csv(bind_rows(score_audit), file.path(OUT, "05_module_scoring_feature_audit.csv"))
score_df <- bind_rows(score_rows)
write_csv(score_df, file.path(OUT, "05_module_scores_by_nucleus.csv.gz"))
score_plot <- bind_rows(score_df |> mutate(view = component),
                        score_df |> mutate(view = "C0+C2 proxy")) |>
  mutate(view = factor(view, levels = c("C0+C2 proxy", "C0", "C2")),
         diet = factor(diet, levels = c("CON", "HFD")))
p05 <- ggplot(score_plot, aes(diet, module_score, fill = diet)) +
  geom_violin(scale = "width", trim = TRUE, color = "#303030", linewidth = .2) +
  geom_boxplot(width = .16, outlier.shape = NA, fill = "white", alpha = .55) +
  facet_grid(program ~ view, scales = "free_y") +
  scale_fill_manual(values = COL, guide = "none") +
  labs(title = "Gourab's Figure 2 programs in VilB RNA",
    subtitle = "C0+C2 proxy and its two component populations, shown separately",
    x = NULL, y = "Mean standardized log-normalized expression per nucleus",
    caption = "Scores use assayed, variable genes in each user-supplied list. One pooled library per diet; violins are descriptive.") +
  theme_slide()
save_pdf(p05, "05_program_module_violins_proxy_and_components.pdf", 14, 8)

effect_table <- function(mat, features, selected_cells, group, modality) {
  have <- features[features %in% rownames(mat)]
  if (!length(have)) return(tibble())
  ix_con <- match(selected_cells[group == "CON"], colnames(mat))
  ix_hfd <- match(selected_cells[group == "HFD"], colnames(mat))
  if (anyNA(ix_con) || anyNA(ix_hfd) || !length(ix_con) || !length(ix_hfd))
    stop("Effect matrix or diet labels are misaligned")
  a <- Matrix::rowSums(mat[have, ix_con, drop = FALSE])
  b <- Matrix::rowSums(mat[have, ix_hfd, drop = FALSE])
  depth <- Matrix::colSums(mat)
  pa <- 100 * Matrix::rowMeans(mat[have, ix_con, drop = FALSE] > 0)
  pb <- 100 * Matrix::rowMeans(mat[have, ix_hfd, drop = FALSE] > 0)
  con_cpm <- 1e6 * a / sum(depth[ix_con])
  hfd_cpm <- 1e6 * b / sum(depth[ix_hfd])
  min_counts <- if (modality == "RNA") 20L else 30L
  min_percent <- if (modality == "RNA") 5 else 1
  tibble(feature = have, modality = modality, n_CON = length(ix_con),
    n_HFD = length(ix_hfd), CON_counts = a, HFD_counts = b,
    CON_CPM = con_cpm, HFD_CPM = hfd_cpm, CON_percent = pa,
    HFD_percent = pb, detection_delta_pp = pb - pa,
    log2FC_HFD_vs_CON = log2((hfd_cpm + .1)/(con_cpm + .1))) |>
    mutate(evaluable = CON_counts + HFD_counts >= min_counts &
             pmax(CON_percent, HFD_percent) >= min_percent,
           clear_effect = evaluable & abs(log2FC_HFD_vs_CON) >= .5 &
             (modality == "ATAC" | abs(detection_delta_pp) >= 5))
}

all_effects <- bind_rows(lapply(c("C0+C2 proxy", "C0", "C2"), function(group_name) {
  group_cells <- if (group_name == "C0+C2 proxy") cells else
    cells[md$focused_cluster == substring(group_name, 2L)]
  e <- effect_table(rna_counts, unique(gene_list$feature), group_cells,
                    md[group_cells, "diet"], "RNA")
  e$population <- group_name
  e
}))
gene_effects <- bind_rows(lapply(c("C0+C2 proxy", "C0", "C2"), function(p) {
  gene_list |> select(program, gene, feature, expected) |>
    mutate(population = p)
})) |>
  left_join(all_effects, by = c("population", "feature")) |>
  mutate(expected_direction = case_when(
    is.na(evaluable) ~ NA, !evaluable ~ NA,
    expected == "up" ~ log2FC_HFD_vs_CON > 0,
    TRUE ~ log2FC_HFD_vs_CON < 0),
    clear_expected = !is.na(expected_direction) & expected_direction & clear_effect)
write_csv(gene_effects, file.path(OUT, "06_RNA_pooled_effects_proxy_C0_C2.csv"))
direction_data <- gene_effects |>
  filter(population == "C0+C2 proxy", evaluable %in% TRUE) |>
  mutate(direction = ifelse(expected_direction, "Same direction", "Opposite direction"),
         gene_order = reorder(gene, log2FC_HFD_vs_CON))
if (nrow(direction_data)) {
  p06b <- ggplot(direction_data, aes(log2FC_HFD_vs_CON, gene_order, color = direction)) +
    geom_vline(xintercept = 0, color = "#777777") +
    geom_point(aes(size = pmax(CON_percent, HFD_percent)), alpha = .85) +
    facet_wrap(~program, scales = "free_y", ncol = 2L) +
    scale_color_manual(values = c("Same direction" = "#007F86",
                                  "Opposite direction" = "#B34A6A")) +
    labs(title = "Which assayed genes point in Gourab's Figure 2 direction?",
      subtitle = "All evaluable genes in the C0+C2 proxy; HFD relative to CON",
      x = "Pooled RNA log2 fold change", y = NULL,
      size = "Maximum detection (%)", color = "Direction",
      caption = "Evaluable: >=20 pooled RNA counts and >=5% detection in either diet. Points are descriptive, not significant.") +
    theme_slide() + theme(axis.text.y = element_text(size = 8))
  save_pdf(p06b, "06b_evaluable_gene_directions.pdf", 14, 12)
}

# Every supplied gene is indexed, including any absent from the RNA assay.
for (program in names(programs)) {
  genes <- gene_list |> filter(.data$program == .env$program)
  stats <- gene_effects |> filter(population == "C0+C2 proxy",
                                  .data$program == .env$program)
  for (page_no in seq_len(ceiling(nrow(genes) / 12L))) {
    page <- genes[seq.int((page_no - 1L) * 12L + 1L,
                          min(page_no * 12L, nrow(genes))), , drop = FALSE]
    detail <- left_join(page |> select(gene, feature, in_RNA),
                        stats |> select(gene, log2FC_HFD_vs_CON,
                                       CON_percent, HFD_percent, evaluable), by = "gene")
    detail$display <- ifelse(!detail$in_RNA,
      paste0(detail$gene, "\nnot assayed"),
      paste0(detail$gene, "  log2FC ", sprintf("%+.2f", detail$log2FC_HFD_vs_CON),
             "\nDetected ", sprintf("%.1f", detail$CON_percent), "% / ",
             sprintf("%.1f", detail$HFD_percent), "%"))
    rows <- bind_rows(lapply(seq_len(nrow(detail)), function(i) {
      vals <- if (detail$in_RNA[i])
        as.numeric(rna_data[detail$feature[i], cells]) else rep(NA_real_, length(cells))
      tibble(gene = detail$gene[i], display = detail$display[i],
             diet = factor(md$diet, levels = c("CON", "HFD")), expression = vals)
    }))
    rows$display <- factor(rows$display, levels = detail$display)
    absent <- rows |> filter(is.na(expression)) |>
      distinct(display) |> mutate(x = 1.5, y = 0, label = "Not assayed")
    p <- ggplot(rows, aes(diet, expression, fill = diet)) +
      geom_violin(scale = "width", trim = TRUE, color = "#303030", linewidth = .2,
                  na.rm = TRUE) +
      geom_boxplot(width = .16, outlier.shape = NA, fill = "white", alpha = .55,
                   na.rm = TRUE) +
      geom_text(data = absent, aes(x, y, label = label), inherit.aes = FALSE,
                size = 3.3, color = "#555555") +
      facet_wrap(~display, ncol = 4L, scales = "free_y", drop = FALSE) +
      scale_fill_manual(values = COL, guide = "none") +
      labs(title = paste(program, "genes in the C0+C2 crypt-core proxy"),
        subtitle = paste("Page", page_no, "of", ceiling(nrow(genes)/12L),
                         "| CON versus HFD pooled VilB nuclei"),
        x = NULL, y = "Saved log-normalized RNA per nucleus",
        caption = "Facet: pooled HFD-versus-CON log2FC; detected CON% / HFD%. Sparse/unevaluable genes remain shown. No P values.") +
      theme_slide() + theme(strip.text = element_text(size = 9),
                            axis.text.x = element_text(size = 9))
    save_pdf(p, sprintf("06_%s_gene_violins_page_%02d.pdf", program, page_no), 15, 10)
  }
}

# Descriptive promoter peak overlap is computed specifically for C0+C2;
# the previous 08E C0/C2 peak effects cannot be substituted for this mixture.
peak_map <- read_csv(promoter_path, show_col_types = FALSE)
require_cols(peak_map, c("gene", "peak"), "Step 08E promoter map")
peak_map <- peak_map |> filter(gene %in% gene_list$feature,
                               peak %in% rownames(atac_counts)) |>
  distinct(gene, peak)
if (!nrow(peak_map)) stop("No supplied program genes have assayed promoter peaks")
atac_effects <- effect_table(atac_counts, unique(peak_map$peak), cells,
                             md$diet, "ATAC")
promoter <- peak_map |> left_join(atac_effects, by = c("peak" = "feature")) |>
  left_join(gene_list |>
              transmute(feature, source_gene = gene, program, expected),
            by = c("gene" = "feature")) |>
  mutate(expected_ATAC_direction = evaluable &
    ifelse(expected == "up", log2FC_HFD_vs_CON > 0,
           log2FC_HFD_vs_CON < 0))
write_csv(promoter, file.path(OUT, "07_promoter_peak_effects_proxy.csv.gz"))

rna_proxy <- gene_effects |> filter(population == "C0+C2 proxy") |>
  select(program, gene, feature, expected, RNA_evaluable = evaluable,
         RNA_clear = clear_effect, RNA_log2FC = log2FC_HFD_vs_CON,
         RNA_CON_percent = CON_percent, RNA_HFD_percent = HFD_percent,
         RNA_expected = expected_direction)
peak_summary <- promoter |> group_by(gene) |>
  summarise(n_promoter_peaks = n_distinct(peak),
            n_evaluable_peaks = sum(evaluable, na.rm = TRUE),
            n_same_direction_peaks = sum(expected_ATAC_direction, na.rm = TRUE),
            n_clear_same_direction_peaks = sum(expected_ATAC_direction & clear_effect,
                                               na.rm = TRUE),
            .groups = "drop")
priority <- rna_proxy |> left_join(peak_summary, by = c("feature" = "gene")) |>
  mutate(n_promoter_peaks = coalesce(n_promoter_peaks, 0L),
         n_evaluable_peaks = coalesce(n_evaluable_peaks, 0L),
         n_same_direction_peaks = coalesce(n_same_direction_peaks, 0L),
         n_clear_same_direction_peaks = coalesce(n_clear_same_direction_peaks, 0L),
         directional_candidate = !is.na(RNA_evaluable) & RNA_evaluable &
           RNA_expected & abs(RNA_log2FC) >= .3 & n_same_direction_peaks > 0,
         priority_score = abs(RNA_log2FC) * pmax(RNA_CON_percent,
                                                 RNA_HFD_percent) / 100)
write_csv(priority, file.path(OUT, "07_coverage_candidate_audit.csv"))

if (length(COVERAGE_GENES_OVERRIDE)) {
  selected_coverage <- unique(COVERAGE_GENES_OVERRIDE)
  if (!all(selected_coverage %in% gene_list$gene))
    stop("Coverage override contains a gene outside supplied program lists")
} else {
  selected_coverage <- priority |>
    filter(directional_candidate) |>
    arrange(desc(n_clear_same_direction_peaks), desc(priority_score)) |>
    slice_head(n = MAX_COVERAGE_GENES) |> pull(gene)
}
write_csv(tibble(selected_gene = selected_coverage),
          file.path(OUT, "07_selected_coverage_genes.csv"))
coverage_status <- tibble(gene = selected_coverage, status = "not attempted", detail = "")

if (RUN_COVERAGE && length(selected_coverage)) {
  coverage_obj <- subset(global, cells = cells)
  SeuratObject::DefaultAssay(coverage_obj) <- "ATAC"
  if (!length(Signac::Fragments(coverage_obj[["ATAC"]]))) {
    coverage_status$status <- "skipped"
    coverage_status$detail <- "Global ATAC assay has no accessible fragment reference"
  } else {
    for (g in selected_coverage) {
      message("Coverage: ", g)
      j <- match(g, coverage_status$gene)
      tryCatch({
        p <- Signac::CoveragePlot(coverage_obj, region = g, features = g,
          assay = "ATAC", expression.assay = "RNA", expression.slot = "data",
          annotation = "gene", peaks = TRUE, links = FALSE, tile = FALSE,
          group.by = "diet", extend.upstream = 5000,
          extend.downstream = 5000, ymax = "q95") +
          patchwork::plot_annotation(
            title = paste0(g, ": C0+C2 pooled ATAC coverage"),
            subtitle = "Same nuclei as RNA violins; CON and HFD pooled libraries",
            caption = "Coverage is visual context. A promoter-overlapping peak is not proof of gene regulation.")
        save_pdf(p, paste0("08_coverage_", g, ".pdf"), 13, 8)
        coverage_status$status[j] <- "saved"
      }, error = function(e) {
        coverage_status$status[j] <<- "failed"
        coverage_status$detail[j] <<- conditionMessage(e)
        warning("Coverage failed for ", g, ": ", conditionMessage(e))
      })
    }
  }
  rm(coverage_obj)
} else if (!length(selected_coverage)) {
  message("No RNA + evaluable same-direction promoter candidates; no coverage plots selected.")
}
write_csv(coverage_status, file.path(OUT, "08_coverage_status.csv"))

writeLines(c(
  "Presentation-only VilB crypt-core walkthrough; source lists: user-supplied Figure 2 programs.",
  paste("Generated UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("C0+C2 CON nuclei:", sum(md$diet == "CON")),
  paste("C0+C2 HFD nuclei:", sum(md$diet == "HFD")),
  "G6/G9 source correspondence motivates C0+C2; G7/G8 are not recapitulated.",
  "Program violins are per-nucleus standardized logged RNA. Gene violins are saved logged RNA.",
  "RNA effects reproduce 08E pooled CPM pseudocount/filter definitions for a NEW C0+C2 contrast.",
  "The RNA/ATAC comparison and motif analysis have no biological replication in this WT dataset.",
  "Neither violin shapes nor coverage tracks establish significance or a regulatory link."
), file.path(OUT, "README_interpretation.txt"))
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))
message("Completed: ", OUT)
