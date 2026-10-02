# Shared by 05F and 08E. Release: wt_crypt_followup.
# Working analysis units only; no mutation of Step 05E or Step 06 inputs.
FOLLOWUP_RELEASE <- "wt_crypt_followup"

need <- function(ok, message) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(message, call. = FALSE)
}
require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  need(length(missing) == 0L,
       paste(label, "is missing:", paste(missing, collapse = ", ")))
}
read_csv_checked <- function(path) {
  need(file.exists(path), paste("Missing input:", path))
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}
write_table <- function(x, directory, filename) {
  readr::write_csv(x, file.path(directory, filename), na = "NA")
}
safe_z <- function(x) {
  out <- rep(NA_real_, length(x)); ok <- is.finite(x)
  if (sum(ok) < 2L || stats::sd(x[ok]) == 0) {
    out[ok] <- 0
  } else {
    out[ok] <- (x[ok] - mean(x[ok])) / stats::sd(x[ok])
  }
  out
}
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0)
    return(NA_real_)
  stats::cor(x[ok], y[ok], method = "spearman")
}
align_columns <- function(x, cells, label) {
  need(!is.null(colnames(x)) && !anyDuplicated(colnames(x)),
       paste(label, "has missing or duplicated cell names."))
  need(!anyDuplicated(cells) && setequal(colnames(x), cells),
       paste(label, "and metadata contain different cell sets."))
  audit <- tibble::tibble(matrix = label, matrix_cells = ncol(x),
    metadata_cells = length(cells), same_cell_set = TRUE,
    identical_order_before_alignment = identical(colnames(x), cells))
  list(matrix = x[, cells, drop = FALSE], audit = audit)
}
load_focused <- function(path) {
  need(file.exists(path), paste("Missing saved focused object:", path))
  obj <- readRDS(path)
  meta <- obj[[]]
  columns <- c("sample_id", "diet", "parent_wnn_cluster", "Phase", "S.Score",
    "G2M.Score", "nCount_RNA", "nCount_ATAC", "percent.mt", "TSS.enrichment",
    "crypt_core_wnn_res_0.3", "crypt_core_rna_res_0.3")
  require_columns(meta, columns, "Step 05E metadata")
  need(!anyNA(meta[, c("sample_id", "diet", "crypt_core_wnn_res_0.3",
    "crypt_core_rna_res_0.3")]), "Missing sample/diet/cluster assignments.")
  meta$cell_barcode <- rownames(meta)
  for (nm in c("sample_id", "diet", "parent_wnn_cluster", "Phase"))
    meta[[nm]] <- as.character(meta[[nm]])
  meta$focused_cluster <- as.character(meta[["crypt_core_wnn_res_0.3"]])
  meta$rna_cluster <- as.character(meta[["crypt_core_rna_res_0.3"]])
  need(setequal(unique(meta$focused_cluster), as.character(0:5)),
    "This release expects the reviewed six r0.3 clusters C0-C5. Do not reuse its label map on a different fit.")
  need(setequal(unique(meta$diet), c("CON", "HFD")), "Expected CON and HFD.")
  samples <- unique(meta[c("sample_id", "diet")])
  need(nrow(samples) == 2L && !anyDuplicated(samples$sample_id),
    "This descriptive release expects one pooled library per diet.")
  mats <- list(); audits <- list()
  for (nm in c("RNA_counts", "RNA_data", "ATAC_counts")) {
    assay <- if (nm == "ATAC_counts") "ATAC.crypt" else "RNA.crypt"
    layer <- if (nm == "RNA_data") "data" else "counts"
    need(assay %in% names(obj), paste("Missing assay", assay))
    need(layer %in% SeuratObject::Layers(obj[[assay]]),
      paste("Expected joined", layer, "layer in", assay))
    a <- align_columns(SeuratObject::LayerData(obj[[assay]], layer = layer),
      meta$cell_barcode, nm)
    mats[[nm]] <- a$matrix; audits[[nm]] <- a$audit
    need(!is.null(rownames(a$matrix)) && !anyDuplicated(rownames(a$matrix)),
      paste("Missing/duplicated features in", nm))
  }
  need(identical(rownames(mats$RNA_counts), rownames(mats$RNA_data)),
    "RNA counts/data feature ordering differs.")
  list(object = obj, metadata = meta, matrices = mats,
    alignment = dplyr::bind_rows(audits))
}
working_labels <- function() {
  tibble::tibble(focused_cluster = as.character(0:5),
    analysis_unit = paste0("WTcrypt_C", 0:5),
    working_label = c("ISC-enriched crypt", "Htr4/Fut9 crypt state",
      "Dmbt1/Hmgcs2 progenitor-like crypt", "Crypt state unresolved",
      "Secretory/goblet-like crypt", "Stem-associated crypt unresolved"),
    focused_annotation_confidence = c("moderate", "moderate", "moderate", "low", "moderate", "low"),
    analysis_role = c("primary", "exploratory", "primary", "sensitivity",
      "secretory_comparator", "sensitivity"),
    label_status = "provisional_working", approved_for_global_annotation = FALSE)
}
sentinel_panel <- function() {
  # Selected markers, not complete published gene signatures or a diagnostic classifier.
  tibble::tribble(
    ~panel, ~gene, ~aliases, ~source,
    "revival_sentinel", "Clu", "Clu", "Ayyaz2019",
    "fetal_repair_sentinels", "Ly6a", "Ly6a", "Yui2018",
    "fetal_repair_sentinels", "Anxa1", "Anxa1", "Yui2018",
    "fetal_repair_sentinels", "Tacstd2", "Tacstd2", "Yui2018",
    "YAP_associated_targets", "Ccn1", "Ccn1;Cyr61", "Yui2018",
    "YAP_associated_targets", "Ccn2", "Ccn2;Ctgf", "Yui2018",
    "YAP_associated_targets", "Ankrd1", "Ankrd1", "Yui2018",
    "YAP_associated_targets", "Ereg", "Ereg", "Yui2018",
    "response_context_not_fetal_specific", "Egr1", "Egr1", "context",
    "response_context_not_fetal_specific", "Jun", "Jun", "context",
    "response_context_not_fetal_specific", "Fos", "Fos", "context",
    "response_context_not_fetal_specific", "Atf3", "Atf3", "context",
    "response_context_not_fetal_specific", "Errfi1", "Errfi1", "context",
    "response_context_not_fetal_specific", "Ifrd1", "Ifrd1", "context",
    "ISC_context", "Lgr5", "Lgr5", "Step05E",
    "ISC_context", "Smoc2", "Smoc2", "Step05E",
    "ISC_context", "Lrig1", "Lrig1", "Step05E",
    "ISC_context", "Ascl2", "Ascl2", "Step05E",
    "ISC_context", "Rnf43", "Rnf43", "Step05E",
    "regional_context", "Olfm4", "Olfm4", "Step05E",
    "secretory_context", "Muc2", "Muc2", "Step05E",
    "secretory_context", "Reg4", "Reg4", "Step05E",
    "sex_associated_context", "Xist", "Xist", "context",
    "sex_associated_context", "Ddx3y", "Ddx3y", "context",
    "sex_associated_context", "Kdm5d", "Kdm5d", "context")
}
resolve_panel <- function(panel, available) {
  hits <- lapply(strsplit(panel$aliases, ";", fixed = TRUE),
    function(v) v[v %in% available])
  panel$matched_feature <- vapply(hits, function(v)
    if (length(v)) v[[1]] else NA_character_, character(1))
  panel$all_matching_aliases <- vapply(hits, paste, character(1), collapse = ";")
  panel$multiple_alias_features <- lengths(hits) > 1L
  panel
}
gene_audit <- function(counts, data, metadata, panel_table, group_column) {
  require_columns(panel_table, c("panel", "gene", "matched_feature"),
    "resolved sentinel panel")
  groups <- unique(metadata[c(group_column, "diet")])
  dplyr::bind_rows(lapply(seq_len(nrow(groups)), function(i) {
    group_value <- as.character(groups[[group_column]][i]); diet_value <- groups$diet[i]
    ix <- which(metadata[[group_column]] == group_value & metadata$diet == diet_value)
    dplyr::bind_rows(lapply(seq_len(nrow(panel_table)), function(j) {
      feature <- panel_table$matched_feature[j]; present <- !is.na(feature)
      panel_name <- panel_table$panel[j]
      gene_name <- panel_table$gene[j]
      tibble::tibble(group = group_value, diet = diet_value, nuclei = length(ix),
        panel = panel_name, gene = gene_name, matched_feature = feature,
        feature_available = present,
        detected_n = if (present) sum(counts[feature, ix] > 0) else NA_integer_,
        percent_detected = if (present) 100 * mean(counts[feature, ix] > 0) else NA_real_,
        mean_log_normalized = if (present) mean(data[feature, ix]) else NA_real_)
    }))
  }))
}
save_vector_plot <- function(plot, directory, filename, width = 12, height = 7) {
  print(plot)
  ggplot2::ggsave(file.path(directory, filename), plot = plot,
    device = grDevices::pdf, useDingbats = FALSE, width = width, height = height)
}
capture_manifest <- function(paths, directory, filename = "input_manifest.csv") {
  need(all(file.exists(paths)), "One or more manifest inputs are absent.")
  info <- file.info(paths)
  out <- tibble::tibble(path = normalizePath(paths), bytes = info$size,
    modified = as.character(info$mtime), md5 = unname(tools::md5sum(paths)))
  write_table(out, directory, filename); out
}
finish_run <- function(directory, filename, lines) {
  writeLines(c(paste("Release:", FOLLOWUP_RELEASE),
    paste("Completed:", format(Sys.time(), tz = "UTC", usetz = TRUE)), lines),
    file.path(directory, filename))
  capture.output(sessionInfo(), file = file.path(directory, "sessionInfo.txt"))
}

# Descriptive pooled-library effect calculations. No biological P values.
pooled_effect <- function(counts, metadata, indices, modality, features = rownames(counts),
                          pc = 0.1, total_depth = Matrix::colSums(counts)) {
  con <- indices[metadata$diet[indices] == "CON"]
  hfd <- indices[metadata$diet[indices] == "HFD"]
  need(length(con) > 0L && length(hfd) > 0L, "A contrast lacks one diet.")
  z <- counts[features, , drop = FALSE]
  a <- Matrix::rowSums(z[, con, drop = FALSE]); b <- Matrix::rowSums(z[, hfd, drop = FALSE])
  da <- sum(total_depth[con]); db <- sum(total_depth[hfd])
  need(da > 0 && db > 0, "Zero pooled library depth.")
  pa <- 100 * Matrix::rowMeans(z[, con, drop = FALSE] > 0)
  pb <- 100 * Matrix::rowMeans(z[, hfd, drop = FALSE] > 0)
  ca <- 1e6 * a / da; cb <- 1e6 * b / db
  min_counts <- if (modality == "RNA") 20 else 30
  min_detect <- if (modality == "RNA") 5 else 1
  tibble::tibble(feature = features, modality = modality, CON_n = length(con), HFD_n = length(hfd),
    CON_counts = a, HFD_counts = b, CON_CPM = ca, HFD_CPM = cb,
    CON_percent = pa, HFD_percent = pb, detection_delta_pp = pb - pa,
    log2FC_HFD_vs_CON = log2((cb + pc)/(ca + pc)),
    evaluable = (a + b >= min_counts) & pmax(pa, pb) >= min_detect) |>
    dplyr::mutate(clear_effect = evaluable & abs(log2FC_HFD_vs_CON) >= 0.5 &
      (modality == "ATAC" | abs(detection_delta_pp) >= 5),
      direction = dplyr::case_when(log2FC_HFD_vs_CON > 0 ~ "HFD_higher",
        log2FC_HFD_vs_CON < 0 ~ "CON_higher", TRUE ~ "unchanged"))
}
balance_draw <- function(indices, diet, depth = NULL, bins = 5L, max_cells = 500L) {
  # Match cell counts in shared depth strata, without replacement. NOT replicates.
  strata <- rep("all", length(indices))
  if (!is.null(depth)) {
    breaks <- unique(stats::quantile(log1p(depth[indices]),
      probs = seq(0, 1, length.out = bins + 1L), na.rm = TRUE))
    if (length(breaks) > 1L) strata <- as.character(cut(log1p(depth[indices]),
      breaks = breaks, include.lowest = TRUE))
  }
  pools <- lapply(unique(strata), function(s) {
    a <- indices[strata == s & diet[indices] == "CON"]
    b <- indices[strata == s & diet[indices] == "HFD"]
    list(a = a, b = b, n = min(length(a), length(b)))
  })
  ns <- vapply(pools, function(x) x$n, integer(1)); total <- sum(ns)
  if (total == 0L) return(integer())
  if (total > max_cells) {
    target <- ns * max_cells / total; ns <- floor(target)
    left <- max_cells - sum(ns)
    if (left > 0L) {
      take <- order(target - ns, decreasing = TRUE)[seq_len(left)]
      ns[take] <- ns[take] + 1L
    }
  }
  # sample.int avoids sample(single_integer) interpreting it as 1:n.
  unlist(lapply(seq_along(pools), function(j) {
    if (ns[j] == 0L) return(integer())
    c(pools[[j]]$a[sample.int(length(pools[[j]]$a), ns[j])],
      pools[[j]]$b[sample.int(length(pools[[j]]$b), ns[j])])
  }), use.names = FALSE)
}
resampling_effects <- function(counts, metadata, indices, full, features,
                               modality, iterations = 20L) {
  selected <- intersect(features, rownames(counts)); depth <- Matrix::colSums(counts)
  if (!length(selected)) return(tibble::tibble())
  x <- counts[selected, , drop = FALSE]
  dplyr::bind_rows(lapply(c("equal_cell_number", "depth_matched"), function(method) {
    draws <- matrix(NA_real_, nrow = length(selected), ncol = iterations)
    ns <- integer(iterations)
    for (i in seq_len(iterations)) {
      ix <- balance_draw(indices, metadata$diet,
        depth = if (method == "depth_matched") depth else NULL)
      ns[i] <- length(ix) %/% 2L
      if (ns[i] < 30L) next
      q <- pooled_effect(x, metadata, ix, modality, total_depth = depth)
      draws[, i] <- q$log2FC_HFD_vs_CON
    }
    anchor <- full$log2FC_HFD_vs_CON[match(selected, full$feature)]
    valid <- rowSums(is.finite(draws))
    stat <- function(fun) apply(draws, 1, function(v)
      if (all(!is.finite(v))) NA_real_ else fun(v[is.finite(v)]))
    tibble::tibble(feature = selected, method = method, iterations = iterations,
      valid_iterations = valid, sampled_n_per_diet_min = min(ns),
      sampled_n_per_diet_max = max(ns),
      median_log2FC = stat(stats::median),
      q10_log2FC = stat(function(v) stats::quantile(v, .1, names = FALSE)),
      q90_log2FC = stat(function(v) stats::quantile(v, .9, names = FALSE)),
      sign_agreement_percent = ifelse(valid > 0,
        100 * rowSums(sign(draws) == sign(anchor), na.rm = TRUE)/pmax(valid, 1), NA_real_))
  }))
}
standardize_strata <- function(counts, metadata, indices, strata, modality,
                               features, min_per_diet = 10L) {
  # Equal stratum weights define a sensitivity estimand, not an adjusted causal effect.
  depth <- Matrix::colSums(counts); strata <- as.character(strata)
  lev <- sort(unique(strata[indices][!is.na(strata[indices])]))
  audit <- dplyr::bind_rows(lapply(lev, function(s) tibble::tibble(stratum = s,
    CON_n = sum(strata[indices] == s & metadata$diet[indices] == "CON", na.rm = TRUE),
    HFD_n = sum(strata[indices] == s & metadata$diet[indices] == "HFD", na.rm = TRUE))))
  if (!nrow(audit)) return(list(effects = tibble::tibble(), audit = audit))
  audit$retained <- pmin(audit$CON_n, audit$HFD_n) >= min_per_diet
  use <- audit$stratum[audit$retained]
  audit$weight <- ifelse(audit$retained, 1/max(length(use), 1), 0)
  if (!length(use)) return(list(effects = tibble::tibble(), audit = audit))
  eff <- dplyr::bind_rows(lapply(use, function(s) {
    ix <- indices[which(strata[indices] == s)]
    pooled_effect(counts, metadata, ix, modality, features, total_depth = depth) |>
      dplyr::mutate(stratum = s)
  })) |>
    dplyr::group_by(feature) |>
    dplyr::summarise(CON_CPM = mean(CON_CPM), HFD_CPM = mean(HFD_CPM),
      CON_percent = mean(CON_percent), HFD_percent = mean(HFD_percent),
      strata_used = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(log2FC_HFD_vs_CON = log2((HFD_CPM + .1)/(CON_CPM + .1)),
      detection_delta_pp = HFD_percent - CON_percent)
  list(effects = eff, audit = audit)
}
promoter_peak_map <- function(object, window = 2000L) {
  # Gene-level annotation only: never infer a TSS from arbitrary exon rows.
  ann <- as.data.frame(Signac::Annotation(object[["ATAC.crypt"]]))
  require_columns(ann, c("seqnames", "start", "end", "strand", "gene_name", "type"),
    "ATAC gene annotation")
  ann <- ann[!is.na(ann$type) & ann$type == "gene" &
    !is.na(ann$gene_name) & nzchar(ann$gene_name) & ann$strand %in% c("+", "-"), ]
  need(nrow(ann) > 0L, "No gene-level TSS annotations; supply a gene-level GRCm39 annotation before promoter integration.")
  ann$tss <- ifelse(ann$strand == "-", ann$end, ann$start)
  ann <- unique(ann[c("seqnames", "strand", "gene_name", "tss")])
  pr <- GenomicRanges::GRanges(as.character(ann$seqnames),
    IRanges::IRanges(pmax(1, ann$tss - window), ann$tss + window), strand = ann$strand)
  peaks <- GenomicRanges::granges(object[["ATAC.crypt"]])
  features <- rownames(object[["ATAC.crypt"]])
  need(length(peaks) == length(features), "ATAC feature/range count mismatch.")
  hyphen_labels <- paste(as.character(GenomicRanges::seqnames(peaks)),
    BiocGenerics::start(peaks), BiocGenerics::end(peaks), sep = "-")
  colon_labels <- paste0(as.character(GenomicRanges::seqnames(peaks)), ":",
    BiocGenerics::start(peaks), "-", BiocGenerics::end(peaks))
  need(identical(features, hyphen_labels) || identical(features, colon_labels),
    "ATAC ranges and count feature names are not identically ordered.")
  hits <- GenomicRanges::findOverlaps(peaks, pr, ignore.strand = TRUE)
  a <- S4Vectors::queryHits(hits); b <- S4Vectors::subjectHits(hits)
  center <- (BiocGenerics::start(peaks)[a] + BiocGenerics::end(peaks)[a])/2
  map <- tibble::tibble(gene = as.character(ann$gene_name[b]), peak = features[a],
    chromosome = as.character(GenomicRanges::seqnames(peaks)[a]),
    peak_start = BiocGenerics::start(peaks)[a], peak_end = BiocGenerics::end(peaks)[a],
    strand = as.character(ann$strand[b]), tss = ann$tss[b],
    promoter_start = pmax(1, ann$tss[b] - window), promoter_end = ann$tss[b] + window,
    peak_center_minus_tss = center - ann$tss[b],
    strand_oriented_center_distance = (center - ann$tss[b]) * ifelse(ann$strand[b] == "-", -1, 1),
    overlap_bp = pmin(BiocGenerics::end(peaks)[a], ann$tss[b] + window) -
      pmax(BiocGenerics::start(peaks)[a], pmax(1, ann$tss[b] - window)) + 1) |>
    dplyr::distinct()
  list(map = map, annotated_genes = unique(as.character(ann$gene_name)))
}
