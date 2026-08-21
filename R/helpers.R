# Shared functions for the stepwise Mana Lab Multiome workflow.
#
# Notebooks contain the biological explanation and decision points. Functions
# here keep repeated implementation details in one place so that a bug fix is
# applied consistently to every step.

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Signac",
  "rtracklayer",
  "sctransform",
  "hdf5r",
  "future",
  "knitr",
  "rmarkdown",
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "Matrix",
  "GenomicRanges",
  "GenomeInfoDb",
  "EnsDb.Mmusculus.v79",
  "ensembldb",
  "biovizBase",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "scDblFinder"
)

minimum_package_versions <- c(
  Seurat = "5.0.0",
  Signac = "1.14.0"
)

qc_numeric_columns <- c(
  "nCount_RNA_min",
  "nCount_RNA_max",
  "nFeature_RNA_min",
  "nFeature_RNA_max",
  "nCount_ATAC_min",
  "nCount_ATAC_max",
  "percent_mt_max",
  "TSS_enrichment_min",
  "nucleosome_signal_max",
  "frip_min"
)

message_step <- function(...) {
  message("\n", paste0(...), "\n", paste(rep("=", 78), collapse = ""))
}

assert_file <- function(path, label = basename(path)) {
  if (!file.exists(path)) {
    visible <- if (dir.exists(dirname(path))) {
      paste(list.files(dirname(path)), collapse = ", ")
    } else {
      "<parent directory does not exist>"
    }
    stop(
      "Missing ", label, ":\n  ", path,
      "\nFiles visible in that directory:\n  ", visible,
      call. = FALSE
    )
  }
  invisible(path)
}

initialize_output_dirs <- function() {
  paths <- file.path(OUTPUT_DIR, RESULT_DIRS)
  invisible(lapply(
    paths,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
  invisible(paths)
}

check_r_version <- function() {
  current <- as.character(getRversion())
  if (!identical(current, EXPECTED_R_VERSION)) {
    warning(
      "This project was tested with R ", EXPECTED_R_VERSION,
      ", but this session is using R ", current,
      ". Package libraries are version-specific on SOL.",
      call. = FALSE
    )
  }
  invisible(current)
}

package_preflight <- function(packages = required_packages, stop_on_failure = TRUE) {
  status <- lapply(packages, function(package_name) {
    tryCatch(
      {
        available <- requireNamespace(package_name, quietly = TRUE)

        if (!available) {
          stop("Package is not installed or its namespace could not be loaded.")
        }

        data.frame(
          package = package_name,
          version = as.character(utils::packageVersion(package_name)),
          status = "OK",
          error = NA_character_,
          stringsAsFactors = FALSE
        )
      },
      error = function(error_condition) {
        data.frame(
          package = package_name,
          version = NA_character_,
          status = "FAILED",
          error = conditionMessage(error_condition),
          stringsAsFactors = FALSE
        )
      }
    )
  })

  status <- do.call(rbind, status)

  failures <- status[status$status != "OK", , drop = FALSE]

  if (stop_on_failure && nrow(failures) > 0L) {
    failure_text <- paste0(
      failures$package,
      ": ",
      failures$error,
      collapse = "\n  "
    )

    stop(
      "Package preflight failed. Resolve these namespace errors before ",
      "continuing:\n  ", failure_text,
      "\n\nRun setup/install_packages.R in the same R 4.4.2 library.",
      call. = FALSE
    )
  }

  if (nrow(failures) == 0L) {
    too_old <- names(minimum_package_versions)[vapply(
      names(minimum_package_versions),
      function(package_name) {
        utils::packageVersion(package_name) <
          minimum_package_versions[[package_name]]
      },
      logical(1)
    )]

    if (length(too_old) > 0L) {
      stop(
        "Package version is below the supported minimum: ",
        paste(
          paste0(
            too_old,
            " >= ",
            minimum_package_versions[too_old]
          ),
          collapse = ", "
        ),
        call. = FALSE
      )
    }
  }

  status
}

configure_parallel_plan <- function() {
  options(
    future.globals.maxSize = FUTURE_GLOBALS_GB * 1024^3,
    Seurat.object.assay.version = "v5"
  )

  if (N_WORKERS > 1L) {
    # Forking is unsafe inside RStudio; multisession launches independent R
    # processes and respects the CPU allocation made in Open OnDemand.
    future::plan(future::multisession, workers = N_WORKERS)
  } else {
    future::plan(future::sequential)
  }
  invisible(future::plan())
}

save_run_information <- function(package_status = NULL) {
  initialize_output_dirs()

  writeLines(
    capture.output(sessionInfo()),
    file.path(OUTPUT_DIR, "00_run_info", "sessionInfo.txt")
  )

  allocation <- data.frame(
    item = c(
      "timestamp",
      "hostname",
      "R_version",
      "R_executable",
      "R_LIBS_USER",
      "SLURM_JOB_ID",
      "SLURM_CPUS_PER_TASK",
      "configured_workers"
    ),
    value = c(
      format(Sys.time(), tz = "America/Phoenix"),
      Sys.info()[["nodename"]],
      R.version.string,
      Sys.which("R"),
      Sys.getenv("R_LIBS_USER"),
      Sys.getenv("SLURM_JOB_ID"),
      Sys.getenv("SLURM_CPUS_PER_TASK"),
      as.character(N_WORKERS)
    ),
    stringsAsFactors = FALSE
  )
  readr::write_csv(
    allocation,
    file.path(OUTPUT_DIR, "00_run_info", "compute_environment.csv")
  )

  if (!is.null(package_status)) {
    readr::write_csv(
      package_status,
      file.path(OUTPUT_DIR, "00_run_info", "package_preflight.csv")
    )
  }
  invisible(allocation)
}

read_sample_sheet <- function(path = SAMPLE_SHEET_PATH) {
  assert_file(path, "sample sheet")
  sample_sheet <- suppressMessages(readr::read_csv(
    path,
    show_col_types = FALSE,
    progress = FALSE
  ))

  required_columns <- c(
    "sample_id",
    "pool_id",
    "genotype",
    "diet",
    "outs_dir"
  )
  missing_columns <- setdiff(required_columns, colnames(sample_sheet))
  if (length(missing_columns) > 0L) {
    stop(
      "The sample sheet is missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  sample_sheet <- as.data.frame(sample_sheet, stringsAsFactors = FALSE)
  if (anyDuplicated(sample_sheet$sample_id)) {
    stop("Every sample_id must be unique.", call. = FALSE)
  }
  if (anyDuplicated(sample_sheet$pool_id)) {
    warning(
      "Repeated pool_id values indicate multiple libraries derived from ",
      "the same biological pool.",
      call. = FALSE
    )
  }

  sample_sheet$h5 <- file.path(
    sample_sheet$outs_dir,
    "filtered_feature_bc_matrix.h5"
  )
  sample_sheet$fragments <- file.path(
    sample_sheet$outs_dir,
    "atac_fragments.tsv.gz"
  )
  sample_sheet$fragment_index <- paste0(sample_sheet$fragments, ".tbi")
  sample_sheet$peaks <- file.path(sample_sheet$outs_dir, "atac_peaks.bed")
  sample_sheet$barcode_metrics <- file.path(
    sample_sheet$outs_dir,
    "per_barcode_metrics.csv"
  )
  sample_sheet$web_summary <- file.path(sample_sheet$outs_dir, "web_summary.html")
  sample_sheet
}

validate_cellranger_inputs <- function(sample_sheet) {
  input_columns <- c(
    "h5",
    "fragments",
    "fragment_index",
    "peaks",
    "barcode_metrics"
  )

  report <- do.call(
    rbind,
    lapply(seq_len(nrow(sample_sheet)), function(row_index) {
      do.call(
        rbind,
        lapply(input_columns, function(column_name) {
          path <- sample_sheet[[column_name]][row_index]
          data.frame(
            sample_id = sample_sheet$sample_id[row_index],
            file_type = column_name,
            path = path,
            exists = file.exists(path),
            size_gb = if (file.exists(path)) {
              unname(file.info(path)$size / 1024^3)
            } else {
              NA_real_
            },
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )

  missing <- report[!report$exists, , drop = FALSE]
  if (nrow(missing) > 0L) {
    stop(
      "Required Cell Ranger ARC inputs are missing:\n  ",
      paste(
        paste0(missing$sample_id, " [", missing$file_type, "]: ", missing$path),
        collapse = "\n  "
      ),
      call. = FALSE
    )
  }
  report
}

build_grcm39_annotation <- function(gtf_path = GRCM39_GTF) {
  assert_file(gtf_path, "GRCm39 GENCODE vM33 GTF")

  annotation <- rtracklayer::import(gtf_path)

  if (!inherits(annotation, "GRanges")) {
    stop(
      "Imported GTF did not produce a GRanges object.",
      call. = FALSE
    )
  }

  # GENCODE uses `gene_type`, while Signac expects `gene_biotype`.
  if (
    !"gene_biotype" %in% colnames(S4Vectors::mcols(annotation)) &&
    "gene_type" %in% colnames(S4Vectors::mcols(annotation))
  ) {
    annotation$gene_biotype <- annotation$gene_type
  }

  required_annotation_columns <- c(
    "gene_name",
    "gene_id",
    "gene_biotype",
    "type"
  )

  missing_columns <- setdiff(
    required_annotation_columns,
    colnames(S4Vectors::mcols(annotation))
  )

  has_transcript_id <- any(
    c("tx_id", "transcript_id") %in%
      colnames(S4Vectors::mcols(annotation))
  )

  if (length(missing_columns) > 0L || !has_transcript_id) {
    stop(
      "GTF annotation is missing columns required by Signac. ",
      "Missing: ",
      paste(missing_columns, collapse = ", "),
      if (!has_transcript_id) {
        "; no tx_id or transcript_id column was found."
      } else {
        ""
      },
      call. = FALSE
    )
  }

  GenomeInfoDb::genome(annotation) <- GRCM39_GENOME_LABEL

  annotation
}

extract_10x_modalities <- function(h5_path) {
  counts <- Seurat::Read10X_h5(
    filename = h5_path,
    use.names = TRUE,
    unique.features = TRUE
  )
  if (!is.list(counts)) {
    stop(
      "Expected a multi-feature H5 file but received one matrix: ", h5_path,
      call. = FALSE
    )
  }

  expected <- c("Gene Expression", "Peaks")
  missing <- setdiff(expected, names(counts))
  if (length(missing) > 0L) {
    stop(
      "Missing feature type(s) in ", h5_path, ": ",
      paste(missing, collapse = ", "),
      "\nAvailable feature types: ", paste(names(counts), collapse = ", "),
      call. = FALSE
    )
  }
  counts
}

add_cellranger_metrics <- function(object, metrics_path) {
  metrics <- suppressMessages(readr::read_csv(
    metrics_path,
    show_col_types = FALSE,
    progress = FALSE
  ))
  if (!"barcode" %in% colnames(metrics)) {
    stop("No barcode column was found in ", metrics_path, call. = FALSE)
  }

  metrics <- as.data.frame(metrics, stringsAsFactors = FALSE)
  rownames(metrics) <- metrics$barcode
  metrics$barcode <- NULL
  colnames(metrics) <- paste0("cr_", colnames(metrics))

  absent_barcodes <- setdiff(colnames(object), rownames(metrics))
  if (length(absent_barcodes) > 0L) {
    stop(
      length(absent_barcodes),
      " filtered-matrix barcodes were absent from per_barcode_metrics.csv.",
      call. = FALSE
    )
  }

  object <- Seurat::AddMetaData(
    object,
    metadata = metrics[colnames(object), , drop = FALSE]
  )

  required_frip_columns <- c(
    "cr_atac_peak_region_fragments",
    "cr_atac_fragments"
  )
  if (all(required_frip_columns %in% colnames(object[[]]))) {
    object$frip_cr <- with(
      object[[]],
      ifelse(
        cr_atac_fragments > 0,
        cr_atac_peak_region_fragments / cr_atac_fragments,
        NA_real_
      )
    )
  } else {
    stop(
      "Cell Ranger fragment columns required for frip_cr were not found. ",
      "Available cr_ columns: ",
      paste(grep("^cr_", colnames(object[[]]), value = TRUE), collapse = ", "),
      call. = FALSE
    )
  }
  object
}

run_scdblfinder <- function(object) {
  rna_counts <- SeuratObject::LayerData(
    object,
    assay = "RNA",
    layer = "counts"
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = rna_counts)
  )

  set.seed(20260727)
  sce <- scDblFinder::scDblFinder(sce, clusters = FALSE)

  doublet_columns <- grep(
    "^scDblFinder",
    colnames(SummarizedExperiment::colData(sce)),
    value = TRUE
  )
  metadata <- as.data.frame(
    SummarizedExperiment::colData(sce)[, doublet_columns, drop = FALSE]
  )
  Seurat::AddMetaData(object, metadata = metadata)
}

robust_log_lower <- function(x, hard_floor, n_mad = 4) {
  transformed <- log10(x[is.finite(x) & x >= 0] + 1)
  if (length(transformed) < 20L) return(hard_floor)
  spread <- stats::mad(
    transformed,
    center = stats::median(transformed),
    constant = 1
  )
  if (!is.finite(spread) || spread == 0) return(hard_floor)
  max(hard_floor, 10^(stats::median(transformed) - n_mad * spread) - 1)
}

robust_log_upper <- function(x, hard_ceiling, n_mad = 4) {
  transformed <- log10(x[is.finite(x) & x >= 0] + 1)
  if (length(transformed) < 20L) return(hard_ceiling)
  spread <- stats::mad(
    transformed,
    center = stats::median(transformed),
    constant = 1
  )
  if (!is.finite(spread) || spread == 0) return(hard_ceiling)
  min(hard_ceiling, 10^(stats::median(transformed) + n_mad * spread) - 1)
}

make_qc_threshold_suggestion <- function(object, sample_id) {
  metadata <- object[[]]
  data.frame(
    sample_id = sample_id,
    nCount_RNA_min = ceiling(robust_log_lower(
      metadata$nCount_RNA,
      QC_LIMITS$nCount_RNA_min
    )),
    nCount_RNA_max = floor(robust_log_upper(
      metadata$nCount_RNA,
      QC_LIMITS$nCount_RNA_max
    )),
    nFeature_RNA_min = ceiling(robust_log_lower(
      metadata$nFeature_RNA,
      QC_LIMITS$nFeature_RNA_min
    )),
    nFeature_RNA_max = floor(robust_log_upper(
      metadata$nFeature_RNA,
      QC_LIMITS$nFeature_RNA_max
    )),
    nCount_ATAC_min = ceiling(robust_log_lower(
      metadata$nCount_ATAC,
      QC_LIMITS$nCount_ATAC_min
    )),
    nCount_ATAC_max = floor(robust_log_upper(
      metadata$nCount_ATAC,
      QC_LIMITS$nCount_ATAC_max
    )),
    percent_mt_max = QC_LIMITS$percent_mt_max,
    TSS_enrichment_min = QC_LIMITS$TSS_enrichment_min,
    nucleosome_signal_max = QC_LIMITS$nucleosome_signal_max,
    frip_min = QC_LIMITS$frip_min,
    stringsAsFactors = FALSE
  )
}

populate_qc_decision_draft <- function(suggestions, decision_path = QC_DECISION_PATH) {
  assert_file(decision_path, "QC decision table")
  decisions <- suppressMessages(readr::read_csv(
    decision_path,
    show_col_types = FALSE,
    progress = FALSE
  ))
  decisions <- as.data.frame(decisions, stringsAsFactors = FALSE)

  missing_samples <- setdiff(suggestions$sample_id, decisions$sample_id)
  extra_samples <- setdiff(decisions$sample_id, suggestions$sample_id)
  if (length(missing_samples) > 0L || length(extra_samples) > 0L) {
    stop(
      "Sample IDs differ between suggestions and the active QC decision table: ",
      decision_path,
      ". Missing in decision file: ",
      paste(missing_samples, collapse = ", "),
      "; extra in decision file: ",
      paste(extra_samples, collapse = ", "),
      call. = FALSE
    )
  }

  suggestions <- suggestions[match(decisions$sample_id, suggestions$sample_id), ]
  for (column_name in qc_numeric_columns) {
    decisions[[column_name]] <- suppressWarnings(as.numeric(decisions[[column_name]]))
    replace <- !is.finite(decisions[[column_name]])
    decisions[[column_name]][replace] <- suggestions[[column_name]][replace]
  }

  readr::write_csv(decisions, decision_path, na = "")
  decisions
}

read_and_validate_qc_decisions <- function(sample_sheet) {
  decisions <- suppressMessages(readr::read_csv(
    QC_DECISION_PATH,
    show_col_types = FALSE,
    progress = FALSE
  ))
  decisions <- as.data.frame(decisions, stringsAsFactors = FALSE)

  required <- c(
    "sample_id",
    qc_numeric_columns,
    "approved",
    "reviewer",
    "review_date",
    "decision_notes"
  )
  missing_columns <- setdiff(required, colnames(decisions))
  if (length(missing_columns) > 0L) {
    stop(
      "QC decision table is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      ". File: ",
      QC_DECISION_PATH,
      call. = FALSE
    )
  }

  missing_samples <- setdiff(sample_sheet$sample_id, decisions$sample_id)
  extra_samples <- setdiff(decisions$sample_id, sample_sheet$sample_id)
  if (length(missing_samples) > 0L || length(extra_samples) > 0L) {
    stop(
      "Sample IDs in the QC decision table do not match the active sample sheet: ",
      SAMPLE_SHEET_PATH,
      call. = FALSE
    )
  }

  for (column_name in qc_numeric_columns) {
    decisions[[column_name]] <- suppressWarnings(as.numeric(decisions[[column_name]]))
  }
  if (any(!is.finite(as.matrix(decisions[, qc_numeric_columns, drop = FALSE])))) {
    stop(
      "Every QC threshold must contain a finite numeric value.",
      call. = FALSE
    )
  }

  approved <- tolower(trimws(as.character(decisions$approved))) %in%
    c("true", "t", "1", "yes", "y")
  if (!all(approved)) {
    stop(
      "QC is paused intentionally. Review the plots, then set approved=TRUE ",
      "for every sample in config/qc_thresholds.csv and rerun this chunk.",
      call. = FALSE
    )
  }
  required_audit_fields <- c(
    "reviewer",
    "review_date",
    "decision_notes"
  )

  missing_audit <- vapply(
    required_audit_fields,
    function(column_name) {
      values <- trimws(as.character(decisions[[column_name]]))
      any(is.na(values) | values == "")
    },
    logical(1)
  )

  if (any(missing_audit)) {
    stop(
      "Approved QC decisions require completed audit fields: ",
      paste(
        required_audit_fields[missing_audit],
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  decisions$approved <- approved
  decisions[match(sample_sheet$sample_id, decisions$sample_id), , drop = FALSE]
}

add_joint_qc_flags <- function(object, thresholds) {
  metadata <- object[[]]

  object$qc_pass_rna <- with(
    metadata,
    nCount_RNA >= thresholds$nCount_RNA_min &
      nCount_RNA <= thresholds$nCount_RNA_max &
      nFeature_RNA >= thresholds$nFeature_RNA_min &
      nFeature_RNA <= thresholds$nFeature_RNA_max &
      percent.mt <= thresholds$percent_mt_max
  )

  object$qc_pass_atac <- with(
    metadata,
    nCount_ATAC >= thresholds$nCount_ATAC_min &
      nCount_ATAC <= thresholds$nCount_ATAC_max &
      TSS.enrichment >= thresholds$TSS_enrichment_min &
      nucleosome_signal <= thresholds$nucleosome_signal_max &
      frip_cr >= thresholds$frip_min
  )

  object$qc_pass_doublet <- metadata$scDblFinder.class == "singlet"
  object$qc_pass <- (
    object$qc_pass_rna &
      object$qc_pass_atac &
      object$qc_pass_doublet
  )

  for (column_name in c(
    "qc_pass_rna",
    "qc_pass_atac",
    "qc_pass_doublet",
    "qc_pass"
  )) {
    value <- object[[column_name, drop = TRUE]]
    value[is.na(value)] <- FALSE
    object[[column_name]] <- value
  }
  object
}

summarize_qc_decision <- function(object, sample_id) {
  metadata <- object[[]]
  data.frame(
    sample_id = sample_id,
    evaluated_barcodes = nrow(metadata),
    pass_rna = sum(metadata$qc_pass_rna),
    pass_atac = sum(metadata$qc_pass_atac),
    pass_doublet = sum(metadata$qc_pass_doublet),
    pass_joint = sum(metadata$qc_pass),
    percent_joint_retained = 100 * mean(metadata$qc_pass),
    fail_rna_only = sum(
      !metadata$qc_pass_rna & metadata$qc_pass_atac & metadata$qc_pass_doublet
    ),
    fail_atac_only = sum(
      metadata$qc_pass_rna & !metadata$qc_pass_atac & metadata$qc_pass_doublet
    ),
    fail_multiple_or_doublet = sum(
      !metadata$qc_pass & !(
        xor(!metadata$qc_pass_rna, !metadata$qc_pass_atac) &
          metadata$qc_pass_doublet
      )
    ),
    stringsAsFactors = FALSE
  )
}

save_sample_qc_plots <- function(
    metadata,
    sample_id,
    output_prefix,
    thresholds = NULL) {
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

  metric_order <- c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "nCount_ATAC",
    "TSS.enrichment",
    "nucleosome_signal",
    "frip_cr"
  )
  long <- tidyr::pivot_longer(
    metadata,
    cols = dplyr::all_of(metric_order),
    names_to = "metric",
    values_to = "value"
  )
  long$metric <- factor(long$metric, levels = metric_order)

  violin <- ggplot2::ggplot(long, ggplot2::aes(x = sample_id, y = value)) +
    ggplot2::geom_violin(fill = "#4C78A8", color = "grey25", scale = "width") +
    ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 4) +
    ggplot2::labs(
      title = paste(sample_id, "QC distributions"),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())

  rna <- ggplot2::ggplot(
    metadata,
    ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)
  ) +
    ggplot2::geom_point(alpha = 0.25, size = 0.35) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(
      title = paste(sample_id, "RNA complexity"),
      color = "% mitochondrial"
    ) +
    ggplot2::theme_classic(base_size = 11)

  atac <- ggplot2::ggplot(
    metadata,
    ggplot2::aes(
      x = nCount_ATAC,
      y = TSS.enrichment,
      color = nucleosome_signal
    )
  ) +
    ggplot2::geom_point(alpha = 0.25, size = 0.35) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(
      title = paste(sample_id, "ATAC quality"),
      color = "nucleosome signal"
    ) +
    ggplot2::theme_classic(base_size = 11)

  if (!is.null(thresholds)) {
    rna <- rna +
      ggplot2::geom_vline(
        xintercept = c(thresholds$nCount_RNA_min, thresholds$nCount_RNA_max),
        linetype = "dashed"
      ) +
      ggplot2::geom_hline(
        yintercept = c(
          thresholds$nFeature_RNA_min,
          thresholds$nFeature_RNA_max
        ),
        linetype = "dashed"
      )
    atac <- atac +
      ggplot2::geom_vline(
        xintercept = c(thresholds$nCount_ATAC_min, thresholds$nCount_ATAC_max),
        linetype = "dashed"
      ) +
      ggplot2::geom_hline(
        yintercept = thresholds$TSS_enrichment_min,
        linetype = "dashed"
      )
  }

  ggplot2::ggsave(
    paste0(output_prefix, "_violin.pdf"),
    violin,
    width = 13,
    height = 7,
    limitsize = FALSE
  )
  ggplot2::ggsave(
    paste0(output_prefix, "_rna_scatter.pdf"),
    rna,
    width = 7,
    height = 5.5,
    limitsize = FALSE
  )
  ggplot2::ggsave(
    paste0(output_prefix, "_atac_scatter.pdf"),
    atac,
    width = 7,
    height = 5.5,
    limitsize = FALSE
  )
  invisible(list(violin = violin, rna = rna, atac = atac))
}