# Lightweight regression checks for dataset/cohort selection.
# Run from the repository root with:
#   Rscript tests/test_cohort_configuration.R

if (!requireNamespace("readr", quietly = TRUE)) {
  stop("The readr package is required for this configuration test.")
}

run_cohort_configuration_tests <- function() {
original_dataset_version <- Sys.getenv(
  "MULTIOME_DATASET_VERSION",
  unset = NA_character_
)
original_cohort_id <- Sys.getenv(
  "MULTIOME_COHORT_ID",
  unset = NA_character_
)

restore_environment_variable <- function(name, value) {
  if (is.na(value)) {
    Sys.unsetenv(name)
  } else {
    do.call(Sys.setenv, stats::setNames(list(value), name))
  }
}

on.exit({
  restore_environment_variable(
    "MULTIOME_DATASET_VERSION",
    original_dataset_version
  )
  restore_environment_variable(
    "MULTIOME_COHORT_ID",
    original_cohort_id
  )
}, add = TRUE)

REPO_ROOT <- if (file.exists("multiome_scRNA-ATAC.Rproj")) {
  normalizePath(".")
} else if (file.exists(file.path("..", "multiome_scRNA-ATAC.Rproj"))) {
  normalizePath("..")
} else {
  stop("Run this test from the repository root or tests directory.")
}

load_case <- function(dataset_version, cohort_id) {
  Sys.setenv(
    MULTIOME_DATASET_VERSION = dataset_version,
    MULTIOME_COHORT_ID = cohort_id
  )

  test_environment <- new.env(parent = globalenv())
  test_environment$REPO_ROOT <- REPO_ROOT
  sys.source(
    file.path(REPO_ROOT, "config", "project_config.R"),
    envir = test_environment
  )
  sys.source(
    file.path(REPO_ROOT, "R", "helpers.R"),
    envir = test_environment
  )

  sample_sheet <- test_environment$read_sample_sheet()
  qc_decisions <- utils::read.csv(
    test_environment$QC_DECISION_PATH,
    na.strings = c("", "NA"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  stopifnot(
    nrow(sample_sheet) == 2L,
    identical(sample_sheet$diet, c("CON", "HFD")),
    identical(unique(sample_sheet$cohort_id), cohort_id),
    identical(
      sort(sample_sheet$sample_id),
      sort(qc_decisions$sample_id)
    ),
    file.exists(test_environment$QC_DECISION_PATH)
  )

  c(
    dataset_version = dataset_version,
    cohort_id = cohort_id,
    group = unique(sample_sheet$group),
    output_dir = test_environment$OUTPUT_DIR
  )
}

valid_cases <- list(
  c("v1_initial_depth", "ppar"),
  c("v2_resequenced", "ppar"),
  c("v2_resequenced", "wt"),
  c("v2_resequenced", "il17")
)

results <- lapply(valid_cases, function(case) {
  load_case(case[[1]], case[[2]])
})

output_directories <- vapply(results, `[[`, character(1), "output_dir")
stopifnot(anyDuplicated(output_directories) == 0L)

expect_configuration_error <- function(dataset_version, cohort_id) {
  error_seen <- FALSE
  tryCatch(
    load_case(dataset_version, cohort_id),
    error = function(condition) {
      error_seen <<- TRUE
    }
  )
  stopifnot(error_seen)
}

expect_configuration_error("v1_initial_depth", "wt")
expect_configuration_error("v2_resequenced", "unknown")

result_table <- as.data.frame(
  do.call(rbind, results),
  stringsAsFactors = FALSE
)
print(result_table, row.names = FALSE)
cat("\nAll cohort configuration checks passed.\n")
}

run_cohort_configuration_tests()
