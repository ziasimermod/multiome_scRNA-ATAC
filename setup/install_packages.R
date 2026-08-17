#!/usr/bin/env Rscript

# One-time package installation/repair for the ASU SOL R 4.4.2 environment.
#
# Run this interactively from the RStudio Console with:
#   source("setup/install_packages.R")
#
# Do not place BiocManager::install() calls inside the analysis notebooks. An
# analysis should stop when its environment is incomplete rather than silently
# altering the package library while data are in memory.

expected_r_version <- "4.4.2"
if (!identical(as.character(getRversion()), expected_r_version)) {
  stop(
    "This setup is for R ", expected_r_version,
    ". Restart the SOL RStudio session with R 4.4.2 before installing.",
    call. = FALSE
  )
}

target_library <- Sys.getenv("R_LIBS_USER")
if (!nzchar(target_library)) target_library <- .libPaths()[1]
dir.create(target_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(target_library, .libPaths())))

lock_directories <- list.files(
  target_library,
  pattern = "^00LOCK",
  full.names = TRUE
)
if (length(lock_directories) > 0L) {
  stop(
    "An incomplete package-installation lock exists:\n  ",
    paste(lock_directories, collapse = "\n  "),
    "\nConfirm that no other R installation is running, remove only the stale ",
    "00LOCK directory, and rerun this setup script.",
    call. = FALSE
  )
}

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_packages <- c(
  "Seurat",
  "SeuratObject",
  "Signac",
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
  # Explicit repair targets for dependency failures previously encountered on
  # the cluster's user library.
  "RANN",
  "spam",
  "cluster",
  "igraph",
  "fastmatch"
)

bioconductor_packages <- c(
  "rtracklayer",
  "EnsDb.Mmusculus.v79",
  "ensembldb",
  "GenomeInfoDb",
  "GenomicRanges",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "scDblFinder",
  "biovizBase"
)

is_loadable <- function(package_name) {
  isTRUE(
    tryCatch(
      requireNamespace(package_name, quietly = TRUE),
      error = function(error_condition) FALSE
    )
  )
}

cran_to_install <- cran_packages[
  !vapply(cran_packages, is_loadable, logical(1))
]
if (length(cran_to_install) > 0L) {
  message("Installing/repairing CRAN packages: ", paste(cran_to_install, collapse = ", "))
  for (package_name in cran_to_install) {
    install.packages(
      package_name,
      lib = target_library,
      dependencies = NA,
      Ncpus = 1L
    )
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    lib = target_library,
    dependencies = NA,
    Ncpus = 1L
  )
}

bioc_to_install <- bioconductor_packages[
  !vapply(bioconductor_packages, is_loadable, logical(1))
]
if (length(bioc_to_install) > 0L) {
  message(
    "Installing/repairing Bioconductor packages: ",
    paste(bioc_to_install, collapse = ", ")
  )
  BiocManager::install(
    bioc_to_install,
    lib = target_library,
    ask = FALSE,
    update = FALSE,
    dependencies = TRUE,
    Ncpus = 1L
  )
}

all_packages <- c(cran_packages, bioconductor_packages)
status <- vapply(all_packages, is_loadable, logical(1))

cat("\nFinal namespace preflight\n")
for (package_name in all_packages) {
  if (status[[package_name]]) {
    cat(
      package_name,
      as.character(utils::packageVersion(package_name)),
      "OK\n"
    )
  } else {
    cat(package_name, "FAILED\n")
  }
}

if (!all(status)) {
  stop(
    "At least one namespace is still not loadable. Copy the first complete ",
    "error from loadNamespace(<package>) before attempting another install.",
    call. = FALSE
  )
}

cat("\nPackage setup completed successfully in:\n", target_library, "\n")
