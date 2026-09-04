# Project-specific settings for the Mana Lab colon Multiome datasets.
#
# Dataset versions and biological cohorts are kept explicit so that the
# original sequencing-depth analysis remains reproducible while PPAR, WT,
# and IL17 cohorts can be analyzed with the same pipeline.

set.seed(20260727)

PROJECT_DIR <- "/scratch/dsaiz/Yesenia_scData2026"

# -------------------------------------------------------------------------
# Dataset version and biological cohort
# -------------------------------------------------------------------------

# v1_initial_depth:
#   Original Cell Ranger ARC analysis before additional ATAC sequencing.
#
# v2_resequenced:
#   Cell Ranger ARC rerun using deeper ATAC sequencing and the original
#   GEX/GEM sequencing. This is the current authoritative dataset.

DATASET_VERSION <- Sys.getenv(
  "MULTIOME_DATASET_VERSION",
  unset = "v2_resequenced"
)

VALID_DATASET_VERSIONS <- c(
  "v1_initial_depth",
  "v2_resequenced"
)

if (!DATASET_VERSION %in% VALID_DATASET_VERSIONS) {
  stop(
    "Unknown DATASET_VERSION: ", DATASET_VERSION,
    "\nValid versions: ",
    paste(VALID_DATASET_VERSIONS, collapse = ", ")
  )
}

COHORT_GROUPS <- c(
  ppar = "PPAR",
  wt = "WT",
  il17 = "IL17"
)

VALID_COHORTS_BY_DATASET <- list(
  v1_initial_depth = "ppar",
  v2_resequenced = names(COHORT_GROUPS)
)

COHORT_ID <- tolower(trimws(Sys.getenv(
  "MULTIOME_COHORT_ID",
  unset = "ppar"
)))

VALID_COHORT_IDS <- VALID_COHORTS_BY_DATASET[[DATASET_VERSION]]

if (!COHORT_ID %in% VALID_COHORT_IDS) {
  stop(
    "Unknown COHORT_ID for ", DATASET_VERSION, ": ", COHORT_ID,
    "\nValid cohorts: ",
    paste(VALID_COHORT_IDS, collapse = ", ")
  )
}

COHORT_GROUP <- unname(COHORT_GROUPS[[COHORT_ID]])

DATASET_CONFIG_DIR <- file.path(
  REPO_ROOT,
  "config",
  "datasets",
  DATASET_VERSION
)

SAMPLE_SHEET_PATH <- file.path(
  DATASET_CONFIG_DIR,
  "samples.csv"
)

COHORT_CONFIG_DIR <- file.path(
  DATASET_CONFIG_DIR,
  "cohorts",
  COHORT_ID
)

QC_DECISION_PATH <- file.path(
  COHORT_CONFIG_DIR,
  "qc_thresholds.csv"
)

CLUSTERING_DECISION_PATH <- file.path(
  COHORT_CONFIG_DIR,
  "clustering_decision.csv"
)

ANNOTATION_DECISION_PATH <- file.path(
  COHORT_CONFIG_DIR,
  "cell_state_annotations.csv"
)

OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  DATASET_VERSION,
  "independent",
  COHORT_ID
)

message(
  "Active Multiome run: dataset=", DATASET_VERSION,
  "; cohort=", COHORT_ID,
  " (", COHORT_GROUP, ")",
  "\nOutput directory: ", OUTPUT_DIR
)

# Expected R version on the ASU SOL Open OnDemand RStudio Server.
EXPECTED_R_VERSION <- "4.4.2"

# Use no more workers than were requested for the RStudio session. RStudio uses
# separate background R processes rather than forked processes.
N_WORKERS <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
if (!is.finite(N_WORKERS) || N_WORKERS < 1L) N_WORKERS <- 1L
N_WORKERS <- min(N_WORKERS, 8L)
FUTURE_GLOBALS_GB <- 50

# Only near-empty barcodes are removed before doublet detection. The more
# consequential biological QC decision is delayed until Step 2.
ULTRA_LOW_RNA_UMI <- 200
ULTRA_LOW_RNA_FEATURES <- 100

# Broad safety limits and automated starting suggestions. These are not
# universal biological truths. Step 1 combines these bounds with sample-aware
# median absolute deviation calculations, and Step 2 requires manual review.
QC_LIMITS <- list(
  nCount_RNA_min = 500,
  nCount_RNA_max = 50000,
  nFeature_RNA_min = 250,
  nFeature_RNA_max = 8000,
  nCount_ATAC_min = 1000,
  nCount_ATAC_max = 100000,
  percent_mt_max = 25,
  TSS_enrichment_min = 1,
  nucleosome_signal_max = 4,
  frip_min = 0.10
)

GRCM39_GENOME_LABEL <- "GRCm39"

GRCM39_GTF <- file.path(
  PROJECT_DIR,
  "reference",
  "GRCm39_2024-A",
  "gencode.vM33.primary_assembly.annotation.gtf.gz"
)

RESULT_DIRS <- c(
  "00_run_info",
  "01_qc/checkpoints",
  "01_qc/plots",
  "01_qc/metadata",
  "02_objects",
  "03_common_peaks",
  "03_objects",
  "04_reduction",
  "05_wnn",
  "06_annotation",
  "07_composition",
  "07_pseudobulk"
)
