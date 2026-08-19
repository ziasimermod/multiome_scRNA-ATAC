# Project-specific settings for the first Mana Lab PPAR iKO Multiome dataset.
#
# Keep data paths here rather than scattering them across notebooks. A future
# universal branch can replace these explicit values with command-line or YAML
# configuration without changing the analysis logic.

set.seed(20260727)

PROJECT_DIR <- "/scratch/dsaiz/Yesenia_scData2026"
OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "Ppar_iKO_multiome_stepwise"
)

SAMPLE_SHEET_PATH <- file.path(REPO_ROOT, "config", "samples.csv")
QC_DECISION_PATH <- file.path(REPO_ROOT, "config", "qc_thresholds.csv")

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

