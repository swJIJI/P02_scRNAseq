# ============================================================
# P02 (HN00262036)
# Dataset-specific configuration
# ============================================================


# ------------------------------------------------------------
# Project information
# ------------------------------------------------------------

PROJECT_ID <- "P02"
DATASET_ID <- "HN00262036"


# ------------------------------------------------------------
# Reproducibility
# ------------------------------------------------------------

RANDOM_SEED <- 42


# ------------------------------------------------------------
# File paths
# ------------------------------------------------------------

# Raw data are stored locally on the Linux analysis workstation.
# `~` refers to the current user's home directory.
RAW_DATA_DIR <- path.expand(
  "~/data/HN00262036/HN00262036_result_10X"
)

# Repository-internal files
METADATA_FILE <- "metadata/sample_metadata.csv"

FIGURE_DIR <- "figures"
TABLE_DIR  <- "tables"
REPORT_DIR <- "report"


# ------------------------------------------------------------
# QC thresholds
# QC 결과 확인 후 dataset-specific하게 결정
# ------------------------------------------------------------

QC_MIN_FEATURES <- NA
QC_MAX_FEATURES <- NA

QC_MIN_COUNTS <- NA
QC_MAX_COUNTS <- NA

QC_MAX_MT <- NA


# ------------------------------------------------------------
# Doublet detection
# ------------------------------------------------------------

DO_DOUBLET_DETECTION <- TRUE

DOUBLET_METHOD <- NA
EXPECTED_DOUBLET_RATE <- NA


# ------------------------------------------------------------
# Variable features
# ------------------------------------------------------------

N_VARIABLE_FEATURES <- NA


# ------------------------------------------------------------
# PCA
# ------------------------------------------------------------

N_PCS <- NA
PCA_DIMS <- NA


# ------------------------------------------------------------
# UMAP / clustering
# ------------------------------------------------------------

UMAP_DIMS <- NA
CLUSTER_RESOLUTION <- NA


# ------------------------------------------------------------
# Marker analysis
# ------------------------------------------------------------

MARKER_LOGFC_THRESHOLD <- NA
MARKER_MIN_PCT <- NA

