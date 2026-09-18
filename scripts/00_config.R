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



# session information저장경로_코드 재현성을 위한 장치

SESSION_INFO_DIR <- file.path(
  REPORT_DIR,
  "session_info"
)


# ------------------------------------------------------------
# QC thresholds
# ------------------------------------------------------------
# Dataset-specific QC thresholds determined after QC review.
# QC 결과를 바탕으로 승인된 P02-specific filtering 기준을 기록한다.

# Minimum detected genes per cell
# Extreme low-complexity cells를 제거하기 위한 conservative lower bound
QC_MIN_FEATURES <- 500

# No upper feature cutoff at this stage
# High-feature cells는 별도의 abnormal population으로 확인되지 않았으며,
# potential doublets는 dedicated doublet-detection step에서 평가한다.
QC_MAX_FEATURES <- NA

# nCount_RNA is used as a diagnostic metric only.
# D7-associated low-RNA population의 biological vs technical origin이
# unresolved 상태이므로 lower count cutoff는 적용하지 않는다.
QC_MIN_COUNTS <- NA #검토 후 의도적으로 hard cutoff를 사용하지 않기로 결정

# No upper count cutoff at this stage
QC_MAX_COUNTS <- NA

# Maximum mitochondrial transcript percentage
# High-mitochondrial cells를 제거하기 위한 conservative upper threshold
QC_MAX_MT <- 20
 
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

