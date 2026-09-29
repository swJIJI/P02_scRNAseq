# ============================================================
# 04_qc_filtering.R
#
# Purpose / 목적:
#   - Apply the approved P02-specific hard QC filters.
#     승인된 P02-specific QC criteria를 적용한다.
#   - Compare per-sample cell numbers before and after filtering.
#
# Input:
#   - scripts/00_config.R
#   - p02 object with percent.mt, created in scripts/03_qc_metrics.R
#
# Output:
#   - p02_qc : QC-filtered P02 Seurat object
#   - qc_filter_summary : per-sample cell retention table
#   - report/session_info/04_qc_filtering_sessionInfo.txt
#
# Hard filters:
#   nFeature_RNA >= QC_MIN_FEATURES
#   percent.mt   <= QC_MAX_MT
#
# nCount_RNA is retained as a diagnostic metric only.
# Upper nFeature/nCount cutoffs are not applied at this stage.
# ============================================================


# ------------------------------------------------------------
# 1. Load dataset-specific configuration
# ------------------------------------------------------------

source("scripts/00_config.R")

set.seed(RANDOM_SEED)


# ------------------------------------------------------------
# 2. Load required packages
# ------------------------------------------------------------

library(Seurat)
library(dplyr)


# ------------------------------------------------------------
# 3. Ensure QC metrics are available
# ------------------------------------------------------------

# 04_qc_filtering.R can be run after 03_qc_metrics.R in the same session.
# If p02 or percent.mt is absent, run 03_qc_metrics.R automatically.
if (
  !exists("p02") ||
  !"percent.mt" %in% colnames(p02@meta.data)
) {
  message(
    "p02 object with percent.mt not found in the current R session. ",
    "Running scripts/03_qc_metrics.R first."
  )

  source("scripts/03_qc_metrics.R")
}

# Hard filter로 사용하는 threshold는 반드시 config에서 확정되어 있어야 한다.
if (is.na(QC_MIN_FEATURES) || is.na(QC_MAX_MT)) {
  stop(
    "QC_MIN_FEATURES and QC_MAX_MT must be set in scripts/00_config.R."
  )
}


# ------------------------------------------------------------
# 4. Apply approved hard QC filters
# ------------------------------------------------------------

qc_keep <- with(
  p02@meta.data,
  nFeature_RNA >= QC_MIN_FEATURES &
    percent.mt <= QC_MAX_MT
)

p02_qc <- subset(
  p02,
  cells = rownames(p02@meta.data)[qc_keep]
)


# ------------------------------------------------------------
# 5. Validate QC-filtered object
# ------------------------------------------------------------

if (ncol(p02_qc) != sum(qc_keep)) {
  stop("Number of retained cells does not match the QC filter result.")
}

if (nrow(p02_qc) != nrow(p02)) {
  stop("Feature count changed during cell-level QC filtering.")
}

# Sample별 filtering 전후 cell 수 비교
qc_filter_summary <- p02@meta.data |>
  dplyr::mutate(qc_keep = qc_keep) |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(
    n_before = dplyr::n(),
    n_after = sum(qc_keep),
    n_removed = n_before - n_after,
    pct_retained = round(100 * n_after / n_before, 2),
    .groups = "drop"
  )

cat(
  "\n===== P02 QC-filtered object =====\n"
)

cat(
  "Number of features:",
  nrow(p02_qc),
  "\n"
)

cat(
  "Number of cells:",
  ncol(p02),
  "->",
  ncol(p02_qc),
  "\n\n"
)

print(qc_filter_summary)


# ------------------------------------------------------------
# 6. Record R / package session information
# ------------------------------------------------------------

dir.create(
  SESSION_INFO_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    SESSION_INFO_DIR,
    "04_qc_filtering_sessionInfo.txt"
  )
)

message("04_qc_filtering.R completed successfully.")
