# ============================================================
# 01_load_data.R
#
# Purpose / 목적:
#   - Load P02 dataset configuration and sample metadata.
#   - Validate sample-to-path mapping.
#   - Read each sample-level 10X filtered count matrix.
#   - Confirm matrix class and dimensions before Seurat object creation.
#
# Input:
#   - scripts/00_config.R
#   - metadata/sample_metadata.csv
#   - sample-level 10X filtered_feature_bc_matrix directories
#
# Output:
#   - sample_metadata : metadata with validated full data paths
#   - count_matrices  : named list of sample-level 10X count matrices
#   - report/session_info/01_load_data_sessionInfo.txt
#
# Important:
#   This script DOES NOT create Seurat objects and DOES NOT apply QC filtering.
# ============================================================


# ------------------------------------------------------------
# 1. Load dataset-specific configuration
# ------------------------------------------------------------

source("scripts/00_config.R")

set.seed(RANDOM_SEED)


# ------------------------------------------------------------
# 2. Load required package
# ------------------------------------------------------------

library(Seurat)


# ------------------------------------------------------------
# 3. Load sample metadata
# ------------------------------------------------------------

sample_metadata <- read.csv(
  METADATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 4. Validate metadata structure
# ------------------------------------------------------------

required_metadata_columns <- c(
  "sample_id",
  "pdo_source",
  "timepoint",
  "biological_replicate",
  "culture_batch",
  "data_path"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(sample_metadata)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    "Missing required metadata column(s): ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

if (anyDuplicated(sample_metadata$sample_id) > 0) {
  stop("Duplicated sample_id values were found in sample_metadata.csv.")
}


# ------------------------------------------------------------
# 5. Construct and validate full 10X data paths
# ------------------------------------------------------------

sample_metadata$full_data_path <- file.path(
  RAW_DATA_DIR,
  sample_metadata$data_path
)

sample_metadata$data_path_exists <- dir.exists(
  sample_metadata$full_data_path
)

print(
  sample_metadata[
    ,
    c(
      "sample_id",
      "full_data_path",
      "data_path_exists"
    )
  ]
)

if (!all(sample_metadata$data_path_exists)) {
  
  missing_samples <- sample_metadata$sample_id[
    !sample_metadata$data_path_exists
  ]
  
  stop(
    "10X data directory not found for sample(s): ",
    paste(missing_samples, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 6. Validate required 10X matrix files
# ------------------------------------------------------------

required_10x_files <- c(
  "barcodes.tsv.gz",
  "features.tsv.gz",
  "matrix.mtx.gz"
)

missing_10x_files <- lapply(
  seq_len(nrow(sample_metadata)),
  function(i) {
    
    sample_id <- sample_metadata$sample_id[i]
    sample_dir <- sample_metadata$full_data_path[i]
    
    required_paths <- file.path(
      sample_dir,
      required_10x_files
    )
    
    missing_files <- required_10x_files[
      !file.exists(required_paths)
    ]
    
    if (length(missing_files) == 0) {
      return(NULL)
    }
    
    paste0(
      sample_id,
      ": ",
      paste(missing_files, collapse = ", ")
    )
  }
)

missing_10x_files <- unlist(missing_10x_files)

if (length(missing_10x_files) > 0) {
  stop(
    "Required 10X file(s) missing:\n",
    paste(missing_10x_files, collapse = "\n")
  )
}

message("All sample data directories and required 10X files were found.")


# ------------------------------------------------------------
# 7. Load sample-level 10X count matrices
# ------------------------------------------------------------

count_matrices <- lapply(
  sample_metadata$full_data_path,
  Read10X
)

names(count_matrices) <- sample_metadata$sample_id


# ------------------------------------------------------------
# 8. Validate loaded matrices
# ------------------------------------------------------------

for (sample_id in names(count_matrices)) {
  
  cat(
    "\n===== ",
    sample_id,
    " =====\n",
    sep = ""
  )
  
  x <- count_matrices[[sample_id]]
  
  # Read10X() may return a list when multiple feature types exist.
  # P02 was verified to return one sparse Gene Expression matrix per sample.
  if (is.list(x)) {
    stop(
      "Read10X returned multiple feature types for sample: ",
      sample_id,
      ". Inspect the feature types before continuing."
    )
  }
  
  print(class(x))
  print(dim(x))
}


# ------------------------------------------------------------
# 9. Confirm feature dimensions are consistent across samples
# ------------------------------------------------------------

n_features_per_sample <- vapply(
  count_matrices,
  nrow,
  integer(1)
)

if (length(unique(n_features_per_sample)) != 1) {
  stop(
    "Feature counts differ across samples: ",
    paste(
      names(n_features_per_sample),
      n_features_per_sample,
      sep = "=",
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 10. Print loading summary
# ------------------------------------------------------------

loading_summary <- data.frame(
  sample_id = names(count_matrices),
  n_features = vapply(count_matrices, nrow, integer(1)),
  n_barcodes = vapply(count_matrices, ncol, integer(1)),
  stringsAsFactors = FALSE
)

print(loading_summary)

cat(
  "\nTotal Cell Ranger-filtered barcodes:",
  sum(loading_summary$n_barcodes),
  "\n"
)


# ------------------------------------------------------------
# 11. Record R / package session information
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
    "01_load_data_sessionInfo.txt"
  )
)

message("01_load_data.R completed successfully.")
