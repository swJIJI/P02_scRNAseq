# ============================================================
# 01_load_data.R
#
# Purpose:
#   P02 sample metadata를 불러오고,
#   각 sample의 10X input directory가 실제로 존재하는지 확인한다.
#
# Input:
#   - scripts/00_config.R
#   - metadata/sample_metadata.csv
#
# This script does NOT create Seurat objects yet.
# ============================================================


# ------------------------------------------------------------
# 1. Load dataset-specific configuration
# ------------------------------------------------------------

source("scripts/00_config.R")


# ------------------------------------------------------------
# 2. Set random seed
# ------------------------------------------------------------

set.seed(RANDOM_SEED)


# ------------------------------------------------------------
# 3. Load sample metadata
# ------------------------------------------------------------

sample_metadata <- read.csv(
  METADATA_FILE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Construct full 10X data paths
# ------------------------------------------------------------

sample_metadata$full_data_path <- file.path(
  RAW_DATA_DIR,
  sample_metadata$data_path
)


# ------------------------------------------------------------
# 5. Check whether directories exist
# ------------------------------------------------------------

sample_metadata$data_path_exists <- dir.exists(
  sample_metadata$full_data_path
)


# Print path validation result
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


# ------------------------------------------------------------
# 6. Stop if any sample directory is missing
# ------------------------------------------------------------

if (!all(sample_metadata$data_path_exists)) {
  
  missing_samples <- sample_metadata$sample_id[
    !sample_metadata$data_path_exists
  ]
  
  stop(
    "10X data directory not found for sample(s): ",
    paste(missing_samples, collapse = ", ")
  )
}


message(
  "All sample data directories were found successfully."
)

# ------------------------------------------------------------
# 7. Load 10X count matrices
# ------------------------------------------------------------

library(Seurat)

count_matrices <- lapply(
  sample_metadata$full_data_path,
  Read10X
)

names(count_matrices) <- sample_metadata$sample_id


# ------------------------------------------------------------
# 8. Inspect loaded matrices
# ------------------------------------------------------------

for (sample_id in names(count_matrices)) {
  
  cat("\n===== ", sample_id, " =====\n", sep = "")
  
  print(class(count_matrices[[sample_id]]))
  print(dim(count_matrices[[sample_id]]))
}

# ------------------------------------------------------------
# 9. Create Seurat objects
# ------------------------------------------------------------

seurat_list <- lapply(
  seq_len(nrow(sample_metadata)),
  function(i) {
    
    sample_id <- sample_metadata$sample_id[i]
    
    obj <- CreateSeuratObject(
      counts = count_matrices[[sample_id]],
      project = PROJECT_ID,
      assay = "RNA",
      
      # QC filtering은 아직 수행하지 않는다.
      # Keep all Cell Ranger-filtered barcodes at this stage.
      min.cells = 0,
      min.features = 0
    )
    
    # sample 간 동일한 10X barcode가 존재할 수 있으므로
    # sample ID를 cell barcode 앞에 추가한다.
    obj <- RenameCells(
      obj,
      add.cell.id = sample_id
    )
    
    # Dataset-specific sample metadata 추가
    obj$sample_id <- sample_id
    obj$pdo_source <- sample_metadata$pdo_source[i]
    obj$timepoint <- sample_metadata$timepoint[i]
    obj$biological_replicate <- sample_metadata$biological_replicate[i]
    obj$culture_batch <- sample_metadata$culture_batch[i]
    
    return(obj)
  }nrow(p02)
ncol(p02)
table(p02$sample_id)
anyDuplicated(colnames(p02))
)

names(seurat_list) <- sample_metadata$sample_id

