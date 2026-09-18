# ============================================================
# 02_create_object.R
#
# Purpose / 목적:
#   - Create one Seurat object per sample from validated 10X matrices.
#   - Add sample-level metadata to every cell.
#   - Add sample ID prefixes to cell barcodes to prevent name collisions.
#   - Merge the four sample-level objects into the initial P02 Seurat object.
#   - Validate dimensions, sample counts, and cell-name uniqueness.
#
# Input:
#   - scripts/00_config.R
#   - objects created in scripts/01_load_data.R:
#       sample_metadata
#       count_matrices
#
# Output:
#   - seurat_list : named list of sample-level Seurat objects
#   - p02         : merged P02 Seurat object
#   - report/session_info/02_create_object_sessionInfo.txt
#
# Important:
#   No dataset-specific QC threshold is applied in this script.
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
# 3. Ensure validated input objects are available
# ------------------------------------------------------------

# 02_create_object.R can be run after 01_load_data.R in the same session.
# If the required objects are absent, run 01_load_data.R automatically.
if (
  !exists("sample_metadata") ||
  !exists("count_matrices")
) {
  message(
    "Validated input objects not found in the current R session. ",
    "Running scripts/01_load_data.R first."
  )
  
  source("scripts/01_load_data.R")
}


# ------------------------------------------------------------
# 4. Create sample-level Seurat objects
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
    
    # Different samples can contain identical raw 10X barcode strings.
    # Prefix each cell with sample_id before merging.
    obj <- RenameCells(
      obj,
      add.cell.id = sample_id
    )
    
    # Add dataset-specific sample metadata to every cell.
    obj$sample_id <- sample_id
    obj$pdo_source <- sample_metadata$pdo_source[i]
    obj$timepoint <- sample_metadata$timepoint[i]
    obj$biological_replicate <- sample_metadata$biological_replicate[i]
    obj$culture_batch <- sample_metadata$culture_batch[i]
    
    return(obj)
  }
)

names(seurat_list) <- sample_metadata$sample_id


# ------------------------------------------------------------
# 5. Validate sample-level Seurat objects
# ------------------------------------------------------------

for (sample_id in names(seurat_list)) {
  
  obj <- seurat_list[[sample_id]]
  x <- count_matrices[[sample_id]]
  
  cat(
    "\n===== ",
    sample_id,
    " =====\n",
    sep = ""
  )
  
  cat("Features:", nrow(obj), "\n")
  cat("Cells:", ncol(obj), "\n")
  
  if (nrow(obj) != nrow(x)) {
    stop(
      "Feature count changed during Seurat object creation for sample: ",
      sample_id
    )
  }
  
  if (ncol(obj) != ncol(x)) {
    stop(
      "Cell count changed during Seurat object creation for sample: ",
      sample_id
    )
  }
}


# ------------------------------------------------------------
# 6. Merge sample-level objects
# ------------------------------------------------------------

p02 <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  project = PROJECT_ID
)


# ------------------------------------------------------------
# 7. Validate merged P02 object
# ------------------------------------------------------------

expected_features <- unique(
  vapply(
    count_matrices,
    nrow,
    integer(1)
  )
)

expected_cells <- sum(
  vapply(
    count_matrices,
    ncol,
    integer(1)
  )
)

if (length(expected_features) != 1) {
  stop(
    "Cannot validate merged feature count because input feature counts differ."
  )
}

if (nrow(p02) != expected_features) {
  stop(
    "Merged feature count does not match the validated input matrices."
  )
}

if (ncol(p02) != expected_cells) {
  stop(
    "Merged cell count does not match the sum of validated input matrices."
  )
}

if (anyDuplicated(colnames(p02)) != 0) {
  stop("Duplicated cell names were found after merging.")
}


# ------------------------------------------------------------
# 8. Print final object summary
# ------------------------------------------------------------

cat(
  "\n===== P02 merged object =====\n"
)

cat(
  "Number of features:",
  nrow(p02),
  "\n"
)

cat(
  "Number of cells:",
  ncol(p02),
  "\n\n"
)

print(
  table(p02$sample_id)
)

cat(
  "\nDuplicated cell names:",
  anyDuplicated(colnames(p02)),
  "\n"
)


# ------------------------------------------------------------
# 9. Record R / package session information
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
    "02_create_object_sessionInfo.txt"
  )
)

message("02_create_object.R completed successfully.")
