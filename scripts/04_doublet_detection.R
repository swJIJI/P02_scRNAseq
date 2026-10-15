# ============================================================
# 04_doublet_detection.R
# Project: P02 (HN00262036) scRNA-seq
#
# Purpose / 목적
# - QC-filtered Seurat object에서 doublet 후보를 탐지한다.
# - scDblFinder score/class를 cell metadata에 연결한다.
# - Sample별 doublet 분포를 확인하여 제거 여부를 검토한다.
#
# Input / 입력
# - p02_qc: 38,606 features x 30,358 cells
# - Approved QC: nFeature_RNA >= 500, percent.mt <= 20
# - nCount_RNA: diagnostic only
#
# Planned workflow / 예정 처리
# 1. Packages, input object, RNA counts layers 확인
# 2. Raw counts 기반 SingleCellExperiment 준비
# 3. scDblFinder 실행 및 sample별 결과 검토
# 4. 승인된 doublet 제거 기준 적용
#
# Planned outputs / 예정 출력
# - Sample별 doublet summary table 및 diagnostic figures
# - Doublet score/class를 포함한 Seurat object
# - 승인된 singlet object
#
# Checkpoint policy / RDS 저장 원칙
# - 의미 있는 stage 완료 및 결과 확인 후 objects/에 저장한다.
# - 중간 점검마다 RDS를 생성하지 않는다.
# - 기존 checkpoint를 임의로 덮어쓰지 않는다.
# ============================================================


# ------------------------------------------------------------
# 1. Preflight checks / 실행 전 조건 확인
# ------------------------------------------------------------

# 필요한 package의 설치 여부와 version 확인
required_packages <- c(
  "Seurat",
  "SeuratObject",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "scDblFinder"
)

package_check <- data.frame(
  package = required_packages,
  installed = vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  ),
  version = vapply(
    required_packages,
    function(pkg) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(utils::packageVersion(pkg))
      } else {
        NA_character_
      }
    },
    character(1)
  )
)

print(package_check)

if (!all(package_check$installed)) {
  stop("필요한 package가 없습니다. package_check 결과를 확인하세요.")
}


# 입력 객체와 예상 크기 확인
stopifnot(exists("p02_qc"))
stopifnot(inherits(p02_qc, "Seurat"))

print(dim(p02_qc))

stopifnot(
  nrow(p02_qc) == 38606L,
  ncol(p02_qc) == 30358L
)


# RNA assay 및 counts layer 확인
stopifnot("RNA" %in% SeuratObject::Assays(p02_qc))

print(SeuratObject::DefaultAssay(p02_qc))
print(class(p02_qc[["RNA"]]))

rna_layers <- SeuratObject::Layers(p02_qc[["RNA"]])
print(rna_layers)

stopifnot(any(grepl("^counts($|\\.)", rna_layers)))


# Metadata 및 sample 구분 확인
qc_metadata <- p02_qc[[]]

stopifnot(
  all(c("sample_id", "nFeature_RNA", "percent.mt") %in%
        colnames(qc_metadata)),
  identical(rownames(qc_metadata), colnames(p02_qc))
)

sample_ids <- as.character(qc_metadata$sample_id)

stopifnot(
  !anyNA(sample_ids),
  all(nzchar(trimws(sample_ids))),
  setequal(
    unique(sample_ids),
    c("FR_D1", "FR_D7", "naive_D1", "naive_D7")
  )
)

print(table(sample_ids))


# 승인된 QC 조건 충족 및 barcode 중복 확인
# nCount_RNA는 filtering 조건에 포함하지 않는다.
input_check <- c(
  no_duplicated_barcodes =
    anyDuplicated(colnames(p02_qc)) == 0L,
  feature_threshold_pass =
    all(qc_metadata$nFeature_RNA >= 500),
  mitochondrial_threshold_pass =
    all(qc_metadata$percent.mt <= 20)
)

print(input_check)
stopifnot(all(input_check))

# ------------------------------------------------------------
# 2. Prepare doublet-analysis object / 분석용 객체 준비
# ------------------------------------------------------------

# 원본 QC object를 보존하고 분석용 복사본 생성
p02_doublet <- p02_qc

# 분리된 raw counts layers를 하나의 "counts" layer로 합치기
p02_doublet <- SeuratObject::JoinLayers(
  object = p02_doublet,
  assay = "RNA",
  layers = "counts",
  new = "counts"
)

# Layer 구조와 객체 크기 확인
print(SeuratObject::Layers(p02_doublet[["RNA"]]))
print(dim(p02_doublet))
print(table(p02_doublet$sample_id))

# Feature 및 cell 순서, metadata 보존 확인
join_check <- c(
  same_features =
    identical(rownames(p02_doublet), rownames(p02_qc)),
  same_cells =
    identical(colnames(p02_doublet), colnames(p02_qc)),
  same_metadata =
    identical(p02_doublet[[]], p02_qc[[]])
)

print(join_check)
stopifnot(all(join_check))


# ------------------------------------------------------------
# 3. Create SingleCellExperiment / scDblFinder 입력 객체 생성
# ------------------------------------------------------------

# Raw counts 추출 — sparse matrix 구조 유지
counts_doublet <- SeuratObject::LayerData(
  object = p02_doublet,
  assay = "RNA",
  layer = "counts"
)

print(class(counts_doublet))
print(dim(counts_doublet))

# Counts와 Seurat object의 feature/cell 순서 확인
stopifnot(
  identical(rownames(counts_doublet), rownames(p02_doublet)),
  identical(colnames(counts_doublet), colnames(p02_doublet))
)

# Cell metadata 추출 및 순서 확인
doublet_metadata <- p02_doublet[[]]

stopifnot(
  identical(rownames(doublet_metadata), colnames(counts_doublet))
)

# Counts와 metadata를 함께 담은 SCE object 생성
sce_doublet <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts_doublet),
  colData = S4Vectors::DataFrame(doublet_metadata)
)

# 입력 객체 구조 확인
print(dim(sce_doublet))
print(SummarizedExperiment::assayNames(sce_doublet))
print(table(sce_doublet$sample_id))

sce_check <- c(
  same_features =
    identical(rownames(sce_doublet), rownames(p02_doublet)),
  same_cells =
    identical(colnames(sce_doublet), colnames(p02_doublet)),
  same_samples =
    identical(
      as.character(sce_doublet$sample_id),
      as.character(p02_doublet$sample_id)
    )
)

print(sce_check)
stopifnot(all(sce_check))

# ------------------------------------------------------------
# 4. Check execution settings / 실행 조건 점검
# ------------------------------------------------------------

# 설치된 버전의 주요 인자와 기본값 확인
scdbl_formals <- formals(scDblFinder::scDblFinder)

print(
  scdbl_formals[
    intersect(
      c("samples", "clusters", "dbr", "dbr.sd",
        "dbr.per1k", "BPPARAM"),
      names(scdbl_formals)
    )
  ]
)

# 이전 scDblFinder 결과가 이미 있는지 확인
# 새로 만든 입력 객체라면 두 항목 모두 NULL이 예상됨
print(sce_doublet$scDblFinder.score)
print(sce_doublet$scDblFinder.class)

# 순차 실행에 사용할 package 확인
print(requireNamespace("BiocParallel", quietly = TRUE))


# ------------------------------------------------------------
# 5. Run scDblFinder / Sample별 doublet 탐지
# ------------------------------------------------------------
# Approved settings / 승인된 조건
# - samples: sample_id (sample별 OCM 포획 단위)
# - clusters: FALSE (cluster에 의존하지 않는 무작위 조합)
# - dbr: NULL (설치된 버전의 자동 설정 사용)
# - seed: 42
# - execution: serial / 순차 실행

set.seed(42)

sce_doublet <- scDblFinder::scDblFinder(
  sce = sce_doublet,
  samples = "sample_id",
  clusters = FALSE,
  dbr = NULL,
  BPPARAM = BiocParallel::SerialParam(),
  verbose = TRUE
)


# ------------------------------------------------------------
# 6. Validate results / 실행 결과 확인
# ------------------------------------------------------------

# 모든 입력 세포가 유지되고 결과가 생성됐는지 확인
stopifnot(
  identical(colnames(sce_doublet), colnames(p02_doublet)),
  !is.null(sce_doublet$scDblFinder.score),
  !is.null(sce_doublet$scDblFinder.class),
  !anyNA(sce_doublet$scDblFinder.score),
  !anyNA(sce_doublet$scDblFinder.class)
)

# Sample별 singlet / doublet 개수
doublet_class_table <- table(
  sample_id = sce_doublet$sample_id,
  doublet_class = sce_doublet$scDblFinder.class
)

print(doublet_class_table)

# Sample별 분류 비율 (%)
print(round(
  100 * prop.table(doublet_class_table, margin = 1),
  digits = 2
))

# 전체 doublet score 분포
print(summary(sce_doublet$scDblFinder.score))

# ------------------------------------------------------------
# 7. Attach doublet results / Seurat metadata에 결과 연결
# ------------------------------------------------------------

# Cell 이름과 순서가 일치하는지 먼저 확인
stopifnot(
  identical(colnames(sce_doublet), colnames(p02_doublet))
)

# 세포별 결과를 barcode가 포함된 data.frame으로 준비
doublet_results <- data.frame(
  scDblFinder.score = sce_doublet$scDblFinder.score,
  scDblFinder.class = as.character(sce_doublet$scDblFinder.class),
  row.names = colnames(sce_doublet)
)

p02_doublet <- SeuratObject::AddMetaData(
  object = p02_doublet,
  metadata = doublet_results
)

# 결과 연결 검증 / Validate barcode order and metadata values
stopifnot(
  ncol(p02_doublet) == 30358L,
  
  identical(
    rownames(p02_doublet[[]]),
    rownames(doublet_results)
  ),
  
  identical(
    unname(p02_doublet$scDblFinder.score),
    unname(doublet_results$scDblFinder.score)
  ),
  
  identical(
    unname(as.character(p02_doublet$scDblFinder.class)),
    unname(as.character(doublet_results$scDblFinder.class))
  )
)

# ------------------------------------------------------------
# 8. Compare QC metrics / 분류별 QC 특성 비교
# ------------------------------------------------------------

doublet_qc_summary <- p02_doublet[[]] |>
  dplyr::group_by(sample_id, scDblFinder.class) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    median_nCount_RNA = median(nCount_RNA),
    median_nFeature_RNA = median(nFeature_RNA),
    median_percent_mt = median(percent.mt),
    median_doublet_score = median(scDblFinder.score),
    .groups = "drop"
  )

print(doublet_qc_summary, n = Inf, width = Inf)

# ------------------------------------------------------------
# 9. Visualize QC distributions / 분류별 QC 분포 확인
# ------------------------------------------------------------

doublet_plot_df <- p02_doublet[[]]

doublet_plot_df$scDblFinder.class <- factor(
  doublet_plot_df$scDblFinder.class,
  levels = c("singlet", "doublet")
)

# Raw UMI counts 분포
doublet_count_plot <- ggplot2::ggplot(
  doublet_plot_df,
  ggplot2::aes(
    x = scDblFinder.class,
    y = nCount_RNA,
    fill = scDblFinder.class
  )
) +
  ggplot2::geom_violin(trim = FALSE) +
  ggplot2::geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white"
  ) +
  ggplot2::facet_wrap(~ sample_id, ncol = 2) +
  ggplot2::scale_y_log10() +
  ggplot2::scale_fill_manual(
    values = c(singlet = "#4C78A8", doublet = "#E45756")
  ) +
  ggplot2::labs(
    x = NULL,
    y = "UMI counts (log10 scale)",
    title = "UMI distribution by doublet classification"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none")

# 검출 유전자 수 분포 — 같은 형식으로 비교
doublet_feature_plot <- doublet_count_plot +
  ggplot2::aes(y = nFeature_RNA) +
  ggplot2::labs(
    y = "Detected genes (log10 scale)",
    title = "Gene distribution by doublet classification"
  )

print(doublet_count_plot)
print(doublet_feature_plot)

# ------------------------------------------------------------
# 10. Save detection checkpoint / 탐지 완료 객체 저장
# ------------------------------------------------------------

# 저장 전 최종 검증 — doublet 후보를 포함한 전체 객체
stopifnot(
  nrow(p02_doublet) == 38606L,
  ncol(p02_doublet) == 30358L,
  !anyNA(p02_doublet$scDblFinder.score),
  !anyNA(p02_doublet$scDblFinder.class),
  identical(colnames(p02_doublet), colnames(p02_qc))
)

# 기존 checkpoint를 덮어쓰지 않도록 timestamp 사용
dir.create("objects", recursive = TRUE, showWarnings = FALSE)

doublet_rds_path <- file.path(
  "objects",
  paste0(
    "P02_04_doublet_annotated_",
    format(Sys.time(), "%Y%m%dT%H%M%S"),
    ".rds"
  )
)

stopifnot(!file.exists(doublet_rds_path))

saveRDS(
  object = p02_doublet,
  file = doublet_rds_path,
  compress = TRUE
)

# 파일 생성 및 크기 확인
stopifnot(
  file.exists(doublet_rds_path),
  file.info(doublet_rds_path)$size > 0
)

print(normalizePath(doublet_rds_path))
print(file.info(doublet_rds_path)[, c("size", "mtime")])

# ------------------------------------------------------------
# 11. Create singlet object / 승인된 doublet 후보 제외
# ------------------------------------------------------------

# 분류값에 결측치나 예상하지 못한 값이 없는지 확인
doublet_class <- as.character(p02_doublet$scDblFinder.class)

stopifnot(
  !anyNA(doublet_class),
  all(doublet_class %in% c("singlet", "doublet"))
)

# 유지할 singlet barcode 선택
singlet_cells <- colnames(p02_doublet)[
  doublet_class == "singlet"
]

# 별도 객체 생성 — 제거 전 p02_doublet는 보존
p02_singlet <- subset(
  x = p02_doublet,
  cells = singlet_cells
)

# 객체 크기, barcode 순서 및 분류 검증
stopifnot(
  nrow(p02_singlet) == nrow(p02_doublet),
  ncol(p02_singlet) == length(singlet_cells),
  ncol(p02_singlet) == 27726L,
  identical(colnames(p02_singlet), singlet_cells),
  all(p02_singlet$scDblFinder.class == "singlet")
)

print(dim(p02_singlet))
print(table(p02_singlet$sample_id))


# Sample별 제거 전후 세포 수 / Before–after summary
sample_levels <- sort(unique(as.character(p02_doublet$sample_id)))

singlet_filter_summary <- data.frame(
  sample_id = sample_levels,
  n_before = as.integer(table(factor(
    p02_doublet$sample_id,
    levels = sample_levels
  ))),
  n_singlet = as.integer(table(factor(
    p02_singlet$sample_id,
    levels = sample_levels
  )))
)

singlet_filter_summary$n_removed <- with(
  singlet_filter_summary,
  n_before - n_singlet
)

singlet_filter_summary$removed_pct <- with(
  singlet_filter_summary,
  round(100 * n_removed / n_before, 2)
)

print(singlet_filter_summary)

# ------------------------------------------------------------
# 12. Save final outputs / Doublet 단계 최종 결과 저장
# ------------------------------------------------------------

# 최종 singlet 객체 검증
stopifnot(
  nrow(p02_singlet) == 38606L,
  ncol(p02_singlet) == 27726L,
  all(p02_singlet$scDblFinder.class == "singlet"),
  sum(singlet_filter_summary$n_removed) == 2632L
)

# 같은 실행에서 생성된 파일에 공통 timestamp 사용
stage_stamp <- format(Sys.time(), "%Y%m%dT%H%M%S")

stage_dirs <- c(
  "objects", "tables", "figures",
  file.path("report", "session_info")
)

for (stage_dir in stage_dirs) {
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
}

# Output paths / 저장 경로
stage_files <- c(
  singlet_rds = file.path(
    "objects", paste0("P02_04_singlet_", stage_stamp, ".rds")
  ),
  filter_summary = file.path(
    "tables", paste0("P02_04_singlet_filter_summary_", stage_stamp, ".csv")
  ),
  qc_summary = file.path(
    "tables", paste0("P02_04_doublet_qc_summary_", stage_stamp, ".csv")
  ),
  cell_results = file.path(
    "tables", paste0("P02_04_doublet_cell_results_", stage_stamp, ".csv")
  ),
  umi_plot = file.path(
    "figures", paste0("P02_04_doublet_UMI_", stage_stamp, ".png")
  ),
  gene_plot = file.path(
    "figures", paste0("P02_04_doublet_genes_", stage_stamp, ".png")
  ),
  session_info = file.path(
    "report", "session_info",
    paste0("P02_04_sessionInfo_", stage_stamp, ".txt")
  )
)

# 기존 파일 덮어쓰기 방지
stopifnot(!any(file.exists(stage_files)))

# Final checkpoint / 다음 분석 단계의 입력 객체
saveRDS(
  object = p02_singlet,
  file = stage_files[["singlet_rds"]],
  compress = TRUE
)

# Summary tables
utils::write.csv(
  singlet_filter_summary,
  stage_files[["filter_summary"]],
  row.names = FALSE
)

utils::write.csv(
  doublet_qc_summary,
  stage_files[["qc_summary"]],
  row.names = FALSE
)

# 제거 후보까지 포함한 전체 세포별 예측 결과 보존
cell_results_export <- data.frame(
  cell_id = colnames(p02_doublet),
  sample_id = p02_doublet$sample_id,
  scDblFinder.score = p02_doublet$scDblFinder.score,
  scDblFinder.class = as.character(p02_doublet$scDblFinder.class)
)

utils::write.csv(
  cell_results_export,
  stage_files[["cell_results"]],
  row.names = FALSE
)

# Diagnostic figures
ggplot2::ggsave(
  filename = stage_files[["umi_plot"]],
  plot = doublet_count_plot,
  width = 9, height = 9, units = "in",
  dpi = 300, bg = "white"
)

ggplot2::ggsave(
  filename = stage_files[["gene_plot"]],
  plot = doublet_feature_plot,
  width = 9, height = 9, units = "in",
  dpi = 300, bg = "white"
)

# 실행 환경 기록
writeLines(
  capture.output(sessionInfo()),
  con = stage_files[["session_info"]]
)

# 모든 파일 생성 및 크기 확인
stopifnot(
  all(file.exists(stage_files)),
  all(file.info(stage_files)$size > 0)
)

print(data.frame(
  output = names(stage_files),
  path = unname(stage_files),
  size_MB = round(file.info(stage_files)$size / 1024^2, 2),
  row.names = NULL
))
