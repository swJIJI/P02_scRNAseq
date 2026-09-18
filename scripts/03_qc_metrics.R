# ==============================================================================
# 03_qc_metrics.R
# P02 (HN00262036) - QC metric calculation and filtering
# QC 지표 계산, 분포 확인 및 승인된 QC criteria 적용
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Load packages / 패키지 불러오기
# ------------------------------------------------------------------------------
source("00_config.R")

library(Seurat)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Object 확인 / Inspect Seurat object
# ------------------------------------------------------------------------------

# 현재 Seurat object의 크기 확인
dim(p02)

# Default assay 확인
DefaultAssay(p02)

# 현재 metadata column 확인
colnames(p02@meta.data)


# ------------------------------------------------------------------------------
# 2. Feature naming 확인 / Inspect feature names
# ------------------------------------------------------------------------------

# 앞부분 gene names 확인
head(rownames(p02), 20)

# Human mitochondrial gene prefix 확인
sum(grepl("^MT-", rownames(p02)))

# 혹시 lowercase naming인지 함께 확인
sum(grepl("^mt-", rownames(p02)))

# ------------------------------------------------------------------------------
# 3. Mitochondrial percentage 계산
#    Calculate mitochondrial RNA percentage
# ------------------------------------------------------------------------------

p02[["percent.mt"]] <- PercentageFeatureSet(
  p02,
  pattern = "^MT-"
)

# 계산 결과 확인
head(
  p02@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt")]
)

summary(
  p02@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt")]
)

# ------------------------------------------------------------------------------
# 4. Sample별 QC summary
#    Summarize QC metrics by sample
# ------------------------------------------------------------------------------

qc_summary_by_sample <- p02@meta.data |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    
    median_nCount_RNA = median(nCount_RNA),
    mean_nCount_RNA   = mean(nCount_RNA),
    
    median_nFeature_RNA = median(nFeature_RNA),
    mean_nFeature_RNA   = mean(nFeature_RNA),
    
    median_percent_mt = median(percent.mt),
    mean_percent_mt   = mean(percent.mt),
    
    .groups = "drop"
  )

qc_summary_by_sample

qc_summary_by_sample |>
  dplyr::select(
    sample_id,
    n_cells,
    median_percent_mt,
    mean_percent_mt
  )

# ------------------------------------------------------------------------------
# 5. QC metric distributions by sample
#    Sample별 주요 QC metric 분포 확인
# ------------------------------------------------------------------------------

qc_violin <- VlnPlot(
  object = p02,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "sample_id",
  pt.size = 0,
  ncol = 3
)

qc_violin

# ------------------------------------------------------------------------------
# 6. Prepare QC plotting data
#    QC plot용 metadata 추출 / Extract metadata for QC plots
# ------------------------------------------------------------------------------

qc_plot_df <- FetchData(
  object = p02,
  vars = c(
    "sample_id",
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt"
  )
)

head(qc_plot_df)

# ------------------------------------------------------------------------------
# 7. Scatter plot: nCount_RNA vs nFeature_RNA
#    count-feature 관계 확인 / Inspect count-feature relationship
# ------------------------------------------------------------------------------

qc_scatter_count_feature <- ggplot(
  qc_plot_df,
  aes(x = nCount_RNA, y = nFeature_RNA)
) +
  geom_point(alpha = 0.25, size = 0.3) +
  facet_wrap(~ sample_id, ncol = 2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "QC scatter plot: nCount_RNA vs nFeature_RNA",
    x = "nCount_RNA (log10 scale)",
    y = "nFeature_RNA (log10 scale)"
  ) +
  theme_bw()

qc_scatter_count_feature

# ------------------------------------------------------------------------------
# 8. Scatter plot: nCount_RNA vs percent.mt
#    count-mito 관계 확인 / Inspect count-mito relationship
# ------------------------------------------------------------------------------

qc_scatter_count_mt <- ggplot(
  qc_plot_df,
  aes(x = nCount_RNA, y = percent.mt)
) +
  geom_point(alpha = 0.25, size = 0.3) +
  facet_wrap(~ sample_id, ncol = 2) +
  scale_x_log10() +
  labs(
    title = "QC scatter plot: nCount_RNA vs percent.mt",
    x = "nCount_RNA (log10 scale)",
    y = "percent.mt"
  ) +
  theme_bw()

qc_scatter_count_mt

# ------------------------------------------------------------------------------
# 9. Evaluate candidate mitochondrial thresholds
#    mitochondrial QC cutoff 후보별 cell 수 확인
# ------------------------------------------------------------------------------

qc_mt_threshold_summary <- p02@meta.data |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    
    mt_gt_10 = sum(percent.mt > 10),
    mt_gt_15 = sum(percent.mt > 15),
    mt_gt_20 = sum(percent.mt > 20),
    
    pct_mt_gt_10 = 100 * mean(percent.mt > 10),
    pct_mt_gt_15 = 100 * mean(percent.mt > 15),
    pct_mt_gt_20 = 100 * mean(percent.mt > 20),
    
    .groups = "drop"
  )

qc_mt_threshold_summary

# ------------------------------------------------------------------------------
# 10. Evaluate candidate lower nFeature thresholds
#     low-complexity cell 후보 규모 확인
# ------------------------------------------------------------------------------

qc_feature_threshold_summary <- p02@meta.data |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    
    feature_lt_500  = sum(nFeature_RNA < 500),
    feature_lt_750  = sum(nFeature_RNA < 750),
    feature_lt_1000 = sum(nFeature_RNA < 1000),
    
    pct_feature_lt_500 =
      100 * mean(nFeature_RNA < 500),
    
    pct_feature_lt_750 =
      100 * mean(nFeature_RNA < 750),
    
    pct_feature_lt_1000 =
      100 * mean(nFeature_RNA < 1000),
    
    .groups = "drop"
  )

qc_feature_threshold_summary

# ------------------------------------------------------------------------------
# 11. Explore overlap between QC criteria
#     QC 조건 간 overlap 확인
#
# NOTE:
# nFeature_RNA < 750 and percent.mt > 15 are exploratory candidates only.
# 아직 최종 QC cutoff가 아님.
# ------------------------------------------------------------------------------

qc_overlap_750_15 <- p02@meta.data |>
  dplyr::mutate(
    low_feature = nFeature_RNA < 750,
    high_mt = percent.mt > 15
  ) |>
  dplyr::mutate(
    qc_category = dplyr::case_when(
      low_feature & high_mt  ~ "both",
      low_feature & !high_mt ~ "low_feature_only",
      !low_feature & high_mt ~ "high_mt_only",
      TRUE                   ~ "pass_both"
    )
  ) |>
  dplyr::count(
    sample_id,
    qc_category,
    name = "n_cells"
  )

qc_overlap_750_15

# ------------------------------------------------------------------------------
# 12. Characterize low-feature cells
#     low-feature population의 특성 확인
#
# NOTE:
# nFeature_RNA < 750 is still an exploratory definition.
# 아직 filtering cutoff로 확정한 값이 아님.
# ------------------------------------------------------------------------------

qc_low_feature_profile <- p02@meta.data |>
  dplyr::mutate(
    feature_group = dplyr::if_else(
      nFeature_RNA < 750,
      "feature_lt_750",
      "feature_ge_750"
    )
  ) |>
  dplyr::group_by(
    sample_id,
    feature_group
  ) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    
    median_nFeature_RNA = median(nFeature_RNA),
    median_nCount_RNA = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    
    q25_nCount_RNA = quantile(nCount_RNA, 0.25),
    q75_nCount_RNA = quantile(nCount_RNA, 0.75),
    
    q25_percent_mt = quantile(percent.mt, 0.25),
    q75_percent_mt = quantile(percent.mt, 0.75),
    
    .groups = "drop"
  )

qc_low_feature_profile

# ------------------------------------------------------------------------------
# 13. Detailed nFeature_RNA distribution
#     low-feature 영역의 세부 분포 확인
# ------------------------------------------------------------------------------

qc_feature_hist <- ggplot(
  qc_plot_df,
  aes(x = nFeature_RNA)
) +
  geom_histogram(
    binwidth = 100,
    boundary = 0
  ) +
  facet_wrap(
    ~ sample_id,
    ncol = 2,
    scales = "free_y"
  ) +
  coord_cartesian(
    xlim = c(0, 2500)
  ) +
  geom_vline(
    xintercept = 750,
    linetype = "dashed"
  ) +
  labs(
    title = "Low nFeature_RNA distribution by sample",
    x = "nFeature_RNA",
    y = "Number of cells"
  ) +
  theme_bw()

qc_feature_hist


# ------------------------------------------------------------------------------
# 14. Detailed low nCount_RNA distribution
#     low-count 영역의 세부 분포 확인
# ------------------------------------------------------------------------------

qc_count_hist <- ggplot(
  qc_plot_df,
  aes(x = nCount_RNA)
) +
  geom_histogram(
    binwidth = 100,
    boundary = 0
  ) +
  facet_wrap(
    ~ sample_id,
    ncol = 2,
    scales = "free_y"
  ) +
  coord_cartesian(
    xlim = c(500, 5000)
  ) +
  labs(
    title = "Low nCount_RNA distribution by sample",
    x = "nCount_RNA",
    y = "Number of cells"
  ) +
  theme_bw()

qc_count_hist


# ------------------------------------------------------------------------------
# 15. Explore candidate lower nCount thresholds
#     low-count cutoff 후보별 cell 수 확인
#
# NOTE:
# These are exploratory reference values only.
# 아직 최종 QC cutoff가 아님.
# ------------------------------------------------------------------------------

qc_count_threshold_summary <- p02@meta.data |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    
    count_lt_750  = sum(nCount_RNA < 750),
    count_lt_1000 = sum(nCount_RNA < 1000),
    count_lt_1500 = sum(nCount_RNA < 1500),
    count_lt_2000 = sum(nCount_RNA < 2000),
    
    pct_count_lt_750 =
      100 * mean(nCount_RNA < 750),
    
    pct_count_lt_1000 =
      100 * mean(nCount_RNA < 1000),
    
    pct_count_lt_1500 =
      100 * mean(nCount_RNA < 1500),
    
    pct_count_lt_2000 =
      100 * mean(nCount_RNA < 2000),
    
    .groups = "drop"
  )

qc_count_threshold_summary

# ------------------------------------------------------------------------------
# 16. Overlap between low-count and low-feature populations
#     low-count와 low-feature population의 overlap 확인
#
# NOTE:
# nCount_RNA < 1000 and nFeature_RNA < 750 are exploratory boundaries only.
# 아직 최종 QC filtering 기준이 아님.
# ------------------------------------------------------------------------------

qc_lowend_overlap <- p02@meta.data |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    
    # 각 조건에 해당하는 cell 수
    count_lt_1000 =
      sum(nCount_RNA < 1000),
    
    feature_lt_750 =
      sum(nFeature_RNA < 750),
    
    # 두 조건을 동시에 만족
    both_low =
      sum(
        nCount_RNA < 1000 &
          nFeature_RNA < 750
      ),
    
    # count만 낮음
    low_count_only =
      sum(
        nCount_RNA < 1000 &
          nFeature_RNA >= 750
      ),
    
    # feature만 낮음
    low_feature_only =
      sum(
        nCount_RNA >= 1000 &
          nFeature_RNA < 750
      ),
    
    # low-count 중 high-MT도 동시에 해당하는 cell
    low_count_high_mt =
      sum(
        nCount_RNA < 1000 &
          percent.mt > 15
      ),
    
    .groups = "drop"
  )

qc_lowend_overlap

# ------------------------------------------------------------------------------
# 17. Characterize exploratory low-count population
#     nCount_RNA < 1000 population의 QC 특성 확인
# ------------------------------------------------------------------------------

qc_low_count_profile <- p02@meta.data |>
  dplyr::mutate(
    count_group = dplyr::if_else(
      nCount_RNA < 1000,
      "count_lt_1000",
      "count_ge_1000"
    )
  ) |>
  dplyr::group_by(
    sample_id,
    count_group
  ) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    
    median_nCount_RNA =
      median(nCount_RNA),
    
    median_nFeature_RNA =
      median(nFeature_RNA),
    
    median_percent_mt =
      median(percent.mt),
    
    q25_nFeature_RNA =
      quantile(nFeature_RNA, 0.25),
    
    q75_nFeature_RNA =
      quantile(nFeature_RNA, 0.75),
    
    q25_percent_mt =
      quantile(percent.mt, 0.25),
    
    q75_percent_mt =
      quantile(percent.mt, 0.75),
    
    .groups = "drop"
  )

qc_low_count_profile


qc_strategy_summary <- p02@meta.data %>%
       mutate(
             # Candidate A: conservative QC
               candidate_A = nFeature_RNA >= 500 & percent.mt <= 20,
             
               # Candidate B: stricter mitochondrial cutoff
               candidate_B = nFeature_RNA >= 500 & percent.mt <= 15,
             
               # Candidate C: stricter low-feature cutoff
               candidate_C = nFeature_RNA >= 750 & percent.mt <= 20,
             
               # Candidate D: nCount-based rule for comparison only
               candidate_D = nCount_RNA >= 1000 & percent.mt <= 20
         ) %>%
       group_by(sample_id) %>%
       summarise(
             n_before = n(),
             
               A_keep = sum(candidate_A),
             A_keep_pct = round(mean(candidate_A) * 100, 2),
           
                B_keep = sum(candidate_B),
             B_keep_pct = round(mean(candidate_B) * 100, 2),
             
               C_keep = sum(candidate_C),
             C_keep_pct = round(mean(candidate_C) * 100, 2),
             
               D_keep = sum(candidate_D),
             D_keep_pct = round(mean(candidate_D) * 100, 2),
             
               .groups = "drop"
         )

qc_strategy_summary

# ------------------------------------------------------------------------------
# 18. Compare candidate QC strategies
#     Candidate QC criteria별 sample-wise cell retention 비교
# ------------------------------------------------------------------------------

qc_strategy_summary <- p02@meta.data %>%
  mutate(
    # Candidate A: conservative QC
    candidate_A = nFeature_RNA >= 500 & percent.mt <= 20,
    
    # Candidate B: stricter mitochondrial cutoff
    candidate_B = nFeature_RNA >= 500 & percent.mt <= 15,
    
    # Candidate C: stricter low-feature cutoff
    candidate_C = nFeature_RNA >= 750 & percent.mt <= 20,
    
    # Candidate D: nCount-based rule for comparison only
    candidate_D = nCount_RNA >= 1000 & percent.mt <= 20
  ) %>%
  group_by(sample_id) %>%
  summarise(
    n_before = n(),
    
    A_keep = sum(candidate_A),
    A_keep_pct = round(mean(candidate_A) * 100, 2),
    
    B_keep = sum(candidate_B),
    B_keep_pct = round(mean(candidate_B) * 100, 2),
    
    C_keep = sum(candidate_C),
    C_keep_pct = round(mean(candidate_C) * 100, 2),
    
    D_keep = sum(candidate_D),
    D_keep_pct = round(mean(candidate_D) * 100, 2),
    
    .groups = "drop"
  )

qc_strategy_summary


# ------------------------------------------------------------------------------
# 19. Characterize removal reasons for Candidate A
#     Conservative candidate에서 cell removal reason 확인
# ------------------------------------------------------------------------------

qc_A_reason <- p02@meta.data %>%
  mutate(
    low_feature = nFeature_RNA < 500,
    high_mt = percent.mt > 20,
    
    removal_reason = case_when(
      low_feature & high_mt ~ "both",
      low_feature           ~ "low_feature_only",
      high_mt               ~ "high_mt_only",
      TRUE                  ~ "pass"
    )
  ) %>%
  count(sample_id, removal_reason) %>%
  group_by(sample_id) %>%
  mutate(
    pct = round(n / sum(n) * 100, 2)
  ) %>%
  ungroup()

qc_A_reason


# ------------------------------------------------------------------------------
# 20. Apply approved hard QC filters
#     승인된 P02-specific QC criteria 적용
#
# Hard filters:
#   nFeature_RNA >= QC_MIN_FEATURES
#   percent.mt <= QC_MAX_MT
#
# nCount_RNA is retained as a diagnostic metric only.
# Upper nFeature/nCount cutoffs are not applied at this stage.
# ------------------------------------------------------------------------------

qc_keep <- with(
  p02@meta.data,
  nFeature_RNA >= QC_MIN_FEATURES &
    percent.mt <= QC_MAX_MT
)

p02_qc <- subset(
  p02,
  cells = rownames(p02@meta.data)[qc_keep]
)


# ------------------------------------------------------------------------------
# 21. Validate QC-filtered object
#     QC filtering 결과 확인
# ------------------------------------------------------------------------------

dim(p02_qc)

# Number of retained cells per sample
table(p02_qc$sample_id)