# DEC02_create_object.md

## 목적

P02 (HN00262036) scRNA-seq 데이터에서 sample별 10X count matrix를 어떤 형태로 읽고,  
초기 Seurat object를 어떤 원칙으로 생성 및 merge할지 결정한다.

대응 script:

- `02_create_object.R`
    

---

## 검토한 데이터 또는 결과

대상 sample:

|sample_id|pdo_source|timepoint|biological_replicate|culture_batch|
|---|---|---|---|---|
|naive_D1|naive|1day|rep1|D1_setup|
|naive_D7|naive|7day|rep1|D7_setup|
|FR_D1|FR|1day|rep1|D1_setup|
|FR_D7|FR|7day|rep1|D7_setup|

Raw data base path:

```
~/data/HN00262036/HN00262036_result_10X
```

각 sample에서 다음 두 종류의 matrix가 존재함을 확인했다.

```
sample_raw_feature_bc_matrix/
sample_filtered_feature_bc_matrix/
```

초기 분석 입력으로 사용할 sample별 filtered matrix 경로:

```
OCM/per_sample_outs/<sample_id>/count/sample_filtered_feature_bc_matrix
```

각 filtered matrix에 다음 10X 파일이 모두 존재함을 확인했다.

```
barcodes.tsv.gz
features.tsv.gz
matrix.mtx.gz
```

### Read10X() 로딩 결과

|sample_id|class|features|filtered barcodes|
|---|---|---|---|
|naive_D1|dgCMatrix|38,606|7,751|
|naive_D7|dgCMatrix|38,606|7,369|
|FR_D1|dgCMatrix|38,606|6,997|
|FR_D7|dgCMatrix|38,606|9,261|

모든 sample에서 feature 수가 38,606으로 동일했다.

### Sample-level Seurat object 확인

각 sample에서 `CreateSeuratObject()` 후 feature 수와 cell 수가 `Read10X()` 결과와 동일함을 확인했다.

### Merge validation 결과

4개 sample-level Seurat object를 merge한 뒤 다음을 확인했다.

```
nrow(p02)
# 38606

ncol(p02)
# 31378

table(p02$sample_id)
#    FR_D1    FR_D7 naive_D1 naive_D7
#     6997     9261     7751     7369

anyDuplicated(colnames(p02))
# 0
```

따라서 merge 후:

- features: 38,606
    
- cells: 31,378
    
- sample별 cell 수 유지
    
- duplicated cell names 없음
    

을 확인했다.

---

## 검토한 후보

### 1. `sample_raw_feature_bc_matrix` 사용

장점:

- Cell Ranger cell calling 이전 barcode까지 포함한다.
    
- ambient RNA correction 등 일부 후속 분석에서 활용 가능하다.
    

주의점:

- 현재 일반적인 Seurat 초기 분석에는 불필요한 barcode까지 포함될 수 있다.
    
- sample-level QC workflow가 복잡해질 수 있다.
    

### 2. `sample_filtered_feature_bc_matrix` 사용

장점:

- Cell Ranger cell calling을 통과한 barcode를 sample별로 바로 사용할 수 있다.
    
- 일반적인 Seurat 분석 시작점으로 적절하다.
    
- 이후 `nFeature_RNA`, `nCount_RNA`, mitochondrial percentage, doublet detection 등 dataset-specific QC를 별도로 수행할 수 있다.
    

### 3. `OCM/multi/count/raw_feature_bc_matrix` 사용

주의점:

- 현재 분석은 sample별 metadata와 sample-specific QC 추적이 중요하다.
    
- `per_sample_outs` 아래에 네 sample의 filtered matrix가 명확히 존재한다.
    

---

## ChatGPT 제안

초기 분석 입력은 각 sample의:

```
sample_filtered_feature_bc_matrix
```

를 사용한다.

`CreateSeuratObject()` 시점에는 추가 QC filtering을 적용하지 않고:

```
min.cells = 0
min.features = 0
```

으로 설정한다.

QC metric을 실제로 계산하고 분포를 검토한 뒤 `04_qc_filtering.R` 단계에서 dataset-specific threshold를 결정한다.

sample 간 동일한 10X barcode 문자열 충돌을 방지하기 위해 cell name 앞에 `sample_id`를 추가한다.

---

## 사용자 결정

승인된 사항:

- sample별 `sample_filtered_feature_bc_matrix`를 초기 입력으로 사용
    
- `Read10X()`로 각 sample matrix 로딩
    
- `CreateSeuratObject()` 단계에서는 추가 QC filtering을 적용하지 않음
    
- `min.cells = 0`, `min.features = 0` 사용
    
- `RenameCells(..., add.cell.id = sample_id)`로 sample ID prefix 추가
    
- 다음 sample-level metadata를 각 cell에 추가
    
    - `sample_id`
        
    - `pdo_source`
        
    - `timepoint`
        
    - `biological_replicate`
        
    - `culture_batch`
        
- 4개 sample-level Seurat object를 하나의 P02 object로 merge
    

---

## 최종 설정

```
CreateSeuratObject(
  counts = count_matrices[[sample_id]],
  project = PROJECT_ID,
  assay = "RNA",
  min.cells = 0,
  min.features = 0
)
```

Cell naming:

```
RenameCells(
  obj,
  add.cell.id = sample_id
)
```

Merge:

```
p02 <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  project = PROJECT_ID
)
```

---

## 결정 근거

QC workflow 자체는 통일하되 QC threshold는 P02 데이터 분포를 확인한 뒤 결정한다는 프로젝트 원칙에 따른다.

따라서 object 생성 단계에서 임의의 `min.features` 또는 `min.cells` threshold로 barcode를 먼저 제거하지 않는다.

sample별 filtered matrix를 사용함으로써 sample metadata와 QC 결과를 명확하게 추적할 수 있다.

merge 후 sample별 cell 수가 원래 count matrix와 정확히 일치하고, cell name duplication이 없음을 확인했으므로 초기 P02 Seurat object 생성이 정상적으로 완료된 것으로 판단한다.

---

## 검증 상태

**Completed**

확인 완료:

- 4개 sample filtered 10X directory 존재
    
- 각 directory의 `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` 존재
    
- `Read10X()` 정상 로딩
    
- 모든 sample에서 38,606 features 확인
    
- sample별 Seurat object 생성 후 cell / feature 수 유지
    
- merged object: 38,606 features × 31,378 cells
    
- sample별 cell count 유지
    
- duplicated cell names = 0