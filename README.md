# P02 (HN00262036) scRNA-seq Analysis

이 repository는 자체 생산한 **P02 (HN00262036) scRNA-seq 데이터**의 분석 코드, 데이터셋별 설정, 분석 의사결정 기록 및 결과물을 관리합니다.

분석의 목표는 단순히 결과를 생성하는 데 그치지 않고, 각 단계에서 **무엇을 실행했는지**, **어떤 값을 사용했는지**, **왜 해당 분석 방법을 선택했는지**를 재현 가능한 형태로 남기는 것입니다.

> 핵심 원칙: **QC 절차는 통일하고, QC 조건은 데이터셋별로 결정합니다.**

## Repository 역할

이 Git repository를 P02 분석 자산의 **Source of Truth**로 사용합니다. 현재 repository의 코드와 문서는 이전 대화, 별도 코드 사본 또는 참고 자료보다 우선합니다.

주요 구성 요소의 역할은 다음과 같습니다.

|구성 요소|역할|
|---|---|
|`[번호]_[annotation].R`|실제 분석 workflow와 실행 코드|
|`config.R`|P02에서 사용할 dataset-specific parameter|
|`DECxx_[annotation].md`|중요한 분석 결정의 검토 과정과 근거|
|`analysis_decisions.md`|사용자가 최종 승인한 결정의 index|
|`sample_metadata.csv`|분석에 사용하는 sample metadata|
|`shared/`|여러 데이터셋에서 재사용 가능한 공통 코드|
|`figures/`|분석 figure|
|`tables/`|분석 table|
|`report/`|분석 report|

## 권장 Repository 구조

아래 구조를 기본으로 사용하되, 아직 필요하지 않은 파일이나 디렉터리는 미리 만들지 않을 수 있습니다.

```
P02_scRNAseq/
├── README.md
├── config.R
├── sample_metadata.csv
├── analysis_decisions.md
├── 01_load_data.R
├── 02_create_object.R
├── 03_qc_metrics.R
├── 04_qc_filtering.R
├── DEC04_qc_filtering.md
├── shared/
│   ├── utils.R
│   ├── qc_functions.R
│   ├── plotting.R
│   ├── seurat_helpers.R
│   └── marker_functions.R
├── figures/
├── tables/
└── report/
```

실제 workflow가 진행되면서 script와 Decision 문서는 순차적으로 추가합니다. 위 파일 목록은 구조를 설명하기 위한 예시이며, 존재하지 않는 분석 단계나 파일을 임의로 생성하지 않습니다.

## Script 명명 규칙

분석 script는 다음 형식을 사용합니다.

```
[번호]_[annotation].R
```

예:

```
01_load_data.R
02_create_object.R
03_qc_metrics.R
04_qc_filtering.R
05_doublet_detection.R
06_normalization.R
07_clustering.R
```

- 번호는 workflow의 실행 순서를 나타냅니다.
    
- annotation은 해당 script의 역할을 짧고 명확한 영어 `snake_case`로 작성합니다.
    
- 기존 번호 체계는 특별한 이유 없이 변경하지 않습니다.
    
- 한 script가 지나치게 많은 분석 단계를 담당하지 않도록 논리적 단위로 분리합니다.
    

## Decision 문서 명명 규칙

분석적 판단이 중요한 단계에는 다음 형식의 Decision 문서를 사용할 수 있습니다.

```
DEC[script 번호]_[annotation].md
```

예:

```
04_qc_filtering.R
DEC04_qc_filtering.md

07_clustering.R
DEC07_clustering.md
```

Decision 문서를 생성하는 경우 `DEC` 번호는 대응하는 script 번호와 일치해야 합니다.

모든 script에 Decision 문서가 필요한 것은 아닙니다. 파일 로딩, 단순 입출력 또는 별도의 분석 판단이 없는 반복 처리에는 Decision 문서를 만들지 않을 수 있습니다.

### Decision 문서에 기록할 내용

각 `DECxx_*.md`에는 가능한 범위에서 다음 내용을 기록합니다.

1. 결정 목적
    
2. 검토한 데이터 또는 결과
    
3. 검토한 후보와 각 장단점
    
4. 제안과 그 근거
    
5. 사용자 결정
    
6. 최종 설정
    
7. 최종 결정 근거
    

즉, Decision 문서는 **제안 → 검토 → 사용자 선택 → 근거**의 과정을 보존합니다.

## `analysis_decisions.md`

`analysis_decisions.md`는 상세한 논의 문서가 아니라 프로젝트의 Decision index 또는 ledger입니다.

- 사용자가 **최종 승인한 결정만** 기록합니다.
    
- 아직 검토 중인 제안은 기록하지 않습니다.
    
- 상세 근거는 대응하는 `DECxx_*.md`에 기록합니다.
    

권장 형식:

```
| Decision | Topic | Final decision | Status |
|---|---|---|---|
| DEC04 | QC filtering | 승인된 QC 설정 요약 | Approved |
| DEC07 | Clustering | 승인된 clustering 설정 요약 | Approved |
```

## `config.R`

공통 workflow는 유지하면서 데이터셋마다 달라지는 값은 가능한 한 `config.R`에서 관리합니다.

예:

- QC thresholds
    
- doublet detection parameter
    
- normalization 또는 integration 관련 승인 설정
    
- PCA 및 사용할 dimension
    
- clustering resolution
    
- marker analysis threshold
    
- figure size, resolution 및 저장 형식
    
- random seed
    

중요한 parameter가 아직 결정되지 않았다면 임의로 확정하지 않고 `NA` 또는 미설정 상태로 둘 수 있습니다. 데이터와 중간 결과를 검토하고 사용자가 승인한 뒤 반영합니다.

## `shared/`

`shared/`에는 특정 데이터셋에 종속되지 않고 여러 프로젝트에서 동일하게 재사용할 수 있는 코드를 둡니다.

|파일|권장 역할|
|---|---|
|`utils.R`|범용 보조 함수, validation, 안전한 파일 저장|
|`qc_functions.R`|공통 QC 계산 및 plotting 함수|
|`plotting.R`|공통 theme, palette, figure 저장 함수|
|`seurat_helpers.R`|반복되는 Seurat 처리 helper|
|`marker_functions.R`|marker 결과 정리 등 공통 함수|

판단 기준은 다음과 같습니다.

> P02에서만 필요한가? → `shared/`에 넣지 않습니다.  
> 다른 데이터셋에서도 같은 코드로 사용할 수 있는가? → `shared/` 후보로 검토합니다.

한 데이터셋에서 한 번 사용한 코드를 바로 공통 함수로 만들지 않습니다. 특정 데이터셋에서 문제가 발생한 경우 먼저 `config.R` 또는 해당 dataset script에서 해결할 문제인지 확인하고, 공통 함수 자체의 문제일 때만 `shared/` 수정을 검토합니다.

## 결과물 관리

분석 결과는 용도에 따라 다음 디렉터리에 저장합니다.

```
figures/   # plots and publication figures
tables/    # result tables and summaries
report/    # HTML, PDF or other analysis reports
```

별도의 `outputs/` 디렉터리는 사용하지 않습니다. 결과물은 가능한 한 생성한 script와 연결해 추적할 수 있도록 파일명을 정합니다.

대용량 raw data, FASTQ 파일, 중간 객체 및 개인 환경에 종속된 파일은 Git에 직접 포함하지 않는 것을 원칙으로 하며, 실제 추적 범위는 `.gitignore`와 데이터 관리 정책에 따릅니다.

## 분석 및 코드 작성 원칙

1. 전체 workflow를 한 번에 확정하지 않고 작은 단계로 작성하고 검증합니다.
    
2. 이전 단계의 실행 결과를 확인한 뒤 다음 단계로 진행합니다.
    
3. 오류가 발생하면 원인을 확인하기 전에 다음 단계로 넘어가지 않습니다.
    
4. QC workflow는 가능한 한 통일하되 threshold는 데이터셋별 결과를 근거로 결정합니다.
    
5. normalization, integration, doublet detection, clustering 및 differential expression처럼 결과에 큰 영향을 주는 방법은 사용자의 승인 없이 확정하지 않습니다.
    
6. 실제 실행하지 않은 코드를 실행했다고 기록하지 않습니다.
    
7. 존재하지 않는 결과, 수치, table 또는 figure를 생성하거나 추정하지 않습니다.
    
8. 데이터에서 관찰한 사실, 가능한 해석, 추가 검증이 필요한 가설을 구분합니다.
    
9. 최종 script에는 목적, 입력, 주요 처리 과정, parameter 및 출력물을 이해할 수 있는 한국어·영어 혼용 주석을 포함합니다.
    

## Git 운영 원칙

- commit, push, merge는 사용자가 직접 결정하고 수행합니다.
    
- 분석 보조 과정에서는 변경이 하나의 논리적 단위로 완료될 때 commit 대상 파일과 commit message를 제안할 수 있습니다.
    
- 사용자 승인 없이 기존 파일을 삭제하거나 덮어쓰지 않습니다.
    
- 분석 코드, 설정 및 문서 변경은 가능한 한 하나의 목적이 분명한 commit 단위로 나눕니다.
    

Commit 예시:

```
feat(qc): define P02 QC filtering workflow and thresholds
```

Commit body 예시:

```
- add dataset-specific QC thresholds to config.R
- document QC threshold rationale in DEC04
- implement filtering workflow in 04_qc_filtering.R
```

## 분석 실행 전 확인 사항

- 현재 작업 중인 branch와 repository 상태를 확인합니다.
    
- `sample_metadata.csv`의 sample 정보가 실제 입력 데이터와 일치하는지 확인합니다.
    
- `config.R`에서 아직 승인되지 않은 필수 parameter가 있는지 확인합니다.
    
- script 번호 순서와 필요한 입력 파일을 확인합니다.
    
- 결과 파일을 저장할 디렉터리가 repository 규칙과 일치하는지 확인합니다.
    

## 재현성과 변경 기록

분석 결과를 재현할 수 있도록 다음 정보를 가능한 범위에서 보존합니다.

- 사용한 package와 version
    
- random seed
    
- 입력 데이터 또는 객체의 식별 정보
    
- 주요 parameter
    
- 실행한 script 순서
    
- 승인된 분석 결정과 근거
    

분석 조건이 변경되면 관련 `config.R`, script, Decision 문서 및 `analysis_decisions.md` 사이의 일관성을 함께 확인합니다.
