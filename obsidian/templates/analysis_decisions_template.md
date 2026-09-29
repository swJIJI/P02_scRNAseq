---
project: P02
dataset: HN00262036
type: dashboard
status: active
tags:
  - scRNAseq
---

# analysis_decisions.md

> Project dashboard / 분석 진행 상황과 승인된 주요 결정을 빠르게 확인하는 문서

## Project overview
- **Project ID:** P02
- **Dataset ID:** HN00262036
- **Analysis type:** scRNA-seq
- **Primary framework:** Seurat
- **Random seed:** 42
- **Repository:** `P02_scRNAseq`

## Analysis purpose
이 dataset의 분석 목적을 요약한다.

## Dataset information
| sample_id | pdo_source | timepoint | biological_replicate | culture_batch |
|---|---|---|---|---|
| ... | ... | ... | ... | ... |

## Current status

### 01. Step name
- Script: `scripts/01_step_name.R`
- Checkpoint: [[CP01_step_name]]
- Decision: [[DEC01_step_name]]
- Status: completed

확인된 핵심 결과를 짧게 기록한다.

## Decision summary
| Decision ID | Topic | Final decision | Status |
|---|---|---|---|
| [[DECXX_name|DECXX]] | ... | ... | approved |

## Parameters not yet finalized
- QC thresholds
- doublet detection method
- PCA dimensions
- clustering resolution

## Reproducibility
```text
report/session_info/
```

## AI-assisted analysis
- **AI system:** OpenAI ChatGPT
- **Model:** GPT-5.6 Sol
- **Role:** Analysis assistant / code and documentation support
- **Final analytical decisions:** User-approved
- **Source of Truth:** Git repository and verified execution results

## Related notes
- [[2026-09-16]]
- [[CPXX_step]]
- [[DECXX_step]]

## Next step
다음 분석 단계와 확인할 내용을 짧게 기록한다.
