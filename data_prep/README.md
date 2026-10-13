# data_prep/ — MOVER 上游数据管线（00–07）

从已准备的 MOVER Epic 输入生成 `analysis_571.csv` 所需的上游中间文件。
脚本按编号顺序运行；每一步的输入/输出与原始脚本对应关系见下表。
原始脚本为 2026-10-10 的版本，此处仅做了路径环境变量化，逻辑逐行一致，
且每一步都已验证重跑输出与现有 CSV 逐格一致（见下方验证表）。

本目录尚不能从原始压缩包完整运行：`cohort_info.json` 的生成、
`pi_pain_meas.txt` 的提取及 B1/B2 的 procedure events / MAR 提取代码尚未收录。
这些输入需由原分析工作目录提供。原始代码正在等待转交，由作者后续上传；
下方既有步骤的验证记录不代表缺失步骤已在本仓库重跑。

## 路径环境变量（另见 ../config.yaml）

| 变量 | 含义 | 默认值 |
|---|---|---|
| `MOVER_DIR` | 原始数据：`patient_information.csv`、`patient_medications.csv`、`cohort_info.json`、`EPIC_EMR.tar.gz`、`pi_pain_meas.txt` | 无（需设置） |
| `WORK_DIR` | 中间文件：`block_cohort.csv`、`meds_*.csv`、`results_logid/`、`results_logid_fix/` 等 | 无（需设置） |
| `OUT_DIR` | 最终输出：`analysis_571.csv` 与模型结果 | 无（需设置） |
| `MOVER_TAR` | CPT 筛选 tar 包（可选；默认 `$MOVER_DIR/EPIC_EMR.tar.gz`） | — |
| `PAIN_MEAS` | 疼痛测量提取文件（可选；默认 `$MOVER_DIR/pi_pain_meas.txt`） | — |

## 步骤

| # | 脚本 | 输入 | 输出 |
|---|---|---|---|
| 00 | `00_cpt_screen.py` | `EPIC_EMR.tar.gz` 内 `patient_coding.csv` | `block_mrns.txt`（CPT 64415/64416 去重 MRN，693） |
| 01 | `01_med_pass.py` | `patient_medications.csv`（7.2GB）、`block_mrns.txt` | `meds_la_given.csv`（全部 Given 局麻药）、`cohort_meds_opioid_given.csv`、`cohort_meds_nonopioid_given.csv`（693 队列）+ 分布统计 json |
| 02 | `02_build_cohort.py` | `block_mrns.txt`、`cohort_info.json`、`meds_la_given.csv` | `block_cohort.csv`（693；matched_la=1 者 644） |
| 03 | `03_cohort_logid.py` | `block_cohort.csv`、`meds_la_given.csv`、`patient_information.csv` | `results_logid/cohort_logid.csv`（600；含 `reconciliation.json`） |
| 04 | `04_multidose_detail.py` | `cohort_logid.csv`、`meds_la_given.csv` | `results_logid/step3_multidose_detail.csv`（775 行，T0 后 60min 逐次给药） |
| 05 | `05_pain_outcomes.py` | `cohort_logid.csv`、`pi_pain_meas.txt`（11.6GB）、`cohort_meds_*.csv` | `results_logid/pain_cohort.csv`、`results_logid/outcomes_logid.csv`（初版结局） |
| 06 | `06_outcomes_fix.py` | `cohort_logid.csv`、`pain_cohort.csv`、`outcomes_logid.csv`、`cohort_meds_*.csv` | `results_logid_fix/outcomes_fix.csv`（600；终版主结局，80 事件） |
| 07 | `07_dex_extract.py` | `cohort_logid.csv`、`patient_medications.csv` | `results_logid/dex_cohort.csv`、`dex_preend_flag.csv` |

之后 `../analysis/build_571_dataset.py` 从上游文件组装 `analysis_571.csv`（571；
规则见该脚本头注）。

## 关键口径（备查）

- **03**：长效局麻 = `MED_ROUTE_NM ∈ {PeriNEURAL, Injection}` 且 `DISPLAY_NAME`
  含 ropivacaine|bupivacaine；按 LOG_ID 左连 `patient_information`（仅 MRN 一致）；
  解剖合理性 YES_KEYS/NO_KEYS（NO 优先）；每人取最早通过给药为索引阻滞，
  T0=该给药时间；`age`=`BIRTH_DATE`（该列即年龄整数）；`surgery_type` 关键词
  （fracture: FRACTURE|ORIF|DISLOCATION；av_fistula: AV FISTULA|FISTULOGRAM|
  AV GRAFT|ARTERIOVENOUS|THROMBECTOMY|SHUNT；其余 other_ue）。
- **05**：疼痛提取 `MEAS == "UC R PAIN INTENSITY/SCORE"`；初版结局：疼痛≥4
  （≥t_end 且 ≥T0）后 4h 内有镇痛药。
- **06（终版主结局）**：触发疼痛 ≥4（T0..T0+96h，按时间依次尝试）；
  候选药须为术后给药（`t > an_stop`；an_stop 缺失时 `RECORD_TYPE=='POST-OP'`）；
  `tp <= t <= tp+4h`（同分钟计入）；删失 `min(出院, T0+96h)`。

## 验证（2026-10-13 重跑，输出与现有 CSV 比对）

| 步骤 | 输出文件 | 行数 | 比对结果 |
|---|---|---|---|
| 00 | `block_mrns.txt` | 693 | 与原始输出完全一致 |
| 01 | `meds_la_given.csv` 等 3 个 | — | 27,961,524 行；la_given=179,468；恢复后 build 零差异 ✓ |
| 02 | `block_cohort.csv` | 693×27 | 零格差异，行序一致 |
| 03 | `cohort_logid.csv` | 600×35 | 零格差异，行序一致 |
| 04 | `step3_multidose_detail.csv` | 775×6 | 零格差异，行序一致 |
| 05 | `pain_cohort.csv` / `outcomes_logid.csv` | 231241×3 / 600×6 | 排序后完全一致 |
| 06 | `outcomes_fix.csv` | 600×6 | 零格差异，行序一致（80 事件） |
| 07 | `dex_cohort.csv` / `dex_preend_flag.csv` | 877×5 / 334×2 | 零格差异，行序一致 |
