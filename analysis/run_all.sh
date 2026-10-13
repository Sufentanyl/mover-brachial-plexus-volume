#!/bin/bash
# MOVER 臂丛阻滞分析 - 端到端运行脚本（571 版主队列）
# 用法： MOVER_DIR=/path/to/MOVER WORK_DIR=/path/to/work OUT_DIR=/path/to/out bash run_all.sh
#   MOVER_DIR: MOVER v2 原始数据目录（含 EPIC_EMR 等）
#   WORK_DIR:  中间文件目录
#   OUT_DIR:   输出目录（analysis_571.csv、结果表、图）
set -euo pipefail

: "${MOVER_DIR:?请设置 MOVER_DIR（MOVER v2 原始数据目录）}"
: "${WORK_DIR:?请设置 WORK_DIR（中间文件目录）}"
: "${OUT_DIR:?请设置 OUT_DIR（输出目录）}"
mkdir -p "$WORK_DIR" "$OUT_DIR"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_PREP="$SCRIPT_DIR/../data_prep"

echo "=== 0. 数据准备（上游） ==="
python3 "$DATA_PREP/00_cpt_screen.py"
python3 "$DATA_PREP/01_med_pass.py"
python3 "$DATA_PREP/02_build_cohort.py"
python3 "$DATA_PREP/03_cohort_logid.py"
python3 "$DATA_PREP/04_multidose_detail.py"
python3 "$DATA_PREP/05_pain_outcomes.py"
python3 "$DATA_PREP/06_outcomes_fix.py"
python3 "$DATA_PREP/07_dex_extract.py"

echo "=== 1. 构建 571 分析数据集 ==="
python3 "$SCRIPT_DIR/build_571_dataset.py"
python3 "$SCRIPT_DIR/verify_trigger_condition.py"

echo "=== 2. 主分析（R） ==="
Rscript "$SCRIPT_DIR/run_v26.R"

echo "=== 3. 敏感性分析 ==="
Rscript "$SCRIPT_DIR/01_mice_extended.R"
Rscript "$SCRIPT_DIR/02_no_shoulder.R"
python3 "$SCRIPT_DIR/03_dedup.py"
Rscript "$SCRIPT_DIR/03_dedup_cox.R"
Rscript "$SCRIPT_DIR/sens_t0_early_excl.R"

echo "=== 3.5 B1/B2 输入准备 ==="
# B1/B2 需要队列级 procedure events 和 MAR 提取文件：
#   $WORK_DIR/rerun5/inputs/procedure_events_cohort.csv、mar_cohort.csv
# 其提取脚本尚未加入本仓库（见 README）。
test -f "$WORK_DIR/rerun5/inputs/procedure_events_cohort.csv" || { echo "缺少 B1/B2 输入文件，见 README"; exit 1; }
echo "=== 4. SAP 补充分析 ==="
python3 "$SCRIPT_DIR/c1_abc.py"
python3 "$SCRIPT_DIR/c2_riskdiff.py"
python3 "$SCRIPT_DIR/d1_ipw_291.py"
python3 "$SCRIPT_DIR/e3_bias_v2.py"
python3 "$SCRIPT_DIR/04_moduleB_revised.py"

echo "=== 完成 ==="
echo "输出目录：$OUT_DIR"
