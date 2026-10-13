#!/usr/bin/env python3
"""extract_pain_meas.py -- 从 EPIC flowsheet 测量表提取疼痛相关测量
(data_prep/05_pain_outcomes.py 的输入 pi_pain_meas.txt).

输入:
  EPIC_patient_measurments.tar.gz (MOVER_TAR 环境变量,
      或 $MOVER_DIR/EPIC_patient_measurments.tar.gz,
      或 /srv/data/gdrive/MOVER/EPIC_patient_measurments.tar.gz)
      内 srv/disk00/EPIC_flowsheets/ 下全部 19 个 flowsheet 成员
      (27 字段, MRN=第4列f3, FLO_MEAS_NAME=第20列f19,
       RECORDED_TIME=第23列f22, MEAS_VALUE=第24列f23; 行以 CRLF 结尾,
       字段内引号用反斜杠转义)
输出:
  $WORK_DIR/pi_pain_meas.txt (33,277,493 行; 原始行原样写出, 含 CRLF,
      无表头; 下游 05 只取 FLO_MEAS_NAME == "UC R PAIN INTENSITY/SCORE" 的行)

提取逻辑 (由现有 pi_pain_meas.txt 反推, 2026-10-13):
  - 不是全量导出, 是过滤导出: 按 tar 内成员顺序, 对每个成员跳过表头后,
    只保留行文本含 "pain" (大小写不敏感, ASCII 子串匹配) 的原始行,
    行内容与顺序保持原表一致 (逐字节)。
  - 成员顺序即 tar 内顺序: 1_2018-4633_flowsheets_20211006.csv (171GB),
    2_2015-1924_flowsheets_201905221700.csv,
    3_2018-4633_FLOWSHEETS_20221014.csv,
    19_2018-4633_FLOWSHEETSB_part15.csv … 4_2018-4633_FLOWSHEETSB_part0.csv。
  - 重要: 现有 pi_pain_meas.txt (33,277,493 行) 似为 90 分钟超时中断的
    部分提取。全成员重跑得 46,398,060 pain 行 (1.44B 行扫描);
    现有文件尾行对应全量输出的第 33,271,416 行 (成员约 part7 处),
    与 2026-10-09 15:45→17:15 (90 分钟) 的生成窗口吻合。
    本脚本做全成员完整提取 (46.4M 行), 为超集; 下游 05 只取
    UC R PAIN INTENSITY/SCORE 行, 功能等价。
  - 另有 6,076 行 (0.018%) 不含 "pain" 的非疼痛测量行 (f19 为
    PATIENT OBSERVATIONS / REMOVAL DATE 等 18 种), 其入选规则未能反推;
    全为非疼痛测量, 对下游无影响。
  - 注意: tar 共 ~17.8GB (压缩) / ~502GB (raw), 流式读取耗时约 1.5-2 小时;
    输出约 16GB (全量)。
  - 行内重复: 原表同一测量行可重复出现 (最高 40 次),
    本脚本不去重, 与现有文件一致。

原始脚本: /tmp/pi_pain_scan (2026-10-09, 已佚);
  提取逻辑由 pi_pain_meas.txt 反推 (输出前 5195 行与成员 pain 行逐字节一致)。
"""
import os
import tarfile

MOVER_DIR = os.environ.get("MOVER_DIR", "/srv/data/mover_analysis")
WORK_DIR = os.environ.get("WORK_DIR", "/srv/data/mover_analysis")

TAR_CANDIDATES = [
    os.environ.get("MOVER_TAR", ""),
    os.path.join(MOVER_DIR, "EPIC_patient_measurments.tar.gz"),
    "/srv/data/gdrive/MOVER/EPIC_patient_measurments.tar.gz",
]
tar_path = next((p for p in TAR_CANDIDATES if p and os.path.exists(p)), None)
if tar_path is None:
    raise SystemExit("EPIC_patient_measurments.tar.gz not found; set MOVER_TAR. "
                     "tried: %s" % TAR_CANDIDATES)

MEMBER_PREFIX = "srv/disk00/EPIC_flowsheets/"

out_path = os.path.join(WORK_DIR, "pi_pain_meas.txt")
print("reading", tar_path, flush=True)

n = 0
kept = 0
n_members = 0
# 注意: "r|gz" 为流模式, 成员须在迭代中即时处理, 不可先建表再回头读取
with tarfile.open(tar_path, "r|gz") as tf:
    with open(out_path, "wb") as out:
        for mm in tf:
            if not (mm.isfile() and mm.name.startswith(MEMBER_PREFIX)):
                continue
            n_members += 1
            print("member %d: %s" % (n_members, mm.name), flush=True)
            fh = tf.extractfile(mm)
            fh.readline()  # 跳过表头 (表头不含 "pain", 不写入)
            for line in fh:
                n += 1
                if b"pain" in line.lower():
                    out.write(line)
                    kept += 1
                if n % 50000000 == 0:
                    print("scanned %d, kept %d" % (n, kept), flush=True)

if n_members == 0:
    raise SystemExit("no flowsheet members found in %s" % tar_path)
print("done: members %d, scanned %d, kept %d -> %s"
      % (n_members, n, kept, out_path), flush=True)
