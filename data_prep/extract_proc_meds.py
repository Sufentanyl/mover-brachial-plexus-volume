#!/usr/bin/env python3
"""extract_proc_meds.py -- 提取 procedure events / MAR 给药 / 队列时间窗
(04_moduleB_revised.py 的输入, rerun5/inputs/).

输入:
  EPIC_EMR.tar.gz (MOVER_TAR 环境变量, 或 $MOVER_DIR/EPIC_EMR.tar.gz,
      或 /srv/data/gdrive/MOVER/EPIC_EMR.tar.gz)
      成员 "EPIC_EMR/EMR/patient_procedure events.csv" (procedure events 全表)
  $MOVER_DIR/patient_medications.csv (MAR 给药全表, 7.7GB 解压副本)
  $MOVER_DIR/patient_information.csv (病人基本信息, 解压副本)
  $WORK_DIR/results_logid/cohort_logid.csv (600 例索引阻滞, 含 T0)
输出 ($WORK_DIR/rerun5/inputs/):
  procedure_events_cohort.csv
      (LOG_ID,EVENT_DISPLAY_NAME,EVENT_TIME; 5452 行; 队列=600 例全集)
  mar_cohort.csv
      (LOG_ID,MRN,DISPLAY_NAME,MEDICATION_NM,MAR_ACTION_NM,MED_ACTION_TIME,
       DOSE_UNIT_NM,MED_ROUTE_NM,ADMIN_SIG,ORDER_CLASS_NM;
       25579 行; 队列=下述 573 例; 不限 MAR 动作)
  cohort_times.csv
      (mrn,LOG_ID_x,T0,IN_OR_DTTM,OUT_OR_DTTM,AN_START_DATETIME,
       AN_STOP_DATETIME; 573 行; 按 mrn 升序)

队列规则:
  - cohort600 = cohort_logid.csv 全部 600 例 (LOG_ID_x)
  - cohort573 = cohort600 去掉 T0 >= 麻醉结束时间的 27 例术后阻滞
    (麻醉结束 = AN_STOP_DATETIME, 缺失则用 OUT_OR_DTTM)
  - procedure_events_cohort 用 cohort600 限定;
    mar_cohort / cohort_times 用 cohort573 限定
  - 行序: procedure_events / mar 保持原表行序;
    cohort_times 去重 (patient_information.csv 同一 LOG_ID 有完全相同的
    重复行, 取首行) 后按 mrn 升序; T0 取自 cohort_logid.csv

注意: EPIC_EMR.tar.gz 内 procedure events 是最后一个成员, 流式读取
  需顺序经过前面成员, 耗时较长属正常。

原始脚本: /tmp/sap/ 下诸脚本 (2026-10-10/11, 已佚);
  提取逻辑由 rerun5/inputs 现有文件反推, 重跑与现有三文件逐行一致。
"""
import csv
import io
import os
import sys
import tarfile
from datetime import datetime

csv.field_size_limit(sys.maxsize)

MOVER_DIR = os.environ.get("MOVER_DIR", "/srv/data/mover_analysis")
WORK_DIR = os.environ.get("WORK_DIR", "/srv/data/mover_analysis")

TAR_CANDIDATES = [
    os.environ.get("MOVER_TAR", ""),
    os.path.join(MOVER_DIR, "EPIC_EMR.tar.gz"),
    "/srv/data/gdrive/MOVER/EPIC_EMR.tar.gz",
]
tar_path = next((p for p in TAR_CANDIDATES if p and os.path.exists(p)), None)
if tar_path is None:
    raise SystemExit("EPIC_EMR.tar.gz not found; set MOVER_TAR. tried: %s" % TAR_CANDIDATES)

MED_PATH = os.path.join(MOVER_DIR, "patient_medications.csv")
INFO_PATH = os.path.join(MOVER_DIR, "patient_information.csv")
for p in (MED_PATH, INFO_PATH):
    if not os.path.exists(p):
        raise SystemExit("%s not found; set MOVER_DIR" % p)

COHORT_PATH = os.path.join(WORK_DIR, "results_logid", "cohort_logid.csv")
if not os.path.exists(COHORT_PATH):
    raise SystemExit("%s not found" % COHORT_PATH)

OUT_DIR = os.path.join(WORK_DIR, "rerun5", "inputs")
os.makedirs(OUT_DIR, exist_ok=True)

PE_MEMBER = "EPIC_EMR/EMR/patient_procedure events.csv"
PE_COLS = ["LOG_ID", "EVENT_DISPLAY_NAME", "EVENT_TIME"]
MAR_COLS = ["LOG_ID", "MRN", "DISPLAY_NAME", "MEDICATION_NM", "MAR_ACTION_NM",
            "MED_ACTION_TIME", "DOSE_UNIT_NM", "MED_ROUTE_NM", "ADMIN_SIG",
            "ORDER_CLASS_NM"]
CT_COLS = ["mrn", "LOG_ID_x", "T0", "IN_OR_DTTM", "OUT_OR_DTTM",
           "AN_START_DATETIME", "AN_STOP_DATETIME"]


def pdt(s):
    s = (s or "").strip()
    if not s:
        return None
    for f in ("%Y-%m-%d %H:%M:%S", "%m/%d/%y %H:%M", "%m/%d/%Y %H:%M"):
        try:
            return datetime.strptime(s, f)
        except ValueError:
            pass
    return None


# ---------- 队列 ----------
coh600 = []          # (LOG_ID_x, MRN, T0, AN_STOP, OUT_OR)
t0map = {}
with open(COHORT_PATH, newline="") as f:
    for r in csv.DictReader(f):
        lid = r["LOG_ID_x"].strip()
        coh600.append((lid, r["MRN"].strip(), r["T0"],
                       r["AN_STOP_DATETIME"], r["OUT_OR_DTTM"]))
        t0map[r["MRN"].strip()] = r["T0"]
coh600_ids = set(lid for lid, _, _, _, _ in coh600)
coh573_ids = set(
    lid for lid, _, t0, anstop, outor in coh600
    if not (pdt(t0) and (pdt(anstop) or pdt(outor))
            and pdt(t0) >= (pdt(anstop) or pdt(outor)))
)
print("cohort600: %d, cohort573: %d (excluded postop: %d)" %
      (len(coh600_ids), len(coh573_ids), len(coh600_ids) - len(coh573_ids)),
      flush=True)


# ---------- 1. procedure events (cohort600) ----------
print("reading procedure events from", tar_path, flush=True)
n_pe = 0
with tarfile.open(tar_path, "r:gz") as tf, \
        open(os.path.join(OUT_DIR, "procedure_events_cohort.csv"),
             "w", newline="") as out:
    m = tf.getmember(PE_MEMBER)
    w = csv.writer(out, lineterminator="\n")  # 原文件为 LF 换行
    w.writerow(PE_COLS)
    fh = io.TextIOWrapper(tf.extractfile(m), encoding="utf-8", errors="replace")
    for r in csv.DictReader(fh):
        if r["LOG_ID"] in coh600_ids:
            w.writerow([r[c] for c in PE_COLS])
            n_pe += 1
print("procedure_events_cohort.csv rows:", n_pe, flush=True)


# ---------- 2. MAR (cohort573) ----------
print("reading", MED_PATH, flush=True)
n_mar = 0
with open(MED_PATH, newline="") as fh, \
        open(os.path.join(OUT_DIR, "mar_cohort.csv"), "w", newline="") as out:
    w = csv.writer(out)
    w.writerow(MAR_COLS)
    for i, r in enumerate(csv.DictReader(fh)):
        if r["LOG_ID"] in coh573_ids:
            w.writerow([r[c] for c in MAR_COLS])
            n_mar += 1
        if (i + 1) % 5000000 == 0:
            print("  meds rows:", i + 1, flush=True)
print("mar_cohort.csv rows:", n_mar, flush=True)


# ---------- 3. cohort times (cohort573, 去重, 按 mrn 排序) ----------
print("reading", INFO_PATH, flush=True)
seen = set()
rows = []
with open(INFO_PATH, newline="") as fh:
    for r in csv.DictReader(fh):
        lid = r["LOG_ID"]
        if lid in coh573_ids and lid not in seen:
            seen.add(lid)
            rows.append([r["MRN"], lid, t0map.get(r["MRN"], ""),
                         r["IN_OR_DTTM"], r["OUT_OR_DTTM"],
                         r["AN_START_DATETIME"], r["AN_STOP_DATETIME"]])
rows.sort(key=lambda x: x[0])
with open(os.path.join(OUT_DIR, "cohort_times.csv"), "w", newline="") as out:
    w = csv.writer(out, lineterminator="\n")  # 原文件为 LF 换行
    w.writerow(CT_COLS)
    w.writerows(rows)
print("cohort_times.csv rows:", len(rows), flush=True)
print("done ->", OUT_DIR, flush=True)
