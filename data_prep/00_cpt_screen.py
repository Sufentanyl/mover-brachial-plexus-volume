#!/usr/bin/env python3
"""00_cpt_screen.py -- CPT 筛选: 从 EPIC_EMR.tar.gz 提取计费编码 64415/64416 的去重 MRN.

输入:  EPIC_EMR.tar.gz (MOVER_TAR 环境变量, 或 $MOVER_DIR/EPIC_EMR.tar.gz)
输出:  $WORK_DIR/block_mrns.txt (排序去重, 693 行)

原始脚本: med_scan1.py (2026-10-10)
"""
import csv
import io
import os
import tarfile

MOVER_DIR = os.environ.get("MOVER_DIR", ".")
WORK_DIR = os.environ.get("WORK_DIR", ".")

CANDIDATES = [
    os.environ.get("MOVER_TAR", ""),
    os.path.join(MOVER_DIR, "EPIC_EMR.tar.gz"),
]
tar_path = next((p for p in CANDIDATES if p and os.path.exists(p)), None)
if tar_path is None:
    raise SystemExit("EPIC_EMR.tar.gz not found; set MOVER_TAR. tried: %s" % CANDIDATES)
print("using", tar_path, flush=True)

tf = tarfile.open(tar_path, "r:gz")
mrns = set()
n = 0
m = tf.getmember("EPIC_EMR/EMR/patient_coding.csv")
f = tf.extractfile(m)
text = io.TextIOWrapper(f, encoding="utf-8", errors="replace")
r = csv.DictReader(text)
for row in r:
    n += 1
    if (row.get("REF_BILL_CODE") or "").strip() in ("64415", "64416"):
        mrns.add((row.get("MRN") or "").strip())
print("coding rows:", n, "block MRNs:", len(mrns), flush=True)

out_path = os.path.join(WORK_DIR, "block_mrns.txt")
with open(out_path, "w") as out:
    for x in sorted(mrns):
        out.write(x + "\n")
print("saved", out_path, flush=True)
