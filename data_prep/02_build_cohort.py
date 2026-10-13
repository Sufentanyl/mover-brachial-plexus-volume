#!/usr/bin/env python3
"""02_build_cohort.py -- 阻滞队列初筛: CPT 队列 693 人, 匹配局麻给药与手术信息.

输入:  $WORK_DIR/block_mrns.txt, $MOVER_DIR/cohort_info.json,
       $WORK_DIR/meds_la_given.csv
输出:  $WORK_DIR/block_cohort.csv (693 行; matched_la=1 者 644)

逻辑:
  - 局麻给药判定: MED_ROUTE_NM in {perineural, nerve catheter, injection}
    (注: 'Injection' 是阻滞记录的另一种写法, 经核查 1131/1188 行 ENC_TYPE='Anesthesia'
     且 1182/1188 与手术同日、主力为 ROPIVACAINE 0.5% IJ SOLN, 可判定为阻滞给药)
  - 按给药日期匹配手术日期(cohort_info.json 的 SURGERY_DATE), 取最早手术日;
    当日最早一条给药为 LA_FIRST_TIME, 手术信息取自该日对应的 info 行。
  - matched_la=0: 无匹配给药 (49 人无 LA 记录 + 44 人给药不在手术日)。

原始脚本: build_cohort.py (2026-10-10)
"""
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime

csv.field_size_limit(sys.maxsize)
MOVER_DIR = os.environ.get("MOVER_DIR", ".")
WORK_DIR = os.environ.get("WORK_DIR", ".")


def parse_conc(name):
    m = re.search(r"(\d+\.?\d*)\s*%", name)
    if m:
        return m.group(1) + "%"
    m = re.search(r"(\d+\.?\d*)\s*MG/ML", name, re.I)
    if m:
        return m.group(1) + "mg/ml"
    return ""


def clean_weight(w):
    try:
        v = float(w)
        kg = v / 16 / 2.20462
        return round(kg, 1) if 30 <= kg <= 180 else ""
    except Exception:
        return ""


def clean_height(h):
    m = re.match(r"""\s*(\d+)\s*'\s*(\d+)""", h or "")
    if m:
        return round(int(m.group(1)) * 30.48 + int(m.group(2)) * 2.54, 1)
    return ""


def norm_sdate(s):
    s = (s or "").strip()
    for fmt in ("%m/%d/%y %H:%M", "%m/%d/%Y %H:%M", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d"):
        try:
            return datetime.strptime(s, fmt).strftime("%Y-%m-%d")
        except ValueError:
            pass
    return ""


BLOCK_ROUTES = {"perineural", "nerve catheter", "injection"}


def is_block_la(row):
    return (row.get("MED_ROUTE_NM") or "").strip().lower() in BLOCK_ROUTES


with open(os.path.join(WORK_DIR, "block_mrns.txt")) as f:
    BLOCK = [x.strip() for x in f if x.strip()]
with open(os.path.join(MOVER_DIR, "cohort_info.json")) as f:
    cohort_info = json.load(f)

la_by_mrn = defaultdict(list)
n_la_total = 0
n_la_regional = 0
n_no_time = 0
with open(os.path.join(WORK_DIR, "meds_la_given.csv"), newline="") as fh:
    r = csv.DictReader(fh)
    for row in r:
        n_la_total += 1
        if not is_block_la(row):
            continue
        n_la_regional += 1
        if not (row.get("MED_ACTION_TIME") or "").strip():
            n_no_time += 1
            continue
        la_by_mrn[row["MRN"].strip()].append(row)

print(f"la_given_total={n_la_total} regional_la={n_la_regional} no_action_time={n_no_time}")

out_rows = []
match_stats = Counter()
for mrn in BLOCK:
    infos = cohort_info.get(mrn, [])
    surg_dates = {}
    for info in infos:
        sd = norm_sdate(info.get("SURGERY_DATE"))
        if sd:
            surg_dates.setdefault(sd, info)
    base = {"MRN": mrn, "matched_la": 0, "LOG_ID": "", "SURGERY_DATE": "", "AGE": "",
            "SEX": "", "ASA": "", "ANES_TYPE": "", "PROCEDURE": "", "LA_DRUG": "",
            "LA_CONC": "", "LA_DOSE_TEXT": "", "LA_DOSE_UNIT": "", "LA_ROUTE": "",
            "LA_FIRST_TIME": "", "AN_START": "", "AN_STOP": "", "IN_OR": "", "OUT_OR": "",
            "ADMSN": "", "DISCH": "", "LOS": "", "WEIGHT_KG": "", "HEIGHT_CM": "",
            "N_LA_RECORDS": 0, "N_SURGERIES": len(infos)}
    la_rows = la_by_mrn.get(mrn, [])
    matched = [rw for rw in la_rows if (rw.get("MED_ACTION_TIME") or "")[:10] in surg_dates]
    if not matched:
        out_rows.append(base)
        match_stats["unmatched"] += 1
        continue
    match_stats["matched"] += 1
    by_date = defaultdict(list)
    for rw in matched:
        by_date[(rw.get("MED_ACTION_TIME") or "")[:10]].append(rw)
    sdate = sorted(by_date)[0]
    recs = sorted(by_date[sdate], key=lambda x: x.get("MED_ACTION_TIME") or "")
    first = recs[0]
    info = surg_dates[sdate]
    base.update({
        "matched_la": 1, "LOG_ID": info["LOG_ID"], "SURGERY_DATE": sdate,
        "AGE": info["AGE"], "SEX": info["SEX"], "ASA": info["ASA"],
        "ANES_TYPE": info["ANES_TYPE"], "PROCEDURE": info["PROCEDURE"],
        "LA_DRUG": (first.get("MEDICATION_NM") or "").strip(),
        "LA_CONC": parse_conc((first.get("MEDICATION_NM") or "") + " " + (first.get("DISPLAY_NAME") or "")),
        "LA_DOSE_TEXT": (first.get("ADMIN_SIG") or "").strip(),
        "LA_DOSE_UNIT": (first.get("DOSE_UNIT_NM") or "").strip(),
        "LA_ROUTE": (first.get("MED_ROUTE_NM") or "").strip(),
        "LA_FIRST_TIME": (first.get("MED_ACTION_TIME") or "").strip(),
        "AN_START": info["AN_START"], "AN_STOP": info["AN_STOP"],
        "IN_OR": info["IN_OR"], "OUT_OR": info["OUT_OR"],
        "ADMSN": info["ADMSN"], "DISCH": info["DISCH"], "LOS": info["LOS"],
        "WEIGHT_KG": clean_weight(info["WEIGHT"]), "HEIGHT_CM": clean_height(info["HEIGHT"]),
        "N_LA_RECORDS": len(recs), "N_SURGERIES": len(infos),
    })
    drugs = sorted(set((x.get("MEDICATION_NM") or "").strip() for x in recs))
    base["LA_DRUG_ALL"] = " | ".join(drugs)
    out_rows.append(base)

cols = ["MRN", "matched_la", "LOG_ID", "SURGERY_DATE", "AGE", "SEX", "ASA", "ANES_TYPE",
        "PROCEDURE", "LA_DRUG", "LA_DRUG_ALL", "LA_CONC", "LA_DOSE_TEXT", "LA_DOSE_UNIT",
        "LA_ROUTE", "LA_FIRST_TIME", "AN_START", "AN_STOP", "IN_OR", "OUT_OR", "ADMSN",
        "DISCH", "LOS", "WEIGHT_KG", "HEIGHT_CM", "N_LA_RECORDS", "N_SURGERIES"]
out_path = os.path.join(WORK_DIR, "block_cohort.csv")
with open(out_path, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
    w.writeheader()
    w.writerows(out_rows)
print(f"cohort rows={len(out_rows)} match_stats={dict(match_stats)}")
print("saved", out_path)
