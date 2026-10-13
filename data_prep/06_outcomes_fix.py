#!/usr/bin/env python3
"""06_outcomes_fix.py -- 主结局终版 (修正时间戳规则).

输入:  $WORK_DIR/results_logid/cohort_logid.csv,
       $WORK_DIR/results_logid/pain_cohort.csv,
       $WORK_DIR/results_logid/outcomes_logid.csv (仅用于新旧对比),
       $WORK_DIR/cohort_meds_opioid_given.csv, $WORK_DIR/cohort_meds_nonopioid_given.csv
输出:  $WORK_DIR/results_logid_fix/outcomes_fix.csv (600 行)
       $WORK_DIR/results_logid_fix/added_mrns.txt

终版定义 (FIXED timestamp rule: pain_time <= med_time, 同分钟计入):
  - 疼痛窗: T0..T0+96h; 触发疼痛: val>=4 (按时间排序, 去重后依次尝试)
  - 候选镇痛药: 术后给药 (t > an_stop; an_stop 缺失时 RECORD_TYPE=='POST-OP')
    阿片类全部计入; 非阿片类按 NSAID/对乙酰氨基酚模式计入 (排除外用/凝胶)
  - 事件: 首个触发疼痛后 4h 内 (tp <= t <= tp+4h) 有候选药, 事件时间=给药时间
  - 删失: min(出院, T0+96h)
  - event_sens/time_sens_h: 术后首个阿片类给药

原始脚本: fix_outcomes_final.py (2026-10-10)
"""
import csv
import os
import re
from datetime import datetime, timedelta
from collections import defaultdict

MOVER_DIR = os.environ.get("MOVER_DIR", ".")
WORK_DIR = os.environ.get("WORK_DIR", ".")
LOGID_DIR = os.path.join(WORK_DIR, "results_logid")
RES = os.path.join(WORK_DIR, "results_logid_fix")
os.makedirs(RES, exist_ok=True)


def parse_dt(s):
    s = (s or "").strip()
    if not s:
        return None
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            pass
    return None


cohort = {}
with open(os.path.join(LOGID_DIR, "cohort_logid.csv")) as f:
    for row in csv.DictReader(f):
        cohort[row["MRN"]] = {
            "t0": parse_dt(row["T0"]),
            "disch": parse_dt(row["t_disch"]),
            "an_stop": parse_dt(row["t_anstop"]),
        }

pain = defaultdict(list)
with open(os.path.join(LOGID_DIR, "pain_cohort.csv")) as f:
    for row in csv.DictReader(f):
        mrn = row["mrn"]
        if mrn not in cohort:
            continue
        t = parse_dt(row["t_pain"])
        if t is None:
            continue
        t0 = cohort[mrn]["t0"]
        if not (t0 <= t <= t0 + timedelta(hours=96)):
            continue
        try:
            v = float(row["val"])
        except (ValueError, TypeError):
            continue
        pain[mrn].append((t, v))

NSAID_PAT = re.compile(r"ketorolac|ibuprofen|diclofenac|celecoxib|meloxicam|naproxen|ketoprofen|piroxicam", re.I)
APAP_PAT = re.compile(r"acetaminophen|tylenol|ofirmev", re.I)


def is_analgesic_nonop(row):
    nm = (row.get("MEDICATION_NM") or "") + " " + (row.get("DISPLAY_NAME") or "")
    route = (row.get("MED_ROUTE_NM") or "").lower()
    if "topical" in route or "gel" in nm.lower():
        return False
    return bool(NSAID_PAT.search(nm) or APAP_PAT.search(nm))


analgesics = defaultdict(list)
for fn, is_op in [("cohort_meds_opioid_given.csv", True), ("cohort_meds_nonopioid_given.csv", False)]:
    with open(os.path.join(WORK_DIR, fn)) as f:
        for row in csv.DictReader(f):
            mrn = row["MRN"]
            if mrn not in cohort:
                continue
            if not is_op and not is_analgesic_nonop(row):
                continue
            t = parse_dt(row["MED_ACTION_TIME"])
            if t is None:
                continue
            t0 = cohort[mrn]["t0"]
            if t0 <= t <= t0 + timedelta(hours=96):
                analgesics[mrn].append((t, is_op, (row.get("RECORD_TYPE") or "").strip()))

out_rows = []
for mrn, c in cohort.items():
    t0 = c["t0"]
    win_end = t0 + timedelta(hours=96)
    censor_t = min([x for x in [c["disch"], win_end] if x is not None])
    if censor_t < t0:
        censor_t = t0

    post_meds = []
    for t, is_op, rt in analgesics.get(mrn, []):
        an = c["an_stop"]
        if (t > an) if an is not None else (rt == "POST-OP"):
            post_meds.append(t)
    post_meds = sorted(post_meds)

    prows = sorted(pain.get(mrn, []))
    ge4 = sorted(set(t for t, v in prows if v >= 4))
    event_p = 0
    time_p = (censor_t - t0).total_seconds() / 3600.0
    reason_p = "censored"
    for tp in ge4:
        cands = [t for t in post_meds if tp <= t <= tp + timedelta(hours=4)]
        if cands and cands[0] <= censor_t:
            event_p = 1
            time_p = (cands[0] - t0).total_seconds() / 3600.0
            reason_p = "event"
            break

    ev2, time_s = 0, (censor_t - t0).total_seconds() / 3600.0
    ocands = []
    for t, is_op, rt in analgesics.get(mrn, []):
        if not is_op:
            continue
        if c["an_stop"] is not None:
            if t > c["an_stop"]:
                ocands.append(t)
        elif rt == "POST-OP":
            ocands.append(t)
    if ocands:
        ocands = sorted(ocands)
        if ocands[0] <= censor_t:
            ev2 = 1
            time_s = (ocands[0] - t0).total_seconds() / 3600.0

    out_rows.append({"mrn": mrn, "event_primary": event_p,
                     "time_primary_h": round(time_p, 2), "reason": reason_p,
                     "event_sens": ev2, "time_sens_h": round(time_s, 2)})

with open(os.path.join(RES, "outcomes_fix.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
    w.writeheader()
    w.writerows(out_rows)

old = {}
with open(os.path.join(LOGID_DIR, "outcomes_logid.csv")) as f:
    for row in csv.DictReader(f):
        old[row["mrn"]] = int(row["event_primary"])
new_ev = {r["mrn"]: r["event_primary"] for r in out_rows}
added = [m for m in new_ev if new_ev[m] == 1 and old.get(m, 0) == 0]
lost = [m for m in new_ev if new_ev[m] == 0 and old.get(m, 0) == 1]
print("OLD events: %d | NEW events: %d | added: %d | lost: %d" % (
    sum(old.values()), sum(new_ev.values()), len(added), len(lost)), flush=True)

with open(os.path.join(RES, "added_mrns.txt"), "w") as f:
    f.write("\n".join(sorted(added)))

import numpy as np
ev_times = np.array([r["time_primary_h"] for r in out_rows if r["event_primary"] == 1])
print("event times: n=%d median=%.2f IQR=[%.2f, %.2f] min=%.2f max=%.2f" % (
    len(ev_times), np.median(ev_times),
    np.percentile(ev_times, 25), np.percentile(ev_times, 75),
    ev_times.min(), ev_times.max()), flush=True)
print("WROTE outcomes_fix.csv", flush=True)
