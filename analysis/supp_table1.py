#!/usr/bin/env python3
"""supp_table1.py -- Supplementary Table I: 定义补救镇痛事件的镇痛药清单.

输入 (WORK_DIR; 路径均来自环境变量):
  $WORK_DIR/results_logid/cohort_logid.csv          (600; MRN, T0, t_disch, t_anstop, LOG_ID_x)
  $WORK_DIR/results_logid_fix/outcomes_fix.csv      (主结局 event_primary / time_primary_h)
  $WORK_DIR/results_logid/pain_cohort.csv           (mrn, t_pain, val)
  $WORK_DIR/cohort_meds_opioid_given.csv            (MAR 'Given' 阿片类, 全部计入)
  $WORK_DIR/cohort_meds_nonopioid_given.csv         (MAR 'Given' 非阿片类, 按下述模式计入)
  $MOVER_DIR/patient_information.csv                (AN_START_DATETIME; 24h 排除规则)

输出 ($OUT_DIR):
  $OUT_DIR/supp_table1.csv
    列: idx, trigger_pain_time, trigger_drug_time, drug, route, class, first_event
    去标识化: 不含 MRN。日期列保留 (与 outputs/trigger_drugs_73.csv 对账用)。
    drug = "MEDICATION_NM | DISPLAY_NAME"; route = MED_ROUTE_NM 原样;
    class = opioid (阿片文件) / nonop (非阿片文件)。

定义 (与 data_prep/06_outcomes_fix.py 主结局完全一致):
  - 队列: 600 -> 573 (剔除 T0 > t_anstop) -> 571 (再剔除 T0 早于麻醉开始 >24h 的 2 例)
  - 触发疼痛: [T0, T0+96h] 内首个评分 >=4 的疼痛
  - 候选镇痛药: 术后给药 (t > t_anstop; 缺失时 RECORD_TYPE=='POST-OP),
    且在 [T0, T0+96h] 内; 阿片类全部计入; 非阿片类按 NSAID/对乙酰氨基酚
    名称模式计入 (排除 route 含 topical / 药名含 gel 的外用)
  - 事件: 触发疼痛后 4h 内 (tp <= t <= tp+4h) 有候选药, 事件药物 = 首个候选药
  - 删失: min(出院, T0+96h); 事件须 <= 删失

行含义:
  - first_event=1: 73 个主结局事件的触发药物 (与 outcomes_fix.csv 事件时间逐一对齐核验,
    容差 120s; 药物本身须满足术后 + [T0, T0+96h] 窗内)
  - first_event=0: 触发疼痛同一 4h 窗内、除事件药物之外的其他候选药
    (即"与触发疼痛同时、但不是定义事件的那一剂"药物)。
    说明: 原 outputs/trigger_drugs_73.csv 的 21 行候选中 19 行是窗外其他就诊记录;
    本脚本对候选药加 [T0, T0+96h] 下界过滤, 只保留窗内候选, 故候选行数会变化。
    另发现原表 73 行 first_event=1 中有 5 行同样是窗外其他就诊的给药记录
    (与事件时间相差 1-2 年), 本脚本按 [T0, T0+96h] 过滤后予以纠正
    (见本文件末尾核验输出)。

排序: first_event=1 在前 (按 trigger_drug_time), 随后 first_event=0 (按 trigger_drug_time)。
"""
import csv
import os
import re
from datetime import datetime, timedelta
from collections import defaultdict

WORK_DIR = os.environ.get("WORK_DIR", "/srv/data/mover_analysis")
MOVER_DIR = os.environ.get("MOVER_DIR", "/srv/data/mover_analysis")
OUT_DIR = os.environ.get("OUT_DIR", "/srv/data/mover_analysis/results_v26")
OUT = os.path.join(OUT_DIR, "supp_table1.csv")

WIN_H = 96


def p(s):
    s = (s or "").strip()
    if not s:
        return None
    # 格式集与 build_571_dataset.py 一致 (AN_START_DATETIME 为美式 %m/%d/%Y 格式)
    for f in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d",
              "%m/%d/%y %H:%M", "%m/%d/%Y %H:%M"):
        try:
            return datetime.strptime(s, f)
        except ValueError:
            pass
    return None


NSAID_PAT = re.compile(r"ketorolac|ibuprofen|diclofenac|celecoxib|meloxicam|naproxen|ketoprofen|piroxicam", re.I)
APAP_PAT = re.compile(r"acetaminophen|tylenol|ofirmev", re.I)


def is_analgesic_nonop(row):
    # 与 06_outcomes_fix.py 完全一致: 排除外用/凝胶
    nm = (row.get("MEDICATION_NM") or "") + " " + (row.get("DISPLAY_NAME") or "")
    route = (row.get("MED_ROUTE_NM") or "").lower()
    if "topical" in route or "gel" in nm.lower():
        return False
    return bool(NSAID_PAT.search(nm) or APAP_PAT.search(nm))


def drug_label(row):
    return "%s | %s" % ((row.get("MEDICATION_NM") or "").strip(),
                        (row.get("DISPLAY_NAME") or "").strip())


# ---------- 1. 队列: 600 -> 573 -> 571 ----------
cohort = {}
with open(os.path.join(WORK_DIR, "results_logid/cohort_logid.csv")) as f:
    for r in csv.DictReader(f):
        cohort[r["MRN"]] = {"t0": p(r["T0"]), "disch": p(r["t_disch"]),
                            "an_stop": p(r["t_anstop"]), "logid": r["LOG_ID_x"]}
print("cohort_logid rows: %d" % len(cohort), flush=True)
assert len(cohort) == 600, len(cohort)

kept = {m: c for m, c in cohort.items()
        if not (c["t0"] and c["an_stop"] and c["t0"] > c["an_stop"])}
print("after T0>t_anstop filter (573): %d" % len(kept), flush=True)
assert len(kept) == 573, len(kept)

an_start = {}
with open(os.path.join(MOVER_DIR, "patient_information.csv")) as f:
    for r in csv.DictReader(f):
        if r["LOG_ID"] not in an_start:
            an_start[r["LOG_ID"]] = p(r["AN_START_DATETIME"])
kept2, excl = {}, []
for m, c in kept.items():
    a = an_start.get(c["logid"])
    if c["t0"] is not None and a is not None and (c["t0"] - a).total_seconds() / 3600 < -24:
        excl.append(m)
    else:
        kept2[m] = c
kept = kept2
print("after an_start_h<-24 filter (571): %d; excluded: %s" % (len(kept), excl), flush=True)
assert len(kept) == 571, len(kept)

# ---------- 2. 疼痛与候选镇痛药 ----------
pain = defaultdict(list)
with open(os.path.join(WORK_DIR, "results_logid/pain_cohort.csv")) as f:
    for r in csv.DictReader(f):
        m = r["mrn"]
        if m not in kept:
            continue
        t = p(r["t_pain"])
        if t is None:
            continue
        t0 = kept[m]["t0"]
        if not (t0 <= t <= t0 + timedelta(hours=WIN_H)):
            continue
        try:
            v = float(r["val"])
        except (ValueError, TypeError):
            continue
        pain[m].append((t, v))

meds = defaultdict(list)  # mrn -> [(t, is_op, row_dict)]
for fn, is_op in [("cohort_meds_opioid_given.csv", True),
                  ("cohort_meds_nonopioid_given.csv", False)]:
    with open(os.path.join(WORK_DIR, fn)) as f:
        for row in csv.DictReader(f):
            m = row["MRN"]
            if m not in kept:
                continue
            if not is_op and not is_analgesic_nonop(row):
                continue
            t = p(row["MED_ACTION_TIME"])
            if t is None:
                continue
            t0 = kept[m]["t0"]
            if t0 <= t <= t0 + timedelta(hours=WIN_H):  # [T0, T0+96h] 过滤
                meds[m].append((t, is_op, row))
for m in meds:
    meds[m].sort(key=lambda x: x[0])


def postop_meds(m):
    """术后候选药 (与 06 一致): t > an_stop; an_stop 缺失时 RECORD_TYPE=='POST-OP'."""
    out = []
    an = kept[m]["an_stop"]
    for t, is_op, row in meds.get(m, []):
        if an is not None:
            if t > an:
                out.append((t, is_op, row))
        elif (row.get("RECORD_TYPE") or "").strip() == "POST-OP":
            out.append((t, is_op, row))
    return out


# ---------- 3. 事件核验: 复现 06 逻辑, 与 outcomes_fix 对齐 ----------
fix = {}
with open(os.path.join(WORK_DIR, "results_logid_fix/outcomes_fix.csv")) as f:
    for r in csv.DictReader(f):
        if r["mrn"] in kept:
            fix[r["mrn"]] = (int(r["event_primary"]), float(r["time_primary_h"]))

rows = []          # (trigger_pain_t, drug_t, drug, route, class, first_event)
n_event = 0
for m, c in kept.items():
    t0 = c["t0"]
    censor_t = min([x for x in [c["disch"], t0 + timedelta(hours=WIN_H)] if x is not None])
    pm = postop_meds(m)
    ge4 = sorted(set(t for t, v in pain.get(m, []) if v >= 4))

    ev_flag, ev_h = fix.get(m, (0, 0.0))
    if ev_flag != 1:
        continue
    ev_time = t0 + timedelta(hours=ev_h)

    # 找到定义事件的 (触发疼痛, 触发药物)
    found = None
    for tp in ge4:
        cands = [(t, is_op, row) for (t, is_op, row) in pm if tp <= t <= tp + timedelta(hours=4)]
        if cands and cands[0][0] <= censor_t:
            found = (tp, cands[0])
            break
    assert found is not None, "no trigger found for event %s" % m
    tp, (dt, is_op, row) = found
    assert abs((dt - ev_time).total_seconds()) < 120, \
        "event time mismatch %s: %s vs %s" % (m, dt, ev_time)
    n_event += 1
    rows.append((tp, dt, drug_label(row), (row.get("MED_ROUTE_NM") or "").strip(),
                 "opioid" if is_op else "nonop", 1, m))

    # 窗内候选触发: 触发疼痛同一 4h 窗内、除事件药物之外的候选药
    trig_cands = [(t2, is_op2, row2) for (t2, is_op2, row2) in pm
                  if tp <= t2 <= tp + timedelta(hours=4)]
    assert trig_cands and trig_cands[0][0] == dt  # 首个即事件药物
    for (t2, is_op2, row2) in trig_cands[1:]:
        rows.append((tp, t2, drug_label(row2), (row2.get("MED_ROUTE_NM") or "").strip(),
                     "opioid" if is_op2 else "nonop", 0, m))

print("events reproduced in 571: %d" % n_event, flush=True)
assert n_event == 73, n_event
n_ev_fix = sum(1 for e, _ in fix.values() if e == 1)
assert n_ev_fix == 73, n_ev_fix

# ---------- 4. 输出 (去标识化: 去掉 MRN) ----------
rows.sort(key=lambda r: (1 - r[5], r[1]))  # first_event=1 在前, 各按药物时间
out_rows = []
for i, (tp, dt, drug, route, cls, fe, _m) in enumerate(rows, start=1):
    out_rows.append({"idx": i,
                     "trigger_pain_time": tp.strftime("%Y-%m-%d %H:%M:%S"),
                     "trigger_drug_time": dt.strftime("%Y-%m-%d %H:%M:%S"),
                     "drug": drug, "route": route, "class": cls, "first_event": fe})

os.makedirs(OUT_DIR, exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
    w.writeheader()
    w.writerows(out_rows)

n1 = sum(1 for r in out_rows if r["first_event"] == 1)
n0 = sum(1 for r in out_rows if r["first_event"] == 0)
print("WROTE %s: %d rows (first_event=1: %d, first_event=0: %d)" % (OUT, len(out_rows), n1, n0),
      flush=True)

# 与 outputs/trigger_drugs_73.csv 的药物清单一致性核验 (需本地仓库可读时)
ref_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "outputs",
                        "trigger_drugs_73.csv")
if os.path.exists(ref_path):
    ref = list(csv.DictReader(open(ref_path)))
    ref1 = sorted((r["drug"], r["route"], r["class"]) for r in ref if r["first_event"] == "1")
    new1 = sorted((r["drug"], r["route"], r["class"]) for r in out_rows if r["first_event"] == 1)
    print("reference trigger_drugs_73.csv: %d rows (fe=1: %d, fe=0: %d)"
          % (len(ref), len(ref1), len(ref) - len(ref1)), flush=True)
    print("drug-list match (fe=1, drug|route|class multiset): %s" % (ref1 == new1), flush=True)
    if ref1 != new1:
        s1, s2 = set(ref1) - set(new1), set(new1) - set(ref1)
        print("  only in reference: %s" % sorted(s1)[:5], flush=True)
        print("  only in new: %s" % sorted(s2)[:5], flush=True)
else:
    print("reference trigger_drugs_73.csv not found (skipped list comparison)", flush=True)
