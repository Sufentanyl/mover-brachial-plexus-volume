#!/usr/bin/env python3
"""05_pain_outcomes.py -- 疼痛评分提取 + LOG_ID 队列结局 (初版).

输入:  $WORK_DIR/results_logid/cohort_logid.csv,
       疼痛测量提取文件 (PAIN_MEAS 环境变量, 或 $MOVER_DIR/pi_pain_meas.txt,
       33M 行, 27 字段, MRN=f3, MEAS=f19, time=f22, value=f23),
       $WORK_DIR/cohort_meds_opioid_given.csv, $WORK_DIR/cohort_meds_nonopioid_given.csv
输出:  $WORK_DIR/results_logid/pain_cohort.csv   (mrn, t_pain, val)
       $WORK_DIR/results_logid/outcomes_logid.csv (初版结局, 600 行; 终版见 06)

疼痛提取: pi_pain_meas.txt 中 MEAS == "UC R PAIN INTENSITY/SCORE" 且 MRN 在队列内的行,
  取 (mrn, t_pain=row[22], val=row[23])。
初版结局定义:
  - 疼痛窗: surgery_date 00:00 ~ +96h; 合格疼痛: val>=4 且 t_pain>=t_end(AN_STOP/OUT_OR/T0)
    且 t_pain>=T0
  - 事件: 首个合格疼痛后 4h 内有镇痛药 (任一), 事件时间=给药时间
  - 删失: min(出院, T0+96h)
  - event_sens/time_sens_h: 术后首个阿片类给药 (次要结局)

原始脚本: rebuild3.py (2026-10-10)
"""
import csv
import os
from datetime import timedelta

import pandas as pd

MOVER_DIR = os.environ.get("MOVER_DIR", ".")
WORK_DIR = os.environ.get("WORK_DIR", ".")
OUT = os.path.join(WORK_DIR, "results_logid")

PAIN_CANDIDATES = [
    os.environ.get("PAIN_MEAS", ""),
    os.path.join(MOVER_DIR, "pi_pain_meas.txt"),
]
pain_path = next((p for p in PAIN_CANDIDATES if p and os.path.exists(p)), None)
if pain_path is None:
    raise SystemExit("pi_pain_meas.txt not found; set PAIN_MEAS. tried: %s" % PAIN_CANDIDATES)
print("pain source:", pain_path, flush=True)

coh = pd.read_csv(os.path.join(OUT, "cohort_logid.csv"), dtype={"MRN": str})
coh["T0"] = pd.to_datetime(coh["T0"])
coh["surg_date"] = pd.to_datetime(coh["surg_date"])
coh["t_disch"] = pd.to_datetime(coh["t_disch"])
coh["t_anstop"] = pd.to_datetime(coh["t_anstop"])
coh["t_outor"] = pd.to_datetime(coh["t_outor"])
# end of surgery: AN_STOP else OUT_OR else T0
coh["t_end"] = coh["t_anstop"].fillna(coh["t_outor"]).fillna(coh["T0"])
mrns = set(coh["MRN"])
pinfo = {r["MRN"]: r for _, r in coh.iterrows()}
print("cohort:", len(coh))

# ---------- pain scores (chunked) ----------
pain_rows = []  # (mrn, time, value)
with open(pain_path) as fh:
    reader = csv.reader(fh)
    for row in reader:
        if len(row) < 24:
            continue
        if row[19] != "UC R PAIN INTENSITY/SCORE":
            continue
        mrn = row[3]
        if mrn not in mrns:
            continue
        try:
            t = pd.to_datetime(row[22])
            v = float(row[23])
        except Exception:
            continue
        pain_rows.append((mrn, t, v))
print("pain records for cohort:", len(pain_rows))
pain = pd.DataFrame(pain_rows, columns=["mrn", "t_pain", "val"])
pain.to_csv(os.path.join(OUT, "pain_cohort.csv"), index=False)

# ---------- analgesics ----------
meds_list = []
for f in ["cohort_meds_opioid_given.csv", "cohort_meds_nonopioid_given.csv"]:
    d = pd.read_csv(os.path.join(WORK_DIR, f), usecols=["MRN", "MED_ACTION_TIME", "DISPLAY_NAME"],
                    dtype={"MRN": str})
    d = d[d["MRN"].isin(mrns)].copy()
    d["t_med"] = pd.to_datetime(d["MED_ACTION_TIME"], errors="coerce")
    d = d[d["t_med"].notna()]
    d["kind"] = "opioid" if "opioid" in f else "nonopioid"
    meds_list.append(d[["MRN", "t_med", "kind", "DISPLAY_NAME"]])
meds = pd.concat(meds_list, ignore_index=True).sort_values(["MRN", "t_med"])
print("analgesic records:", len(meds))
opioid = meds[meds["kind"] == "opioid"].copy()

# ---------- outcome per patient ----------
res = []
for mrn, r in pinfo.items():
    T0 = r["T0"]
    sd = r["surg_date"]
    t_end = r["t_end"]
    w0 = sd.normalize()
    w1 = sd.normalize() + timedelta(hours=96)
    censor = T0 + timedelta(hours=96)
    if pd.notna(r["t_disch"]) and r["t_disch"] > T0:
        censor = min(censor, r["t_disch"])
    pp = pain[(pain["mrn"] == mrn) & (pain["t_pain"] >= w0) & (pain["t_pain"] <= w1)
              & (pain["t_pain"] >= t_end) & (pain["t_pain"] >= T0) & (pain["val"] >= 4)]
    pm = meds[(meds["MRN"] == mrn) & (meds["t_med"] >= w0) & (meds["t_med"] <= w1) & (meds["t_med"] >= T0)]
    event_t = None
    for _, prow in pp.sort_values("t_pain").iterrows():
        cand = pm[(pm["t_med"] > prow["t_pain"]) & (pm["t_med"] <= prow["t_pain"] + timedelta(hours=4))]
        if len(cand):
            event_t = cand["t_med"].min()
            break
    if event_t is not None and event_t <= censor:
        res.append((mrn, 1, (event_t - T0).total_seconds() / 3600, "event"))
    else:
        res.append((mrn, 0, (censor - T0).total_seconds() / 3600, "censored"))
res = pd.DataFrame(res, columns=["mrn", "event_primary", "time_primary_h", "reason"])

sens_rows = []
for mrn, r in pinfo.items():
    T0 = r["T0"]
    sd = r["surg_date"]
    t_end = r["t_end"]
    w0 = sd.normalize()
    w1 = w0 + timedelta(hours=96)
    censor = T0 + timedelta(hours=96)
    if pd.notna(r["t_disch"]) and r["t_disch"] > T0:
        censor = min(censor, r["t_disch"])
    po = opioid[(opioid["MRN"] == mrn) & (opioid["t_med"] >= t_end) & (opioid["t_med"] >= T0)
                & (opioid["t_med"] >= w0) & (opioid["t_med"] <= w1)]
    st = po["t_med"].min() if len(po) else None
    if st is not None and st <= censor:
        sens_rows.append((mrn, 1, (st - T0).total_seconds() / 3600))
    else:
        sens_rows.append((mrn, 0, (censor - T0).total_seconds() / 3600))
sens = pd.DataFrame(sens_rows, columns=["mrn", "event_sens", "time_sens_h"])
out = res.merge(sens, on="mrn")
out.to_csv(os.path.join(OUT, "outcomes_logid.csv"), index=False)
print("\nprimary events:", out["event_primary"].sum(), "/", len(out))
print("sens events:", out["event_sens"].sum())
ev = out[out["event_primary"] == 1]["time_primary_h"]
print("event time median: %.1f IQR [%.1f, %.1f]" % (ev.median(), ev.quantile(.25), ev.quantile(.75)))
print("within 24h:", (ev <= 24).sum(), "/", len(ev))
