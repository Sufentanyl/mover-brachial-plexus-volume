#!/usr/bin/env python3
"""07_dex_extract.py -- 地塞米松给药提取 (600 队列).

输入:  $WORK_DIR/results_logid/cohort_logid.csv, $MOVER_DIR/patient_medications.csv
输出:  $WORK_DIR/results_logid/dex_cohort.csv      (mrn, logid, drug, route, t_med)
       $WORK_DIR/results_logid/dex_preend_flag.csv (mrn, pre_end: 每人是否有
              t_med <= AN_STOP 的地塞米松给药; pre_end 取各条 max)

逻辑: DISPLAY_NAME 含 dexamethasone (不分大小写), MRN 在 600 队列,
  MAR_ACTION_NM=='Given'。pre_end 判定用 AN_STOP_DATETIME (非 t_anstop 规范化列,
  二者一致); t_med 或 AN_STOP 缺失 -> False。
注: 队列中无神经周围地塞米松 (perineural dex 记录数为 0)。

原始脚本: dex_extract.py (2026-10-10)
"""
import os

import pandas as pd

MOVER_DIR = os.environ.get("MOVER_DIR", ".")
WORK_DIR = os.environ.get("WORK_DIR", ".")
OUT = os.path.join(WORK_DIR, "results_logid")

coh = pd.read_csv(os.path.join(OUT, "cohort_logid.csv"), dtype={"MRN": str})
mrns = set(coh["MRN"])
anstop = dict(zip(coh["MRN"], pd.to_datetime(coh["AN_STOP_DATETIME"], errors="coerce")))

rows = []
chunks = pd.read_csv(os.path.join(MOVER_DIR, "patient_medications.csv"), dtype=str, chunksize=500000,
                     usecols=["LOG_ID", "MRN", "DISPLAY_NAME", "MEDICATION_NM", "MAR_ACTION_NM",
                              "MED_ACTION_TIME", "MED_ROUTE_NM", "DOSE_UNIT_NM", "ADMIN_SIG"])
n = 0
for ch in chunks:
    n += 1
    f = ch[ch["DISPLAY_NAME"].str.contains("dexamethasone", case=False, na=False)]
    f = f[f["MRN"].isin(mrns) & (f["MAR_ACTION_NM"] == "Given")]
    for _, r in f.iterrows():
        rows.append({"mrn": r["MRN"], "logid": r["LOG_ID"], "drug": r["DISPLAY_NAME"],
                     "route": r["MED_ROUTE_NM"], "t_med": r["MED_ACTION_TIME"]})
    if n % 10 == 0:
        print("chunks:", n, "found:", len(rows), flush=True)
dex = pd.DataFrame(rows)
dex["t_med"] = pd.to_datetime(dex["t_med"], errors="coerce")
dex.to_csv(os.path.join(OUT, "dex_cohort.csv"), index=False)
print("\ndex records:", len(dex))
print("route:", dex["route"].value_counts().to_dict())
print("patients:", dex["mrn"].nunique())
pn = dex[dex["route"].str.contains("PeriNEURAL", na=False)]
print("perineural dex records:", len(pn))


def before_end(r):
    ae = anstop.get(r["mrn"])
    tm = r["t_med"]
    if pd.isna(ae) or pd.isna(tm):
        return False
    return tm <= ae


dex["pre_end"] = dex.apply(before_end, axis=1)
iv = dex[~dex["route"].str.contains("PeriNEURAL", na=False)]
print("IV pre-end patients:", iv[iv["pre_end"]]["mrn"].nunique())
print("IV any-time patients:", iv["mrn"].nunique())
dex[["mrn", "pre_end"]].groupby("mrn")["pre_end"].max().to_csv(
    os.path.join(OUT, "dex_preend_flag.csv"))
print("saved.")
