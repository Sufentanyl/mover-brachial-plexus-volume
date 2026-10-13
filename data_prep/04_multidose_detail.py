#!/usr/bin/env python3
"""04_multidose_detail.py -- T0 后 60min 内逐次给药明细 (累计容量计算基础).

输入:  $WORK_DIR/results_logid/cohort_logid.csv, $WORK_DIR/meds_la_given.csv
输出:  $WORK_DIR/results_logid/step3_multidose_detail.csv  (775 行)
       $WORK_DIR/results_logid/step3_group_compare.csv

逻辑:
  - 长效局麻神经周围给药: MED_ROUTE_NM in {PeriNEURAL, Injection},
    DISPLAY_NAME 含 ropivacaine|bupivacaine, 剂量/时间非空。
  - 窗口 [T0, T0+60min]; 列: mrn, drug, dose, route, min_from_T0(分钟,1位小数),
    same_logid (给药 LOG_ID == T0 给药的 LOG_ID, 即 LOG_ID_x)。
  - 附带输出每人的 max_single_all/max_single_60/cum60 分组对照表。

原始脚本: step3d.py (2026-10-10)
"""
import os

import numpy as np
import pandas as pd

MOVER_DIR = os.environ.get("MOVER_DIR", ".")
WORK_DIR = os.environ.get("WORK_DIR", ".")
OUT = os.path.join(WORK_DIR, "results_logid")

coh = pd.read_csv(os.path.join(OUT, "cohort_logid.csv"), dtype={"MRN": str})
coh["T0"] = pd.to_datetime(coh["T0"])
coh = coh.rename(columns={"LOG_ID_x": "LOG_ID"})

meds = pd.read_csv(os.path.join(WORK_DIR, "meds_la_given.csv"), dtype={"MRN": str, "LOG_ID": str})
pn = meds[meds["MED_ROUTE_NM"].isin(["PeriNEURAL", "Injection"])].copy()
la = pn[pn["DISPLAY_NAME"].str.contains("ropivacaine|bupivacaine", case=False, na=False)].copy()
la["dose_ml"] = pd.to_numeric(la["ADMIN_SIG"], errors="coerce")
la["t_med"] = pd.to_datetime(la["MED_ACTION_TIME"], errors="coerce")
la = la[la["dose_ml"].notna() & la["t_med"].notna()].copy()

t0map = dict(zip(coh["MRN"], coh["T0"]))
logmap = dict(zip(coh["MRN"], coh["LOG_ID"].astype(str)))
recs = []
det_rows = []
for mrn, t0 in t0map.items():
    sub = la[la["MRN"] == mrn].sort_values("t_med")
    w60 = sub[(sub["t_med"] >= t0) & (sub["t_med"] <= t0 + pd.Timedelta(minutes=60))]
    recs.append({"MRN": mrn,
                 "max_single_all": sub["dose_ml"].max(),
                 "max_single_60": w60["dose_ml"].max() if len(w60) else np.nan,
                 "cum60": w60["dose_ml"].sum(),
                 "n60": len(w60), "nall": len(sub)})
    for _, a in w60.iterrows():
        dt = (a["t_med"] - t0).total_seconds() / 60
        det_rows.append({"mrn": mrn, "drug": a["DISPLAY_NAME"], "dose": a["dose_ml"],
                         "route": a["MED_ROUTE_NM"], "min_from_T0": round(dt, 1),
                         "same_logid": str(a["LOG_ID"]) == logmap[mrn]})
rc = pd.DataFrame(recs)
coh = coh.merge(rc, on="MRN", how="left")
det = pd.DataFrame(det_rows)


def vc(v):
    if pd.isna(v):
        return "NA"
    if v <= 15:
        return "<=15"
    if v < 25:
        return ">15-<25"
    return ">=25"


coh["g_maxall"] = coh["max_single_all"].apply(vc)
coh["g_max60"] = coh["max_single_60"].apply(vc)
coh["g_cum"] = coh["cum60"].apply(vc)

print("=== A: max_single_ALL vs cum60 ===")
print(pd.crosstab(coh["g_maxall"], coh["g_cum"]).to_string())
print("\n=== B: max_single_60min vs cum60 (fair) ===")
print(pd.crosstab(coh["g_max60"], coh["g_cum"]).to_string())
print("\nchanged (fair):", (coh["g_max60"] != coh["g_cum"]).sum())
print("\n60-min administrations:", len(det))

det.to_csv(os.path.join(OUT, "step3_multidose_detail.csv"), index=False)
coh[["MRN", "max_single_all", "max_single_60", "cum60", "g_maxall", "g_max60", "g_cum"]].to_csv(
    os.path.join(OUT, "step3_group_compare.csv"), index=False)
print("saved.")
