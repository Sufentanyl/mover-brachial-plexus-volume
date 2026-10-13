#!/usr/bin/env python3
"""03_cohort_logid.py -- LOG_ID 连接重建臂丛阻滞队列 (600 人).

输入:  $WORK_DIR/block_cohort.csv, $WORK_DIR/meds_la_given.csv,
       $MOVER_DIR/patient_information.csv,
       $WORK_DIR/results/outcomes_ue.csv (仅用于新旧队列对账)
输出:  $WORK_DIR/results_logid/cohort_logid.csv (600 行)
       $WORK_DIR/results_logid/reconciliation.json

逻辑 (两步, 原为 rebuild1.py + patch1.py, 此处合并为单脚本直接产出终版):
  1. 取 block_cohort 693 人; meds_la_given 中筛选长效局麻神经周围给药
     (MED_ROUTE_NM in {PeriNEURAL, Injection}, DISPLAY_NAME 含 ropivacaine|bupivacaine,
      剂量/时间非空)。
  2. 按 LOG_ID 左连接 patient_information, 仅保留 MRN 一致的连接。
  3. 解剖合理性 (PRIMARY_PROCEDURE_NM): NO_KEYS 优先排除, 再看 YES_KEYS;
     未知 -> 不通过 (保守)。
  4. 每人取最早的一条通过给药为索引阻滞; T0 = 该给药时间;
     max_single_ml = 最大单次剂量; cum60_ml = T0 后 60min 内累计;
     n_admin = 通过给药总条数; LOG_ID = T0 给药的 LOG_ID (与手术信息同源,
     故 LOG_ID_x == LOG_ID_y)。
  5. 协变量: age = BIRTH_DATE (该列即年龄整数); asa3 按 ASA_RATING 映射
     (Severe/Incapacitating/Moribund/Brain Dead=1; Healthy/Mild=0);
     bmi = 体重(oz->kg)/身高(ft'in"->m)^2; surgery_type 按关键词
     (fracture: FRACTURE|ORIF|DISLOCATION; av_fistula: AV FISTULA|FISTULOGRAM|
      AV GRAFT|ARTERIOVENOUS|THROMBECTOMY|SHUNT; 其余 other_ue);
     shoulder = 术名含 SHOULDER; inpatient = PATIENT_CLASS_GROUP==INPATIENT;
     t_disch/t_anstop/t_outor/surg_date 为规范化时间; 体积分组 <=15 / >15-<25 / >=25。

原始脚本: rebuild1.py + patch1.py (2026-10-10)
"""
import json
import os
import re

import numpy as np
import pandas as pd

MOVER_DIR = os.environ.get("MOVER_DIR", ".")
WORK_DIR = os.environ.get("WORK_DIR", ".")
OUT = os.path.join(WORK_DIR, "results_logid")
os.makedirs(OUT, exist_ok=True)

# ---------- load ----------
bc = pd.read_csv(os.path.join(WORK_DIR, "block_cohort.csv"), dtype=str)
block_mrns = set(bc["MRN"].astype(str))
print("block MRNs (CPT 64415/64416):", len(block_mrns))

meds = pd.read_csv(os.path.join(WORK_DIR, "meds_la_given.csv"), dtype={"MRN": str, "LOG_ID": str})
pn = meds[meds["MED_ROUTE_NM"].isin(["PeriNEURAL", "Injection"])].copy()
la = pn[pn["DISPLAY_NAME"].str.contains("ropivacaine|bupivacaine", case=False, na=False)].copy()
la = la[la["MRN"].isin(block_mrns)].copy()
la["dose_ml"] = pd.to_numeric(la["ADMIN_SIG"], errors="coerce")
la["t_med"] = pd.to_datetime(la["MED_ACTION_TIME"], errors="coerce")
la = la[la["dose_ml"].notna() & la["t_med"].notna()].copy()
print("long-acting perineural records w/ dose+time:", len(la), "patients:", la["MRN"].nunique())

pi = pd.read_csv(os.path.join(MOVER_DIR, "patient_information.csv"), dtype={"MRN": str, "LOG_ID": str})

j = la.merge(pi, on="LOG_ID", how="left", suffixes=("", "_pi"))
j = j[j["MRN"] == j["MRN_pi"]].copy()  # keep only MRN-consistent joins
print("after LOG_ID join w/ MRN consistency:", len(j))

# ---------- anatomical plausibility ----------
YES_KEYS = ["UPPER EXTREMITY", "UPPER EXT", "AV FISTULA", "FISTULOGRAM", "AV GRAFT",
            "ARTERIOVENOUS", "THROMBECTOMY", "HAND", "WRIST", "FOREARM", "ELBOW", "BRACHIAL",
            "FINGER", "THUMB", "CARPAL", "CUBITAL", "ROTATOR CUFF", "SHOULDER", "CLAVICLE",
            "ACROMIO", "HUMERUS", "RADIUS", "ULNA", "PALM", "OLECRANON", "MCP JOINT",
            "SCAPHOID", "GANGLION CYST", "CMC JOINT", "HYPOTHENAR", "SYNDACTYLY",
            "TENOSYNOVECTOMY", "CARPECTOMY", "AMPUTATION, RAY", "PINNING, FINGER",
            "GROIN TO HAND"]
NO_KEYS = ["LOWER EXTREMITY", "KNEE", "ANKLE", "FOOT", "TOE", "FEMUR", "TIBIA",
           "TIBIAL", "FIBULA", "HIP", "KIDNEY", "HEART", "HYSTERECTOMY", "COLECTOMY",
           "CHOLECYSTECTOMY", "APPENDECTOMY", "LAPAROSCOP", "LAPAROTOMY", "HERNIA",
           "NEPHRECTOMY", "PROSTATECTOMY", "CATARACT", "SEPTOPLASTY", "ANORECTAL",
           "PERIANAL", "BLADDER", "URETER", "LUNG", "MITRAL VALVE", "STERNOTOMY", "VATS",
           "THORACOTOMY", "HAMSTRING", "NECK", "ABDOMINAL", "VISCERAL", "COLONOSCOPY",
           "EGD", "ESOPHAGO", "PERITONEAL", "VENA CAVA", "CHEST TUBE", "GASTRECTOMY",
           "ILEOSTOMY", "COLOSTOMY", "HEMODIALYSIS CATHETER", "TUNNELED CATH",
           "CENTRAL VENOUS", "VENOGRAM", "ANGIOGRAM", "THROMBOLYSIS", "EMBOLIZATION",
           "ECHOCARDIOGRAM"]


def plausible(proc):
    u = str(proc).upper()
    if any(k in u for k in NO_KEYS):
        return False
    if any(k in u for k in YES_KEYS):
        return True
    return False  # conservative: unknown -> not plausible


j["plausible"] = j["PRIMARY_PROCEDURE_NM"].apply(plausible)
print("plausible administrations:", j["plausible"].sum(), "/", len(j))

# ---------- per-person first plausible block ----------
jp = j[j["plausible"]].sort_values(["MRN", "t_med"]).copy()
first = jp.groupby("MRN").first().reset_index()
print("patients with >=1 plausible block administration:", len(first))


def exposures(g):
    g = g.sort_values("t_med")
    t0 = g["t_med"].iloc[0]
    maxdose = g["dose_ml"].max()
    cum60 = g[g["t_med"] <= t0 + pd.Timedelta(minutes=60)]["dose_ml"].sum()
    return pd.Series({"T0": t0, "max_single_ml": maxdose, "cum60_ml": cum60,
                      "n_admin": len(g), "LOG_ID": g["LOG_ID"].iloc[0]})


exp = jp.groupby("MRN").apply(exposures).reset_index()

cohort = first[["MRN"]].merge(exp, on="MRN", how="left")
surg = jp.drop_duplicates("MRN").set_index("MRN")[
    ["LOG_ID", "PRIMARY_PROCEDURE_NM", "SURGERY_DATE", "PRIMARY_ANES_TYPE_NM",
     "ASA_RATING", "PATIENT_CLASS_GROUP", "HOSP_ADMSN_TIME", "HOSP_DISCH_TIME",
     "AN_STOP_DATETIME", "OUT_OR_DTTM", "SEX", "BIRTH_DATE", "HEIGHT", "WEIGHT"]]
cohort = cohort.merge(surg, left_on="MRN", right_index=True, how="left")

# ---------- covariates ----------
# 历史列顺序: age, sex, asa, bmi (rebuild1) ... vol_cum60_cat, vol_max_cat, asa3 (patch1 追加)
# 取值以 patch1 为准: age=BIRTH_DATE(该列即年龄整数); asa3 按 ASA_RATING 映射;
# bmi/sex 按 LOG_ID 直连 patient_information。
pimap = pi.drop_duplicates("LOG_ID").set_index("LOG_ID")
logids = cohort["LOG_ID_x"].astype(str)
cohort["age"] = [float(pimap.loc[lg, "BIRTH_DATE"]) if lg in pimap.index else np.nan
                 for lg in logids]
cohort["sex"] = [pimap.loc[lg, "SEX"] if lg in pimap.index else np.nan for lg in logids]


def _asa3(x):
    x = str(x)
    if x in ("Severe Systemic Disease", "Incapacitating Disease", "Moribund", "Brain Dead"):
        return 1
    if x in ("Healthy", "Mild Systemic Disease"):
        return 0
    return np.nan


_asa3num = [_asa3(pimap.loc[lg, "ASA_RATING"]) if lg in pimap.index else np.nan
            for lg in logids]
cohort["asa"] = pd.Series(_asa3num).map({1: "III+", 0: "I-II"})


def parse_bmi(h, w):
    try:
        w = float(w)
        kg = w / 16 / 2.20462
        m = re.match(r"\s*(\d+)\s*'\s*(\d+)", str(h))
        if not m:
            return np.nan
        inch = int(m.group(1)) * 12 + int(m.group(2))
        return kg / (inch * 0.0254) ** 2
    except Exception:
        return np.nan


cohort["bmi"] = [parse_bmi(pimap.loc[lg, "HEIGHT"], pimap.loc[lg, "WEIGHT"])
                 if lg in pimap.index else np.nan for lg in logids]


def surg_type(p):
    u = str(p).upper()
    if any(k in u for k in ["FRACTURE", "ORIF", "DISLOCATION"]):
        return "fracture"
    if any(k in u for k in ["AV FISTULA", "FISTULOGRAM", "AV GRAFT", "ARTERIOVENOUS",
                            "THROMBECTOMY", "SHUNT"]):
        return "av_fistula"
    return "other_ue"


cohort["surgery_type"] = cohort["PRIMARY_PROCEDURE_NM"].apply(surg_type)
cohort["shoulder"] = cohort["PRIMARY_PROCEDURE_NM"].str.upper().str.contains(
    "SHOULDER", na=False).astype(int)
cohort["anes_type"] = cohort["PRIMARY_ANES_TYPE_NM"]
cohort["inpatient"] = (cohort["PATIENT_CLASS_GROUP"].str.upper() == "INPATIENT").astype(int)
cohort["t_disch"] = pd.to_datetime(cohort["HOSP_DISCH_TIME"], errors="coerce")
cohort["t_anstop"] = pd.to_datetime(cohort["AN_STOP_DATETIME"], errors="coerce")
cohort["t_outor"] = pd.to_datetime(cohort["OUT_OR_DTTM"], errors="coerce")
cohort["surg_date"] = pd.to_datetime(cohort["SURGERY_DATE"], errors="coerce")


def volcat(v):
    if v <= 15:
        return "<=15mL"
    if v < 25:
        return ">15-<25mL"
    return ">=25mL"


cohort["vol_cum60_cat"] = cohort["cum60_ml"].apply(volcat)
cohort["vol_max_cat"] = cohort["max_single_ml"].apply(volcat)
cohort["asa3"] = _asa3num  # patch1: 追加在末列

out_path = os.path.join(OUT, "cohort_logid.csv")
cohort.to_csv(out_path, index=False)
print("\nFINAL N:", len(cohort))
print(cohort["vol_cum60_cat"].value_counts().to_string())
print("\nshoulder n:", cohort["shoulder"].sum())
print("inpatient n:", cohort["inpatient"].sum())
print("surgery_type:", cohort["surgery_type"].value_counts().to_dict())
print("age median:", cohort["age"].median(), "missing:", cohort["age"].isna().sum())
print("asa3 mean:", cohort["asa3"].mean(), "missing:", cohort["asa3"].isna().sum())
print("bmi median:", cohort["bmi"].median(), "missing:", cohort["bmi"].isna().sum())

# ---------- reconciliation vs old 542 ----------
old_path = os.path.join(WORK_DIR, "results", "outcomes_ue.csv")
old = pd.read_csv(old_path, dtype=str) if os.path.exists(old_path) else None
if old is None:
    print("跳过新旧对账：无 outcomes_ue.csv")
old = old[(old["ue_plausible"] == "1") & (old["la_dose_corr_ml"].notna())]
old_mrns = set(old["mrn"].astype(str))
new_mrns = set(cohort["MRN"].astype(str))
print("\n=== reconciliation ===")
print("old N:", len(old_mrns), "new N:", len(new_mrns))
print("retained:", len(old_mrns & new_mrns))
print("dropped (in old, not new):", len(old_mrns - new_mrns))
print("added (in new, not old):", len(new_mrns - old_mrns))
old_st = dict(zip(old["mrn"].astype(str), old["surgery_type"]))
new_st = dict(zip(cohort["MRN"].astype(str), cohort["surgery_type"]))
both = old_mrns & new_mrns
agree = sum(1 for m in both if old_st.get(m) == new_st.get(m))
print(f"surgery_type agreement on retained: {agree}/{len(both)} = {agree/len(both)*100:.1f}%")

json.dump({"old_n": len(old_mrns), "new_n": len(new_mrns),
           "retained": len(old_mrns & new_mrns), "dropped": len(old_mrns - new_mrns),
           "added": len(new_mrns - old_mrns), "surgtype_agree": agree,
           "surgtype_n": len(both)},
          open(os.path.join(OUT, "reconciliation.json"), "w"), indent=1)
print("\nsaved", out_path)
