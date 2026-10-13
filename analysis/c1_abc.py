"""C1 术后 6h 阿片使用（关键次要终点）：方案 A/B/C。

定义（与 SAP 一致）：
  观察窗：[AN_STOP（麻醉结束；缺失时用 T0 代替），min(T0+6h，出院)]；
  窗口<=0 者排除。端点含闭区间。
  MME：窗内所有途径阿片 -> IV MME（Table S1 预设系数）：
    静脉：吗啡 1mg=1；芬太尼 10mcg=1；氢吗啡酮 0.15mg=1；
          哌替啶 7.5mg=1；瑞芬太尼不换算（仅标记）。
    口服：羟考酮（单方）1mg=0.5；羟考酮复方片按片剂中 mg 数×0.5×片数；
          氢可酮复方片按片剂中 mg 数×1/3×片数；
          口服氢吗啡酮 mg×4/3；口服吗啡 mg×1/3。其余零星药物记 0
          并单独计数报告。
模型（均为结局"窗内任用阿片"的修正 Poisson 回归，稳健标准误 HC1）：
  A（主估计）：不加 offset，log(观察时长，小时) 作协变量。
  B（敏感性）：同 A，但不调整观察时长。
  C（敏感性）：log(观察时长) 作 offset，协变量同旧 C1（去麻醉方式）；
     结果报告为"率比"（rate ratio）。
协变量：
  A/B：容量组 + 主模型 10 协变量（vol_le15, vol_ge25, age, asa3,
     surg_fracture, surg_other, shoulder, anes_general, anes_mac, inpatient）。
  C：主模型去麻醉方式（vol_le15, vol_ge25, age, asa3, surg_fracture,
     surg_other, shoulder, inpatient）。
完整病例以各自协变量集为准。

输入（$WORK_DIR 下）：
  results_v26/analysis_571.csv、results_logid/cohort_logid.csv、
  cohort_meds_opioid_given.csv
输出（$OUT_DIR）：
  moduleC_postop.csv（C2 行保留；旧 C1-A/B/C 六行被替换为本轮结果）、
  c1_abc.json（本轮三套估计值的完整记录）
"""
import json
import os
import re
from collections import defaultdict
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import statsmodels.api as sm

WORK_DIR = os.environ["WORK_DIR"]
OUT_DIR = os.environ["OUT_DIR"]
SEED = 20261013
np.random.seed(SEED)

COVS10 = ["vol_le15", "vol_ge25", "age", "asa3", "surg_fracture",
          "surg_other", "shoulder", "anes_general", "anes_mac", "inpatient"]
COVS8 = ["vol_le15", "vol_ge25", "age", "asa3", "surg_fracture",
         "surg_other", "shoulder", "inpatient"]


def p(s):
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return None
    s = str(s).strip()
    if not s:
        return None
    for f in ("%Y-%m-%d %H:%M:%S", "%m/%d/%y %H:%M", "%m/%d/%Y %H:%M"):
        try:
            return datetime.strptime(s, f)
        except ValueError:
            pass
    return None


# ---------------- MME 换算 ----------------
unmapped = defaultdict(int)


def to_mme_iv(name, dose, unit):
    n = name.lower()
    u = unit.lower()
    if "remifentanil" in n:
        return None  # 不换算
    if "fentanyl" in n:
        mg = dose / 1000.0 if u == "mcg" else dose
        return mg / 0.01
    if "hydromorphone" in n:
        mg = dose if u == "mg" else dose / 1000.0
        return mg / 0.15
    if "morphine" in n:
        mg = dose if u == "mg" else dose / 1000.0
        return mg
    if "meperidine" in n or "pethidine" in n:
        mg = dose if u == "mg" else dose / 1000.0
        return mg / 7.5
    unmapped[("IV", name[:50])] += 1
    return 0.0


def to_mme_oral(name, dose, unit):
    n = name.lower()
    u = unit.lower()
    if "remifentanil" in n:
        return None
    # 复方片：羟考酮 ×0.5；氢可酮 ×1/3（口服 MME ÷ 3）
    if "acetaminophen" in n and u == "tablet":
        m = re.search(r"(\d+)-325", name)
        if m:
            mg_per_tab = float(m.group(1))
            if "hydrocodone" in n:
                return mg_per_tab * (1.0 / 3.0) * dose
            return mg_per_tab * 0.5 * dose
        unmapped[("Oral-combo", name[:50])] += 1
        return 0.0
    if "oxycodone" in n and u == "mg":
        return dose * 0.5
    if "hydromorphone" in n and u == "mg":
        return dose * 4.0 / 3.0
    if "morphine" in n and u == "mg":
        return dose / 3.0
    unmapped[("Oral", name[:50])] += 1
    return 0.0


def to_mme(name, dose, unit, route):
    if route == "IntraVENOUS":
        return to_mme_iv(name, dose, unit)
    if route == "Oral":
        return to_mme_oral(name, dose, unit)
    # 其他途径（鼻喷等）：按静脉换算（C1 为所有途径）
    v = to_mme_iv(name, dose, unit)
    if v is None:
        return None
    return v


# ---------------- 人群与窗口 ----------------
cohort = pd.read_csv(os.path.join(os.getenv("OUT_DIR", WORK_DIR), "analysis_571.csv"))
cohort["asa3"] = pd.to_numeric(cohort["asa3"], errors="coerce")
co = pd.read_csv(f"{WORK_DIR}/results_logid/cohort_logid.csv",
                 usecols=["MRN", "T0", "t_anstop", "t_disch"])
comap = {r["MRN"]: r for _, r in co.iterrows()}

ops = pd.read_csv(f"{WORK_DIR}/cohort_meds_opioid_given.csv")
ops = ops[ops["MAR_ACTION_NM"] == "Given"].copy()
ops["mt"] = pd.to_datetime(ops["MED_ACTION_TIME"], errors="coerce")
ops_bypt = dict(tuple(ops.groupby("MRN")))

rows = []
n_win0 = 0
n_anstop_missing = 0
for _, r in cohort.iterrows():
    m = r["mrn"]
    c = comap[m]
    t0 = p(c["T0"])
    an = p(c["t_anstop"])
    disch = p(c["t_disch"])
    if an is None:
        n_anstop_missing += 1
        start = t0
    else:
        start = an
    end = t0 + timedelta(hours=6)
    if disch is not None:
        end = min(end, disch)
    win_hr = (end - start).total_seconds() / 3600.0
    if win_hr <= 0:
        n_win0 += 1
        continue
    mme = 0.0
    n_used = 0
    remi = False
    g = ops_bypt.get(m)
    if g is not None:
        for _, o in g.iterrows():
            mt = o["mt"]
            if pd.isna(mt) or not (start <= mt <= end):
                continue
            v = to_mme(o["DISPLAY_NAME"], float(o["ADMIN_SIG"]),
                       o["DOSE_UNIT_NM"], o["MED_ROUTE_NM"])
            if v is None:
                remi = True
                continue
            if v > 0:
                n_used += 1
            mme += v
    rows.append({"mrn": m, "win_hr": win_hr, "mme": mme,
                 "any_opi": int(mme > 0), "n_rows": n_used,
                 "remi_any": bool(remi), "anstop_missing": an is None})

d = pd.DataFrame(rows).merge(cohort, on="mrn", how="inner")
print(f"C1 eligible: {len(d)} (window<=0 excluded: {n_win0}, "
      f"AN_STOP missing->T0: {n_anstop_missing})")
print("any opioid users:", int(d["any_opi"].sum()))
print("unmapped drugs:", dict(unmapped))


def est(fit, v):
    b, se = fit.params[v], fit.bse[v]
    return {"estimate": float(np.exp(b)),
            "lo": float(np.exp(b - 1.96 * se)),
            "hi": float(np.exp(b + 1.96 * se)),
            "p": float(fit.pvalues[v])}


def run_model(dd, covs, with_logtime, use_offset):
    """修正 Poisson 回归（HC1）。with_logtime: log(win_hr) 作协变量；
    use_offset: log(win_hr) 作 offset（二者互斥）。"""
    dd = dd.copy()
    X = dd[covs].copy()
    if with_logtime:
        X["log_win_hr"] = np.log(dd["win_hr"])
    X = sm.add_constant(X, has_constant="add")
    if use_offset:
        fit = sm.GLM(dd["any_opi"], X, family=sm.families.Poisson(),
                      offset=np.log(dd["win_hr"])).fit(cov_type="HC1")
    else:
        fit = sm.GLM(dd["any_opi"], X,
                     family=sm.families.Poisson()).fit(cov_type="HC1")
    return fit


schemes = {}
# A（主估计）：10 协变量，log(观察时长) 作协变量
ddA = d.dropna(subset=COVS10).copy()
fitA = run_model(ddA, COVS10, with_logtime=True, use_offset=False)
schemes["A"] = {"n": len(ddA), "users": int(ddA["any_opi"].sum()),
                "covariates": "10 (vol_le15, vol_ge25 + main 8)",
                "time_handling": "log(win_hr) as covariate, no offset",
                "log_win_hr_coef": float(fitA.params["log_win_hr"]),
                "log_win_hr_p": float(fitA.pvalues["log_win_hr"]),
                "le15": est(fitA, "vol_le15"),
                "ge25": est(fitA, "vol_ge25")}
print(f"scheme A complete cases: {len(ddA)}, users: {int(ddA['any_opi'].sum())}")

# B（敏感性）：同 A，不调整观察时长
fitB = run_model(ddA, COVS10, with_logtime=False, use_offset=False)
schemes["B"] = {"n": len(ddA), "users": int(ddA["any_opi"].sum()),
                "covariates": "10 (vol_le15, vol_ge25 + main 8)",
                "time_handling": "no observation-time adjustment",
                "le15": est(fitB, "vol_le15"),
                "ge25": est(fitB, "vol_ge25")}

# C（敏感性）：原 offset 模型，8 协变量（去麻醉方式），报告率比
ddC = d.dropna(subset=COVS8).copy()
fitC = run_model(ddC, COVS8, with_logtime=False, use_offset=True)
schemes["C"] = {"n": len(ddC), "users": int(ddC["any_opi"].sum()),
                "covariates": "8 (vol_le15, vol_ge25 + main 6, anesthesia type excluded)",
                "time_handling": "log(win_hr) as offset",
                "le15": est(fitC, "vol_le15"),
                "ge25": est(fitC, "vol_ge25")}
print(f"scheme C complete cases: {len(ddC)}, users: {int(ddC['any_opi'].sum())}")

for k, s in schemes.items():
    for grp in ("le15", "ge25"):
        e = s[grp]
        print(f"  {k}/{grp}: {e['estimate']:.2f} "
              f"({e['lo']:.2f}-{e['hi']:.2f}), p={e['p']:.4f}")

# ---------------- 写 moduleC_postop.csv ----------------
csv_path = f"{OUT_DIR}/moduleC_postop.csv"
old = pd.read_csv(csv_path) if os.path.exists(csv_path) else pd.DataFrame(
    columns=["analysis", "population", "contrast", "estimate", "lo", "hi",
             "method", "note"])
keep = old[~((old["analysis"] == "C1") &
              old["note"].str.contains("scheme A|sensitivity B|sensitivity C",
                                       na=False))]
new_rows = [
    ("C1", f"n={schemes['A']['n']}, MME>0 {schemes['A']['users']}",
     "<=15mL vs ref", "PRIMARY estimate (scheme A)",
     "Modified Poisson, robust HC1, log(obs_h) as covariate (no offset); 10 covariates",
     schemes["A"]["le15"]),
    ("C1", f"n={schemes['A']['n']}, MME>0 {schemes['A']['users']}",
     ">=25mL vs ref", "PRIMARY estimate (scheme A)",
     "Modified Poisson, robust HC1, log(obs_h) as covariate (no offset); 10 covariates",
     schemes["A"]["ge25"]),
    ("C1", f"n={schemes['B']['n']}, MME>0 {schemes['B']['users']}",
     "<=15mL vs ref", "sensitivity B (no observation-time adjustment)",
     "Modified Poisson, robust HC1, no time adjustment; 10 covariates",
     schemes["B"]["le15"]),
    ("C1", f"n={schemes['B']['n']}, MME>0 {schemes['B']['users']}",
     ">=25mL vs ref", "sensitivity B (no observation-time adjustment)",
     "Modified Poisson, robust HC1, no time adjustment; 10 covariates",
     schemes["B"]["ge25"]),
    ("C1", f"n={schemes['C']['n']}, MME>0 {schemes['C']['users']}",
     "<=15mL vs ref", "sensitivity C (rate ratio, supplementary)",
     "Modified Poisson, robust HC1, offset log(obs_h); 8 covariates",
     schemes["C"]["le15"]),
    ("C1", f"n={schemes['C']['n']}, MME>0 {schemes['C']['users']}",
     ">=25mL vs ref", "sensitivity C (rate ratio, supplementary)",
     "Modified Poisson, robust HC1, offset log(obs_h); 8 covariates",
     schemes["C"]["ge25"]),
]
new_df = pd.DataFrame([{
    "analysis": a, "population": pop, "contrast": ct,
    "estimate": round(e["estimate"], 2), "lo": round(e["lo"], 2),
    "hi": round(e["hi"], 2), "method": m, "note": nt}
    for a, pop, ct, nt, m, e in new_rows])
out = pd.concat([keep, new_df], ignore_index=True)
out.to_csv(csv_path, index=False)

res = {
    "method": ("any-opioid-use within [AN_STOP (T0 if missing), min(T0+6h, discharge)] "
               "analyzed with modified Poisson regression (robust HC1). "
               "A (primary): 10 main covariates + log(obs_h) as covariate, no offset. "
               "B: same as A without obs-time adjustment. "
               "C: log(obs_h) offset, 8 covariates (anesthesia type excluded); "
               "reported as rate ratio (supplementary)."),
    "seed": SEED,
    "n_eligible": len(d),
    "n_window_le0_excluded": n_win0,
    "n_anstop_missing_T0_used": n_anstop_missing,
    "remi_any_n": int(d["remi_any"].sum()),
    "unmapped_drugs": {f"{k[0]}|{k[1]}": v for k, v in unmapped.items()},
    "schemes": schemes,
    "script": "mover-sap/analysis/c1_abc.py",
}
with open(f"{OUT_DIR}/c1_abc.json", "w") as f:
    json.dump(res, f, indent=1)
print(f"wrote {csv_path} and {OUT_DIR}/c1_abc.json")
