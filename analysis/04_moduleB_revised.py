"""rerun5-04: B1/B2 修订版重跑。

B1（稿中 n=420，39 事件，<=15mL RR 0.43 (0.07-2.75)）：
  人群：anes_type=='Regional' 或 anes_mac==1（n=426，asa3 完整 420）。
  结局：T0 之后出现 Intubation / "LMA  Placed"（procedure events 带时间戳，
        EVENT_TIME > T0，严格大于）。
  模型：修正泊松回归（robust HC1 SE）估计 RR；logistic 估计 OR；
        协变量为主模型去掉麻醉方式（age, asa3, surg_fracture, surg_other,
        shoulder, inpatient）；<=15mL 仅 1 事件，加报 Firth logistic
        （手工实现 Jeffreys 先验惩罚，profile-likelihood CI）。
B2（稿中 n=272，<=15mL 均值比 0.90 (0.61-1.33)）：
  人群：区域/MAC 且 procedure events 有 Anesthesia Start/Stop（起<止）。
  结局：麻醉起止间静脉阿片 IV MME（吗啡 1mg=1；芬太尼 10mcg=1；
        氢吗啡酮 0.15mg=1；哌替啶 7.5mg=1；瑞芬太尼不换算、仅标记）。
  模型：两部分模型——Part1 任用阿片修正泊松；Part2 用药者 log(MME) OLS；
        另报中位数差（分位数回归）；协变量 = B1 协变量去 shoulder（全0）
        + 手术时长 surg_hr。

输入（已固化）：rerun5/inputs/procedure_events_cohort.csv、
  rerun5/inputs/mar_cohort.csv、$OUT_DIR/analysis_571.csv、
  $WORK_DIR/results_logid/cohort_logid.csv（T0 与 mrn-LOG_ID 映射）。
输出：$OUT_DIR/rerun5/moduleB_rerun5.json
"""
import json
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from datetime import datetime

WORK_DIR = os.environ.get("WORK_DIR", ".")
OUT_DIR = os.environ.get("OUT_DIR", ".")
INP = os.path.join(WORK_DIR, "rerun5", "inputs")
OUT = os.path.join(OUT_DIR, "rerun5")
os.makedirs(OUT, exist_ok=True)
SEED = 20261011
np.random.seed(SEED)

# ---------------- 通用 ----------------
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

def rr_table(fit, vars_):
    o = {}
    for v in vars_:
        b, se = fit.params[v], fit.bse[v]
        o[v] = {"est": float(np.exp(b)), "lo": float(np.exp(b - 1.96 * se)),
                "hi": float(np.exp(b + 1.96 * se)), "p": float(fit.pvalues[v])}
    return o

co = pd.read_csv(os.path.join(OUT_DIR, "analysis_571.csv"))
co["asa3"] = pd.to_numeric(co["asa3"], errors="coerce")
ct = pd.read_csv(os.path.join(WORK_DIR, "results_logid", "cohort_logid.csv"), usecols=["MRN", "LOG_ID_x", "T0"])
t0map = dict(zip(ct["MRN"], pd.to_datetime(ct["T0"])))
logid2mrn = dict(zip(ct["LOG_ID_x"].astype(str), ct["MRN"]))

sub = co[(co["anes_type"] == "Regional") | (co["anes_mac"] == 1)].copy()
print("Regional/MAC:", len(sub))

pe = pd.read_csv(f"{INP}/procedure_events_cohort.csv",
                 usecols=["LOG_ID", "EVENT_DISPLAY_NAME", "EVENT_TIME"])
pe["LOG_ID"] = pe["LOG_ID"].astype(str)
pe["mrn"] = pe["LOG_ID"].map(logid2mrn)
pe = pe[pe["mrn"].notna()].copy()
pe["et"] = pd.to_datetime(pe["EVENT_TIME"], errors="coerce")

# ---------------- B1 ----------------
AIR = ["Intubation", "LMA  Placed"]
air = pe[pe["EVENT_DISPLAY_NAME"].isin(AIR) & pe["et"].notna()].copy()
air = air[air["mrn"].isin(set(sub["mrn"]))]
air["t0"] = air["mrn"].map(t0map)
b1_mrns = set(air[air["et"] > air["t0"]]["mrn"].unique())
sub["b1_event"] = sub["mrn"].isin(b1_mrns).astype(int)
print("B1 events:", int(sub["b1_event"].sum()))
print(sub.assign(vgrp=np.where(sub["vol_le15"] == 1, "le15",
       np.where(sub["vol_ge25"] == 1, "ge25", "ref")))
      .groupby("vgrp")["b1_event"].agg(["sum", "count"]).to_string())

COVS_B1 = ["age", "asa3", "surg_fracture", "surg_other", "shoulder", "inpatient"]
d1 = sub.dropna(subset=COVS_B1).copy()
X1 = sm.add_constant(d1[["vol_le15", "vol_ge25"] + COVS_B1], has_constant="add")
y1 = d1["b1_event"]

mod_pois = sm.GLM(y1, X1, family=sm.families.Poisson()).fit(cov_type="HC1")
mod_logit = sm.Logit(y1, X1).fit(disp=0)

res = {"B1_modpoisson": {"n": len(d1), "events": int(y1.sum()),
                        **{k: {"RR": v["est"], "lo": v["lo"], "hi": v["hi"], "p": v["p"]}
                             for k, v in rr_table(mod_pois, ["vol_le15", "vol_ge25"]).items()}},
       "B1_logistic": {"n": len(d1), "events": int(y1.sum()),
                       **{k: {"OR": v["est"], "lo": v["lo"], "hi": v["hi"], "p": v["p"]}
                            for k, v in rr_table(mod_logit, ["vol_le15", "vol_ge25"]).items()}}}

# ---- Firth logistic（Jeffreys 惩罚，手工 Newton-Raphson） ----
def firth_logit(X, y, maxit=200, tol=1e-10):
    Xn = np.asarray(X, dtype=float)
    yn = np.asarray(y, dtype=float)
    n, pp = Xn.shape
    beta = np.zeros(pp)
    for _ in range(maxit):
        eta = Xn @ beta
        mu = 1.0 / (1.0 + np.exp(-eta))
        W = mu * (1 - mu)
        XtWX = (Xn.T * W) @ Xn
        try:
            XtWX_inv = np.linalg.inv(XtWX)
        except np.linalg.LinAlgError:
            XtWX_inv = np.linalg.pinv(XtWX)
        # hat 对角线
        Z = Xn * np.sqrt(W)[:, None]
        h = np.einsum("ij,jk,ik->i", Z, XtWX_inv, Z)
        # Firth 修正 score
        s = Xn.T @ (yn - mu + h * (0.5 - mu))
        step = XtWX_inv @ s
        beta_new = beta + step
        if np.max(np.abs(step)) < tol:
            beta = beta_new
            break
        beta = beta_new
    return beta, XtWX_inv

def penalized_loglik(beta, Xn, yn):
    eta = Xn @ beta
    mu = 1.0 / (1.0 + np.exp(-eta))
    mu = np.clip(mu, 1e-12, 1 - 1e-12)
    ll = np.sum(yn * np.log(mu) + (1 - yn) * np.log(1 - mu))
    W = mu * (1 - mu)
    XtWX = (Xn.T * W) @ Xn
    sign, logdet = np.linalg.slogdet(XtWX)
    return ll + 0.5 * logdet

def firth_profile_ci(Xn, yn, beta_hat, j):
    # 对 beta_j 做 profile penalized likelihood CI（双侧二分）
    ll_max = penalized_loglik(beta_hat, Xn, yn)
    crit = ll_max - 0.5 * 3.8414588  # chisq_1(0.95)/2
    def profile(bj):
        # 固定 beta_j=bj，用 Firth Newton 优化其余参数
        beta = beta_hat.copy()
        beta[j] = bj
        for _ in range(200):
            eta = Xn @ beta
            mu = 1.0 / (1.0 + np.exp(-eta))
            W = mu * (1 - mu)
            XtWX = (Xn.T * W) @ Xn
            try:
                XtWX_inv = np.linalg.inv(XtWX)
            except np.linalg.LinAlgError:
                XtWX_inv = np.linalg.pinv(XtWX)
            Z = Xn * np.sqrt(W)[:, None]
            h = np.einsum("ij,jk,ik->i", Z, XtWX_inv, Z)
            s = Xn.T @ (yn - mu + h * (0.5 - mu))
            s[j] = 0.0
            step = XtWX_inv @ s
            step[j] = 0.0
            beta += step
            if np.max(np.abs(step)) < 1e-10:
                break
        return penalized_loglik(beta, Xn, yn)
    # Wald SE 仅用于定搜索步长
    eta0 = Xn @ beta_hat
    mu0 = 1.0 / (1.0 + np.exp(-eta0))
    W0 = mu0 * (1 - mu0)
    se_wald = float(np.sqrt(np.linalg.pinv((Xn.T * W0) @ Xn)[j, j]))
    out = {}
    for side, sgn in (("lo", -1), ("hi", 1)):
        inner, outer = beta_hat[j], beta_hat[j] + sgn * 10 * se_wald
        # 向外扩展直到 profile(outer) < crit；保持 profile(inner) >= crit
        for _ in range(60):
            if profile(outer) < crit:
                break
            outer += sgn * 10 * se_wald
        for _ in range(80):
            mid = 0.5 * (inner + outer)
            if profile(mid) >= crit:
                inner = mid
            else:
                outer = mid
        out[side] = float(0.5 * (inner + outer))
    return out

Xn1 = np.asarray(X1, dtype=float)
yn1 = np.asarray(y1, dtype=float)
beta_f, _ = firth_logit(Xn1, yn1)
cols1 = list(X1.columns)
firth_out = {}
for v in ["vol_le15", "vol_ge25"]:
    j = cols1.index(v)
    ci = firth_profile_ci(Xn1, yn1, beta_f, j)
    # Wald p（惩罚信息矩阵）
    eta = Xn1 @ beta_f
    mu = 1.0 / (1.0 + np.exp(-eta))
    W = mu * (1 - mu)
    V = np.linalg.pinv((Xn1.T * W) @ Xn1)
    se = float(np.sqrt(V[j, j]))
    z = float(beta_f[j] / se)
    from math import erf, sqrt
    pval = float(2 * (1 - 0.5 * (1 + erf(abs(z) / sqrt(2)))))
    firth_out[v] = {"OR": float(np.exp(beta_f[j])),
                    "lo": float(np.exp(ci["lo"])), "hi": float(np.exp(ci["hi"])),
                    "p": pval}
res["B1_firth"] = {"n": len(d1), "events": int(y1.sum()), **firth_out,
                   "method": "Firth penalized logistic (Jeffreys prior), "
                             "profile penalized-likelihood 95% CI, manual Newton-Raphson"}
print("B1 Firth:", json.dumps(firth_out, indent=1))

# ---------------- B2 ----------------
ev2 = pe[pe["EVENT_DISPLAY_NAME"].isin(["Anesthesia Start", "Anesthesia Stop"])
         & pe["et"].notna()]
t0s = ev2[ev2["EVENT_DISPLAY_NAME"] == "Anesthesia Start"].groupby("mrn")["et"].min()
t1s = ev2[ev2["EVENT_DISPLAY_NAME"] == "Anesthesia Stop"].groupby("mrn")["et"].max()
win = pd.DataFrame({"anstart": t0s, "anstop": t1s}).dropna()
win = win[win["anstop"] > win["anstart"]]
print("anesthesia window available:", len(win))
b2 = sub[sub["mrn"].isin(win.index)].copy()
b2 = b2.merge(win, left_on="mrn", right_index=True, how="inner")
print("B2 n:", len(b2))

mar = pd.read_csv(f"{INP}/mar_cohort.csv",
                  usecols=["LOG_ID", "MRN", "MEDICATION_NM", "MAR_ACTION_NM",
                           "MED_ACTION_TIME", "DOSE_UNIT_NM", "MED_ROUTE_NM", "ADMIN_SIG"])
mar = mar[mar["MAR_ACTION_NM"] == "Given"].copy()
mar["mt"] = pd.to_datetime(mar["MED_ACTION_TIME"], errors="coerce")
opi_pat = ("morphine|fentanyl|hydromorphone|meperidine|oxycodone|remifentanil|"
           "sufentanil|alfentanil")
mar["is_opi"] = mar["MEDICATION_NM"].str.contains(opi_pat, case=False, na=False)
mar["is_iv"] = mar["MED_ROUTE_NM"].str.contains("IntraVENOUS", case=False, na=False)
opi = mar[mar["is_opi"] & mar["is_iv"]].copy()
opi["LOG_ID"] = opi["LOG_ID"].astype(str)
opi["mrn"] = opi["LOG_ID"].map(logid2mrn)
opi = opi[opi["mrn"].isin(set(b2["mrn"]))].copy()
opi = opi.merge(b2[["mrn", "anstart", "anstop"]], on="mrn", how="inner")
opi = opi[(opi["mt"] >= opi["anstart"]) & (opi["mt"] <= opi["anstop"])].copy()
print("IV opioid rows in window:", len(opi))
opi["dose"] = pd.to_numeric(opi["ADMIN_SIG"], errors="coerce")
opi["is_remi"] = opi["MEDICATION_NM"].str.contains("remifentanil", case=False, na=False)
print("remifentanil rows:", int(opi["is_remi"].sum()))

def to_mme(r):
    unit = str(r["DOSE_UNIT_NM"]).lower()
    dose = r["dose"]
    nm = str(r["MEDICATION_NM"]).lower()
    if pd.isna(dose):
        return np.nan
    if "fentanyl" in nm:
        return dose / 10.0 if unit == "mcg" else dose * 1000 / 10.0
    if "hydromorphone" in nm:
        return dose / 0.15 if unit == "mg" else dose / 1000 / 0.15
    if "morphine" in nm:
        return dose if unit == "mg" else dose / 1000
    if "meperidine" in nm or "pethidine" in nm:
        return dose / 7.5 if unit == "mg" else dose / 1000 / 7.5
    return np.nan  # 瑞芬太尼/羟考酮等不换算

opi["mme"] = opi.apply(to_mme, axis=1)
print("mme null (unmapped, excl remi):", int(opi["mme"].isna().sum()))
remi_bypt = opi.groupby("mrn")["is_remi"].max()
pm = opi.groupby("mrn")["mme"].sum().reset_index()
b2 = b2.merge(pm, on="mrn", how="left")
b2["mme"] = b2["mme"].fillna(0)
b2["any_opi"] = (b2["mme"] > 0).astype(int)
b2["remi_any"] = b2["mrn"].map(remi_bypt).fillna(False).astype(bool)
b2["surg_hr"] = (pd.to_datetime(b2["anstop"]) - pd.to_datetime(b2["anstart"])).dt.total_seconds() / 3600
print("B2 any opioid:", int(b2["any_opi"].sum()), "| remi_any:", int(b2["remi_any"].sum()))

assert (b2["shoulder"] == 0).all(), "shoulder not all zero in B2"
COVS_B2 = ["age", "asa3", "surg_fracture", "surg_other", "inpatient", "surg_hr"]
dd = b2.dropna(subset=COVS_B2).copy()
X2 = sm.add_constant(dd[["vol_le15", "vol_ge25"] + COVS_B2], has_constant="add")

m1 = sm.GLM(dd["any_opi"], X2, family=sm.families.Poisson()).fit(cov_type="HC1")
res["B2_any"] = {"n": len(dd), "events": int(dd["any_opi"].sum()),
                 **{k: {"RR": v["est"], "lo": v["lo"], "hi": v["hi"], "p": v["p"]}
                      for k, v in rr_table(m1, ["vol_le15", "vol_ge25"]).items()}}

du = dd[dd["mme"] > 0].copy()
Xu = sm.add_constant(du[["vol_le15", "vol_ge25"] + COVS_B2], has_constant="add")
m2 = sm.OLS(np.log(du["mme"]), Xu).fit(cov_type="HC1")
res["B2_amount"] = {"n": len(du),
                    **{k: {"mean_ratio": v["est"], "lo": v["lo"], "hi": v["hi"], "p": v["p"]}
                         for k, v in rr_table(m2, ["vol_le15", "vol_ge25"]).items()}}

mq = smf.quantreg("mme ~ vol_le15 + vol_ge25 + " + " + ".join(COVS_B2), dd).fit(q=0.5)
res["B2_median"] = {"n": len(dd),
                    **{v: {"median_diff": float(mq.params[v]), "p": float(mq.pvalues[v])}
                         for v in ["vol_le15", "vol_ge25"]}}

res["meta"] = {
    "seed": SEED,
    "B2_covariates": COVS_B2,
    "B2_note": "shoulder all zero in B2 subset, dropped; surg_hr added per SAP v1.1",
    "B1_covariates": COVS_B1,
    "manuscript_target": {
        "B1": {"n": 420, "events": 39, "vol_le15_RR": 0.43, "vol_le15_CI": [0.07, 2.75]},
        "B2": {"n": 272, "vol_le15_mean_ratio": 0.90, "vol_le15_CI": [0.61, 1.33]}},
    "script": os.path.join(OUT, "04_moduleB_revised.py"),
}
with open(f"{OUT}/moduleB_rerun5.json", "w") as f:
    json.dump(res, f, indent=1)
print("wrote moduleB_rerun5.json")
print(json.dumps({k: v for k, v in res.items() if k != "meta"}, indent=1)[:3000])
