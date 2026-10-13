"""D1 (revised, 2026-10-10; 571版 2026-10-13): 缺失疼痛评分的选择偏倚 IPW 评估.
人群: 主队列完整病例 561 例中,有术后疼痛评分者 290 例(has_pain 时间过滤定义).
步骤:
  1) 在 561 例上 logistic 估计 P(has_pain=1 | X),
     X = 容量组(vol_le15/vol_ge25)、年龄、ASA3、手术类别、肩部、麻醉方式、patient class
  2) 290 例中权重 w = p_bar / P (稳定化), 1st/99th 百分位截断
  3) 加权 Cox (稳健 SE) vs 未加权 Cox (同 290 例, 同 10 协变量) 对比
"""
import os
import pandas as pd, numpy as np, json
from sklearn.linear_model import LogisticRegression
from lifelines import CoxPHFitter

WORK_DIR = os.environ.get("WORK_DIR", ".")
OUT_DIR = os.environ.get("OUT_DIR", ".")
OUT = os.path.join(OUT_DIR, "sap")
os.makedirs(OUT, exist_ok=True)
d = pd.read_csv(os.path.join(OUT_DIR, "analysis_571.csv"))
d2 = d[d["asa3"].notna()].copy().reset_index(drop=True)  # 571版: 561, 与主模型一致
print("complete cases:", len(d2), "events:", int(d2["event_primary"].sum()))

XCOLS = ["vol_le15", "vol_ge25", "age", "asa3", "surg_fracture", "surg_other",
         "shoulder", "anes_general", "anes_mac", "inpatient"]
X = d2[XCOLS].astype(float)
y = d2["has_pain"].astype(int)
N_PAIN = int(y.sum())
print("P(has_pain):", round(y.mean(), 4), "| n_pain:", N_PAIN)
assert N_PAIN > 0

lr = LogisticRegression(max_iter=1000).fit(X, y)
ps = lr.predict_proba(X)[:, 1]
print("P range: [%.4f, %.4f]  (positivity check)" % (ps.min(), ps.max()))
p_bar = y.mean()
w_all = p_bar / ps
pain_mask = (y == 1).to_numpy()
w = w_all[pain_mask]
lo, hi = np.percentile(w, [1, 99])
w = np.clip(w, lo, hi)
ess = (w.sum() ** 2) / (w ** 2).sum()
print(f"weights({N_PAIN}): mean={w.mean():.3f} sd={w.std():.3f} min={w.min():.3f} "
      f"max={w.max():.3f} median={np.median(w):.3f} | trunc=[{lo:.3f},{hi:.3f}]")
print(f"ESS: {ess:.1f} / {N_PAIN}")

pain = d2[pain_mask].copy().reset_index(drop=True)
pain["w"] = w
print("pain-cc events:", int(pain["event_primary"].sum()))

def fit_cox(dd, wcol=None, label=""):
    cdf = pd.DataFrame({
        "time": dd["time_primary"], "event": dd["event_primary"].astype(int),
        "vol_le15": dd["vol_le15"], "vol_ge25": dd["vol_ge25"],
        "age": dd["age"].astype(float), "asa3": dd["asa3"].astype(float),
        "surg_fracture": dd["surg_fracture"], "surg_other": dd["surg_other"],
        "shoulder": dd["shoulder"], "anes_general": dd["anes_general"],
        "anes_mac": dd["anes_mac"], "inpatient": dd["inpatient"]})
    cph = CoxPHFitter()
    if wcol:
        cdf["w"] = dd[wcol].to_numpy()
        cph.fit(cdf, "time", "event", weights_col="w", robust=True)
    else:
        cph.fit(cdf, "time", "event")
    s = cph.summary
    res = {"label": label, "n": len(dd), "events": int(dd["event_primary"].sum())}
    for v in ["vol_le15", "vol_ge25"]:
        r = s.loc[v]
        se = r["se(coef)"]
        res[v] = {"HR": round(float(np.exp(r["coef"])), 2),
                  "lo": round(float(np.exp(r["coef"] - 1.96 * se)), 2),
                  "hi": round(float(np.exp(r["coef"] + 1.96 * se)), 2),
                  "p": round(float(r["p"]), 4)}
    hr = res["vol_le15"]["HR"]
    print(f"{label}: le15 HR={res['vol_le15']['HR']} "
          f"({res['vol_le15']['lo']}-{res['vol_le15']['hi']}) p={res['vol_le15']['p']}")
    return res

unw = fit_cox(pain, label=f"unweighted({N_PAIN})")
wtd = fit_cox(pain, wcol="w", label=f"IPW-weighted({N_PAIN})")

out = {
    "n_complete": len(d2), "n_pain": int(y.sum()),
    "events_pain": int(pain["event_primary"].sum()),
    "P_has_pain": round(float(p_bar), 4),
    "ps_range": [round(float(ps.min()), 4), round(float(ps.max()), 4)],
    "weight": {"mean": round(float(w.mean()), 3), "sd": round(float(w.std()), 3),
               "min": round(float(w.min()), 3), "max": round(float(w.max()), 3),
               "median": round(float(np.median(w)), 3),
               "trunc_1_99": [round(float(lo), 3), round(float(hi), 3)]},
    "ESS": round(float(ess), 1),
    "X": XCOLS,
    "unweighted": unw, "weighted": wtd,
    "method": ("logistic P(has_pain|X) on {nc} complete cases; stabilized w=p_bar/P "
               "on {np} pain cases, 1st/99th truncated; weighted cause-specific Cox "
               "(10 covariates, robust SE) vs unweighted Cox on same {np}").format(nc=len(d2), np=N_PAIN)
}
with open(f"{OUT}/d1_results_291.json", "w") as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
print("wrote d1_results_291.json")
