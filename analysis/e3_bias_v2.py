import os
"""E3 (revised 2026-10-10): Lin/Psaty/Kronmal 未测混杂偏倚分析 + 主 HR 抽样误差.
每次 MC 迭代: 先从 log HR ~ N(log(2.03), SE^2) 抽样 (SE=0.3612, 由主分析 95% CI 反推),
再叠加偏倚因子校正. 报告总不确定性区间 (抽样误差 + 偏倚参数不确定性).
"""
import numpy as np, json
np.random.seed(20261010)
hr_obs = 2.03
se = 0.3612  # 主分析 95% CI (1.00-4.12) 反推: SE = (log(4.12)-log(1.00))/(2*1.96)
log_hr_mean = np.log(hr_obs)
print("SE check: (log(4.12)-log(1.00))/(2*1.96) = %.4f" % ((np.log(4.12)-np.log(1.00))/(2*1.96)))
N = 10000

def run(p1_lo, p1_hi, p0_lo, p0_hi, label):
    p1 = np.random.uniform(p1_lo, p1_hi, N)
    p0 = np.random.uniform(p0_lo, p0_hi, N)
    hr_u = np.random.uniform(1.2, 3.0, N)
    bf = (p1*(hr_u-1) + 1) / (p0*(hr_u-1) + 1)
    log_hr_draw = np.random.normal(log_hr_mean, se, N)
    hr_adj = np.exp(log_hr_draw - np.log(bf))
    med = float(np.median(hr_adj))
    lo = float(np.percentile(hr_adj, 2.5))
    hi = float(np.percentile(hr_adj, 97.5))
    pb1 = float((hr_adj < 1).mean())
    print(label + ": median %.2f (%.2f-%.2f), P(HR<1)=%.3f" % (med, lo, hi, pb1))
    return {"median": round(med,2), "p2.5": round(lo,2), "p97.5": round(hi,2),
            "prob_below_1": round(pb1,3)}

out = {"observed_HR": hr_obs, "logHR_SE": se, "N": N, "seed": 20261010,
       "method": ("Lin/Psaty/Kronmal external adjustment, Monte Carlo; "
                  "each iteration draws logHR ~ N(log(2.03), 0.3612^2) then applies bias factor; "
                  "reported interval = total uncertainty (sampling + bias-parameter)")}
out["scenario_original"] = run(0.15, 0.50, 0.02, 0.15, "Original priors")
out["scenario_shoulder"] = run(0.098, 0.098, 0.039, 0.039, "Shoulder reference")
json.dump(out, open(os.path.join(os.getenv("OUT_DIR", "."), "e3_results.json"), "w"), indent=1)
print("saved e3_results.json")
