#!/usr/bin/env python3
"""evalue.py -- E-value 计算: 评估未测混杂需要多强才能解释掉主效应.

公式 (VanderWeele & Ding, Ann Intern Med 2017):
  HR > 1 时:  E = HR + sqrt(HR * (HR - 1))
  对 95% CI 下限同样计算 -> E_lower: 把 CI 下限拉到 1 所需的混杂强度。
  若 CI 下限 <= 1 (区间已包含无效值), 按惯例 E_lower = 1
  (即任意弱混杂都可使区间包含 1)。

输入 (硬编码; 来源见 outputs/number_provenance.csv 与 outputs/SAP_results_summary_571.md):
  - 主估计 (校正 Cox, n=571):              HR = 1.93 (95% CI 0.95-3.94)
  - PS 重叠加权 (n=571):                   HR = 2.12 (95% CI 1.06-4.23)
  - Firth 惩罚似然 (n=571):                HR = 2.00 (95% CI 0.95-3.89)

输出 ($OUT_DIR):
  $OUT_DIR/evalue.json  (各估计的 E-value 及其 CI 下限的 E-value)
  同时打印可读结果。

验证: 公式手工验算见本文件末尾 assert (亦在汇报中给出)。
"""
import json
import math
import os

OUT_DIR = os.environ.get("OUT_DIR", "/srv/data/mover_analysis/results_v26")
OUT = os.path.join(OUT_DIR, "evalue.json")

ESTIMATES = [
    {"label": "adjusted_Cox_main", "HR": 1.93, "CI_lo": 0.95, "CI_hi": 3.94,
     "note": "主估计: 校正 Cox, n=571"},
    {"label": "PS_overlap_weighted", "HR": 2.12, "CI_lo": 1.06, "CI_hi": 4.23,
     "note": "PS 重叠加权, n=571"},
    {"label": "Firth_penalized", "HR": 2.00, "CI_lo": 0.95, "CI_hi": 3.89,
     "note": "Firth 惩罚似然, n=571"},
]


def evalue(hr):
    """HR>1 时的 E-value; HR<=1 (含 CI 下限已跨过无效值) 返回 1.0。"""
    if hr <= 1:
        return 1.0
    return hr + math.sqrt(hr * (hr - 1))


def main():
    res = {}
    for e in ESTIMATES:
        ev = evalue(e["HR"])
        ev_lo = evalue(e["CI_lo"])
        res[e["label"]] = {
            "HR": e["HR"], "CI": [e["CI_lo"], e["CI_hi"]],
            "E_value": round(ev, 2),
            "E_value_CI_lower": round(ev_lo, 2),
            "note": e["note"],
        }
        print("%-18s HR=%.2f (%.2f-%.2f): E-value=%.2f, E-value(CI下限)=%.2f  [%s]"
              % (e["label"], e["HR"], e["CI_lo"], e["CI_hi"], ev, ev_lo, e["note"]),
              flush=True)

    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1, ensure_ascii=False)
    print("WROTE %s" % OUT, flush=True)

    # ---- 手工验算 (公式正确性) ----
    # E(1.93) = 1.93 + sqrt(1.93*0.93) = 1.93 + sqrt(1.7949) ≈ 1.93+1.3397 = 3.2697 -> 3.27
    assert abs(evalue(1.93) - 3.2697) < 1e-3, evalue(1.93)
    # E(2.12) = 2.12 + sqrt(2.12*1.12) = 2.12 + sqrt(2.3744) ≈ 2.12+1.5409 = 3.6609 -> 3.66
    assert abs(evalue(2.12) - 3.6609) < 1e-3, evalue(2.12)
    # E(2.00) = 2.00 + sqrt(2.00) ≈ 2.00+1.4142 = 3.4142 -> 3.41
    assert abs(evalue(2.00) - 3.4142) < 1e-3, evalue(2.00)
    # E(1.06) = 1.06 + sqrt(1.06*0.06) = 1.06 + sqrt(0.0636) ≈ 1.06+0.2522 = 1.3122 -> 1.31
    assert abs(evalue(1.06) - 1.3122) < 1e-3, evalue(1.06)
    # CI 下限 <= 1 -> 1.0
    assert evalue(0.95) == 1.0 and evalue(1.0) == 1.0
    print("manual verification: OK (公式手工验算通过)", flush=True)


if __name__ == "__main__":
    main()
