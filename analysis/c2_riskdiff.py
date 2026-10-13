import os
import pandas as pd, numpy as np, json
from lifelines import AalenJohansenFitter
np.random.seed(20261010)

WORK_DIR = os.environ.get("WORK_DIR", ".")
OUT_DIR = os.environ.get("OUT_DIR", ".")
OUT = os.path.join(OUT_DIR, "sap")
os.makedirs(OUT, exist_ok=True)

d = pd.read_csv(os.path.join(OUT_DIR, "analysis_571.csv"))
d["vol"] = np.where(d["vol_le15"]==1, "<=15mL", np.where(d["vol_ge25"]==1, ">=25mL", ">15-<25mL"))

def cif_at(times, status, t):
    # break ties with tiny jitter to avoid duplicate-label reindex error
    jt = times.values.astype(float) + np.random.RandomState(0).rand(len(times)) * 1e-9
    ajf = AalenJohansenFitter()
    ajf.fit(pd.Series(jt, index=times.index), status, event_of_interest=1)
    cd = ajf.cumulative_density_
    vals = cd[cd.index <= t]
    return float(vals.iloc[-1, 0]) if len(vals) else 0.0

res = {}
for g in ["<=15mL", ">15-<25mL", ">=25mL"]:
    sub = d[d["vol"]==g]
    for t in [6, 24]:
        cf = cif_at(sub["fg_time"], sub["fg_status"], t)
        res[f"{g}_{t}h"] = {"n": len(sub), "cif": round(cf, 4)}
        print(f"{g} {t}h: n={len(sub)} CIF={cf:.4f}", flush=True)

def boot_diff(t, B=1000):
    diffs = []
    le = d[d["vol"]=="<=15mL"]; ref = d[d["vol"]==">15-<25mL"]
    for b in range(B):
        ble = le.sample(n=len(le), replace=True)
        bref = ref.sample(n=len(ref), replace=True)
        diffs.append(cif_at(ble["fg_time"], ble["fg_status"], t) - cif_at(bref["fg_time"], bref["fg_status"], t))
    diffs = np.array(diffs)
    return np.percentile(diffs, 2.5), np.percentile(diffs, 97.5)

out = {"cif": res, "seed": 20261010, "B": 1000}
for t in [6, 24]:
    lo, hi = boot_diff(t)
    pt = (cif_at(d[d["vol"]=="<=15mL"]["fg_time"], d[d["vol"]=="<=15mL"]["fg_status"], t)
          - cif_at(d[d["vol"]==">15-<25mL"]["fg_time"], d[d["vol"]==">15-<25mL"]["fg_status"], t))
    nnh = round(1/pt, 1) if lo > 0 else None
    out[f"riskdiff_{t}h"] = {"point": round(pt,4), "lo": round(lo,4), "hi": round(hi,4), "NNH": nnh}
    print(f"Risk diff {t}h: {pt:.4f} (95% CI {lo:.4f}-{hi:.4f}) NNH={nnh}", flush=True)

json.dump(out, open(os.path.join(OUT, "c2_results.json"), "w"), indent=1)
print("saved")
