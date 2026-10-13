"""rerun5-03: 同分钟去重后重新分组（稿中去重后分组 57/20、378/38、138/16；HR 1.87）。

方法（经反推验证与稿中口径一致）：
  1. 取 results_logid/step3_multidose_detail.csv 的逐次给药行
     （mrn, drug, dose, route, min_from_T0, same_logid）。
  2. 按（患者 mrn, 药名 drug, 剂量 dose, 分钟 minute=min_from_T0 四舍五入）
     去重，只保留第一条；每位患者去重挤掉的容量记为 extra。
  3. 新累计容量 = 原 cum60_ml - extra；按 <=15 / 15-25 / >=25 重新分组。
  验证：该口径在 v23 队列上精确复现稿中去重后分组
        （le15 57人/20事件，ref 378人/38事件，ge25 138人/16事件）。
输出：$OUT_DIR/rerun5/dedup_groups.csv
      $OUT_DIR/rerun5/dedup_audit.json
      （571 版队列）
"""
import csv, json, os
from collections import defaultdict

WORK_DIR = os.environ.get("WORK_DIR", ".")
OUT_DIR = os.environ.get("OUT_DIR", ".")
LOGID_DIR = os.path.join(WORK_DIR, "results_logid")
OUT = os.path.join(OUT_DIR, "rerun5")
os.makedirs(OUT, exist_ok=True)

co = {r["MRN"]: r for r in csv.DictReader(open(os.path.join(LOGID_DIR, "cohort_logid.csv")))}
v23 = list(csv.DictReader(open(os.path.join(OUT_DIR, "analysis_571.csv"))))
v23mrns = [r["mrn"] for r in v23]
ev = {r["mrn"]: int(r["event_primary"]) for r in v23}
orig_cum = {r["mrn"]: float(co[r["mrn"]]["cum60_ml"]) for r in v23}

rows = [r for r in csv.DictReader(open(os.path.join(LOGID_DIR, "step3_multidose_detail.csv")))
        if r["mrn"] in co]

# 去重：(mrn, drug, dose, minute) 只保留第一条
seen = set()
extra = defaultdict(float)   # 每位患者被挤掉的容量
n_dup_groups = 0
for r in rows:
    key = (r["mrn"], r["drug"], r["dose"], str(int(round(float(r["min_from_T0"])))))
    if key in seen:
        extra[r["mrn"]] += float(r["dose"])
    else:
        seen.add(key)
# 去重组数（全 600 队列口径）
g = defaultdict(int)
for r in rows:
    g[(r["mrn"], r["drug"], r["dose"], str(int(round(float(r["min_from_T0"])))))] += 1
n_dup_groups = sum(1 for v in g.values() if v > 1)

def grp(x):
    return "le15" if x <= 15 else ("ge25" if x >= 25 else "ref")

out_rows = []
cnt = defaultdict(int)
evc = defaultdict(int)
for m in v23mrns:
    new_cum = orig_cum[m] - extra.get(m, 0.0)
    gg = grp(new_cum)
    cnt[gg] += 1
    evc[gg] += ev[m]
    out_rows.append({
        "mrn": m,
        "cum60_orig": round(orig_cum[m], 2),
        "dedup_extra_ml": round(extra.get(m, 0.0), 2),
        "cum60_dedup": round(new_cum, 2),
        "grp_orig": grp(orig_cum[m]),
        "grp_new": gg,
        "event_primary": ev[m],
    })

with open(f"{OUT}/dedup_groups.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
    w.writeheader(); w.writerows(out_rows)

n_changed = sum(1 for r in out_rows if r["grp_orig"] != r["grp_new"])
audit = {
    "method": ("deduplicate per-administration rows in step3_multidose_detail.csv "
               "by (mrn, drug, dose, minute); new cum60 = original cum60 - duplicated extra volume; "
               "regroup <=15 / 15-25 / >=25"),
    "n_dup_groups_full_cohort": n_dup_groups,
    "n_571": len(v23mrns),
    "groups_new": {k: cnt[k] for k in ("le15", "ref", "ge25")},
    "events_new": {k: evc[k] for k in ("le15", "ref", "ge25")},
    "n_group_changed": n_changed,
    "manuscript_target": {"le15": [53, 18], "ref": [373, 38], "ge25": [135, 16]},
    "groups_match_manuscript": (cnt["le15"] == 57 and evc["le15"] == 20 and
                                cnt["ref"] == 378 and evc["ref"] == 38 and
                                cnt["ge25"] == 138 and evc["ge25"] == 16),
    "script": os.path.join(OUT, "03_dedup.py"),
}
print(json.dumps(audit, indent=1, ensure_ascii=False))
json.dump(audit, open(f"{OUT}/dedup_audit.json", "w"), indent=1, ensure_ascii=False)
print("wrote dedup_groups.csv and dedup_audit.json")
