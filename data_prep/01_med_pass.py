#!/usr/bin/env python3
"""01_med_pass.py -- 药物分拣: 从 patient_medications.csv (7.2GB) 分拣三类给药记录.

输入:  $MOVER_DIR/patient_medications.csv, $WORK_DIR/block_mrns.txt
输出:  $WORK_DIR/meds_la_given.csv              (全部 Given 的局麻药, 不限 MRN)
       $WORK_DIR/cohort_meds_opioid_given.csv    (BLOCK 队列 Given 阿片类)
       $WORK_DIR/cohort_meds_nonopioid_given.csv (BLOCK 队列 Given 非阿片镇痛药)
       + 若干分布统计 json (order_class_dist, la_route_dist, opioid/nonopioid_top40,
         adjuvant_*)

分类规则(与下游一致):
  LA_RE     = ropivacaine|bupivacaine|mepivacaine|chloroprocaine|levobupivacaine|prilocaine|lidocaine
  OPIOID_RE = morphine|hydromorphone|fentanyl|oxycodone|hydrocodone|tramadol|codeine|
              oxymorphone|methadone|sufentanil|remifentanil|alfentanil|meperidine|
              pethidine|buprenorphine|nalbuphine|butorphanol|tapentadol|dezocine|diamorphine
  NONOP_RE  = ketorolac|acetaminophen|paracetamol|ibuprofen|celecoxib|diclofenac|
              ketamine|gabapentin|pregabalin
只保留 MAR_ACTION_NM == 'Given' 的行; 输出字段做 strip。

原始脚本: med_pass.py (2026-10-10)
"""
import csv
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict

csv.field_size_limit(sys.maxsize)
MOVER_DIR = os.environ.get("MOVER_DIR", ".")
WORK_DIR = os.environ.get("WORK_DIR", ".")
MED = os.path.join(MOVER_DIR, "patient_medications.csv")
LOG = open(os.path.join(WORK_DIR, "med_pass.log"), "a")


def log(msg):
    msg = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(msg, flush=True)
    LOG.write(msg + "\n")
    LOG.flush()


with open(os.path.join(WORK_DIR, "block_mrns.txt")) as f:
    BLOCK = set(x.strip() for x in f if x.strip())
log(f"block MRNs loaded: {len(BLOCK)}")

LA_RE = re.compile(r"ropivacaine|bupivacaine|mepivacaine|chloroprocaine|levobupivacaine|prilocaine|lidocaine", re.I)
OPIOID_RE = re.compile(r"morphine|hydromorphone|fentanyl|oxycodone|hydrocodone|tramadol|codeine|oxymorphone|methadone|sufentanil|remifentanil|alfentanil|meperidine|pethidine|buprenorphine|nalbuphine|butorphanol|tapentadol|dezocine|diamorphine", re.I)
NONOP_RE = re.compile(r"ketorolac|acetaminophen|paracetamol|ibuprofen|celecoxib|diclofenac|ketamine|gabapentin|pregabalin", re.I)
ADJ = [("dexamethasone", re.compile(r"dexamethasone", re.I)),
       ("clonidine", re.compile(r"clonidine", re.I)),
       ("dexmedetomidine", re.compile(r"dexmedetomidine", re.I)),
       ("epinephrine", re.compile(r"epinephrine|adrenaline", re.I))]

KEEP = ["MRN", "LOG_ID", "MEDICATION_ID", "DISPLAY_NAME", "MEDICATION_NM", "MAR_ACTION_NM",
        "MED_ACTION_TIME", "ADMIN_SIG", "DOSE_UNIT_NM", "MED_ROUTE_NM", "RECORD_TYPE",
        "ORDER_CLASS_NM", "ENC_TYPE_NM", "START_DATE", "END_DATE", "ORDERING_DATE"]

f_la = open(os.path.join(WORK_DIR, "meds_la_given.csv"), "w", newline="")
f_op = open(os.path.join(WORK_DIR, "cohort_meds_opioid_given.csv"), "w", newline="")
f_an = open(os.path.join(WORK_DIR, "cohort_meds_nonopioid_given.csv"), "w", newline="")
w_la = csv.DictWriter(f_la, fieldnames=KEEP)
w_la.writeheader()
w_op = csv.DictWriter(f_op, fieldnames=KEEP)
w_op.writeheader()
w_an = csv.DictWriter(f_an, fieldnames=KEEP)
w_an.writeheader()

order_class = Counter()
la_route = Counter()
opioid_top = Counter()
nonop_top = Counter()
adj_counts = defaultdict(Counter)
adj_given = Counter()
adj_samples = defaultdict(list)
stats = Counter()
t0 = time.time()
n = 0

with open(MED, newline="") as fh:
    r = csv.DictReader(fh)
    for row in r:
        n += 1
        stats["total"] += 1
        mar = (row.get("MAR_ACTION_NM") or "").strip()
        oc = (row.get("ORDER_CLASS_NM") or "").strip()
        order_class[(oc, mar)] += 1
        name = (row.get("DISPLAY_NAME") or "") + " " + (row.get("MEDICATION_NM") or "")
        mrn = (row.get("MRN") or "").strip()
        if mar != "Given":
            if n % 4000000 == 0:
                log(f"rows={n} elapsed={int(time.time()-t0)}s")
            continue
        stats["given"] += 1
        out = {k: (row.get(k) or "").strip() for k in KEEP}
        route = out["MED_ROUTE_NM"]
        if LA_RE.search(name):
            stats["la_given"] += 1
            la_route[route] += 1
            w_la.writerow({k: out[k] for k in KEEP})
        if OPIOID_RE.search(name):
            stats["opioid_given"] += 1
            opioid_top[name.strip()[:80]] += 1
            if mrn in BLOCK:
                w_op.writerow(out)
        if NONOP_RE.search(name):
            nonop_top[name.strip()[:80]] += 1
            if mrn in BLOCK:
                w_an.writerow(out)
        for anm, rx in ADJ:
            if rx.search(name):
                adj_given[anm] += 1
                adj_counts[anm][route] += 1
                if len(adj_samples[anm]) < 5:
                    adj_samples[anm].append(out)
        if n % 2000000 == 0:
            log(f"rows={n} given={stats['given']} la={stats['la_given']} elapsed={int(time.time()-t0)}s")

f_la.close()
f_op.close()
f_an.close()
json.dump({"total_rows": stats["total"], "given_rows": stats["given"],
           "la_given_rows": stats["la_given"], "opioid_given_rows": stats["opioid_given"],
           "elapsed_s": int(time.time() - t0)},
          open(os.path.join(WORK_DIR, "meds_pass_stats.json"), "w"), indent=1)
json.dump([[oc, mar, c] for (oc, mar), c in order_class.most_common(60)],
          open(os.path.join(WORK_DIR, "order_class_dist.json"), "w"))
json.dump(la_route.most_common(40), open(os.path.join(WORK_DIR, "la_route_dist.json"), "w"))
json.dump(opioid_top.most_common(40), open(os.path.join(WORK_DIR, "opioid_top40.json"), "w"))
json.dump(nonop_top.most_common(40), open(os.path.join(WORK_DIR, "nonopioid_top40.json"), "w"))
json.dump({k: dict(v) for k, v in adj_counts.items()},
          open(os.path.join(WORK_DIR, "adjuvant_routes.json"), "w"), indent=1)
json.dump(dict(adj_given), open(os.path.join(WORK_DIR, "adjuvant_counts.json"), "w"), indent=1)
json.dump({k: v for k, v in adj_samples.items()},
          open(os.path.join(WORK_DIR, "adjuvant_samples.json"), "w"), indent=1)
log(f"DONE rows={n} stats={dict(stats)}")
