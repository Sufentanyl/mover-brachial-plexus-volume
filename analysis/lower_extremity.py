#!/usr/bin/env python3
"""lower_extremity.py -- 下肢阻滞扩展队列 (lower-extremity cohort) 构建与结局计算.

CPT 64447 (股神经阻滞) / 64445 (坐骨神经阻滞), 按 LOG_ID 重建, 终点
Supplementary Figure 1 (按阻滞类型分层的 CIF) 的数据源.

输入 (原始数据, MOVER_DIR):
  EPIC_EMR.tar.gz 内 EPIC_EMR/EMR/patient_coding.csv  (MOVER_TAR 环境变量,
      或 $MOVER_DIR/EPIC_EMR.tar.gz, 或 /srv/data/gdrive/MOVER/EPIC_EMR.tar.gz)
  $MOVER_DIR/patient_medications.csv   (7.2GB, 全量扫描, 仅取队列 MRN)
  $MOVER_DIR/patient_information.csv   (手术记录, LOG_ID 连接)
  疼痛测量提取文件 (PAIN_MEAS 环境变量, 或 $MOVER_DIR/pi_pain_meas.txt,
      或 /tmp/pi_pain_meas.txt; 11.6GB)

输出:
  $WORK_DIR/lower_extremity/cohort_expanded_logid.csv  (245 行; 中间文件,
      含 mrn/log_id/t0/block/event/time_h/reason/disch_comp;
      time_h 为 round(et,4), 精度高于原脚本的 round(et,2), 以便下游
      精确再生 lower.csv)
  $OUT_DIR/lower.csv  (245 行; 去标识作图输入: time_h, status, block)

队列构建 (与 data_prep/00-03 主队列逻辑平行, 差异见下):
  - CPT 筛选: REF_BILL_CODE in {64447, 64445}, MRN 级: 仅 64447 -> femoral,
    仅 64445 -> sciatic, 两者皆有 -> combined (CPT 账单 n=356)
  - 长效局麻给药: MED_ROUTE_NM in {PeriNEURAL, Injection} 且 DISPLAY_NAME
    含 ropivacaine|bupivacaine, 按 LOG_ID 精确连接 patient_information 手术记录;
    一人多 LOG_ID 时取最早给药的 LOG_ID; T0=该 LOG_ID 下最早给药时间
  - 解剖合理性: PRIMARY_PROCEDURE_NM 关键词过滤 (见 plausible(); 与臂丛队列
    的 YES_KEYS/NO_KEYS 思路一致, 按阻滞部位定制), 不吻合者排除 (n=90);
    无同手术神经周围给药记录者排除 (n=21)
  - 要求 T0 后 60 min 内局麻累计容量 >0 (mL); 术后阻滞不排除
    (手稿方法注: 术后阻滞排除仅用于臂丛主队列, 下肢队列无对应麻醉结束时间)

结局 (与 06_outcomes_fix.py 终版主结局同口径):
  - 触发疼痛: T0..T0+96h 内疼痛评分 >=4, 按时间依次尝试; 候选镇痛药
    (阿片类 + 全身性 NSAID/对乙酰氨基酚, 排除外用/凝胶) 须为术后给药
    (t > an_stop; an_stop 缺失时 RECORD_TYPE=='POST-OP');
    tp <= t <= tp+4h (同分钟计入); 首个满足者为事件, 事件时间=给药时间
  - 删失: min(出院, T0+96h)

lower.csv 列定义:
  - time_h:  从 T0 起的小时数, 截断 96h, round(x,4)
  - status:  1=补救镇痛事件; 2=出院未补救 (删失于出院, 出院<T0+96h);
             0=随访满 96h 删失
  - block:   femoral (股神经, 115/54) / femsci (股-坐骨联合, 119/35) /
             sciatic (坐骨神经, 11/4)   [括号内: n/事件数]

原始脚本: /tmp/expanded_logid.py (2026-10-10, VPS). 本脚本仅移植其中的
队列构建 + 结局计算部分 (原脚本第 6-7 节的 Aalen-Johansen 描述性与数据
声明核查为一次性分析, 不再重复; 其所需的作图输入即本脚本输出的 lower.csv).
与原脚本的唯一有意差异: 中间文件 time_h 由 round(et,2) 改为 round(et,4),
以便 export_cif_inputs.py 从中间文件精确再生 lower.csv.

路径用环境变量 (见 ../config.yaml): MOVER_DIR / WORK_DIR / OUT_DIR.
"""
import csv
import io
import os
import re
import sys
import tarfile
import time
from collections import Counter, defaultdict
from datetime import datetime, timedelta

csv.field_size_limit(sys.maxsize)

MOVER_DIR = os.environ.get("MOVER_DIR", "/srv/data/mover_analysis")
WORK_DIR = os.environ.get("WORK_DIR", "/srv/data/mover_analysis")
OUT_DIR = os.environ.get("OUT_DIR", "/srv/data/mover_analysis/results_v26")

LE_DIR = os.path.join(WORK_DIR, "lower_extremity")
os.makedirs(LE_DIR, exist_ok=True)


def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


def parse_dt(s):
    s = (s or "").strip()
    if not s:
        return None
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d",
                "%m/%d/%y %H:%M", "%m/%d/%Y %H:%M", "%m/%d/%y %H:%M:%S"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            pass
    return None


def parse_num(s):
    try:
        return float(str(s).strip())
    except (ValueError, TypeError):
        return None


# ============ 0. CPT 筛选 ============
TAR_CANDIDATES = [
    os.environ.get("MOVER_TAR", ""),
    os.path.join(MOVER_DIR, "EPIC_EMR.tar.gz"),
    "/srv/data/gdrive/MOVER/EPIC_EMR.tar.gz",
]
tar_path = next((p for p in TAR_CANDIDATES if p and os.path.exists(p)), None)
if tar_path is None:
    raise SystemExit("EPIC_EMR.tar.gz not found; set MOVER_TAR. tried: %s" % TAR_CANDIDATES)
log("tar: %s" % tar_path)

FEM, SCI = "64447", "64445"
mrn_codes = defaultdict(set)
with tarfile.open(tar_path) as tf:
    member = tf.getmember("EPIC_EMR/EMR/patient_coding.csv")
    with tf.extractfile(member) as fh:
        txt = io.TextIOWrapper(fh, encoding="utf-8", errors="replace")
        for row in csv.DictReader(txt):
            code = (row.get("REF_BILL_CODE") or "").strip()
            if code in (FEM, SCI):
                mrn_codes[row["MRN"]].add(code)

block_of = {}
for mrn, codes in mrn_codes.items():
    if codes == {FEM}:
        block_of[mrn] = "femoral"
    elif codes == {SCI}:
        block_of[mrn] = "sciatic"
    else:
        block_of[mrn] = "combined"
log("CPT MRNs: %d -> %s" % (len(block_of), dict(Counter(block_of.values()))))

# ============ 1. 用药扫描 (仅候选 MRN) ============
LA_RE = re.compile(r"ropivacaine|bupivacaine|mepivacaine|chloroprocaine|"
                   r"levobupivacaine|prilocaine|lidocaine", re.I)
OPIOID_RE = re.compile(r"morphine|hydromorphone|fentanyl|oxycodone|hydrocodone|"
                       r"tramadol|codeine|oxymorphone|methadone|sufentanil|"
                       r"remifentanil|alfentanil|meperidine|pethidine|buprenorphine|"
                       r"nalbuphine|butorphanol|tapentadol|dezocine|diamorphine", re.I)
NONOP_RE = re.compile(r"ketorolac|acetaminophen|paracetamol|ibuprofen|celecoxib|"
                      r"diclofenac|ketamine|gabapentin|pregabalin", re.I)

la_recs = defaultdict(list)       # mrn -> list of dicts
opioid_recs = defaultdict(list)   # mrn -> list of (t, record_type)
nonop_recs = defaultdict(list)    # mrn -> list of (t, name, route, record_type)
n = 0
t0s = time.time()
meds_path = os.path.join(MOVER_DIR, "patient_medications.csv")
with open(meds_path, newline="") as fh:
    for row in csv.DictReader(fh):
        n += 1
        if n % 7000000 == 0:
            log("meds scan rows=%d elapsed=%ds" % (n, int(time.time() - t0s)))
        mrn = (row.get("MRN") or "").strip()
        if mrn not in block_of:
            continue
        if (row.get("MAR_ACTION_NM") or "").strip() != "Given":
            continue
        name = (row.get("DISPLAY_NAME") or "") + " " + (row.get("MEDICATION_NM") or "")
        t = parse_dt(row.get("MED_ACTION_TIME"))
        if t is None:
            continue
        if LA_RE.search(name):
            la_recs[mrn].append({
                "log_id": (row.get("LOG_ID") or "").strip(),
                "t": t,
                "display": (row.get("DISPLAY_NAME") or "").strip(),
                "dose": parse_num(row.get("ADMIN_SIG")),
                "unit": (row.get("DOSE_UNIT_NM") or "").strip(),
                "route": (row.get("MED_ROUTE_NM") or "").strip(),
                "record_type": (row.get("RECORD_TYPE") or "").strip(),
            })
        if OPIOID_RE.search(name):
            opioid_recs[mrn].append((t, (row.get("RECORD_TYPE") or "").strip()))
        if NONOP_RE.search(name):
            nonop_recs[mrn].append((t, (row.get("MEDICATION_NM") or "") + " " +
                                    (row.get("DISPLAY_NAME") or ""),
                                    (row.get("MED_ROUTE_NM") or ""),
                                    (row.get("RECORD_TYPE") or "").strip()))
log("meds scan done rows=%d | MRNs with LA recs=%d, opioid=%d, nonop=%d"
    % (n, len(la_recs), len(opioid_recs), len(nonop_recs)))

# ============ 2. 疼痛评分扫描 ============
PAIN_CANDIDATES = [
    os.environ.get("PAIN_MEAS", ""),
    os.path.join(MOVER_DIR, "pi_pain_meas.txt"),
    "/tmp/pi_pain_meas.txt",
]
pain_path = next((p for p in PAIN_CANDIDATES if p and os.path.exists(p)), None)
if pain_path is None:
    raise SystemExit("pi_pain_meas.txt not found; set PAIN_MEAS. tried: %s" % PAIN_CANDIDATES)
log("pain source: %s" % pain_path)

pain = defaultdict(list)
np_ = 0
with open(pain_path, newline="") as fh:
    for row in csv.reader(fh):
        if len(row) < 24:
            continue
        if row[19] != "UC R PAIN INTENSITY/SCORE":
            continue
        mrn = row[3]
        if mrn not in block_of:
            continue
        try:
            t = datetime.strptime(row[22].strip(), "%Y-%m-%d %H:%M:%S")
            v = float(row[23])
        except (ValueError, IndexError):
            continue
        pain[mrn].append((t, v))
        np_ += 1
log("pain rows for candidate MRNs: %d, patients with any pain: %d" % (np_, len(pain)))

# ============ 3. 手术信息 (LOG_ID 连接) ============
surg = {}
with open(os.path.join(MOVER_DIR, "patient_information.csv"), newline="") as fh:
    for row in csv.DictReader(fh):
        surg[row["LOG_ID"]] = {
            "mrn": row["MRN"],
            "surgery_date": parse_dt(row.get("SURGERY_DATE")),
            "procedure": (row.get("PRIMARY_PROCEDURE_NM") or "").strip(),
            "disch": parse_dt(row.get("HOSP_DISCH_TIME")),
            "an_stop": parse_dt(row.get("AN_STOP_DATETIME")),
        }
log("surgery records loaded: %d" % len(surg))

# ============ 4. 神经周围过滤 + LOG_ID 队列 ============
PN_ROUTES = {"PeriNEURAL", "Injection"}
PN_NAME = re.compile(r"ropivacaine|bupivacaine", re.I)

UPPER_EXCL = re.compile(r"SHOULDER|ELBOW|WRIST|HAND|HUMERUS|RADIUS|ULNA|FINGER|"
                        r"FISTULA|ARTERIOVENOUS|CARPAL|BRACHIAL|AXILL|CLAVICLE|"
                        r"SCAPULA|THUMB", re.I)
LE_FOOT = r"TOE|ACHILLES|TARSAL|METATARSAL|HALLUX|CALCANEUS|PLANTAR|CHEILECTOMY"
FEM_OK = re.compile(r"HIP|KNEE|FEMOR|FEMUR|PATELLA|THIGH|TIBIA|FIBULA|ANKLE|FOOT|"
                    r"\bLEG\b|LOWER EXTREMITY|ARTHROPLASTY|" + LE_FOOT, re.I)
SCI_OK = re.compile(r"KNEE|TIBIA|FIBULA|ANKLE|FOOT|POPLITEAL|BELOW KNEE|\bLEG\b|"
                    r"LOWER EXTREMITY|" + LE_FOOT, re.I)
SCI_EXCL = re.compile(r"HIP|THIGH|FEMUR", re.I)
COMB_OK = re.compile(r"HIP|KNEE|FEMOR|FEMUR|PATELLA|THIGH|TIBIA|FIBULA|ANKLE|FOOT|"
                     r"POPLITEAL|BELOW KNEE|\bLEG\b|LOWER EXTREMITY|ARTHROPLASTY|"
                     + LE_FOOT, re.I)


def plausible(block, proc):
    p = proc or ""
    if UPPER_EXCL.search(p):
        return False
    if block == "femoral":
        return bool(FEM_OK.search(p))
    if block == "sciatic":
        return bool(SCI_OK.search(p)) and not SCI_EXCL.search(p)
    if block == "combined":
        return bool(COMB_OK.search(p))
    return False


excluded = Counter()
exc_procs = Counter()
cohort = {}   # mrn -> dict
multi_logid = 0
for mrn in block_of:
    blk = block_of[mrn]
    pn = [r for r in la_recs.get(mrn, [])
          if r["route"] in PN_ROUTES and PN_NAME.search(r["display"])]
    by_case = defaultdict(list)
    n_blank = 0
    for r in pn:
        if r["log_id"] and r["log_id"] in surg:
            by_case[r["log_id"]].append(r)
        else:
            n_blank += 1
    if not by_case:
        excluded["no_perineural_same_case"] += 1
        continue
    if len(by_case) > 1:
        multi_logid += 1
    lid = min(by_case, key=lambda k: min(x["t"] for x in by_case[k]))
    recs = sorted(by_case[lid], key=lambda x: x["t"])
    t0 = recs[0]["t"]
    s = surg[lid]
    if not plausible(blk, s["procedure"]):
        excluded["implausible_%s" % blk] += 1
        exc_procs[(blk, s["procedure"][:70])] += 1
        continue
    cum60 = sum(r["dose"] for r in recs
                if r["dose"] is not None and (r["unit"] or "").lower() == "ml"
                and t0 <= r["t"] <= t0 + timedelta(minutes=60))
    if not cum60 or cum60 <= 0:
        excluded["missing_dose_%s" % blk] += 1
        continue
    cohort[mrn] = {"mrn": mrn, "block": blk, "log_id": lid, "t0": t0,
                   "disch": s["disch"], "an_stop": s["an_stop"],
                   "procedure": s["procedure"],
                   "surgery_date": s["surgery_date"],
                   "cum60_ml": round(cum60, 1), "n_pn": len(recs),
                   "n_cases": len(by_case), "n_pn_blank_logid": n_blank}
log("cohort n=%d | block types: %s"
    % (len(cohort), dict(Counter(c["block"] for c in cohort.values()))))
log("excluded: %s" % dict(excluded))
log("patients with >1 perineural-linked case: %d" % multi_logid)
log("top implausible (block, procedure): %s" % exc_procs.most_common(15))
assert len(cohort) == 245, len(cohort)

# ============ 5. 结局 (终版主结局口径) ============
NSAID_PAT = re.compile(r"ketorolac|ibuprofen|diclofenac|celecoxib|meloxicam|"
                       r"naproxen|ketoprofen|piroxicam", re.I)
APAP_PAT = re.compile(r"acetaminophen|tylenol|ofirmev", re.I)


def nonop_ok(nm, route):
    r = (route or "").lower()
    if "topical" in r or "gel" in nm.lower():
        return False
    return bool(NSAID_PAT.search(nm) or APAP_PAT.search(nm))


out_rows = []
for mrn, c in cohort.items():
    t0 = c["t0"]
    win_end = t0 + timedelta(hours=96)
    censor_t = min([x for x in [c["disch"], win_end] if x is not None])
    if censor_t < t0:
        censor_t = t0
    # 镇痛药 (窗内)
    ana = []
    for t, rt in opioid_recs.get(mrn, []):
        if t0 <= t <= win_end:
            ana.append((t, rt))
    for t, nm, route, rt in nonop_recs.get(mrn, []):
        if t0 <= t <= win_end and nonop_ok(nm, route):
            ana.append((t, rt))
    post_meds = []
    for t, rt in ana:
        an = c["an_stop"]
        if (t > an) if an is not None else (rt == "POST-OP"):
            post_meds.append(t)
    post_meds = sorted(post_meds)
    ge4 = sorted(set(t for t, v in pain.get(mrn, []) if v >= 4 and t0 <= t <= win_end))
    event, et = 0, (censor_t - t0).total_seconds() / 3600.0
    reason = "censored"
    if not pain.get(mrn):
        reason = "censored_no_pain_in_window"
    elif not ge4:
        reason = "censored_no_pain_ge4"
    for tp in ge4:
        cands = [t for t in post_meds if tp <= t <= tp + timedelta(hours=4)]
        if cands and cands[0] <= censor_t:
            event, et, reason = 1, (cands[0] - t0).total_seconds() / 3600.0, "event"
            break
    disch_comp = (event == 0 and c["disch"] is not None and c["disch"] < win_end)
    out_rows.append({"mrn": mrn, "block": c["block"], "log_id": c["log_id"],
                     "t0": t0.isoformat(), "procedure": c["procedure"],
                     "cum60_ml": c["cum60_ml"],
                     "event": event, "time_h": round(et, 4), "reason": reason,
                     "disch_comp": disch_comp})
    c.update({"event": event, "time_h": et, "reason": reason,
              "disch_comp": disch_comp})
n_ev = sum(r["event"] for r in out_rows)
log("events: %d/%d" % (n_ev, len(out_rows)))
assert n_ev == 93, n_ev

# ============ 6. 输出 ============
cohort_path = os.path.join(LE_DIR, "cohort_expanded_logid.csv")
with open(cohort_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["mrn", "block", "log_id", "t0", "procedure",
                                      "cum60_ml", "n_pn", "n_cases",
                                      "n_pn_blank_logid", "event", "time_h",
                                      "reason", "disch_comp"],
                       lineterminator="\n")
    w.writeheader()
    for c in cohort.values():
        r = {k: (v.isoformat() if hasattr(v, "isoformat") else v)
             for k, v in c.items() if k in w.fieldnames}
        r["t0"] = c["t0"].isoformat()
        r["time_h"] = round(c["time_h"], 4)
        w.writerow(r)
log("wrote %s" % cohort_path)

lower_path = os.path.join(OUT_DIR, "lower.csv")
os.makedirs(OUT_DIR, exist_ok=True)
with open(lower_path, "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["time_h", "status", "block"])
    for r in out_rows:
        blk = "femsci" if r["block"] == "combined" else r["block"]
        status = 1 if r["event"] == 1 else (2 if r["disch_comp"] else 0)
        w.writerow([r["time_h"], status, blk])
log("wrote %s (rows=%d)" % (lower_path, len(out_rows)))
