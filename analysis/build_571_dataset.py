#!/usr/bin/env python3
"""Build analysis_571.csv (primary analysis dataset, n=571).

Upstream inputs (from WORK_DIR):
  results_logid/cohort_logid.csv          (600; rebuild1.py + patch1.py)
  results_logid_fix/outcomes_fix.csv      (primary outcome; fix_outcomes_final.py)
  results_logid/dex_preend_flag.csv       (dexamethasone flag; dex_extract.py)
  results_logid/pain_cohort.csv           (pain scores; rebuild3.py)
  cohort_meds_opioid_given.csv            (med_pass.py)
  cohort_meds_nonopioid_given.csv         (med_pass.py)
Raw input (from MOVER_DIR):
  patient_information.csv                 (AN_START_DATETIME for the 24h rule)

Cohort rule (Methods inclusion criterion):
  1. drop blocks after end of anesthesia (T0 > t_anstop)            -> 573
  2. drop blocks >24h before anesthesia start of the index surgery
     (an_start_h < -24, where an_start_h = (T0 - AN_START)/3600;
      missing AN_START is kept)                                     -> 571

Other definitions (unchanged from the 573 build):
  - has_pain: time-filtered pain flag: pain score in [t_anstop, T0+96h]
    (or [T0, T0+96h] when t_anstop missing)
  - primary outcome from outcomes_fix.csv
    (trigger pain >=4 at/after t_anstop, analgesic within 4h)
  - time_opioid/event_opioid: first post-(anesthesia-start) opioid
    within min(discharge, 96h)
  - time_t0end/event_t0end: time origin shifted to end of anesthesia

Paths come from environment variables (see config.yaml):
  MOVER_DIR  raw MOVER data
  WORK_DIR   pipeline intermediates
  OUT_DIR    final outputs (analysis_571.csv)

"""
import csv
import os
from datetime import datetime, timedelta
from collections import Counter

MOVER_DIR = os.environ.get("MOVER_DIR", ".")
WORK_DIR = os.environ.get("WORK_DIR", ".")
OUT_DIR = os.environ.get("OUT_DIR", ".")

LOGID_DIR = os.path.join(WORK_DIR, "results_logid")
FIX_DIR = os.path.join(WORK_DIR, "results_logid_fix")


def p(s):
    if not s or not str(s).strip():
        return None
    for f in ('%Y-%m-%d %H:%M:%S', '%m/%d/%y %H:%M', '%m/%d/%Y %H:%M'):
        try:
            return datetime.strptime(str(s).strip(), f)
        except ValueError:
            pass
    return None


def main():
    rows = list(csv.DictReader(open(os.path.join(LOGID_DIR, "cohort_logid.csv"))))
    print("cohort_logid rows:", len(rows))
    assert len(rows) == 600, len(rows)

    # --- rule 1: drop post-anesthesia blocks (T0 > t_anstop) -> 573
    kept = [r for r in rows
            if not (p(r['T0']) and p(r['t_anstop']) and p(r['T0']) > p(r['t_anstop']))]
    print("after T0>t_anstop filter:", len(kept))
    assert len(kept) == 573, len(kept)

    # --- rule 2: drop blocks >24h before anesthesia start (rule-based) -> 571
    an_start = {}
    with open(os.path.join(MOVER_DIR, "patient_information.csv")) as f:
        for r in csv.DictReader(f):
            # first occurrence wins (matches patch1.py drop_duplicates behaviour)
            if r["LOG_ID"] not in an_start:
                an_start[r["LOG_ID"]] = p(r["AN_START_DATETIME"])
    excl, kept2, n_missing = [], [], 0
    for r in kept:
        t0, an = p(r['T0']), an_start.get(r['LOG_ID_x'])
        if t0 is None or an is None:
            n_missing += 1
            kept2.append(r)
            continue
        h = (t0 - an).total_seconds() / 3600.0  # negative = block before anesthesia
        if h < -24:
            excl.append((r['MRN'], round(h / 24, 2)))
        else:
            kept2.append(r)
    print("excluded by an_start_h<-24:", excl)
    print("missing AN_START (kept):", n_missing)
    assert len(kept2) == 571, len(kept2)
    kept = kept2

    oc = {}
    with open(os.path.join(FIX_DIR, "outcomes_fix.csv")) as f:
        for r in csv.DictReader(f):
            oc[r['mrn']] = (float(r['time_primary_h']), int(r['event_primary']))
    dex = {}
    with open(os.path.join(LOGID_DIR, "dex_preend_flag.csv")) as f:
        for r in csv.DictReader(f):
            dex[r['mrn']] = 1 if r['pre_end'] == 'True' else 0

    # time-filtered pain flag: t_pain in [t_anstop, T0+96h]
    # (or [T0, T0+96h] when t_anstop missing)
    tm = {r['MRN']: (p(r['T0']), p(r['t_anstop'])) for r in kept}
    pain_mrns = set()
    with open(os.path.join(LOGID_DIR, "pain_cohort.csv")) as f:
        for r in csv.DictReader(f):
            m = r['mrn']
            if m not in tm:
                continue
            t0, an = tm[m]
            tp = p(r['t_pain'])
            if tp is None or t0 is None:
                continue
            start = an if an is not None else t0
            if start <= tp <= t0 + timedelta(hours=96):
                pain_mrns.add(m)
    print("pain_mrns time-filtered:", len(pain_mrns))

    oprows = list(csv.DictReader(open(os.path.join(WORK_DIR, "cohort_meds_opioid_given.csv"))))
    op_by_mrn = {}
    for r in oprows:
        if r['MAR_ACTION_NM'] != 'Given':
            continue
        t = p(r['MED_ACTION_TIME'])
        if t is None:
            continue
        op_by_mrn.setdefault(r['MRN'], []).append(t)
    for m in op_by_mrn:
        op_by_mrn[m].sort()

    out = []
    for r in kept:
        m = r['MRN']
        t0 = p(r['T0'])
        an = p(r['t_anstop'])
        dd = p(r['t_disch'])
        tpri, epri = oc[m]
        cum60 = float(r['cum60_ml'])
        d = {
            'mrn': m,
            'vol_le15': 1 if r['vol_cum60_cat'] == '<=15mL' else 0,
            'vol_ge25': 1 if r['vol_cum60_cat'] == '>=25mL' else 0,
            'cum60': cum60,
            'max_le15': 1 if r['vol_max_cat'] == '<=15mL' else 0,
            'max_ge25': 1 if r['vol_max_cat'] == '>=25mL' else 0,
            'age': float(r['age']),
            'asa3': '' if r['asa3'].strip() == '' else int(float(r['asa3'])),
            'surg_fracture': 1 if r['surgery_type'] == 'fracture' else 0,
            'surg_other': 1 if r['surgery_type'] == 'other_ue' else 0,
            'surgery_type': r['surgery_type'],
            'shoulder': int(r['shoulder']),
            'anes_general': 1 if r['anes_type'].strip() == 'General' else 0,
            'anes_mac': 1 if 'MAC' in r['anes_type'] else 0,
            'anes_type': r['anes_type'],
            'inpatient': int(r['inpatient']),
            'sex': r['sex'],
            'bmi': '' if r['bmi'].strip() == '' else float(r['bmi']),
            'dex_preend': dex.get(m, ''),
            'has_pain': 1 if m in pain_mrns else 0,
            'time_primary': tpri, 'event_primary': epri,
            'time_24': min(tpri, 24.0), 'event_24': epri if tpri <= 24.0 else 0,
        }
        if an is not None:
            shift = (t0 - an).total_seconds() / 3600.0
            d['time_t0end'] = tpri + shift
            d['event_t0end'] = epri
            d['t0end_ok'] = 1
        else:
            d['time_t0end'] = ''
            d['event_t0end'] = ''
            d['t0end_ok'] = 0
        disch_h = (dd - t0).total_seconds() / 3600.0 if dd else 999
        if epri == 1:
            d['fg_status'] = 1
            d['fg_time'] = tpri
        elif disch_h < 96:
            d['fg_status'] = 2
            d['fg_time'] = disch_h
        else:
            d['fg_status'] = 0
            d['fg_time'] = 96.0
        start = an if an is not None else t0
        censor = min(disch_h, 96.0)
        ev_t = ''
        ev = 0
        for ot in op_by_mrn.get(m, []):
            oh = (ot - t0).total_seconds() / 3600.0
            if ot > start and oh <= censor:
                ev_t = oh
                ev = 1
                break
        d['time_opioid'] = ev_t if ev else censor
        d['event_opioid'] = ev
        out.append(d)

    print("t0end non-positive:", sum(1 for d in out if d['t0end_ok'] == 1 and float(d['time_t0end']) <= 0))
    print("opioid events:", sum(d['event_opioid'] for d in out))
    print("primary events:", sum(d['event_primary'] for d in out))
    print("fg status:", dict(Counter(d['fg_status'] for d in out)))
    print("anes missing:", sum(1 for d in out if d['anes_type'].strip() == ''),
          "| bmi miss:", sum(1 for d in out if d['bmi'] == ''),
          "| asa3 miss:", sum(1 for d in out if d['asa3'] == ''),
          "| dex miss:", sum(1 for d in out if d['dex_preend'] == ''))
    print("has_pain:", sum(d['has_pain'] for d in out))
    print("max_single groups:", Counter((d['max_le15'], d['max_ge25']) for d in out))
    evmax = [d for d in out if d['event_primary'] == 1]
    print("events by max_single:", Counter(
        ('le15' if d['max_le15'] else 'ge25' if d['max_ge25'] else 'ref') for d in evmax))

    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = os.path.join(OUT_DIR, "analysis_571.csv")
    cols = list(out[0].keys())
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(out)
    print("wrote", out_path, len(out))


if __name__ == "__main__":
    main()
