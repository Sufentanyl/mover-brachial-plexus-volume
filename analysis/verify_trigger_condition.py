import os
import csv, re
from datetime import datetime, timedelta
from collections import defaultdict

BASE = os.getenv('WORK_DIR', '.')
def p(s):
    if not s or not str(s).strip(): return None
    for f in ('%Y-%m-%d %H:%M:%S','%m/%d/%y %H:%M','%m/%d/%Y %H:%M','%Y-%m-%d %H:%M','%Y-%m-%d'):
        try: return datetime.strptime(str(s).strip(), f)
        except: pass
    return None

# 573 cohort
rows = list(csv.DictReader(open(f'{BASE}/results_logid/cohort_logid.csv')))
kept = [r for r in rows if not (p(r['T0']) and p(r['t_anstop']) and p(r['T0']) > p(r['t_anstop']))]
cohort = {r['MRN']: {'t0': p(r['T0']), 'an_stop': p(r['t_anstop'])} for r in kept}

# 74 events from outcomes_fix
events = {}
with open(f'{BASE}/results_logid_fix/outcomes_fix.csv') as f:
    for r in csv.DictReader(f):
        if r['mrn'] in cohort and r['event_primary'] == '1':
            events[r['mrn']] = float(r['time_primary_h'])
print('events in 573:', len(events))
assert len(events) == 74

# pain scores
pain = defaultdict(list)
with open(f'{BASE}/results_logid/pain_cohort.csv') as f:
    for r in csv.DictReader(f):
        m = r['mrn']
        if m not in cohort: continue
        t = p(r['t_pain'])
        try: v = float(r['val'])
        except: continue
        if t: pain[m].append((t, v))

# analgesics (same definition as rebuild2)
NSAID_PAT = re.compile(r'ketorolac|ibuprofen|diclofenac|celecoxib|meloxicam|naproxen|ketoprofen|piroxicam', re.I)
APAP_PAT = re.compile(r'acetaminophen|tylenol|ofirmev', re.I)
KET_PAT = re.compile(r'ketamine', re.I)
OPIOID_IN_NAME = re.compile(r'hydrocodone|oxycodone|codeine|tramadol|morphine|fentanyl|hydromorphone|oxymorphone|buprenorphine|tapentadol', re.I)
def is_nonop(row):
    nm = (row.get('MEDICATION_NM') or '') + ' ' + (row.get('DISPLAY_NAME') or '')
    route = (row.get('MED_ROUTE_NM') or '').lower()
    if 'topical' in route or 'gel' in nm.lower(): return False
    if OPIOID_IN_NAME.search(nm): return False
    return bool(NSAID_PAT.search(nm) or APAP_PAT.search(nm) or KET_PAT.search(nm))
analgesics = defaultdict(list)
with open(f'{BASE}/cohort_meds_opioid_given.csv') as f:
    for row in csv.DictReader(f):
        m = row['MRN']
        if m not in cohort: continue
        if (row.get('MAR_ACTION_NM') or '').strip() != 'Given': continue
        t = p(row['MED_ACTION_TIME'])
        if t: analgesics[m].append(t)
with open(f'{BASE}/cohort_meds_nonopioid_given.csv') as f:
    for row in csv.DictReader(f):
        m = row['MRN']
        if m not in cohort: continue
        if (row.get('MAR_ACTION_NM') or '').strip() != 'Given': continue
        if not is_nonop(row): continue
        t = p(row['MED_ACTION_TIME'])
        if t: analgesics[m].append(t)
for m in analgesics: analgesics[m].sort()

# for each event, find the trigger pain (first >=4 with analgesic in 4h) and check condition
affected = []
for m, ev_h in events.items():
    t0, an = cohort[m]['t0'], cohort[m]['an_stop']
    start = an if an is not None else t0
    ev_time = t0 + timedelta(hours=ev_h)
    # find trigger: pain >=4 such that analgesic in [pain, pain+4h] matches ev_time
    found = None
    for tp, v in sorted(pain.get(m, [])):
        if v < 4: continue
        for at in analgesics.get(m, []):
            if tp <= at <= tp + timedelta(hours=4) and abs((at - ev_time).total_seconds()) < 120:
                found = (tp, v, at)
                break
        if found: break
    if found:
        tp, v, at = found
        ok = tp >= start
        if not ok:
            affected.append((m, tp, an, at))
            print(f'  AFFECTED {m}: pain {tp} (val {v}) < an_stop {an}, analgesic {at}')
    else:
        print(f'  NO TRIGGER FOUND {m}: event at {ev_time}')

print(f'affected by new condition: {len(affected)}')
