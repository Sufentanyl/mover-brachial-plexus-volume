#!/usr/bin/env python3
"""export_cif_inputs.py -- 导出 CIF 作图脚本 (plot_cif_figures.py) 的去标识输入文件.

Figure 3 (主队列按容量分层 CIF) 与 Supplementary Figure 1 (下肢队列按阻滞
类型分层 CIF) 的作图输入; 不含 MRN/LOG_ID/日期, 仅保留作图所需的
(time_h, status, 分层变量).

输入:
  主队列: $OUT_DIR/analysis_571.csv (无则回退 $WORK_DIR/results_v26/analysis_571.csv)
          由 analysis/build_571_dataset.py 生成, 含 fg_status/fg_time 列
  下肢队列: $WORK_DIR/lower_extremity/cohort_expanded_logid.csv
          由 analysis/lower_extremity.py 生成

输出 (去标识, 行序与输入一致):
  $OUT_DIR/primary.csv  (571 行: time_h, status, vol_group)
  $OUT_DIR/lower.csv    (245 行: time_h, status, block)

列定义 (两文件共用):
  - time_h: 从 T0 起的小时数, 截断 96h, round(x,4)
  - status: 1=补救镇痛事件; 2=出院未补救 (删失于出院, 出院<T0+96h);
            0=随访满 96h 删失
  - primary.csv 的 vol_group: le15 (<=15mL) / ref (>15-<25mL) / ge25 (>=25mL),
    取自 analysis_571.csv 的 vol_le15/vol_ge25
  - lower.csv 的 block: femoral (股神经) / femsci (股-坐骨联合) / sciatic (坐骨神经)

主队列 status 口径即 build_571_dataset.py 的 fg_status:
  event_primary=1 -> 1; 否则出院小时数<96 -> 2; 否则 -> 0, time_h=96.

路径用环境变量 (见 ../config.yaml): MOVER_DIR / WORK_DIR / OUT_DIR.
"""
import csv
import os
import sys

MOVER_DIR = os.environ.get("MOVER_DIR", "/srv/data/mover_analysis")
WORK_DIR = os.environ.get("WORK_DIR", "/srv/data/mover_analysis")
OUT_DIR = os.environ.get("OUT_DIR", "/srv/data/mover_analysis/results_v26")


def find_main():
    cands = [
        os.path.join(OUT_DIR, "analysis_571.csv"),
        os.path.join(WORK_DIR, "results_v26", "analysis_571.csv"),
    ]
    hit = next((p for p in cands if os.path.exists(p)), None)
    if hit is None:
        raise SystemExit("analysis_571.csv not found; tried: %s" % cands)
    return hit


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ---- primary.csv: 主队列 ----
    main_path = find_main()
    print("main input:", main_path, flush=True)
    rows = list(csv.DictReader(open(main_path)))
    assert len(rows) == 571, len(rows)
    for c in ("fg_status", "fg_time", "vol_le15", "vol_ge25"):
        if c not in rows[0]:
            raise SystemExit("analysis_571.csv 缺列: %s" % c)

    primary_path = os.path.join(OUT_DIR, "primary.csv")
    with open(primary_path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["time_h", "status", "vol_group"])
        for r in rows:
            vg = ("le15" if r["vol_le15"] == "1"
                  else "ge25" if r["vol_ge25"] == "1" else "ref")
            w.writerow([round(float(r["fg_time"]), 4), r["fg_status"], vg])
    n1 = sum(1 for r in rows if r["fg_status"] == "1")
    print("wrote %s (rows=%d, events=%d)" % (primary_path, len(rows), n1))

    # ---- lower.csv: 下肢队列 ----
    le_path = os.path.join(WORK_DIR, "lower_extremity", "cohort_expanded_logid.csv")
    if not os.path.exists(le_path):
        raise SystemExit("下肢队列中间文件不存在, 请先运行 analysis/lower_extremity.py: %s"
                         % le_path)
    print("lower input:", le_path, flush=True)
    le = list(csv.DictReader(open(le_path)))
    assert len(le) == 245, len(le)

    lower_path = os.path.join(OUT_DIR, "lower.csv")
    with open(lower_path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["time_h", "status", "block"])
        for r in le:
            blk = "femsci" if r["block"] == "combined" else r["block"]
            if r["event"] == "1":
                status = 1
            elif str(r["disch_comp"]).strip() in ("True", "true", "1"):
                status = 2
            else:
                status = 0
            w.writerow([round(float(r["time_h"]), 4), status, blk])
    n1 = sum(1 for r in le if r["event"] == "1")
    print("wrote %s (rows=%d, events=%d)" % (lower_path, len(le), n1))


if __name__ == "__main__":
    main()
