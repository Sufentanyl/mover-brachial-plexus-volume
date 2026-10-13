# -*- coding: utf-8 -*-
"""Figures 2, 3 and S1 of the MOVER brachial plexus block manuscript.

Input files contain no identifiers or dates:
  primary.csv  columns time_h, status, vol_group   (primary cohort, 571 rows)
  lower.csv    columns time_h, status, block       (lower-extremity cohort, 245 rows)

  time_h     hours from the first perineural injection (T0) to the first event or censoring
  status     1 = in-hospital rescue analgesia; 2 = discharge without rescue analgesia;
             0 = censored (death, or still in hospital at 96 h)
  vol_group  le15 / ref / ge25
  block      femoral / femsci / sciatic

Usage:
  python plot_cif_figures.py primary.csv lower.csv out_dir [median_discharge_h]

Writes Figure2, Figure3 and FigureS1 as .png (preview) and .tif (RGB, LZW) at the width
and resolution set by W and DPI below (Minerva Anestesiologica: double column 14.9 cm,
600 dpi for graphs); all text 8 pt. Line styles and markers differ between groups so
that the figures remain legible in black and white. Prints the cumulative incidences,
numbers at risk and log-rank P values so they can be checked against the manuscript.
"""
import csv
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import blended_transform_factory
from PIL import Image

TMAX = 96.0
XTICKS = [0, 6, 12, 24, 48, 72, 96]
W = 5.87          # figure width, inches (14.9 cm)
DPI = 600
L = 1.3           # left edge of the plotting area, inches
AXW = W - L - 0.2  # width of the plotting area, inches

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "legend.frameon": False,
})

VOL = [("le15", "≤15 mL", "#D55E00"),
       ("ref", ">15 to <25 mL", "#0072B2"),
       ("ge25", "≥25 mL", "#009E73")]
BLOCK = [("femoral", "Femoral", "#0072B2"),
         ("femsci", "Femoral–sciatic", "#E69F00"),
         ("sciatic", "Sciatic", "#009E73")]
LINESTYLES = ["-", (0, (5, 2)), (0, (1.2, 1.2))]
MARKERS = ["o", "s", "^"]


def read(path, gcol):
    t, s, g = [], [], []
    with open(path, newline="", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            t.append(float(r["time_h"]))
            s.append(int(float(r["status"])))
            g.append(r[gcol].strip())
    t, s, g = np.array(t), np.array(s), np.array(g)
    late = t > TMAX
    t[late], s[late] = TMAX, 0
    return t, s, g


def aj(t, s, cause):
    """Aalen-Johansen cumulative incidence of `cause`; the other non-zero status competes."""
    xs, ys, surv, cif = [0.0], [0.0], 1.0, 0.0
    for tj in np.unique(t[s > 0]):
        n = np.sum(t >= tj)
        at = t == tj
        cif += surv * np.sum(at & (s == cause)) / n
        surv *= 1 - np.sum(at & (s > 0)) / n
        xs.append(tj)
        ys.append(cif)
    xs.append(max(t.max(), xs[-1]))
    ys.append(cif)
    return np.array(xs), np.array(ys)


def cif_at(xs, ys, x):
    return ys[np.searchsorted(xs, x, side="right") - 1]


def at_risk(t, x):
    return int(np.sum(t >= x))


def chi2_sf(x, df):
    if df == 1:
        return math.erfc(math.sqrt(x / 2))
    if df == 2:
        return math.exp(-x / 2)
    raise ValueError("df must be 1 or 2")


def logrank(t, s, g, keys, cause=1):
    """Log-rank test of cause-specific hazards (other statuses treated as censored)."""
    k = len(keys)
    O, E, V = np.zeros(k), np.zeros(k), np.zeros((k, k))
    for tj in np.unique(t[s == cause]):
        risk = t >= tj
        n = risk.sum()
        hit = (t == tj) & (s == cause)
        d = hit.sum()
        ng = np.array([np.sum(risk & (g == key)) for key in keys], float)
        O += [np.sum(hit & (g == key)) for key in keys]
        E += ng * d / n
        if n > 1:
            V += d * (n - d) / (n - 1) * (np.diag(ng / n) - np.outer(ng, ng) / n ** 2)
    z = (O - E)[:-1]
    chi2 = float(z @ np.linalg.solve(V[:-1, :-1], z))
    return chi2, k - 1, chi2_sf(chi2, k - 1)


def fmt_p(p):
    return "P < 0.001" if p < 0.001 else "P = %.3f" % p


def axes_in(fig, left, bottom, width, height, **kw):
    fh = fig.get_figheight()
    return fig.add_axes([left / W, bottom / fh, width / W, height / fh], **kw)


def time_axis(ax, label=True):
    ax.set_xlim(0, TMAX)
    ax.set_xticks(XTICKS)
    if label:
        ax.set_xlabel("Time from first perineural injection (h)")


def risk_table(fig, ax, rows, left, bottom, row_h):
    """Numbers at risk under each x tick; rows = [(label, color, times)] top to bottom."""
    rt = axes_in(fig, left, bottom, ax.get_position().width * W, row_h * (len(rows) + 1))
    rt.set_xlim(ax.get_xlim())
    rt.set_ylim(len(rows) + 0.5, -0.5)
    rt.axis("off")
    lab = blended_transform_factory(rt.transAxes, rt.transData)
    off = -0.28 / (ax.get_position().width * W)
    rt.text(off, -0.1, "No. at risk", transform=lab, ha="right", va="center", fontweight="bold")
    for i, (name, color, tt) in enumerate(rows, start=1):
        rt.text(off, i, name, transform=lab, ha="right", va="center", color=color)
        for x in XTICKS:
            rt.text(x, i, str(at_risk(tt, x)), ha="center", va="center")
    return rt


def save(fig, out_dir, name):
    png = os.path.join(out_dir, name + ".png")
    fig.savefig(png, dpi=DPI, facecolor="white")
    plt.close(fig)
    Image.open(png).convert("RGB").save(os.path.join(out_dir, name + ".tif"),
                                        compression="tiff_lzw", dpi=(DPI, DPI))


def report(title, t, s, g, groups):
    print("\n== %s ==" % title)
    for key, label, _ in groups:
        m = g == key
        xs, ys = aj(t[m], s[m], 1)
        print("%-16s n=%3d  events=%3d  CIF24=%5.1f%%  CIF96=%5.1f%%  at risk %s" % (
            label, m.sum(), np.sum(s[m] == 1), 100 * cif_at(xs, ys, 24), 100 * cif_at(xs, ys, 96),
            [at_risk(t[m], x) for x in XTICKS]))


def cif_by_group(t, s, g, groups, out_dir, name, legend_anchor):
    keys = [k for k, _, _ in groups]
    chi2, df, p = logrank(t, s, g, keys)
    report(name, t, s, g, groups)
    print("log-rank chi2=%.3f df=%d %s" % (chi2, df, fmt_p(p)))

    fig = plt.figure(figsize=(W, 4.3), dpi=DPI)
    ax = axes_in(fig, L, 1.45, AXW, 2.7)
    top = 0.0
    for (key, label, color), ls in zip(groups, LINESTYLES):
        m = g == key
        xs, ys = aj(t[m], s[m], 1)
        top = max(top, ys.max())
        ax.step(xs, 100 * ys, where="post", color=color, lw=1.3, ls=ls, clip_on=False,
                label="%s (n = %d, %d events)" % (label, m.sum(), np.sum(s[m] == 1)))
    ymax = 5 * math.ceil((100 * top + 4) / 5)
    ax.set_ylim(0, ymax)
    ax.set_yticks(np.arange(0, ymax + 0.1, 5 if ymax <= 40 else 10))
    ax.set_ylabel("Cumulative incidence of rescue analgesia (%)")
    time_axis(ax)
    leg = ax.legend(loc="lower right", bbox_to_anchor=legend_anchor, handlelength=2.6,
                    borderaxespad=0, labelspacing=0.4)
    leg.set_title("Log-rank %s" % fmt_p(p), prop={"size": 8})
    leg._legend_box.align = "left"
    risk_table(fig, ax, [(label, color, t[g == key]) for key, label, color in groups],
               L, 0.12, 0.2)
    save(fig, out_dir, name)


def figure2(t, s, g, out_dir, median_dis):
    print("\n== Figure2 ==")
    xr, yr = aj(t, s, 1)
    xd, yd = aj(t, s, 2)
    for x in (6, 12, 24, 48, 96):
        print("t=%2d h  rescue %5.1f%%  discharge %5.1f%%  at risk %d" % (
            x, 100 * cif_at(xr, yr, x), 100 * cif_at(xd, yd, x), at_risk(t, x)))
    bins = [(0, 4), (4, 6), (6, 12), (12, 24), (24, TMAX + 1)]
    print("rescue events in [0,4) [4,6) [6,12) [12,24) [24,96] h, median time to rescue:")
    for key, label, _ in [("all", "All", "")] + VOL:
        ev = t[(s == 1) & ((g == key) | (key == "all"))]
        print("  %-14s %s  %.1f h" % (label, [int(np.sum((ev >= lo) & (ev < hi))) for lo, hi in bins],
                                       np.median(ev)))
    print("median time among discharged without rescue: %.2f h" % np.median(t[s == 2]))

    fig = plt.figure(figsize=(W, 6.0), dpi=DPI)
    ax = axes_in(fig, L, 3.35, AXW, 2.35)
    ax.step(xd, 100 * yd, where="post", color="#4D4D4D", lw=1.3, ls=(0, (5, 2)), clip_on=False)
    ax.step(xr, 100 * yr, where="post", color="#D55E00", lw=1.3, clip_on=False)
    ax.text(TMAX, 100 * yd[-1] + 3, "Discharge without rescue analgesia", ha="right",
            va="bottom", color="#4D4D4D")
    ax.text(TMAX, 100 * yr[-1] + 3, "Rescue analgesia", ha="right", va="bottom", color="#D55E00")
    ax.set_ylim(0, 100)
    ax.set_yticks(range(0, 101, 20))
    ax.set_ylabel("Cumulative incidence (%)")
    time_axis(ax)
    ax.text(-0.27, 1.04, "A", transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
    risk_table(fig, ax, [("All patients", "black", t)], L, 2.45, 0.2)

    bx = axes_in(fig, L, 0.45, AXW, 1.55)
    rng = np.random.default_rng(2026)
    labels = []
    for i, ((key, label, color), mk) in enumerate(zip(VOL, MARKERS)):
        e = t[(g == key) & (s == 1)]
        bx.scatter(e, i + rng.uniform(-0.22, 0.22, e.size), s=13, color=color, marker=mk,
                   edgecolors="white", linewidths=0.3, zorder=3, clip_on=False)
        labels.append("%s\n%d events" % (label, e.size))
    bx.set_yticks(range(len(VOL)))
    bx.set_yticklabels(labels)
    bx.set_ylim(len(VOL) - 0.5, -0.5)
    bx.tick_params(axis="y", length=0)
    bx.spines["left"].set_visible(False)
    bx.axvline(median_dis, color="black", lw=0.8, ls=(0, (4, 2)), zorder=1)
    bx.text(median_dis + 1.2, -0.62, "Median time to discharge, %.1f h" % median_dis,
            ha="left", va="bottom")
    time_axis(bx)
    bx.text(-0.27, 1.10, "B", transform=bx.transAxes, fontsize=10, fontweight="bold", va="bottom")
    save(fig, out_dir, "Figure2")


def main():
    primary, lower, out_dir = sys.argv[1:4]
    median_dis = float(sys.argv[4]) if len(sys.argv) > 4 else 5.6
    os.makedirs(out_dir, exist_ok=True)
    t, s, g = read(primary, "vol_group")
    figure2(t, s, g, out_dir, median_dis)
    cif_by_group(t, s, g, VOL, out_dir, "Figure3", (1.0, 0.42))
    t, s, g = read(lower, "block")
    cif_by_group(t, s, g, BLOCK, out_dir, "FigureS1", (1.0, 0.04))


if __name__ == "__main__":
    main()
