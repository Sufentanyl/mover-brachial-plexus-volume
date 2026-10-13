import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image

# usage: python fig1_flowchart.py out.png out.tif [width_in] [dpi] [font_pt]
WT = float(sys.argv[3]) if len(sys.argv) > 3 else 5.87   # 14.9 cm, Minerva double column
DPI = int(sys.argv[4]) if len(sys.argv) > 4 else 600
FS = float(sys.argv[5]) if len(sys.argv) > 5 else 7.5
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": FS})
W, H = 6.7, 5.0   # layout units; the width is compressed to WT inches
fig = plt.figure(figsize=(WT, H), dpi=DPI)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

def box(xc, yc, w, h, text, lw=0.8, bold=False):
    ax.add_patch(FancyBboxPatch((xc - w/2, yc - h/2), w, h, boxstyle="square,pad=0", lw=lw, ec="black", fc="white"))
    ax.text(xc, yc, text, ha="center", va="center", fontsize=FS, fontweight="bold" if bold else "normal", linespacing=1.25)
def arrow(x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", lw=0.8, color="black", mutation_scale=8, shrinkA=0, shrinkB=0))

# left column (brachial plexus)
xL, wL = 1.08, 2.05
ys = [4.62, 3.80, 2.98, 2.16, 1.34, 0.42]
hs = [0.46, 0.34, 0.46, 0.34, 0.34, 0.34]
texts = ["Brachial plexus block billing records\n(CPT 64415/64416)\nn = 693",
         "Single-shot brachial plexus blocks\nn = 690",
         "Perineural local anesthetic record\nin the same surgical case\nn = 643",
         "Eligible single-shot blocks\nn = 600",
         "Blocks before the end of anesthesia\nn = 573",
         "Primary cohort\nn = 571"]
for k, (y, h, t) in enumerate(zip(ys, hs, texts)):
    box(xL, y, wL, h, t, lw=1.6 if k == 5 else 0.8, bold=(k == 5))
for k in range(5):
    arrow(xL, ys[k] - hs[k]/2, xL, ys[k+1] + hs[k+1]/2)

# exclusion boxes (middle column)
xE, wE = 3.33, 1.95
excl = ["Continuous catheter (CPT 64416)\nn = 3",
        "No perineural local anesthetic\nrecord in the same surgical case\nn = 47",
        "Anatomically implausible surgery\nn = 43",
        "Block after the end of anesthesia\nn = 27",
        "Block >24 h before the start\nof anesthesia\nn = 2"]
eh = [0.34, 0.46, 0.34, 0.34, 0.46]
for k in range(5):
    ym = (ys[k] - hs[k]/2 + ys[k+1] + hs[k+1]/2) / 2
    box(xE, ym, wE, eh[k], excl[k])
    arrow(xL, ym, xE - wE/2, ym)

# right column (lower extremity)
xR, wR = 5.57, 2.05
box(xR, ys[0], wR, 0.46, "Femoral/sciatic block billing records\n(CPT 64447/64445)\nn = 356")
yRb = 0.55; hRb = 0.62
box(xR, yRb, wR, hRb, "Lower-extremity cohort\nn = 245\nFemoral 115; femoral–sciatic 119;\nsciatic 11", lw=1.6)
xline = 4.75
ax.annotate("", xy=(xline, yRb + hRb/2), xytext=(xline, ys[0] - 0.23), arrowprops=dict(arrowstyle="-|>", lw=0.8, color="black", mutation_scale=8, shrinkA=0, shrinkB=0))
wX = 1.75; xX = xline + 0.15 + wX/2
box(xX, 2.6, wX, 0.76, "Excluded\nNo perineural local\nanesthetic record, n = 21\nAnatomically implausible\nsurgery, n = 90")
arrow(xline, 2.6, xX - wX/2, 2.6)

out_png = sys.argv[1]; out_tif = sys.argv[2]
fig.savefig(out_png, dpi=DPI, facecolor="white")
img = Image.open(out_png).convert("RGB")
img.save(out_tif, compression="tiff_lzw", dpi=(DPI, DPI))
print(img.size)
