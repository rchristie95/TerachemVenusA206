# label_toc.py — add the paper's punchline caption to the raytraced TOC render.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pathlib import Path

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["cmr10", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
    "axes.formatter.use_mathtext": True,
})

ROOT = Path(__file__).resolve().parent
# The retained manuscript asset already contains the previous labels.  Blank
# only the white title/caption bands before drawing the audited values; the
# coordinate-derived molecular render itself is left byte-for-pixel unchanged.
SRC = ROOT / "manuscript" / "toc_steom.png"
OUT = ROOT / "manuscript" / "toc_steom.png"

img = mpimg.imread(SRC)
img[:135, :, :3] = 1.0
img[-145:, :, :3] = 1.0
h, w = img.shape[:2]
fig = plt.figure(figsize=(w / 300, h / 300), dpi=300)
ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(img); ax.axis("off")

ax.text(0.5, 0.925,
        "Excitonic coupling in the Venus fluorescent-protein dimer",
        transform=ax.transAxes, ha="center", va="center",
        color="black", fontsize=13)
ax.text(0.5, 0.075,
        r"Tandem ensemble: $J=32.8\pm1.6\ \mathrm{cm^{-1}}$"
        r"     Davydov splitting: $2|J|=65.6\pm3.1\ \mathrm{cm^{-1}}$",
        transform=ax.transAxes, ha="center", va="center",
        color="black", fontsize=9.5)

fig.savefig(OUT, dpi=300, facecolor="white")
print("wrote", OUT, w, "x", h)
