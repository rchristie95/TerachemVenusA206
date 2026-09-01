"""Relabel the TOC render with the current result.

The shipped figures/toc_steom.png carried J ~ 111 and "Davydov splitting
2|J| ~ 221 cm^-1"; label_toc_white.py was itself one generation behind at
J = 32.8 / 2|J| = 65.6. Both predate the per-frame ensemble.

The wording changes too, not just the numbers. "Davydov splitting" asserts a
resolvable doublet, which is the opposite of what the paper concludes: the dimer
is detuning-dominated and the couplet is spectroscopically unresolved. The TOC
now states the comparison the paper actually makes.

Reads the CLEAN unlabelled render rather than re-labelling an already-labelled
file, which is what label_toc_white.py does (SRC == OUT).
"""
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
})

SRC = Path("/home/robson/PetaChem/plots_steom/toc_steom_white.png")
OUT = Path("/home/robson/PetaChem/UoP_Surrey_Venus_Collaboration/figures/toc_steom.png")

img = mpimg.imread(SRC)
h, w = img.shape[:2]
fig = plt.figure(figsize=(w / 300, h / 300), dpi=300)
ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(img); ax.axis("off")

ax.text(0.5, 0.93, "Excitonic coupling in the Venus fluorescent-protein dimer",
        transform=ax.transAxes, ha="center", va="center",
        color="black", fontsize=13)
ax.text(0.5, 0.075,
        r"$J = 27.1 \pm 2.5\ \mathrm{cm^{-1}}$"
        r"$\qquad \langle|\Delta|\rangle = 609\ \mathrm{cm^{-1}} \gg 2|J| = 54.1\ \mathrm{cm^{-1}}$",
        transform=ax.transAxes, ha="center", va="center",
        color="black", fontsize=9.5)

fig.savefig(OUT, dpi=300, facecolor="white")
print("wrote", OUT, w, "x", h)
