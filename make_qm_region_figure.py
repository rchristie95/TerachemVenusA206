"""Draw the exact 44-atom ORCA QM region as a publication-ready 2D scheme.

Connectivity and bond orders are inferred from the archived XYZ coordinates.
The three link hydrogens added at the QM/MM boundary are highlighted explicitly.
"""

from __future__ import annotations

from io import BytesIO
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import rdDepictor, rdDetermineBonds
from rdkit.Chem.Draw import rdMolDraw2D


ROOT = Path(__file__).resolve().parent
# Use the byte-retained geometry referenced by both production STEOM inputs so
# the published diagram and the electronic-structure model have one canonical
# coordinate source.  (steom_qm.xyz contains the same atoms to sub-microangstrom
# rounding, but is not the file executed by ORCA.)
XYZ = ROOT / "neo_model" / "orca_steom" / "geom_cthrp.xyz"
OUT = ROOT / "manuscript" / "Fig_QM_Region_44.png"

# One-based XYZ indices 42--44 are the QM/MM link hydrogens.
CAP_INDICES = {41, 42, 43}


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        Path("C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf"),
        Path("C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return ImageFont.truetype(str(candidate), size)
    return ImageFont.load_default()


def load_region() -> Chem.Mol:
    mol = Chem.MolFromXYZBlock(XYZ.read_text(encoding="utf-8"))
    if mol is None:
        raise RuntimeError(f"Could not read {XYZ}")
    rdDetermineBonds.DetermineConnectivity(mol)
    rdDetermineBonds.DetermineBondOrders(mol, charge=-1)
    if mol.GetNumAtoms() != 44:
        raise RuntimeError(f"Expected 44 atoms, found {mol.GetNumAtoms()}")
    for atom in mol.GetAtoms():
        atom.SetIntProp("originalIndex", atom.GetIdx())
    return mol


def draw_fragment(mol: Chem.Mol, width: int, height: int) -> Image.Image:
    rdDepictor.Compute2DCoords(mol)
    cap_atoms = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetIntProp("originalIndex") in CAP_INDICES
    ]
    for idx in cap_atoms:
        mol.GetAtomWithIdx(idx).SetProp("atomNote", "link H")

    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    options = drawer.drawOptions()
    options.padding = 0.08
    options.bondLineWidth = 3.0
    options.multipleBondOffset = 0.14
    options.additionalAtomLabelPadding = 0.12
    options.annotationFontScale = 0.62
    options.minFontSize = 22
    options.maxFontSize = 42
    options.explicitMethyl = True
    options.clearBackground = False
    drawer.DrawMolecule(
        mol,
        highlightAtoms=cap_atoms,
        highlightAtomColors={idx: (0.95, 0.55, 0.10) for idx in cap_atoms},
        highlightAtomRadii={idx: 0.38 for idx in cap_atoms},
    )
    drawer.FinishDrawing()
    return Image.open(BytesIO(drawer.GetDrawingText())).convert("RGBA")


def centred(draw: ImageDraw.ImageDraw, text: str, box: tuple[int, int, int, int], fnt) -> None:
    left, top, right, bottom = box
    bbox = draw.textbbox((0, 0), text, font=fnt)
    x = left + (right - left - (bbox[2] - bbox[0])) / 2
    y = top + (bottom - top - (bbox[3] - bbox[1])) / 2
    draw.text((x, y), text, fill=(25, 25, 25), font=fnt)


def main() -> None:
    mol = load_region()
    fragments = list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True))
    if sorted(fragment.GetNumAtoms() for fragment in fragments) != [13, 31]:
        raise RuntimeError("Expected 31-atom CR2 and 13-atom Tyr203 fragments")
    cr2 = next(fragment for fragment in fragments if fragment.GetNumAtoms() == 31)
    tyr = next(fragment for fragment in fragments if fragment.GetNumAtoms() == 13)

    canvas = Image.new("RGB", (2600, 1500), "white")
    draw = ImageDraw.Draw(canvas)
    title_font = font(58, bold=True)
    label_font = font(42, bold=True)
    detail_font = font(34)
    small_font = font(30)

    centred(draw, "44-atom QM region (net charge -1)", (0, 20, 2600, 120), title_font)
    draw.rounded_rectangle((70, 150, 1720, 1285), radius=32, outline=(95, 110, 125), width=4)
    draw.rounded_rectangle((1770, 150, 2530, 1285), radius=32, outline=(95, 110, 125), width=4)

    cr2_image = draw_fragment(cr2, 1570, 930)
    tyr_image = draw_fragment(tyr, 680, 930)
    canvas.paste(cr2_image, (110, 255), cr2_image)
    canvas.paste(tyr_image, (1810, 255), tyr_image)

    centred(draw, "Anionic CR2 chromophore", (70, 160, 1720, 240), label_font)
    centred(draw, "31 atoms; two link-H caps", (70, 1200, 1720, 1270), detail_font)
    centred(draw, "Tyr203 phenol", (1770, 160, 2530, 240), label_font)
    centred(draw, "13 atoms; one link-H cap", (1770, 1200, 2530, 1270), detail_font)

    draw.line((1745, 330, 1745, 1110), fill=(135, 135, 135), width=4)
    for y in range(350, 1110, 42):
        draw.ellipse((1738, y, 1752, y + 14), fill=(135, 135, 135))
    centred(draw, "pi-stacked", (1610, 665, 1880, 730), small_font)

    draw.ellipse((820, 1350, 860, 1390), fill=(242, 140, 26), outline=(140, 75, 5), width=2)
    draw.text(
        (880, 1344),
        "QM/MM boundary hydrogen cap (three total; XYZ atoms 42-44)",
        fill=(45, 45, 45),
        font=detail_font,
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(OUT, dpi=(300, 300), optimize=True)
    print(f"Wrote {OUT} ({mol.GetNumAtoms()} atoms; charge -1)")


if __name__ == "__main__":
    main()
