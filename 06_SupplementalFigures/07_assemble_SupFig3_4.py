#!/usr/bin/env python3
"""Supplementary Figures 3 and 4: assemble the panel SVGs into A4 figure pages.

Purpose:      Lays out the panels written by 03_SupFig3A.py / 03_SupFig3B.py (Supplementary
              Figure 3) and 04_SupFig4AB.py / 04_SupFig4C.py (Supplementary Figure 4) on an A4
              page with boxed panel letters, replacing the manual Inkscape assembly. The panels
              are embedded as vector graphics, so all text stays editable.

Inputs:       figures/SupFig3A.svg, SupFig3B.svg, SupFig4AB.svg, SupFig4C.svg.

Outputs:      figures/SupFig3.svg and SupFig4.svg; PDF and 300-dpi PNG as well when Inkscape is
              available (on PATH, or its path in the INKSCAPE environment variable).

Dependencies: Python + lxml; Inkscape 1.x optional for the PDF and PNG exports.
"""
# REVISION 2026-10-08: new script (Supplementary Figures 3 and 4 previously assembled by hand).
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from lxml import etree

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import FIGURES_DIR

SVG = "http://www.w3.org/2000/svg"
XLINK = "http://www.w3.org/1999/xlink"
PT_MM = 25.4 / 72
PAGE_W, PAGE_H, LEFT, TOP = 210.0, 297.0, 12.0, 6.0   # A4, mm


def _pt(length):
    value = float(re.match(r"[\d.]+", length).group())
    return value * (1 / PT_MM if length.endswith("mm") else 1.0)


def panel(path, prefix, x, y, width):
    """Nested <svg> for one panel at (x, y) mm, `width` mm wide; returns (element, height, mm per pt)."""
    root = etree.parse(str(path)).getroot()
    w_pt, h_pt = _pt(root.get("width")), _pt(root.get("height"))
    scale = width / w_pt
    for el in root.iter():
        if el.get("id"):
            el.set("id", prefix + el.get("id"))
        for attr, value in list(el.attrib.items()):
            if "url(#" in value:
                el.set(attr, value.replace("url(#", f"url(#{prefix}"))
            elif attr in (f"{{{XLINK}}}href", "href") and value.startswith("#"):
                el.set(attr, "#" + prefix + value[1:])
    nested = etree.Element(f"{{{SVG}}}svg", x=f"{x:.3f}", y=f"{y:.3f}", width=f"{width:.3f}",
                           height=f"{h_pt * scale:.3f}", viewBox=root.get("viewBox"))
    for child in root:
        nested.append(child)
    return nested, h_pt * scale, scale


def letter(page, letter_text, x, y, size_pt):
    """Boxed panel letter with its top-left corner at (x, y) mm."""
    side = size_pt * PT_MM * 1.4
    g = etree.SubElement(page, f"{{{SVG}}}g", id=f"panel_letter_{letter_text}")
    etree.SubElement(g, f"{{{SVG}}}rect", x=f"{x:.3f}", y=f"{y:.3f}", width=f"{side:.3f}",
                     height=f"{side:.3f}", style="fill:none;stroke:#000000;stroke-width:0.35")
    t = etree.SubElement(g, f"{{{SVG}}}text", x=f"{x + side / 2:.3f}", y=f"{y + side * 0.74:.3f}",
                         style=f"font-family:Arial;font-size:{size_pt * PT_MM:.3f}px;"
                               "text-anchor:middle;fill:#000000")
    t.text = letter_text


def page():
    return etree.Element(f"{{{SVG}}}svg", nsmap={None: SVG, "xlink": XLINK}, version="1.1",
                         width=f"{PAGE_W}mm", height=f"{PAGE_H}mm", viewBox=f"0 0 {PAGE_W} {PAGE_H}")


def write(fig, name):
    out = FIGURES_DIR / f"{name}.svg"
    etree.ElementTree(fig).write(str(out), xml_declaration=True, encoding="utf-8", pretty_print=True)
    print(f"Saved: {out.name}")
    inkscape = os.environ.get("INKSCAPE") or shutil.which("inkscape")
    if inkscape:
        subprocess.run([inkscape, str(out), "--export-type=pdf", f"--export-filename={out.with_suffix('.pdf')}"], check=True)
        subprocess.run([inkscape, str(out), "--export-type=png", "--export-dpi=300", "--export-background=white",
                        "--export-background-opacity=1", f"--export-filename={out.with_suffix('.png')}"], check=True)
        print(f"Saved: {name}.pdf, {name}.png")
    else:
        print("Inkscape not found: PDF and PNG not exported")


def main():
    width = PAGE_W - LEFT - 6.0
    # Supplementary Figure 3: A (bone-marrow myeloid, Zavidij) above B (shipped samples only).
    fig = page()
    a, h_a, _ = panel(FIGURES_DIR / "SupFig3A.svg", "A_", LEFT, TOP, 120.0)
    b, h_b, _ = panel(FIGURES_DIR / "SupFig3B.svg", "B_", LEFT, TOP + h_a + 8.0, width)
    fig.extend([a, b])
    letter(fig, "A", 2.0, TOP, 12)
    letter(fig, "B", 2.0, TOP + h_a + 8.0, 12)
    write(fig, "SupFig3")

    # Supplementary Figure 4: A and B (lettered inside SupFig4AB) above C (Boiarsky PCs). The C
    # letter takes the size the A and B letters have after scaling, so the three match.
    fig = page()
    ab, h_ab, scale = panel(FIGURES_DIR / "SupFig4AB.svg", "AB_", LEFT, TOP, width)
    c, _, _ = panel(FIGURES_DIR / "SupFig4C.svg", "C_", LEFT, TOP + h_ab + 8.0, 85.0)   # text size matched to A and B
    fig.extend([ab, c])
    letter(fig, "C", 2.0, TOP + h_ab + 8.0, 15 * scale / PT_MM)
    write(fig, "SupFig4")


if __name__ == "__main__":
    main()
