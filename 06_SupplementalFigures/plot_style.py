"""Shared text style for the Supplementary Figure scripts.

Purpose:      One place for the typography of the supplementary figures: Arial (the metrically
              identical Liberation Sans as a fallback, rewritten to Arial in the SVG), editable
              SVG text, italic statistical symbols (n, p, q, r, z, t) through matplotlib
              mathtext exactly as in 05_Figure5/02_Figure5C.py, small values written as
              a×10^b rather than in e-notation, true minus signs, and SVGs saved without a
              background so that they drop into Inkscape without a page-sized white box.

Inputs:       none.

Outputs:      sym(), num(), stat(), panel_letter() and save(), imported by the scripts in this folder.

Dependencies: matplotlib.
"""
# REVISION 2026-10-08: new module (journal request: italic scalar variables, roman numbers and
# functions, "×" for multiplication, no e-notation).
import glob
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager

FONT = "Arial"
MINUS = "−"


def _setup():
    available = {f.name for f in font_manager.fontManager.ttflist}
    if FONT not in available:
        for path in glob.glob("/usr/share/fonts/**/LiberationSans-*.ttf", recursive=True):
            font_manager.fontManager.addfont(path)
        available = {f.name for f in font_manager.fontManager.ttflist}
    font = FONT if FONT in available else ("Liberation Sans" if "Liberation Sans" in available else "sans-serif")
    plt.rcParams.update({
        "font.family": font,
        "svg.fonttype": "none",          # text stays text in the SVG
        "pdf.fonttype": 42,
        "axes.unicode_minus": True,      # tick labels use the minus sign, not a hyphen
        "mathtext.fontset": "custom",    # math in the plot font: $q$ is Arial italic
        "mathtext.rm": font,
        "mathtext.it": f"{font}:italic",
        "mathtext.bf": f"{font}:bold",
        "mathtext.bfit": f"{font}:italic:bold",
        "mathtext.cal": font,
        "mathtext.sf": font,
    })
    return font


PLOT_FONT = _setup()


def sym(s, bold=False):
    """A statistical symbol in italic ($q$); bold italic inside a bold label."""
    return f"$\\mathbfit{{{s}}}$" if bold else f"${s}$"


def num(v, digits=2, sci_below=0.001, sig=2, bold=False):
    """A label number: fixed decimals, or a×10^b below `sci_below`; never e-notation."""
    if v is None or v != v:
        return "NA"
    if v != 0 and abs(v) < sci_below:
        mantissa, exponent = f"{v:.{sig - 1}e}".split("e")
        exponent = str(int(exponent)).replace("-", MINUS)
        style = "\\mathbf" if bold else "\\mathrm"
        return f"{mantissa.replace('-', MINUS)}${style}{{\\times10}}^{{{style}{{{exponent}}}}}$"
    return f"{v:.{digits}f}".replace("-", MINUS)


def stat(s, v, bold=False, **kw):
    """'q=0.042' with an italic q, e.g. stat('q', 0.042) or stat('q', 3.9e-4, digits=3)."""
    return f"{sym(s, bold)}={num(v, bold=bold, **kw)}"


def panel_letter(fig, x, y, letter, fontsize=15):
    """Boxed panel letter, the style of the assembled main and supplementary figures."""
    fig.text(x, y, letter, fontsize=fontsize, fontweight="normal",
             bbox=dict(boxstyle="square,pad=0.35", facecolor="none", edgecolor="black", linewidth=1.0))


def save(fig, basename, figures_dir, dpi=300):
    """PNG (white background), PDF and SVG (transparent background, editable Arial text)."""
    for ext in ("png", "pdf", "svg"):
        out = Path(figures_dir) / f"{basename}.{ext}"
        if ext == "svg":
            fig.savefig(out, bbox_inches="tight", transparent=True)
            if PLOT_FONT != FONT:
                text = out.read_text(encoding="utf-8")
                for quote in ('"', "'"):
                    text = text.replace(f"font-family: {quote}{PLOT_FONT}{quote}", f"font-family: {FONT}")
                text = text.replace(f"font-family: {PLOT_FONT}", f"font-family: {FONT}")
                out.write_text(text, encoding="utf-8")
        else:
            fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
        print(f"Saved: {out.name}")
