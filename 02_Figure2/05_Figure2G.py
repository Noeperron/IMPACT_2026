#!/usr/bin/env python3
# ============================================================================
# Purpose:      Figure 2G: paired pre/post-teclistamab measurements in 10 SMM
#               patients enrolled on the teclistamab arm of Immuno-PRISM
#               (NCT05469893, DFCI 22-154). Four sub-panels ordered L->R by
#               the mechanistic causal chain:
#                 (1) tumor burden (serum M-spike) -- clinical response,
#                     excluding the 1 Kappa-LC patient with no baseline
#                     M-spike (n=9); IFX+/IFX- treated as below detection
#                     (0 g/dL) with connector lines
#                 (2) sBCMA -- soluble decoy; collapses with plasma-cell clearance
#                 (3) APRIL -- rises as tumor clears
#                 (4) BAFF  -- rises concurrently; rules out IVIG-anti-BAFF
#                              compensatory induction (would predict fall).
#               Paired Wilcoxon signed-rank tests with Benjamini-Hochberg
#               correction ACROSS THE FOUR SUB-PANELS (was 3; updated 2026-07-08
#               after Bucket 1 of the Romanos-feedback overhaul).
# Inputs:       data/external/PrePostTEC_olink_deid.csv (deid: P01..P10).
#               data/external/PrePostTEC_cohort_deid.csv (baseline + post M-spike).
# Outputs:      Figure2G_PrePostTEC_4panels.{png,pdf,svg} in figures/.
# Dependencies: Python + pandas, scipy, matplotlib. Sources ../config.py.
# ============================================================================
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# Every figure in this paper is set in Arial. Arial is not redistributable and is usually absent on
# Linux, so fall back to Liberation Sans, which is metrically identical (same glyph advance widths).
# svg.fonttype="none" keeps text as TEXT rather than converting glyphs to outlines, so the panel
# stays editable in the figure-assembly software; the saved SVG is then rewritten to say Arial.
FONT = "Arial"
_avail = {f.name for f in font_manager.fontManager.ttflist}
if FONT not in _avail:
    for _p in glob.glob("/usr/share/fonts/**/LiberationSans-*.ttf", recursive=True):
        font_manager.fontManager.addfont(_p)
    _avail = {f.name for f in font_manager.fontManager.ttflist}
PLOT_FONT = FONT if FONT in _avail else ("Liberation Sans" if "Liberation Sans" in _avail else "sans-serif")
plt.rcParams["font.family"] = PLOT_FONT
plt.rcParams["svg.fonttype"] = "none"
# REVISION 2026-10-08: math text ($q$, $n$) in the plot font, italic.
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = PLOT_FONT
plt.rcParams["mathtext.it"] = f"{PLOT_FONT}:italic"

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
try:
    from config import DATA_DIR, FIGURES_DIR
except Exception:
    DATA_DIR = HERE.parent / "data"
    FIGURES_DIR = HERE.parent / "figures"

IN_OLINK  = Path(DATA_DIR) / "external" / "PrePostTEC_olink_deid.csv"
IN_COHORT = Path(DATA_DIR) / "external" / "PrePostTEC_cohort_deid.csv"
FIGURES_DIR = Path(FIGURES_DIR)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# L->R order per the mechanistic causal chain (Bucket 1 of Romanos overhaul):
# tumor burden DOWN -> soluble decoy DOWN -> ligands (APRIL/BAFF) UP.
SUBPANELS = ["Mspike", "sBCMA", "APRIL", "BAFF"]
GENE_LABELS = {"Mspike": "Serum M-spike (g/dL)",
               "APRIL":  "TNFSF13 (APRIL)",
               "BAFF":   "TNFSF13B (BAFF)",
               "sBCMA":  "TNFRSF17 (sBCMA)"}
BOX_COLORS = ["#A6CEE3", "#FBB4AE"]  # light blue baseline, light red post-TEC

# ---- BH adjustment helper (no external dep) ----
def bh_qvalues(p_array):
    """Benjamini-Hochberg adjusted q-values (FDR)."""
    p_array = np.asarray(p_array, dtype=float)
    n = len(p_array)
    order = np.argsort(p_array)
    ranks_0idx = np.empty(n, dtype=int)
    ranks_0idx[order] = np.arange(n)
    raw_q = p_array * n / (ranks_0idx + 1)
    q = np.empty(n)
    min_so_far = 1.0
    for r in range(n - 1, -1, -1):
        idx = order[r]
        min_so_far = min(min_so_far, raw_q[idx])
        q[idx] = min(min_so_far, 1.0)
    return q

# ---- Load Olink (proteins) ----
df = pd.read_csv(IN_OLINK)
df["NPX"] = pd.to_numeric(df["NPX"], errors="coerce")
print(f"Olink: {len(df)} rows; {df['Patient_ID'].nunique()} patients; proteins: {df['Protein'].unique()}")

# ---- Load M-spike cohort table (baseline + post) ----
coh = pd.read_csv(IN_COHORT)
def _coerce_mspike(x):
    """IFX- (immunofixation negative, no measurable) and IFX+ (trace, not
    quantifiable) both map to 0.0 g/dL for plotting; captured in caption."""
    x = str(x).strip()
    if x in ("IFX-", "IFX+"):
        return 0.0
    try:
        return float(x)
    except ValueError:
        return np.nan
coh["Baseline_M_spike_num"] = coh["Baseline_M_spike_g_per_dL"].map(_coerce_mspike)
coh["PostTEC_M_spike_num"]  = coh["PostTEC_M_spike_g_per_dL"].map(_coerce_mspike)


def _imwg_response(baseline_str, post_str):
    """IMWG response classification by serum M-spike.

    CR   : post-treatment M-spike is immunofixation-negative (IFX-).
    VGPR : post is IFX+ (detectable only by immunofixation) OR quantifiable
           with >=90% reduction from baseline.
    PR   : quantifiable with 50-89% reduction from baseline.
    """
    post = str(post_str).strip()
    if post == "IFX-":
        return "CR"
    if post == "IFX+":
        return "VGPR"
    try:
        b = float(baseline_str)
        p = float(post_str)
    except ValueError:
        return "N/A"
    if b <= 0:
        return "N/A"
    reduction = (b - p) / b
    if reduction >= 0.90:
        return "VGPR"
    if reduction >= 0.50:
        return "PR"
    return "MR/SD"


coh["IMWG"] = [_imwg_response(b, p) for b, p in
               zip(coh["Baseline_M_spike_g_per_dL"], coh["PostTEC_M_spike_g_per_dL"])]

# Exclude Kappa-LC patient (has no measurable M-spike at baseline; light-chain only).
coh_mspike = coh[coh["Isotype"] != "Kappa_LC"].copy()
n_mspike = len(coh_mspike)
imwg_counts = coh_mspike["IMWG"].value_counts().to_dict()
n_cr   = int(imwg_counts.get("CR", 0))
n_vgpr = int(imwg_counts.get("VGPR", 0))
n_pr   = int(imwg_counts.get("PR", 0))
print(f"M-spike: n={n_mspike} (excluded 1 Kappa_LC patient with no baseline M-spike)")
print(f"M-spike IMWG response: {n_cr} CR / {n_vgpr} VGPR / {n_pr} PR")

# ---- Compute paired Wilcoxon stats for the FOUR sub-panels UP FRONT (BH across 4) ----
stats_data = {}

# (1) M-spike: different scale, so track separately but include in the BH family
sub_b = coh_mspike["Baseline_M_spike_num"].values
sub_p = coh_mspike["PostTEC_M_spike_num"].values
delta = sub_p - sub_b
n_up = int(np.sum(delta > 0))
w_ms  = stats.wilcoxon(sub_p, sub_b, alternative="two-sided", zero_method="wilcox")
# For percent reduction: use only baseline > 0 (all 9 remaining have baseline > 0)
pct_red = 100.0 * (sub_b - sub_p) / sub_b
median_pct_red = float(np.median(pct_red))
stats_data["Mspike"] = dict(
    baseline=sub_b, post=sub_p, delta=delta, n=n_mspike,
    p=float(w_ms.pvalue),
    median_delta=float(np.median(delta)),
    median_pct_red=median_pct_red,
    n_up=n_up,
    n_cr=n_cr, n_vgpr=n_vgpr, n_pr=n_pr,
    # IMWG clinical response categories carry both the direction and magnitude
    # of change without depending on an arbitrary LOD imputation (fold-change
    # is undefined when post-treatment values hit IFX-/IFX+ censoring).
    display_effect=f"{n_cr} CR / {n_vgpr} VGPR / {n_pr} PR",
    ylabel="Serum M-spike (g/dL)",
    is_protein=False,
)

# (2-4) Proteins
for protein in ["sBCMA", "APRIL", "BAFF"]:
    sub = df[df["Protein"] == protein].copy()
    wide = sub.pivot(index="Patient_ID", columns="Timepoint", values="NPX").dropna()
    baseline = wide["Baseline"].values
    post = wide["PostTEC"].values
    delta = post - baseline
    n = len(wide)
    w_two = stats.wilcoxon(post, baseline, alternative="two-sided", zero_method="wilcox")
    median_delta = float(np.median(delta))
    median_fc_linear = float(2 ** median_delta)
    n_up = int(np.sum(delta > 0))
    stats_data[protein] = dict(
        baseline=baseline, post=post, delta=delta, n=n,
        p=float(w_two.pvalue), median_delta=median_delta,
        median_fc_linear=median_fc_linear, n_up=n_up,
        display_effect=f"median fold-change = {median_fc_linear:.2f}x",
        ylabel="Plasma levels (Olink NPX)",
        is_protein=True,
    )

# Benjamini-Hochberg across the FOUR sub-panels (family-wise correction; matches
# the manuscript's q<0.1 threshold convention).
q_array = bh_qvalues([stats_data[s]["p"] for s in SUBPANELS])
for i, s in enumerate(SUBPANELS):
    stats_data[s]["q"] = float(q_array[i])

# ---- Four sub-panels (L->R by mechanistic causal chain) ----
fig, axes = plt.subplots(1, 4, figsize=(12.5, 3.9), dpi=200)
positions = [1, 2]
np.random.seed(42)

for ax, key in zip(axes, SUBPANELS):
    s = stats_data[key]
    baseline = s["baseline"]; post = s["post"]; delta = s["delta"]; n = s["n"]
    n_up = s["n_up"]

    # Boxplot
    bp = ax.boxplot([baseline, post], positions=positions, widths=0.55, showfliers=False,
                    patch_artist=True,
                    medianprops=dict(color="black", linewidth=1.4),
                    whiskerprops=dict(color="black", linewidth=1),
                    capprops=dict(color="black", linewidth=1))
    for patch, color in zip(bp["boxes"], BOX_COLORS):
        patch.set_facecolor(color)
        patch.set_edgecolor("black")
        patch.set_linewidth(1)
        patch.set_zorder(1)

    # Connector lines + points (lines OVER boxes)
    jitter = np.random.uniform(-0.06, 0.06, size=n)
    for i in range(n):
        ax.plot([positions[0] + jitter[i], positions[1] + jitter[i]],
                [baseline[i], post[i]],
                color="grey", linewidth=0.6, alpha=0.7, zorder=2.5)
    for i in range(n):
        ax.scatter([positions[0] + jitter[i], positions[1] + jitter[i]],
                   [baseline[i], post[i]],
                   s=22, color="#333333", edgecolor="black", linewidth=0.4, zorder=3)

    # Significance annotation bar
    ymax = max(post.max(), baseline.max())
    ymin = min(post.min(), baseline.min())
    y_span = max(ymax - ymin, 1e-6)
    bar_y = ymax + 0.10 * y_span
    ax.plot([positions[0], positions[1]], [bar_y, bar_y], color="black", linewidth=0.9)
    ax.plot([positions[0], positions[0]], [bar_y - 0.025*y_span, bar_y], color="black", linewidth=0.9)
    ax.plot([positions[1], positions[1]], [bar_y - 0.025*y_span, bar_y], color="black", linewidth=0.9)
    q_val = s["q"]
    q_str = f"$q$ = {q_val:.3f}" if q_val >= 0.001 else f"$q$ = {q_val:.1e}"
    # Only the q-value is annotated. The IMWG response categories were removed from the M-spike
    # panel and the fold-change from the protein panels; magnitude and direction are already
    # carried by the boxplots and the connector slopes.
    ax.text((positions[0] + positions[1]) / 2, bar_y + 0.05 * y_span, q_str,
            ha="center", va="bottom", fontsize=15)

    # Cosmetics
    # REVISION 2026-10-08: larger text (was 14 / 13 / 13 / 12).
    ax.set_title(GENE_LABELS[key], fontsize=16, loc="left", fontweight="bold", pad=4)
    ax.set_xticks(positions)
    ax.set_xticklabels(["Pre-TEC", "Post-TEC"], fontsize=15)
    ax.set_ylabel(s["ylabel"], fontsize=15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=15)
    ax.set_ylim(ymin - 0.10*y_span, ymax + 0.40*y_span)

    if s["is_protein"]:
        print(f"{key:>7}: median delta NPX = {s['median_delta']:+.3f}, "
              f"fold-change = {s['median_fc_linear']:.2f}, "
              f"paired Wilcoxon p = {s['p']:.3g}, "
              f"BH q = {s['q']:.3g}, {n_up}/{n} directional")
    else:
        print(f"{key:>7}: IMWG {s['n_cr']} CR / {s['n_vgpr']} VGPR / {s['n_pr']} PR, "
              f"median delta = {s['median_delta']:+.2f} g/dL, "
              f"median % reduction = {s['median_pct_red']:.1f}%, "
              f"paired Wilcoxon p = {s['p']:.3g}, "
              f"BH q = {s['q']:.3g}")

# REVISION 2026-10-08: title as used in the assembled figure, editable text, italic n.
fig.suptitle("Pre vs post teclistamab plasma proteomics ($n$=10 SMM patients, Immuno-PRISM)",
             fontsize=18, y=0.965)
plt.tight_layout(rect=[0, 0, 1, 0.95])

for ext in ("png", "pdf", "svg"):
    out = FIGURES_DIR / f"Figure2G_PrePostTEC_4panels.{ext}"
    # REVISION 2026-10-08: transparent SVG (no white background box in Inkscape).
    plt.savefig(out, dpi=200, bbox_inches="tight", transparent=(ext == "svg"))
    if ext == "svg" and PLOT_FONT != FONT:
        # matplotlib records the family it actually resolved. Rewrite it to Arial: Liberation Sans
        # is metrically identical, so every glyph position already written stays correct, and the
        # panel then opens in real Arial in the figure-assembly software.
        txt = out.read_text(encoding="utf-8")
        for q in ('"', "'"):
            txt = txt.replace(f"font-family: {q}{PLOT_FONT}{q}", f"font-family: {FONT}")
        txt = txt.replace(f"font-family: {PLOT_FONT}", f"font-family: {FONT}")
        out.write_text(txt, encoding="utf-8")
    print(f"Saved: {out}")
plt.close()

# ---- Write a machine-readable numbers file for downstream text edits ----
import json
numbers = {}
for k in SUBPANELS:
    s = stats_data[k]
    numbers[k] = dict(
        n=int(s["n"]),
        p=float(s["p"]),
        q=float(s["q"]),
        median_delta=float(s["median_delta"]),
        n_up=int(s["n_up"]),
    )
    if s["is_protein"]:
        numbers[k]["median_fc_linear"] = float(s["median_fc_linear"])
    else:
        numbers[k]["median_pct_reduction"] = float(s["median_pct_red"])
numbers_out = FIGURES_DIR / "Figure2G_4panels_numbers.json"
with numbers_out.open("w") as f:
    json.dump(numbers, f, indent=2)
print(f"Numbers written to: {numbers_out}")
