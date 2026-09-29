#!/usr/bin/env python3
"""Figure 3E (manuscript Figure 3E): APRIL-responsive module score in tumor versus non-malignant
plasma cells, SWIFT-seq cohort (Lightbody et al., Nat Cancer 2025).

Purpose:      Within-sample paired comparison of the 15-gene APRIL-responsive module between
              clonal (tumor) and non-malignant plasma cells, in bone marrow and peripheral blood.
              Both populations come from the same specimen, so patient, batch and sample handling
              are held constant by construction.

Inputs:       data/external/swiftseq_april_persample_deid.csv -- de-identified per-specimen
              summary (Sample_ID, Patient_ID, Compartment, Disease_Stage, library count,
              per-population cell counts and mean module scores). One row is one participant and
              one tissue at baseline, with multiple libraries of the same specimen pooled.
              Cell-level scoring was performed upstream on the SWIFT-seq plasma-cell object with
              scanpy sc.tl.score_genes() using the gene list in Supplementary Table 5; source data
              are in dbGaP phs003855.v1.p1 (controlled access).

              Cohort filtering follows the convention of the source study (Lightbody et al.).
              Secondary tumour clones (Tumor2) are excluded upstream, so a specimen's tumour mean
              is always its dominant clone. No covariate adjustment is applied or possible here:
              both populations come from one specimen, so age and sex are identical within a pair
              and cancel in the within-pair difference. The design is controlled for them by
              construction, as the paired Olink panels of Figure 5 are.

Outputs:      figures/Figure3E.png (and PDF + SVG).

Dependencies: Python + pandas, numpy, scipy, statsmodels, matplotlib; reads config.py.
"""
import sys, glob
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import pandas as pd, numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib import font_manager as _fm

# Every figure in this paper is set in Arial; fall back to the metrically identical Liberation
# Sans and rewrite the SVG font-family. svg.fonttype="none" keeps text editable.
FONT = "Arial"
_av = {f.name for f in _fm.fontManager.ttflist}
if FONT not in _av:
    for _p in glob.glob("/usr/share/fonts/**/LiberationSans-*.ttf", recursive=True):
        _fm.fontManager.addfont(_p)
    _av = {f.name for f in _fm.fontManager.ttflist}
PLOT_FONT = FONT if FONT in _av else ("Liberation Sans" if "Liberation Sans" in _av else "sans-serif")
plt.rcParams["font.family"] = PLOT_FONT
plt.rcParams["svg.fonttype"] = "none"


def save_figure(basename, dpi=300):
    for ext in ("png", "pdf", "svg"):
        out = FIGURES_DIR / f"{basename}.{ext}"
        plt.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
        if ext == "svg" and PLOT_FONT != FONT:
            t = Path(out).read_text(encoding="utf-8")
            for q in ('"', "'"):
                t = t.replace(f"font-family: {q}{PLOT_FONT}{q}", f"font-family: {FONT}")
            t = t.replace(f"font-family: {PLOT_FONT}", f"font-family: {FONT}")
            Path(out).write_text(t, encoding="utf-8")
        print(f"Saved: {out.name}")


MIN_CELLS = 10          # per population, matching the threshold used in Figure 3D
NORMAL_C, TUMOR_C = '#4DBBD5', '#E64B35'   # blue = non-malignant, red = malignant (as Figure 3D)

df = pd.read_csv(DATA_DIR / "external" / "swiftseq_april_persample_deid.csv")
df = df[(df.n_normal_PC >= MIN_CELLS) & (df.n_tumor_PC >= MIN_CELLS)]
# One row is already one participant per compartment at baseline, so no participant can enter the
# paired test twice. Assert rather than filter, so a regression in the input surfaces here.
assert not df.duplicated(['Patient_ID', 'Compartment']).any(), \
    "input carries more than one specimen per participant per compartment"
print(f"{len(df)} paired specimens from {df.Patient_ID.nunique()} participants "
      f"({int((df.n_libraries > 1).sum())} pooled from multiple libraries)")

COMPS = [('BM', 'Bone marrow'), ('PB', 'Peripheral blood')]
data, colors, labels, ps, pairs = [], [], [], [], []
for key, _ in COMPS:
    s = df[df.Compartment == key]
    data += [s.APRIL_normal.values, s.APRIL_tumor.values]
    colors += [NORMAL_C, TUMOR_C]
    labels += ['Normal', 'Tumor']
    pairs.append((s.APRIL_normal.values, s.APRIL_tumor.values))
    # Two-sided paired Wilcoxon signed-rank on the within-sample differences.
    ps.append(stats.wilcoxon(s.APRIL_tumor, s.APRIL_normal, alternative='two-sided').pvalue)
qs = multipletests(ps, method='fdr_bh')[1]
for (key, nice), p, q in zip(COMPS, ps, qs):
    s = df[df.Compartment == key]
    d = (s.APRIL_tumor - s.APRIL_normal).values
    print(f"{nice}: n={len(s)} paired samples, tumor higher {int((d>0).sum())}/{len(s)}, "
          f"median diff {np.median(d):+.4f}, paired Wilcoxon p={p:.3g}, q={q:.3g}")

pos = [1, 2, 3.6, 4.6]
# Narrow enough to sit beside Figure 3D on one page width; 3D is 6.0 in, this is 3.6 in.
fig, ax = plt.subplots(figsize=(3.6, 4.6))
bp = ax.boxplot(data, positions=pos, widths=0.6, patch_artist=True, showfliers=False,
                medianprops=dict(color='black', linewidth=1.8), zorder=2)
for patch, c in zip(bp['boxes'], colors):
    patch.set_facecolor(c); patch.set_alpha(0.55)
rng = np.random.default_rng(42)
jit = [rng.normal(pos[i], 0.08, size=len(v)) for i, v in enumerate(data)]
# The design is paired within sample, so each sample's two populations are joined by a light line.
for k, (nv, tv) in enumerate(pairs):
    xn, xt = jit[2 * k], jit[2 * k + 1]
    for i in range(len(nv)):
        # REVISION (reviewer 2, round 2): darkened from #9a9a9a/0.25/0.18, which was too faint
        # to follow. The pairing is the design of this panel, so the connectors must be visible.
        ax.plot([xn[i], xt[i]], [nv[i], tv[i]], color='#666666', linewidth=0.5,
                alpha=0.45, zorder=1)
for i, (vals, c) in enumerate(zip(data, colors)):
    ax.scatter(jit[i], vals, s=28, color=c, edgecolor='black', linewidth=0.5, alpha=0.92, zorder=3)

ax.set_xticks(pos)
ax.set_xticklabels([f"{lab}\nn={len(v)}" for lab, v in zip(labels, data)], fontsize=8.5)
ax.set_ylabel('APRIL-responsive module score\n(per-sample mean, plasma cells)', fontsize=9)
ax.tick_params(axis='y', labelsize=11)

ymin = min(np.min(v) for v in data); yr = max(np.max(v) for v in data) - ymin
YTOP, LVL = 0.50, 0.470
ax.set_ylim(ymin - 0.05 * yr, YTOP)
for (x1, x2), q in zip([(pos[0], pos[1]), (pos[2], pos[3])], qs):
    h = 0.018 * yr
    ax.plot([x1, x1, x2, x2], [LVL - h, LVL, LVL, LVL - h], color='black', linewidth=1.0)
    ax.text((x1 + x2) / 2, LVL + 0.4 * h, f"q={q:.2f}" if q >= 0.01 else f"q={q:.1e}",
            ha='center', va='bottom', fontsize=8)
for side in ('top', 'right'):
    ax.spines[side].set_visible(False)
ax.set_xlim(0.4, 5.2)
ax.set_title('Tumor vs normal plasma cells\n(Lightbody et al., phs003855)', fontsize=10, pad=6)
plt.tight_layout(); plt.subplots_adjust(bottom=0.20)
for (x1, x2), (_, nice) in zip([(pos[0], pos[1]), (pos[2], pos[3])], COMPS):
    ax.text((x1 + x2) / 2, -0.125, nice, ha='center', va='top', fontsize=9,
            transform=ax.get_xaxis_transform())
save_figure("Figure3E")
