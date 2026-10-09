#!/usr/bin/env python3
"""Supplementary Figure 4C: APRIL-responsive gene module score in non-malignant bone-marrow
plasma cells across NBM -> MGUS -> SMM (Boiarsky et al. 2022, GSE193531).

Purpose:      Published-cohort validation of the APRIL-responsive gene module in the bystander
              normal plasma cell compartment, the population the niche-depletion argument
              concerns. Malignant plasma cells, classified by inferCNV by the original authors,
              are excluded.

Statistics:   Jonckheere-Terpstra ordered trend across NBM -> MGUS -> SMM, plus Benjamini-
              Hochberg-corrected cross-stage contrasts versus NBM. Both are shown on the figure.

Inputs:       data/external/GSE193531_umi-count-matrix.csv.gz and
              data/external/GSE193531_cell-level-metadata.csv, the latter carrying the CNV-based
              malignant calls used here only to exclude malignant cells.

Outputs:      figures/SupFig4C.png (and PDF + SVG).

Dependencies: Python + scanpy, numpy, pandas, scipy, matplotlib; reads config.py.
"""
import sys
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import pandas as pd
import numpy as np
import scanpy as sc
from scipy import stats
from scipy.sparse import csr_matrix
from statsmodels.stats.multitest import multipletests
import matplotlib as mpl
import matplotlib.pyplot as plt
# Every figure in this paper is set in Arial. Register the metrically identical Liberation Sans
# as a fallback and rewrite the SVG font-family afterwards; svg.fonttype="none" keeps text editable.
import glob as _glob
from matplotlib import font_manager as _fm
FONT = "Arial"
_av = {f.name for f in _fm.fontManager.ttflist}
if FONT not in _av:
    for _p in _glob.glob("/usr/share/fonts/**/LiberationSans-*.ttf", recursive=True):
        _fm.fontManager.addfont(_p)
    _av = {f.name for f in _fm.fontManager.ttflist}
PLOT_FONT = FONT if FONT in _av else ("Liberation Sans" if "Liberation Sans" in _av else "sans-serif")
plt.rcParams["font.family"] = PLOT_FONT
plt.rcParams["svg.fonttype"] = "none"


def save_figure(basename, dpi=300):
    """Write PNG + PDF + SVG; rewrite the SVG font-family to Arial (see note above)."""
    for ext in ("png", "pdf", "svg"):
        out = FIGURES_DIR / f"{basename}.{ext}"
        # REVISION 2026-10-08: transparent SVG (no page-sized white box in Inkscape).
        plt.savefig(out, dpi=dpi, bbox_inches="tight",
                    **({"transparent": True} if ext == "svg" else {"facecolor": "white"}))
        if ext == "svg" and PLOT_FONT != FONT:
            t = Path(out).read_text(encoding="utf-8")
            for q in ('"', "'"):
                t = t.replace(f"font-family: {q}{PLOT_FONT}{q}", f"font-family: {FONT}")
            t = t.replace(f"font-family: {PLOT_FONT}", f"font-family: {FONT}")
            Path(out).write_text(t, encoding="utf-8")
        print(f"Saved: {out.name}")

import matplotlib.patches as mpatches
# REVISION 2026-10-08: italic statistical symbols through mathtext set in Arial by the shared
# plot_style, which keeps the labels editable text (was: no mathtext, dejavusans fallback).
import plot_style as ps

np.random.seed(42)

# Same 15-gene module as Figure 3D / SupFig 3
APRIL_TARGETS = {
    'Survival': ['MCL1', 'BCL2'],
    'NF-kB_pathway': ['NFKB2', 'RELB', 'NFKB1', 'NFKBIA'],
    'Adhesion': ['CD44', 'ICAM1'],
    'Immunomodulation': ['CCL3', 'CCL4', 'VEGFA', 'IL10', 'CD274', 'TGFB1', 'CXCL8'],
}
ALL_APRIL_GENES = [g for genes in APRIL_TARGETS.values() for g in genes]

EXT_DATA_DIR = DATA_DIR / "external"
MIN_CELLS_PER_SAMPLE = 10   # per-sample mean only included when >=10 PCs of that population in the sample

print("Loading GSE193531 (Boiarsky 2022)...")
metadata = pd.read_csv(EXT_DATA_DIR / "GSE193531_cell-level-metadata.csv", index_col=0)
print(f"  metadata: {len(metadata):,} cells")
counts = pd.read_csv(EXT_DATA_DIR / "GSE193531_umi-count-matrix.csv.gz", index_col=0)
print(f"  count matrix: {counts.shape[0]:,} genes x {counts.shape[1]:,} cells")

counts_T = counts.T
common = counts_T.index.intersection(metadata.index)
counts_T = counts_T.loc[common]; metadata = metadata.loc[common]
adata = sc.AnnData(X=csr_matrix(counts_T.values), obs=metadata,
                   var=pd.DataFrame(index=counts_T.columns))
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
available = [g for g in ALL_APRIL_GENES if g in adata.var_names]
print(f"  APRIL genes present: {len(available)}/{len(ALL_APRIL_GENES)}")
sc.tl.score_genes(adata, gene_list=available, score_name='APRIL_score', random_state=42)

obs = adata.obs[['sample_ID', 'disease_stage', 'normal_or_neoplastic']].copy()
obs['score'] = adata.obs['APRIL_score'].values
obs = obs.dropna(subset=['normal_or_neoplastic'])

# Per-(sample, neoplastic-status) aggregate; require >= MIN_CELLS_PER_SAMPLE
agg = (obs.groupby(['sample_ID', 'disease_stage', 'normal_or_neoplastic'])
          .agg(n_cells=('score', 'size'), mean=('score', 'mean'))
          .reset_index())
agg = agg[agg['n_cells'] >= MIN_CELLS_PER_SAMPLE].copy()

# Box order, left -> right along increasing plasma-cell tumor mass:
# NBM normal | MGUS normal | MGUS neo | SMM normal | SMM neo
# (NDMM dropped to match the cohort scope of Fig 3D, which has no treatment-naive MM. With only
# two neoplastic groups remaining (MGUS, SMM), the Jonckheere-Terpstra ordered-trend test cannot
# be computed for the neoplastic axis -- the test requires >=3 groups -- so only the JT across
# normal PCs (NBM -> MGUS -> SMM, 3 groups) is reported.)
# Box colors encode the population (blue shades = non-malignant normal PCs going light-to-dark
# along NBM -> SMM; red shades = malignant neoplastic PCs going light-to-dark along MGUS -> SMM).
# Same disease palette as Figure 3C and Figure 3D: NBM takes the HD colour, since it is the
# healthy-marrow comparator, so the same disease stage reads the same wherever it appears.
NORMAL_SHADES     = ['#3498db', '#f1c40f', '#e74c3c']   # NBM(=HD), MGUS, SMM
# Only non-malignant plasma cells are shown, which is the population the APRIL-niche argument
# concerns.
SPECS = [
    ('NBM',  'NBM',  'normal', NORMAL_SHADES[0], False),
    ('MGUS', 'MGUS', 'normal', NORMAL_SHADES[1], False),
    ('SMM',  'SMM',  'normal', NORMAL_SHADES[2], False),
]
NORMAL_BRACKET_COLOR = 'black'     # matches the bracket colour in manuscript Figure 3D
data, colors, hatches, labels = [], [], [], []
for lab, stage, status, color, hatch in SPECS:
    vals = agg.loc[(agg.disease_stage == stage) & (agg.normal_or_neoplastic == status), 'mean'].values
    data.append(vals); colors.append(color); hatches.append('///' if hatch else None); labels.append(lab)

# Statistics ----
def jonckheere_terpstra(groups_in_order):
    """Two-sided JT z, p; normal approximation, mid-rank ties."""
    g = [np.asarray(x) for x in groups_in_order]; k = len(g); J = 0.0
    for i in range(k):
        for j in range(i + 1, k):
            xi, xj = g[i], g[j]
            for v in xi:
                J += float(np.sum(v < xj)) + 0.5 * float(np.sum(v == xj))
    n = [len(x) for x in g]; N = sum(n)
    mu = (N * N - sum(ni * ni for ni in n)) / 4.0
    var = (N * N * (2 * N + 3) - sum(ni * ni * (2 * ni + 3) for ni in n)) / 72.0
    z = (J - mu) / np.sqrt(var); p = 2 * (1 - stats.norm.cdf(abs(z)))
    return z, p

norm_groups = [data[0], data[1], data[2]]   # NBM, MGUS, SMM  (3 groups for JT)
z_n, p_n = jonckheere_terpstra(norm_groups)
# Neoplastic JT not reported: only 2 groups remain (MGUS_neo, SMM_neo) after dropping NDMM,
# and the Jonckheere-Terpstra test requires at least 3 ordered groups.

# Cross-stage pairwise comparisons within each population:
# Normal PCs: NBM-vs-MGUS, NBM-vs-SMM (2 comparisons, BH/2).
# Neoplastic PCs: MGUS-vs-SMM only (single comparison; reported as raw p, no BH needed).
# Indices into `data` (NBM_n=0, MGUS_n=1, MGUS_neo=2, SMM_n=3, SMM_neo=4).
norm_cross = [('NBM-MGUS', 0, 1), ('NBM-SMM', 0, 2)]

def cross_qs(pairs):
    ps = []
    for _, i, j in pairs:
        a, b = data[i], data[j]
        ps.append(stats.mannwhitneyu(a, b, alternative='two-sided').pvalue
                  if (len(a) >= 2 and len(b) >= 2) else np.nan)
    valid_p = [p for p in ps if not np.isnan(p)]
    qs_valid = multipletests(valid_p, method='fdr_bh')[1] if valid_p else []
    qmap, k = {}, 0
    for (lbl, _, _), p in zip(pairs, ps):
        qmap[lbl] = (qs_valid[k] if not np.isnan(p) else np.nan)
        if not np.isnan(p): k += 1
    return qmap, ps

norm_cross_q, norm_cross_p = cross_qs(norm_cross)
print("\n  cross-stage NORMAL (BH/2):")
for (lbl, _, _), p in zip(norm_cross, norm_cross_p):
    print(f"    {lbl}: p={p:.4g}, q={norm_cross_q[lbl]:.4g}")

print(f"JT(normal across NBM -> SMM):     z={z_n:.3f}, p={p_n:.4g}  [primary figure statistic]")
print("counts per box (NBM, MGUS, SMM; non-malignant PCs):", [len(d) for d in data])

# Plot ----
# Positions: NBM at 1; MGUS pair at (3, 4.6); SMM pair at (6.8, 8.4). Within-pair gap 1.6,
# between-stage gap ~2.2 -- enough room to keep "MGUS normal" / "MGUS malignant" tick labels
# from colliding at the small panel size.
pos = [1, 2, 3]
# REVISION 2026-10-08: the healthy-marrow donors are labelled HD, as in the legend; italic n.
tick_labels = [f"{'HD' if lab == 'NBM' else lab}\n$n$={len(d)}" for lab, d in zip(labels, data)]

# Compact figure tuned for use as a single panel in a multi-panel figure: small physical size,
# but fonts/strokes sized to survive scale-down. Within-stage gray brackets removed (all ns at
# the available n; their absence reduces clutter without losing information that the JT and
# cross-stage brackets already convey).
fig, ax = plt.subplots(figsize=(5.0, 5.0))
bp = ax.boxplot(data, positions=pos, widths=0.6, patch_artist=True, showfliers=False,
                medianprops=dict(color='black', linewidth=1.8))
for patch, c, h in zip(bp['boxes'], colors, hatches):
    patch.set_facecolor(c); patch.set_alpha(0.55)
    if h: patch.set_hatch(h)
rng = np.random.default_rng(42)
for i, (vals, c) in enumerate(zip(data, colors)):
    if len(vals) == 0: continue
    ax.scatter(rng.normal(pos[i], 0.08, size=len(vals)), vals, s=28, color=c,
               edgecolor='black', linewidth=0.5, alpha=0.92)
ax.set_xticks(pos); ax.set_xticklabels(tick_labels, fontsize=11)
ax.set_ylabel('APRIL-responsive module score\n(per-sample mean, plasma cells)', fontsize=12)
ax.tick_params(axis='y', labelsize=11)

ymax_data = max(np.max(d) if len(d) else 0 for d in data)
ymin_data = min(np.min(d) if len(d) else 0 for d in data)
yr = ymax_data - ymin_data
# Bracket levels above the boxes (staggered so spanning brackets don't collide).
# Normal NBM-MGUS (1-3) and neoplastic MGUS-SMM (4-7) have no horizontal overlap -> same level.
# NBM-SMM (1-6) overlaps both -> one level higher.
norm_lvl  = [ymax_data + 0.10 * yr,  # NBM-MGUS  (short, level 1)
             ymax_data + 0.22 * yr]  # NBM-SMM   (spans the panel, level 2)
jt_norm_y = ymax_data + 0.42 * yr   # JT two-line label (same level as Figure 3D)
ax.set_ylim(ymin_data - 0.05 * yr, jt_norm_y + 0.18 * yr)

def fmt_q(q):
    if np.isnan(q): return "n/a"
    return ps.stat('q', q, digits=2, sci_below=0.01)   # REVISION 2026-10-08: 4.0×10^-3, not 4.0e-03

def draw_bracket(x1, x2, y, label, color, lw=1.0, fontsize=11):
    h = 0.018 * yr
    ax.plot([x1, x1, x2, x2], [y - h, y, y, y - h], color=color, linewidth=lw)
    ax.text((x1 + x2) / 2, y + 0.4 * h, label,
            ha='center', va='bottom', fontsize=fontsize, color=color)

# Cross-stage normal-PC brackets (dark blue; BH-corrected within the 2 normal comparisons)
for (lbl, ia, ib), y in zip(norm_cross, norm_lvl):
    draw_bracket(pos[ia], pos[ib], y, fmt_q(norm_cross_q[lbl]),
                 color=NORMAL_BRACKET_COLOR, lw=1.0, fontsize=11)

center_x = 2    # midpoint of the 3-box layout. Two-line JT label centered here.
# JT normal PC line: prefix in normal weight, stat (z, p) bold on a second line.
# REVISION 2026-10-08: z and p in bold italic, a true minus sign, and p as a×10^b (the legend
# reports p=2.3×10^-3). plot_style sets mathtext in Arial, so the label stays editable text.
ax.text(center_x, jt_norm_y + 0.055 * yr,
        "Jonckheere-Terpstra ordered trend (HD → SMM)",
        ha='center', va='center', fontsize=11.5, fontweight='normal')
ax.text(center_x, jt_norm_y - 0.055 * yr,
        f"{ps.sym('z', bold=True)} = {ps.num(z_n, bold=True)},  "
        f"{ps.sym('p', bold=True)} = {ps.num(p_n, sci_below=0.01, bold=True)}",
        ha='center', va='center', fontsize=11.5, fontweight='bold')

ax.set_title('Non-malignant BM plasma cells (Boiarsky GSE193531)',
             fontsize=12, pad=8, fontweight='bold')

# No legend: a single population (non-malignant PCs) is plotted.
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

plt.tight_layout()
save_figure("SupFig4C")
plt.close()

