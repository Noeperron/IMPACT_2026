#!/usr/bin/env python3
"""Supplementary Figures 4A and 4B: validation of the APRIL-responsive gene signature (GSE173644).

Purpose:      External validation of the 15-gene APRIL-responsive module in an in-vitro APRIL
              stimulation time course of human B cells (Stephenson et al., J Immunol 2022),
              scoring the same module and comparing APRIL-treated with untreated cells.

Inputs:       data/external/GSE173644_timecourse.txt.gz, which ships in the Zenodo deposit.

Outputs:      figures/SupFig4AB.png (and PDF + SVG), heatmap and 120-min fold changes.

Dependencies: Python + pandas, numpy, matplotlib, seaborn, scipy, statsmodels; reads config.py.
"""
import sys
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
from matplotlib.patches import Patch

# REVISION 2026-10-08: Arial with editable SVG text and italic symbols (shared plot_style).
import plot_style as ps
plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.size': 11})

# APRIL-responsive signature (15 genes)
APRIL_GENES = ['MCL1', 'BCL2', 'NFKB1', 'NFKB2', 'RELB', 'NFKBIA', 'CD44', 'ICAM1',
               'CCL3', 'CCL4', 'VEGFA', 'IL10', 'CD274', 'TGFB1', 'CXCL8']

GENE_CATEGORIES = {
    'MCL1': 'Pro-survival', 'BCL2': 'Pro-survival',
    'NFKB1': 'NF-\u03baB', 'NFKB2': 'NF-\u03baB', 'RELB': 'NF-\u03baB', 'NFKBIA': 'NF-\u03baB',
    'CD44': 'Adhesion/Migration', 'ICAM1': 'Adhesion/Migration',
    'CCL3': 'Immunomodulatory', 'CCL4': 'Immunomodulatory', 'VEGFA': 'Immunomodulatory',
    'IL10': 'Immunomodulatory', 'CD274': 'Immunomodulatory', 'TGFB1': 'Immunomodulatory',
    'CXCL8': 'Immunomodulatory',
}

CAT_ORDER = ['Pro-survival', 'NF-\u03baB', 'Adhesion/Migration', 'Immunomodulatory']
CAT_COLORS = {
    'Pro-survival': '#E64B35', 'NF-\u03baB': '#4DBBD5',
    'Adhesion/Migration': '#00A087', 'Immunomodulatory': '#F39B7F',
}

# Gene display order (grouped by category)
GENE_ORDER = [g for cat in CAT_ORDER for g in APRIL_GENES if GENE_CATEGORIES[g] == cat]

# Load time course (GSE173644: 4 donors, APRIL stimulation at 0/30/60/120/360 min)
DATA_FILE = DATA_DIR / "external" / "GSE173644_timecourse.txt.gz"
with gzip.open(DATA_FILE, 'rt') as f:
    tc_raw = pd.read_csv(f, sep='\t', index_col=0, low_memory=False)
tc = tc_raw.drop('OfficialSymbol').drop(columns=[c for c in tc_raw.columns if 'Unnamed' in str(c)])
tc = tc.astype(float)
timepoints = [0]*4 + [30]*4 + [60]*4 + [120]*4 + [360]*4
tc.columns = [f't{tp}_D{d}' for tp, d in zip(timepoints, list(range(1, 5))*5)]
print(f"Time course: {tc.shape[0]:,} genes x {tc.shape[1]} samples")

# Differential expression at each timepoint (paired t-test vs baseline)
baseline_cols = [f't0_D{d}' for d in range(1, 5)]
results_list = []
for tp in [30, 60, 120, 360]:
    tp_cols = [f't{tp}_D{d}' for d in range(1, 5)]
    for gene in APRIL_GENES:
        bl = tc.loc[gene, baseline_cols].values
        tr = tc.loc[gene, tp_cols].values
        diff = tr - bl
        fc = diff.mean()
        sd = diff.std(ddof=1)
        se = sd / np.sqrt(4)
        if sd > 0:
            t_stat = fc / se
            p = 2 * stats.t.sf(abs(t_stat), df=3)
        else:
            t_stat, p, se = 0, 1.0, 0
        results_list.append({'gene': gene, 'timepoint': tp, 'log2FC': fc,
                             'baseline_mean': bl.mean(), 'treated_mean': tr.mean(),
                             't_stat': t_stat, 'p_value': p, 'se': se})

results = pd.DataFrame(results_list)
for tp in [30, 60, 120, 360]:
    mask = results['timepoint'] == tp
    _, padj, _, _ = multipletests(results.loc[mask, 'p_value'].values, method='fdr_bh')
    results.loc[mask, 'padj'] = padj
results['significant'] = (results['padj'] < 0.05) & (results['log2FC'] > 0)

# Primary timepoint: 120 min
res120 = results[results['timepoint'] == 120].copy()
n_sig = res120['significant'].sum()
print(f"120 min: {n_sig}/15 significant (FDR<0.05, positive FC)")

# Panel A: Heatmap (z-scored across timepoints) ---
heatmap_data = pd.DataFrame(index=GENE_ORDER, columns=[0, 30, 60, 120, 360], dtype=float)
for gene in GENE_ORDER:
    for tp in [0, 30, 60, 120, 360]:
        heatmap_data.loc[gene, tp] = tc.loc[gene, [f't{tp}_D{d}' for d in range(1, 5)]].mean()
heatmap_z = heatmap_data.sub(heatmap_data.mean(axis=1), axis=0).div(heatmap_data.std(axis=1), axis=0)

annot_mat = pd.DataFrame('', index=GENE_ORDER, columns=[0, 30, 60, 120, 360])
for _, row in results.iterrows():
    if row['gene'] in GENE_ORDER:
        if row['padj'] < 0.001:
            annot_mat.loc[row['gene'], row['timepoint']] = '***'
        elif row['padj'] < 0.01:
            annot_mat.loc[row['gene'], row['timepoint']] = '**'
        elif row['padj'] < 0.05:
            annot_mat.loc[row['gene'], row['timepoint']] = '*'

# Panel B: Bar chart at 120 min ---
res120_sorted = res120.sort_values('log2FC', ascending=True).copy()
colors_bar = [CAT_COLORS[GENE_CATEGORIES[g]] if sig else '#CCCCCC'
              for g, sig in zip(res120_sorted['gene'], res120_sorted['significant'])]

# Figure ---
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12, 6.5), gridspec_kw={'width_ratios': [1.1, 1]})

sns.heatmap(heatmap_z, annot=annot_mat.values, fmt='', cmap='RdBu_r', center=0,
            xticklabels=['0', '30', '60', '120', '360'], yticklabels=GENE_ORDER,
            linewidths=0.5, linecolor='white',
            cbar_kws={'label': '$z$-score (VST)', 'shrink': 0.7, 'aspect': 30},
            ax=ax_a, annot_kws={'fontsize': 9, 'fontweight': 'bold'})
ax_a.set_xlabel('Minutes post-APRIL stimulation', fontsize=12)
ax_a.set_ylabel('')
# REVISION 2026-10-08: boxed panel letters (the style of the assembled figures) and italic gene symbols.
ax_a.set_title('A', fontsize=15, fontweight='normal', loc='left', pad=10,
               bbox=dict(boxstyle='square,pad=0.35', facecolor='none', edgecolor='black', linewidth=1.0))
ax_a.tick_params(axis='y', labelsize=10)
for _t in ax_a.get_yticklabels():
    _t.set_fontstyle('italic')

ax_b.barh(range(len(res120_sorted)), res120_sorted['log2FC'], color=colors_bar,
          xerr=res120_sorted['se'], capsize=3, edgecolor='white', linewidth=0.5,
          error_kw={'linewidth': 1, 'color': '#555555'})
# REVISION 2026-10-08: each donor's 120-min log2 fold change as a dot on its bar (n=4 donors).
donor_fc = np.array([[tc.loc[g, f't120_D{d}'] - tc.loc[g, f't0_D{d}'] for d in range(1, 5)]
                     for g in res120_sorted['gene']])
for i, vals in enumerate(donor_fc):
    ax_b.scatter(vals, i + np.linspace(-0.18, 0.18, len(vals)), s=14, facecolor='white',
                 edgecolor='black', linewidth=0.7, zorder=4)
ax_b.set_yticks(range(len(res120_sorted)))
ax_b.set_yticklabels(res120_sorted['gene'], fontsize=10, fontstyle='italic')
ax_b.axvline(0, color='black', linewidth=0.8, zorder=0)
ax_b.set_xlabel('Log$_2$ fold change (120 min vs baseline)', fontsize=12, labelpad=6)
ax_b.set_title('B', fontsize=15, fontweight='normal', loc='left', pad=10,
               bbox=dict(boxstyle='square,pad=0.35', facecolor='none', edgecolor='black', linewidth=1.0))
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)

for i, (_, row) in enumerate(res120_sorted.iterrows()):
    padj = row['padj']
    star = '***' if padj < 0.001 else ('**' if padj < 0.01 else ('*' if padj < 0.05 else 'ns'))
    color = 'black' if star != 'ns' else '#999999'
    # REVISION 2026-10-08: stars sit beyond the error bar and the donor dots.
    if row['log2FC'] >= 0:
        xpos = max(row['log2FC'] + row['se'], donor_fc[i].max()) + 0.04
    else:
        xpos = min(row['log2FC'] - row['se'], donor_fc[i].min()) - 0.04
    ha = 'left' if row['log2FC'] >= 0 else 'right'
    ax_b.text(xpos, i, star, va='center', ha=ha, fontsize=9, fontweight='bold', color=color)

elements = [Patch(facecolor=CAT_COLORS[c], edgecolor='black', linewidth=0.5, label=c) for c in CAT_ORDER]
elements.append(Patch(facecolor='#CCCCCC', edgecolor='black', linewidth=0.5, label='Not significant'))
ax_b.legend(handles=elements, fontsize=8, framealpha=0.95, edgecolor='#CCCCCC',
            loc='lower right', handlelength=1.2, handletextpad=0.5)

fig.text(0.5, 0.98, 'Validation of APRIL-responsive gene signature (GSE173644, ' + ps.sym('n', bold=True) + '=4 donors)',
         ha='center', fontsize=13, fontweight='bold')
fig.text(0.5, 0.955, f'{n_sig}/15 genes significantly upregulated at 120 min (paired $t$-test, FDR < 0.05)',
         ha='center', fontsize=10, color='#555555')

plt.tight_layout(rect=[0, 0, 1, 0.94])
ps.save(fig, "SupFig4AB", FIGURES_DIR)   # REVISION 2026-10-08: vector PDF and SVG as well as PNG
plt.close()
