#!/usr/bin/env python3
"""Supplementary Figure 1: Vaccine response by SMM 20/2/20 risk tier.

Design:
  - Panels A-C: tier reproduction of Fig 1B / 1C / 1F. Columns are HD / MGUS / LR-SMM /
    IR-SMM / HR-SMM / MM. HD, MGUS and MM use the FULL eligible ELISA cohort; MM are all
    treated, as in Figure 1. SMM is restricted to TREATMENT-NAIVE participants only.
    Previously treated SMM were removed because the treated fraction rose with risk tier
    (at peak post-2nd dose, 6/8 LR, 5/7 IR and 4/4 HR were treated), so any tier gradient
    was inseparable from a treatment gradient.
  - The 20/2/20 tier is computed IN PLACE from the three IMWG criteria carried in the
    ELISA tables themselves (M_spike >= 2 g/dL, FLC_ratio >= 20, BM_PC >= 20%): 0 factors
    LR, 1 factor IR, >=2 factors HR. A participant is classified only when all three
    components are present. This replaces the previous external bridge derived from a
    genomics reference cohort, which covered only 48 of 241 SMM and carried an assessment
    a median of about two years before vaccination. Computing in place roughly triples the
    treatment-naive n per tier and, for untreated participants, uses a contemporaneous
    rather than a historical assessment.
  - Continuous tumor-burden panels were dropped after the user judged the
    M-spike / BM PC% / FLC ratio scatters within treatment-naive SMM to be
    underpowered and inconclusive across metrics (M-spike rho=-0.21 p=0.28,
    BM PC% rho=+0.15 p=0.48, FLC ratio rho=-0.31 p=0.081); the numerical
    results are reported in the manuscript text and in the response letter only.
  - Cytogenetics panels are NOT included: copy-number calls linked to only 68 of the
    731 serology participants, and n_pos <= 2 per CNA among risk-classed SMM.
  - Fig 1E (waning slope) NOT reproduced -- per-tier serial-sample n = 2-6.
  - Age + sex adjustment via rank-based ANCOVA (Fig 1 main convention). JT ordered-trend
    test on age+sex residuals across HD < MGUS < LR-SMM < IR-SMM < HR-SMM; "Jonckheere-Terpstra"
    spelled out in full per manuscript convention. Statistics drawn IN-PANEL as brackets.

Purpose:      Supplementary Figure 1: vaccine response by SMM 20/2/20 risk tier, in treatment-naive SMM only. Three panels (Fig 1B/C/F equivalents) split treatment-naive SMM into LR/IR/HR using the IMWG 20/2/20 criteria computed from the M-spike, BM plasma cell percentage and free light chain ratio recorded in the ELISA tables; statistics use age+sex-adjusted rank-based ANCOVA + JT ordered-trend on HD<MGUS<LR<IR<HR.

Inputs:       data/elisa/elisa_spike_post2nd.csv, data/elisa/elisa_spike_post3rd.csv (both carry M_spike, BM_PC and FLC_ratio, from which the 20/2/20 tier is computed).

Outputs:      figures/SupFig1.png.

Dependencies: Python + pandas, numpy, matplotlib, scipy, statsmodels; reads config.py.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# REVISION 2026-10-08: italic statistical symbols, a×10^b instead of e-notation, transparent SVG.
import plot_style as ps
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

RNG_SEED = 2026

ELISA_2ND = DATA_DIR / 'elisa' / 'elisa_spike_post2nd.csv'
ELISA_3RD = DATA_DIR / 'elisa' / 'elisa_spike_post3rd.csv'

# IMWG 20/2/20 thresholds
MSPIKE_CUT, FLC_CUT, BMPC_CUT = 2.0, 20.0, 20.0

TIER_ORDER  = ['HD', 'MGUS', 'LRSMM', 'IRSMM', 'HRSMM', 'MM']
TIER_LABELS = {'HD':'HD', 'MGUS':'MGUS', 'LRSMM':'LR-SMM', 'IRSMM':'IR-SMM',
               'HRSMM':'HR-SMM', 'MM':'MM'}
TIER_COLORS = ['steelblue','orange','#FCD0A1','#EE5C42','#8B0000','#458B00']

# Helpers ---
def tier_2_20(df):
    """IMWG 20/2/20 tier from the three criteria carried in the ELISA tables.

    Counts risk factors (M-spike >= 2 g/dL, involved:uninvolved FLC ratio >= 20,
    bone marrow plasma cells >= 20%): 0 -> LRSMM, 1 -> IRSMM, >=2 -> HRSMM.
    Returns NaN unless all three components are present, so a participant is never
    classified on partial information.
    """
    v = df[['M_spike', 'BM_PC', 'FLC_ratio']].apply(pd.to_numeric, errors='coerce')
    n = ((v['M_spike'] >= MSPIKE_CUT).astype(float)
         + (v['FLC_ratio'] >= FLC_CUT).astype(float)
         + (v['BM_PC'] >= BMPC_CUT).astype(float))
    n[v.isna().any(axis=1)] = np.nan
    return n.map({0.0: 'LRSMM', 1.0: 'IRSMM', 2.0: 'HRSMM', 3.0: 'HRSMM'})

def make_panel_df(elisa_file, t_lo, t_hi, day_col):
    """Per-patient first-sample DataFrame.

    NOTE: SMM is restricted to TREATMENT-NAIVE participants (Ever_treated == 'No').
    Previously treated SMM are dropped because the treated fraction rose with risk tier,
    which made any tier gradient inseparable from a treatment gradient. HD, MGUS and MM
    use the full eligible cohort, matching Figure 1.
    """
    d = pd.read_csv(elisa_file)
    d = d[(d[day_col] >= t_lo) & (d[day_col] <= t_hi)].copy()
    d = d.sort_values(['Common_ID', day_col]).drop_duplicates('Common_ID', keep='first')
    d['Group'] = d['Disease'].map({'Healthy':'HD','MGUS':'MGUS','IgM-MGUS':'MGUS','MM':'MM','SMM':None})
    smm_naive = (d['Disease'] == 'SMM') & (d['Ever_treated'] == 'No')
    d.loc[smm_naive, 'Group'] = tier_2_20(d.loc[smm_naive])
    return d[d['Group'].isin(TIER_ORDER)]

def jt_trend(df, value_col, group_col, order, adjust_with):
    sub = df.dropna(subset=[value_col, group_col] + list(adjust_with)).copy()
    sub['rk'] = sub[value_col].rank()
    formula = 'rk ~ ' + ' + '.join([f'C({c})' if sub[c].dtype == object else c for c in adjust_with])
    sub['resid'] = smf.ols(formula, data=sub).fit().resid
    sub = sub[sub[group_col].isin(order)]
    groups = [sub[sub[group_col]==g]['resid'].values for g in order]
    if any(len(g) < 1 for g in groups) or sum(len(g) for g in groups) < 3:
        return float('nan'), float('nan'), 0
    U = 0
    for i in range(len(order)-1):
        for j in range(i+1, len(order)):
            for a in groups[i]:
                for b in groups[j]:
                    if   b > a: U += 1
                    elif b == a: U += 0.5
    n = [len(g) for g in groups]; N = sum(n)
    mu = (N*N - sum(ni*ni for ni in n)) / 4.0
    sigma2 = (N*N*(2*N+3) - sum(ni*ni*(2*ni+3) for ni in n)) / 72.0
    z = (U - mu) / np.sqrt(sigma2) if sigma2 > 0 else float('nan')
    p = 2 * stats.norm.sf(abs(z))
    return z, p, N

def adjusted_q_vs_hd(df, value_col, group_col, ref='HD'):
    sub = df.dropna(subset=[value_col, group_col, 'Age', 'Sex']).copy()
    sub['rk'] = sub[value_col].rank()
    sub[group_col] = pd.Categorical(sub[group_col], categories=[ref] + [g for g in TIER_ORDER if g != ref])
    m = smf.ols(f"rk ~ C({group_col}) + Age + C(Sex)", data=sub).fit()
    ps, names = [], []
    for g in TIER_ORDER:
        if g == ref: continue
        key = f"C({group_col})[T.{g}]"
        if key in m.pvalues.index:
            ps.append(m.pvalues[key]); names.append(g)
    if not ps: return {}
    q = multipletests(ps, method='fdr_bh')[1]
    return dict(zip(names, q))

def fmt_p(p):
    if pd.isna(p): return 'NA'
    # REVISION 2026-10-08: 3.9×10^-4 instead of 3.9e-04.
    return ps.num(p, digits=3, sci_below=0.001, sig=2)

def draw_bracket(ax, x1, x2, y, label, tip_frac=0.012, fontsize=7.2, linewidth=0.6,
                  text_x=None):
    """Draw a square bracket between x1 and x2 at height y with a label above.
    If text_x is provided, the label is positioned at that x-coordinate (default = midpoint).
    """
    ylim = ax.get_ylim(); tip = (ylim[1] - ylim[0]) * tip_frac
    ax.plot([x1, x1, x2, x2], [y - tip, y, y, y - tip],
            color='black', linewidth=linewidth, clip_on=False, solid_capstyle='butt')
    tx = (x1 + x2) / 2 if text_x is None else text_x
    ax.text(tx, y + tip * 0.6, label, ha='center', va='bottom',
            fontsize=fontsize, clip_on=False)

# ============================================================
# Tier panels A/B/C, the Figure 1B / 1C / 1F equivalents (treatment-naive SMM only)
# ============================================================
b = make_panel_df(ELISA_2ND, 14, 60,    'Days_post2nd')
c = make_panel_df(ELISA_2ND, 60.0001, 120, 'Days_post2nd')
f_ = make_panel_df(ELISA_3RD, 0, 100000, 'Days_post3rd')

def tier_panel_data(df):
    d = df.dropna(subset=['ELISA_Titer'])
    counts = {t: (d['Group']==t).sum() for t in TIER_ORDER}
    q = adjusted_q_vs_hd(d, 'ELISA_Titer', 'Group')
    z, p_jt, _ = jt_trend(d, 'ELISA_Titer', 'Group',
                          order=['HD','MGUS','LRSMM','IRSMM','HRSMM'],
                          adjust_with=['Age','Sex'])
    return {'df': d, 'counts': counts, 'q': q, 'jt_z': z, 'jt_p': p_jt}

tier_results = [
    ('A', tier_panel_data(b),  '2 weeks – 2 months\npost-2nd dose (Fig 1B)'),
    ('B', tier_panel_data(c),  '2 – 4 months\npost-2nd dose (Fig 1C)'),
    ('C', tier_panel_data(f_), 'Post-3rd dose (Fig 1F)'),
]

# ============================================================
# Figure: 1 row x 3 cols
# ============================================================
fig = plt.figure(figsize=(15.0, 6.5))
gs = fig.add_gridspec(1, 3, hspace=0.30, wspace=0.28,
                       left=0.06, right=0.985, top=0.85, bottom=0.20)

for col, (lbl, r, title) in enumerate(tier_results):
    ax = fig.add_subplot(gs[0, col])
    d = r['df']
    data = [d[d['Group']==t]['ELISA_Titer'].values for t in TIER_ORDER]
    bp = ax.boxplot(data, positions=range(len(TIER_ORDER)), widths=0.65,
                    patch_artist=True, showfliers=False)
    for patch, c0 in zip(bp['boxes'], TIER_COLORS):
        patch.set_facecolor(c0); patch.set_alpha(0.55)
    for el in ['whiskers','caps','medians']:
        plt.setp(bp[el], color='black', linewidth=1.0)
    for i, (v, c0) in enumerate(zip(data, TIER_COLORS)):
        if len(v):
            ax.scatter(np.random.default_rng(RNG_SEED+i).normal(i, 0.08, len(v)), v,
                       s=18, color=c0, edgecolor='white', linewidth=0.4, alpha=0.85, zorder=3)
    ax.set_xticks(range(len(TIER_ORDER)))
    xlab = []
    for t in TIER_ORDER:
        xlab.append(f"{TIER_LABELS[t]}\n($n$={r['counts'][t]})")
    ax.set_xticklabels(xlab, fontsize=8.5, rotation=30, ha='right')
    if col == 0: ax.set_ylabel(r'Spike IgG titer (OD$_{\mathrm{450–570}}$)', fontsize=10)
    ax.set_title(title, fontsize=10.5, fontweight='bold', pad=6)
    ax.spines[['top','right']].set_visible(False)

    # In-panel brackets: q vs HD (5 brackets, stacked) + JT trend (top bracket)
    flat_all = np.concatenate([v for v in data if len(v)])
    y_top = float(np.max(flat_all)); y_min = float(np.min(flat_all))
    step  = (y_top - y_min) * 0.16
    for k, g in enumerate(['MGUS','LRSMM','IRSMM','HRSMM','MM']):
        if g not in r['q']: continue
        y = y_top + (k + 1) * step
        draw_bracket(ax, 0, k + 1, y, f"$q$={fmt_p(r['q'][g])}", fontsize=9.0)
    jt_y = y_top + 6.4 * step
    draw_bracket(ax, 0, 4, jt_y,
                 f"Jonckheere-Terpstra (HD→HR-SMM, age+sex-adj):\n$z$={ps.num(r['jt_z'])}, $p$={fmt_p(r['jt_p'])}",
                 fontsize=9.0, text_x=2.5)
    ax.set_ylim(top=jt_y + 3.2 * step)

# Panel-letter labels
for lbl, spec in [('A', gs[0,0]),('B', gs[0,1]),('C', gs[0,2])]:
    pos = spec.get_position(fig)
    fig.text(pos.x0 - 0.012, pos.y1 + 0.005, lbl, fontsize=13, fontweight='bold')

fig.suptitle('Vaccine response by SMM 20/2/20 risk tier, treatment-naive SMM only',
             fontsize=11.5, fontweight='bold', y=0.96)

OUT_PNG = FIGURES_DIR / 'SupFig1.png'
ps.save(fig, 'SupFig1', FIGURES_DIR)   # REVISION 2026-10-08: transparent SVG, editable Arial text
print(f"\nSaved: {OUT_PNG}")

# Console dump ---
print("\n=== TIER PANELS (treatment-naive SMM only; 20/2/20 computed from M-spike, BM PC%, FLC ratio) ===")
for lbl, r, title in tier_results:
    print(f"  {lbl} ({title.splitlines()[0]}):")
    print(f"    counts:")
    for t in TIER_ORDER:
        print(f"      {TIER_LABELS[t]:8s} n={r['counts'][t]:3d}")
    q_items = ', '.join([f"{k}={v:.3g}" for k, v in r['q'].items()])
    print(f"    q vs HD (age+sex adj, BH): {{ {q_items} }}")
    print(f"    JT (HD→HR-SMM): z={r['jt_z']:.2f}, p={r['jt_p']:.4g}")
plt.close(fig)
