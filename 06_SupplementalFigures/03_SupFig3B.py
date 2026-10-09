#!/usr/bin/env python3
"""Supplementary Figure 3B: Figure 3C reproduced on shipped samples only.

Purpose:      Shipment control for Figure 3C. Figure 3C reports no
              difference in myeloid TNFSF13 (APRIL) expression between disease groups. Because
              shipment is unevenly distributed across those groups (every healthy donor and all
              but one participant with MGUS was shipped, whereas SMM is mixed), that null could
              in principle be produced by shipment rather than by biology. This panel repeats
              Figure 3C using only samples documented as shipped, holding shipment constant.
              The contrasts remain non-significant, so the negative result is not a shipment
              artifact.

Inputs:       H5AD_NORM (scRNAseq_IMPACT_Zenodo.h5ad) for TNFSF13 expression and cell labels;
              data/metadata/Supplementary_Table_4_scRNAseq_sample_list.csv for Shipment_FedEx
              and the age and sex covariates.

Outputs:      figures/SupFig3B.png (and PDF + SVG); tables/SupFig3B_stats.csv.

Dependencies: Python + scanpy, pandas, numpy, statsmodels, matplotlib; reads config.py.
"""
import sys
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
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

# REVISION 2026-10-08: italic statistical symbols, a×10^b instead of e-notation, transparent SVG.
import plot_style as ps


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


APRIL_GENE = 'TNFSF13'
DISEASE_GROUPS = ['HD', 'MGUS', 'SMM']
PLOT_GROUPS = ['HD', 'MGUS', 'SMM (Untreated)', 'SMM (Treated)']
RNG_SEED = 2026
# Shipment is read from Supplementary Table 4: 1 shipped, 0 not shipped. A blank field is EXCLUDED
# rather than assumed shipped, so a sample with no record never enters the shipped-only analysis.
SHIP_MAP = {1.0: 'Shipped', 0.0: 'Not shipped'}


def per_patient_expression():
    """Per-participant mean myeloid TNFSF13, using the cell selection of 02_Figure3C.py."""
    adata = sc.read_h5ad(H5AD_NORM, backed='r')
    april_idx = list(adata.var_names).index(APRIL_GENE)
    cells = pd.DataFrame({
        'Annotation_Level_1': adata.obs['Annotation_Level_1'].astype(str).values,
        'l2': adata.obs['Annotation_Level_2'].astype(str).values,
        'diagnosis': adata.obs['Diagnosis'].astype(str).values,
        'tp_raw': adata.obs['Timepoint'].astype(str).values,
        'treat': adata.obs['TreatmentStatus'].astype(str).values,
        'patient': adata.obs['Deidentified_Patient_ID'].astype(str).values,
    })
    cells['pos'] = np.arange(len(cells))
    tp_map = {'Pre-Vx': 'Pre-Vx', 'Post-2nd': 'Post-Vx'}
    clean = (~cells['l2'].isin(['QC_removed', 'CLL'])) & (~cells['l2'].str.startswith('db:'))
    keep = (clean
            & cells['Annotation_Level_1'].isin(['Mono', 'DC'])
            & cells['tp_raw'].isin(tp_map.keys())
            & cells['diagnosis'].isin(DISEASE_GROUPS))
    myeloid = cells[keep].copy()
    myeloid['Timepoint'] = myeloid['tp_raw'].map(tp_map)
    myeloid['group'] = np.where(
        myeloid['diagnosis'] != 'SMM', myeloid['diagnosis'],
        np.where(myeloid['treat'] == 'Never_treated', 'SMM (Untreated)', 'SMM (Treated)'))
    positions = myeloid['pos'].to_numpy()
    expr = []
    for i in range(0, len(positions), 100000):
        chunk = adata.X[list(positions[i:i + 100000]), april_idx]
        chunk = chunk.toarray().flatten() if hasattr(chunk, 'toarray') else np.asarray(chunk).flatten()
        expr.extend(chunk)
    myeloid['APRIL_expr'] = expr
    return myeloid.groupby(['patient', 'group', 'Timepoint'])['APRIL_expr'].mean().reset_index()


def contrasts(pat, label):
    """Figure 3C's test: age- and sex-adjusted rank ANCOVA vs HD, BH across the six contrasts."""
    rows = []
    for tp in ['Pre-Vx', 'Post-Vx']:
        t = pat[pat['Timepoint'] == tp].dropna(subset=['Age', 'Sex']).copy()
        t['rk'] = t['APRIL_expr'].rank()
        t['group'] = pd.Categorical(t['group'], categories=PLOT_GROUPS)
        fit = smf.ols('rk ~ C(group, Treatment(reference="HD")) + Age + C(Sex)', data=t).fit()
        for g in PLOT_GROUPS[1:]:
            rows.append(dict(cohort=label, Timepoint=tp, Comparison=f'HD vs {g}',
                             n_HD=int((t['group'] == 'HD').sum()), n=int((t['group'] == g).sum()),
                             p=float(fit.pvalues[f'C(group, Treatment(reference="HD"))[T.{g}]'])))
    out = pd.DataFrame(rows)
    out['q'] = multipletests(out['p'], method='fdr_bh')[1]
    return out


def attach_metadata(pat):
    """Join shipment status and the age and sex covariates from Supplementary Table 4."""
    meta_path = DATA_DIR / "metadata" / "Supplementary_Table_4_scRNAseq_sample_list.csv"
    if not meta_path.exists():
        meta_path = REPO_DIR.parent / "Supplementary_Tables" / "Supplementary_Table_4_scRNAseq_sample_list.csv"
    meta = pd.read_csv(meta_path)
    meta['tp'] = meta['Timepoint'].map({'Pre-Vx': 'Pre-Vx', 'Post-2nd': 'Post-Vx'})
    ship = (meta[meta['tp'].notna()]
            .groupby(['Deidentified_Patient_ID', 'tp'])['Shipment_FedEx']
            .agg(lambda s: s.dropna().unique().tolist()).reset_index())
    assert (ship['Shipment_FedEx'].map(len) <= 1).all(), "conflicting shipment records"
    ship['Shipping'] = ship['Shipment_FedEx'].map(lambda v: SHIP_MAP.get(v[0]) if len(v) == 1 else None)
    ship = ship.rename(columns={'Deidentified_Patient_ID': 'patient', 'tp': 'Timepoint'})
    pat = pat.merge(ship[['patient', 'Timepoint', 'Shipping']], on=['patient', 'Timepoint'], how='left')

    cov = (meta.groupby('Deidentified_Patient_ID')
                .agg(Age=('Age', 'median'),
                     Sex=('Sex', lambda s: s.dropna().mode().iloc[0] if s.dropna().size else None))
                .reset_index().rename(columns={'Deidentified_Patient_ID': 'patient'}))
    return pat.merge(cov, on='patient', how='left')


def main():
    pat = attach_metadata(per_patient_expression())

    shipped = pat[pat['Shipping'] == 'Shipped'].copy()
    print(f"Participant-timepoints: {len(pat)} in total, {len(shipped)} documented as shipped, "
          f"{int(pat['Shipping'].isna().sum())} with no record (excluded)")
    print("\nShipment by group (participant-timepoints):")
    print(pat.groupby(['group', 'Timepoint'])['Shipping'].value_counts(dropna=False)
             .unstack(fill_value=0).to_string())

    # Gate: the full-cohort run must reproduce the published Figure 3C before the restricted
    # run means anything. Figure 3C reports all six contrasts non-significant.
    full = contrasts(pat, 'all samples')
    print("\nFull cohort (reproduces Figure 3C):")
    print(full.round(4).to_string(index=False))
    assert full['q'].min() > 0.1, "full-cohort run does not reproduce the published null"

    rest = contrasts(shipped, 'shipped only')
    print("\nShipped samples only:")
    print(rest.round(4).to_string(index=False))
    print(f"\nminimum q, all samples {full['q'].min():.3f} | shipped only {rest['q'].min():.3f}")

    tables_dir = REPO_DIR / "tables"
    tables_dir.mkdir(exist_ok=True)
    pd.concat([full, rest], ignore_index=True).to_csv(tables_dir / "SupFig3B_stats.csv", index=False)

    # --- plot: Figure 3C's layout and styling, on shipped samples only ----------------------
    diag_colors = {'HD': '#3498db', 'MGUS': '#f1c40f',
                   'SMM (Untreated)': '#e74c3c', 'SMM (Treated)': '#A93226'}
    rng = np.random.default_rng(RNG_SEED)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
    plt.subplots_adjust(wspace=0.12)
    # y-limit is taken from the FULL cohort so this panel is on the same scale as Figure 3C.
    ymax_all = pat['APRIL_expr'].max()

    for ax_idx, (ax, tp) in enumerate(zip(axes, ['Pre-Vx', 'Post-Vx'])):
        tp_data = shipped[shipped['Timepoint'] == tp]
        box_data = [tp_data.loc[tp_data['group'] == g, 'APRIL_expr'].values for g in PLOT_GROUPS]
        bp = ax.boxplot(box_data, positions=range(len(PLOT_GROUPS)), widths=0.65,
                        patch_artist=True, showfliers=False)
        for patch, diag in zip(bp['boxes'], PLOT_GROUPS):
            patch.set_facecolor(diag_colors[diag])
            patch.set_alpha(0.7)
        for element in ['whiskers', 'caps', 'medians']:
            plt.setp(bp[element], color='black', linewidth=1.2)
        for i, (data, diag) in enumerate(zip(box_data, PLOT_GROUPS)):
            ax.scatter(rng.normal(i, 0.1, size=len(data)), data, alpha=0.6, color=diag_colors[diag],
                       edgecolor='white', s=45, zorder=3, linewidth=0.5)

        # Both p and q are printed. Benjamini-Hochberg is a step-up procedure, so the largest
        # p in the family sets a ceiling that every smaller q inherits: here all six q values
        # equal 0.9998. Printing q alone would show six identical numbers and hide the fact that
        # six different tests ran, and rounding 0.9998 to "1.00" would claim a precision the
        # value does not have, so a q at or above 0.995 is reported as ">0.99".
        def _one(sym, v):
            # A value of 0.9998 must not be printed as "1.00": that claims a precision it does
            # not have and reads as a placeholder. Anything that would round to 1.00 is reported
            # as ">0.99" instead.
            if v >= 0.995:
                return f'${sym}$>0.99'
            # REVISION 2026-10-08: italic p and q.
            return f'${sym}$={v:.2f}' if v >= 0.01 else f'${sym}$={v:.3f}'
        def _fmt(pv, qv):
            return f'{_one("p", pv)}, {_one("q", qv)}'
        for _j, _grp in enumerate(PLOT_GROUPS[1:], start=1):
            _row = rest[(rest['Timepoint'] == tp) & (rest['Comparison'] == f'HD vs {_grp}')].iloc[0]
            _y = ymax_all + 0.08 + 0.15 * (_j - 1)
            ax.plot([0, _j], [_y, _y], 'k-', linewidth=1)
            ax.text(_j / 2, _y + 0.02, _fmt(_row['p'], _row['q']), ha='center', va='bottom', fontsize=9.5)

        ax.set_title(tp, fontsize=14, fontweight='bold', pad=6)
        ax.set_xticks(range(len(PLOT_GROUPS)))
        _disp = {'SMM (Untreated)': 'SMM\nUntreated', 'SMM (Treated)': 'SMM\nTreated'}
        ax.set_xticklabels([f'{_disp.get(d, d)}\n($n$={len(v)})' for d, v in zip(PLOT_GROUPS, box_data)],
                           fontsize=9.5)
        ax.tick_params(axis='y', labelsize=11)
        if ax_idx == 0:
            # REVISION 2026-10-08: the expression values are ln(counts per 10,000 + 1).
            ax.set_ylabel('APRIL Expression\n[ln(CP10k+1)]', fontsize=12, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(-0.5, len(PLOT_GROUPS) - 0.5)
        ax.set_ylim(0, ymax_all + 0.70)

    plt.suptitle('APRIL Expression (Myeloid Cells), shipped samples only',
                 fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    save_figure("SupFig3B")


if __name__ == "__main__":
    main()
