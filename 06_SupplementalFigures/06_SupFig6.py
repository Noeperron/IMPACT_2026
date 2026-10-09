#!/usr/bin/env python3
"""Supplementary Figure 6: sample shipment and the IL-1beta response signature.

Purpose:      Supplementary Figure 6: shipment control for Figure 5C. Figure 5C is a
              within-participant paired test, so it is protected by design wherever a
              participant's two samples share a shipment status, which is the case for
              45 of the 52 paired participants. Both panels are restricted to treatment-naive
              SMM, the only Figure 5C group containing both shipped and non-shipped
              participants: every paired HD and every paired MGUS participant was shipped.
              (A) repeats the Figure 5C paired test separately in each shipment stratum;
              (B) asks directly whether the pre-to-post change differs by shipment. The
              underlying cohort matches Figure 5C exactly, with previously treated SMM and
              MM excluded.

Inputs:       H5AD_IL1B (scRNAseq_IMPACT_Zenodo.h5ad), data/il1b_response_genes_human.csv,
              data/hvg_2678_genes.txt, and data/metadata/Supplementary_Table_4_scRNAseq_sample_list.csv
              for Shipment_FedEx.

Outputs:      figures/SupFig6.png (and PDF + SVG); tables/SupFig6_stats.csv.

Dependencies: Python + scanpy, pandas, numpy, scipy, matplotlib; reads config.py and reuses the
              cell filter, effect size and statistics helpers of 05_Figure5/02_Figure5C.py so the
              two panels are computed identically.
"""
import sys
import gc
import warnings
warnings.filterwarnings('ignore')
import importlib.util
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
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

# Reuse the canonical Figure 5C helpers rather than restating them, so the cell-type filter,
# the untreated-SMM rule and the paired effect size cannot drift between the two figures.
_spec = importlib.util.spec_from_file_location(
    "fig5c", Path(__file__).parent.parent / "05_Figure5" / "02_Figure5C.py")
fig5c = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fig5c)

GROUPS = ['HD', 'MGUS', 'SMM']            # exactly Figure 5C's groups (the canonical gate)
PLOT_GROUP = 'SMM'                        # the only group with both shipment strata
ARMS = ['Shipped', 'Not shipped']
SHIP_MAP = {1.0: 'Shipped', 0.0: 'Not shipped'}
MIN_PAIRS = 3
RNG_SEED = 2026
# Published Figure 5C paired p-values. The pipeline below must reproduce them before any
# stratified result is trusted; a mismatch means the two figures no longer describe the same data.
CANONICAL_P = {'HD': 0.0208, 'MGUS': 0.0833, 'SMM': 0.6705}


def score_cells():
    """Reproduce the Figure 5C per-participant scores. Mirrors 02_Figure5C.py main()."""
    il1b_genes = pd.read_csv(fig5c.IL1B_GENES_PATH)['gene'].tolist()
    adata = sc.read_h5ad(H5AD_IL1B)
    import scipy.sparse as sp
    adata.X = sp.csr_matrix(adata.X.shape, dtype='float32')
    gc.collect()
    adata.layers['counts'] = adata.layers['counts'].astype('float32')
    gc.collect()
    hvg = [l.strip() for l in open(fig5c.HVG_GENES_PATH) if l.strip()]
    hvg = [g for g in hvg if g in adata.var_names]
    adata = adata[:, hvg].copy()
    adata.X = adata.layers['counts'].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    available = [g for g in il1b_genes if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=available, score_name='IL1B_response',
                      ctrl_size=min(100, len(available)), use_raw=False)

    tp_map = {'Pre-Vx': 'Pre-Vx', 'Post-2nd': 'Post-Vx', 'Post-1st': 'Post-1st'}
    cell_df = pd.DataFrame({
        'PatientID': adata.obs['Deidentified_Patient_ID'].astype(str).values,
        'Annotation_Level_2': adata.obs['Annotation_Level_2'].astype(str).values,
        'IL1B_response': adata.obs['IL1B_response'].values,
        'Disease': adata.obs['Diagnosis'].astype(str).values,
        'Timepoint_Clean': pd.Series(adata.obs['Timepoint'].astype(str).values).map(tp_map).values,
        'Treatment_Status': np.where(
            adata.obs['TreatmentStatus'].astype(str).values == 'Never_treated', 0.0, 1.0),
    })
    cell_df = cell_df[cell_df['Timepoint_Clean'].notna()].copy()
    cell_df = cell_df[~cell_df['Annotation_Level_2'].map(fig5c.should_exclude_celltype)].copy()
    del adata
    gc.collect()

    per = (cell_df.groupby(['PatientID', 'Timepoint_Clean', 'Disease', 'Treatment_Status'])
                  ['IL1B_response'].mean().reset_index())
    per.columns = ['Patient_ID', 'Timepoint_Clean', 'Disease', 'Treatment_Status', 'IL1B_score']
    per = fig5c.filter_untreated_smm(per)
    per = per[per['Disease'].isin(GROUPS) & per['Timepoint_Clean'].isin(['Pre-Vx', 'Post-Vx'])]
    return per


def attach_shipment(paired):
    meta_path = DATA_DIR / "metadata" / "Supplementary_Table_4_scRNAseq_sample_list.csv"
    if not meta_path.exists():
        meta_path = REPO_DIR.parent / "Supplementary_Tables" / "Supplementary_Table_4_scRNAseq_sample_list.csv"
    meta = pd.read_csv(meta_path)
    meta['tp'] = meta['Timepoint'].map({'Pre-Vx': 'Pre-Vx', 'Post-2nd': 'Post-Vx'})
    sh = (meta[meta['tp'].notna()]
          .groupby(['Deidentified_Patient_ID', 'tp'])['Shipment_FedEx']
          .agg(lambda s: s.dropna().unique().tolist()).reset_index())
    assert (sh['Shipment_FedEx'].map(len) <= 1).all(), "conflicting shipment records"
    sh['Shipping'] = sh['Shipment_FedEx'].map(lambda v: SHIP_MAP.get(v[0]) if len(v) == 1 else None)
    w = sh.pivot_table(index='Deidentified_Patient_ID', columns='tp', values='Shipping',
                       aggfunc='first').reset_index()
    w.columns.name = None
    w = w.rename(columns={'Deidentified_Patient_ID': 'Patient_ID',
                          'Pre-Vx': 'ship_pre', 'Post-Vx': 'ship_post'})
    out = paired.merge(w[['Patient_ID', 'ship_pre', 'ship_post']], on='Patient_ID', how='left')
    # A participant whose two samples differ in shipment status belongs to neither stratum.
    out['stratum'] = np.where(out['ship_pre'].notna() & out['ship_pre'].eq(out['ship_post']),
                              out['ship_pre'], 'Mixed or not recorded')
    return out


def main():
    per = score_cells()
    paired = per.pivot_table(index=['Patient_ID', 'Disease'], columns='Timepoint_Clean',
                             values='IL1B_score').reset_index().dropna(subset=['Pre-Vx', 'Post-Vx'])
    print(f"Paired participants: {len(paired)}")

    # Gate: this must be the Figure 5C cohort and the Figure 5C result.
    canon = fig5c.compute_statistics(paired, GROUPS)
    print("\nCanonical Figure 5C reproduction:")
    for _, r in canon.iterrows():
        print(f"  {r['Disease']} (n={int(r['n'])}): p={r['p_value']:.4f} q={r['q_value']:.4f}")
        assert abs(r['p_value'] - CANONICAL_P[r['Disease']]) < 5e-4, \
            f"{r['Disease']} p={r['p_value']:.4f} does not reproduce Figure 5C ({CANONICAL_P[r['Disease']]})"
    print("  gate passed: reproduces the published Figure 5C")

    d = attach_shipment(paired)
    d['delta'] = d['Post-Vx'] - d['Pre-Vx']
    print("\nParticipants by group and shipment stratum:")
    print(d.groupby(['Disease', 'stratum']).size().unstack(fill_value=0).to_string())

    # --- panel A: the paired test within each stratum ---------------------------------------
    rows = []
    for g in GROUPS:
        for arm in ARMS:
            s = d[(d['Disease'] == g) & (d['stratum'] == arm)]
            rec = dict(panel='A', Disease=g, stratum=arm, n=len(s), p=np.nan, effsize=np.nan)
            if len(s) >= MIN_PAIRS:
                rec['p'] = float(stats.wilcoxon(s['Post-Vx'], s['Pre-Vx']).pvalue)
                rec['effsize'] = fig5c.paired_effect_size(s['Post-Vx'].values, s['Pre-Vx'].values)
            rows.append(rec)
    stA = pd.DataFrame(rows)
    # Every group is computed for the record, but the correction family is the two contrasts the
    # panel actually shows; HD and MGUS have no non-shipped stratum, so they carry no contrast.
    plotted = stA['p'].notna() & (stA['Disease'] == PLOT_GROUP)
    stA.loc[plotted, 'q'] = multipletests(stA.loc[plotted, 'p'], method='fdr_bh')[1]

    # --- panel B: does the change itself depend on shipment? --------------------------------
    smm = d[(d['Disease'] == 'SMM') & (d['stratum'].isin(ARMS))]
    a = smm.loc[smm['stratum'] == 'Shipped', 'delta'].values
    b = smm.loc[smm['stratum'] == 'Not shipped', 'delta'].values
    u = stats.mannwhitneyu(a, b, alternative='two-sided')
    n1, n2 = len(a), len(b)
    z = (u.statistic - n1 * n2 / 2) / np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    rB = dict(panel='B', Disease='SMM', stratum='Shipped vs Not shipped', n=n1 + n2,
              p=float(u.pvalue), effsize=abs(z) / np.sqrt(n1 + n2), q=np.nan)
    st = pd.concat([stA, pd.DataFrame([rB])], ignore_index=True)
    print("\nStatistics:")
    print(st.round(4).to_string(index=False))
    tables_dir = REPO_DIR / "tables"
    tables_dir.mkdir(exist_ok=True)
    st.to_csv(tables_dir / "SupFig6_stats.csv", index=False)

    # --- plot -------------------------------------------------------------------------------
    PRE, POST = "#4682B4", "#EE5C42"      # the Figure 5C timepoint colours
    BOX_LINE = "#333333"
    rng = np.random.default_rng(RNG_SEED)
    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.5, wspace=0.3)

    vals = np.r_[d['Pre-Vx'].values, d['Post-Vx'].values]
    ymin, ymax = vals.min(), vals.max()
    span = max(ymax - ymin, 1e-9)

    for ci, arm in enumerate(ARMS):
            g = PLOT_GROUP
            ax = fig.add_subplot(gs[0, ci])
            s = d[(d['Disease'] == g) & (d['stratum'] == arm)]
            row = stA[(stA['Disease'] == g) & (stA['stratum'] == arm)].iloc[0]
            if len(s) == 0:
                ax.text(0.5, 0.5, 'no participants', ha='center', va='center',
                        fontsize=10, style='italic', color='#555555', transform=ax.transAxes)
                ax.set_xticks([]); ax.set_yticks([])
            else:
                pre, post = s['Pre-Vx'].values, s['Post-Vx'].values
                bp = ax.boxplot([pre, post], positions=[1, 2], widths=0.75, showfliers=False,
                                showcaps=False, patch_artist=True,
                                medianprops=dict(color=BOX_LINE, linewidth=2.14),
                                whiskerprops=dict(color=BOX_LINE, linewidth=1.07))
                for patch, colour in zip(bp["boxes"], [PRE, POST]):
                    patch.set_facecolor(colour); patch.set_alpha(0.5)
                    patch.set_edgecolor(BOX_LINE); patch.set_linewidth(1.07); patch.set_zorder(1)
                jp = rng.uniform(-0.15, 0.15, size=len(pre))
                jq = rng.uniform(-0.15, 0.15, size=len(pre))
                for i in range(len(pre)):
                    ax.plot([1 + jp[i], 2 + jq[i]], [pre[i], post[i]], color="black",
                            linewidth=0.3, alpha=0.5, zorder=2)
                ax.scatter(1 + jp, pre, s=26, facecolor=PRE, edgecolor="black", linewidth=0.5, zorder=3)
                ax.scatter(2 + jq, post, s=26, facecolor=POST, edgecolor="black", linewidth=0.5, zorder=3)
                if pd.notna(row['q']):
                    lab = ps.stat('q', row['q'], digits=3, sci_below=0.001)   # REVISION 2026-10-08
                    if row['q'] < 0.1:
                        lab += f", $r$={abs(row['effsize']):.2f}"
                    bar = ymax + 0.08 * span
                    ax.plot([1, 1, 2, 2], [bar - 0.02 * span, bar, bar, bar - 0.02 * span],
                            color="black", linewidth=0.9)
                    ax.text(1.5, bar + 0.03 * span, lab, ha="center", va="bottom", fontsize=11)
                else:
                    ax.text(1.5, ymax + 0.10 * span, f'$n$={len(s)}, too few\nto test', ha='center',
                            va='bottom', fontsize=8.5, style='italic', color='#555555')
                ax.set_xticks([1, 2]); ax.set_xticklabels(["Pre", "Post"], fontsize=12)
                ax.tick_params(axis="y", labelsize=11)
                ax.set_xlim(0.5, 2.5)
                ax.set_ylim(ymin - 0.05 * span, ymax + 0.32 * span)
            for side in ("top", "right", "bottom", "left"):
                ax.spines[side].set_visible(True); ax.spines[side].set_color("black")
            ax.set_title(f"Treatment-naive SMM, {arm.lower()}\n($n$={len(s)})", fontsize=12)
            if ci == 0 and len(s):
                ax.set_ylabel("IL-1β response score", fontsize=12)

    axB = fig.add_subplot(gs[1, 0])
    ARM_COLORS = {'Shipped': '#4682B4', 'Not shipped': '#EE5C42'}
    dv = [a, b]
    bp = axB.boxplot(dv, positions=[1, 2], widths=0.7, showfliers=False, patch_artist=True,
                     medianprops=dict(color=BOX_LINE, linewidth=2.0),
                     whiskerprops=dict(color=BOX_LINE, linewidth=1.07))
    for patch, arm in zip(bp["boxes"], ARMS):
        patch.set_facecolor(ARM_COLORS[arm]); patch.set_alpha(0.5)
        patch.set_edgecolor(BOX_LINE); patch.set_linewidth(1.07)
    for i, (v, arm) in enumerate(zip(dv, ARMS), start=1):
        axB.scatter(rng.normal(i, 0.08, size=len(v)), v, s=26, facecolor=ARM_COLORS[arm],
                    edgecolor="black", linewidth=0.5, zorder=3)
    dspan = max(np.r_[a, b].max() - np.r_[a, b].min(), 1e-9)
    ytop = np.r_[a, b].max() + 0.10 * dspan
    axB.plot([1, 1, 2, 2], [ytop - 0.02 * dspan, ytop, ytop, ytop - 0.02 * dspan],
             color="black", linewidth=0.9)
    axB.text(1.5, ytop + 0.03 * dspan, f"$p$={rB['p']:.2f}", ha="center", va="bottom", fontsize=11)
    axB.set_xticks([1, 2])
    axB.set_xticklabels([f"Shipped\n($n$={n1})", f"Not shipped\n($n$={n2})"], fontsize=11)
    axB.tick_params(axis="y", labelsize=11)
    axB.set_ylabel("Change in IL-1β response\n(post − pre)", fontsize=12)
    axB.set_title("Treatment-naive SMM", fontsize=12)
    axB.set_ylim(np.r_[a, b].min() - 0.08 * dspan, ytop + 0.16 * dspan)
    for side in ("top", "right", "bottom", "left"):
        axB.spines[side].set_visible(True); axB.spines[side].set_color("black")

    _tag = dict(fontsize=15, fontweight="normal",
                bbox=dict(boxstyle="square,pad=0.35", facecolor="white",
                          edgecolor="black", linewidth=1.0))
    fig.text(0.02, 0.975, "A", **_tag)
    fig.text(0.02, 0.475, "B", **_tag)
    for ext in ("png", "pdf", "svg"):
        out = FIGURES_DIR / f"SupFig6.{ext}"
        # REVISION 2026-10-08: transparent SVG (no page-sized white box in Inkscape).
        plt.savefig(out, dpi=300, bbox_inches="tight",
                    **({"transparent": True} if ext == "svg" else {"facecolor": "white"}))
        if ext == "svg" and PLOT_FONT != FONT:
            t = Path(out).read_text(encoding="utf-8")
            for q in ('"', "'"):
                t = t.replace(f"font-family: {q}{PLOT_FONT}{q}", f"font-family: {FONT}")
            t = t.replace(f"font-family: {PLOT_FONT}", f"font-family: {FONT}")
            Path(out).write_text(t, encoding="utf-8")
        print(f"Saved: {out.name}")


if __name__ == "__main__":
    main()
